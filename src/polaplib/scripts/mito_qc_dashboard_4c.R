#!/usr/bin/env Rscript
# scripts/mito_qc_dashboard_4c.R
# Version: v0.3.1
#
# 4C-style QC dashboard with splicing visualization:
#   - Radar plot of normalized 4C metrics
#   - Bar panel of selected raw metrics
#   - Contig scatter
#   - Splicing panel (cis vs trans, canonical vs non-canonical)
#
# Usage:
#   mito_qc_dashboard_4c.R 4c_metrics.tsv contig_table.tsv spliced_genes.tsv output.pdf
#

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: mito_qc_dashboard_4c.R <4c_metrics.tsv> <contig_table.tsv> <spliced_genes.tsv> <output.pdf>\n",
    file = stderr()
  )
  quit(status = 1)
}

metrics_file <- args[1]
contig_file <- args[2]
spliced_file <- args[3]
output_file <- args[4]

# ------------------ 4C metrics ------------------
m4c <- read_tsv(
  metrics_file,
  col_types = cols(
    category    = col_character(),
    metric      = col_character(),
    value       = col_double(),
    description = col_character()
  )
)

get_m <- function(cat, name) {
  m4c %>%
    filter(category == cat, metric == name) %>%
    pull(value) %>%
    {
      if (length(.) == 0) NA_real_ else .[1]
    }
}

geneset_comp <- get_m("completeness", "geneset_completeness_prop")
high_int_prop <- get_m("completeness", "high_integrity_gene_prop")
num_contigs <- get_m("contiguity", "num_contigs")
total_len <- get_m("contiguity", "total_length")
N50 <- get_m("contiguity", "N50")
frac_multi_copy <- get_m("consistency", "frac_multi_copy_genes")
frac_spliced <- get_m("consistency", "frac_spliced_genes")
frac_canon_spliced <- get_m("consistency", "frac_canonical_trans_spliced_detected")
circular_ratio <- get_m("consistency", "circular_contig_ratio")

# ------------------ Radar ------------------
radar <- tibble::tibble(
  metric_label = c(
    "Gene completeness",
    "High-integrity genes",
    "N50 / total length",
    "1 / # contigs",
    "Circular contigs ratio",
    "1 - frac multi-copy",
    "Canonical trans-spliced\nrecovery"
  ),
  score = c(
    geneset_comp,
    high_int_prop,
    ifelse(is.na(N50) || is.na(total_len) || total_len == 0,
      NA_real_, N50 / total_len
    ),
    ifelse(is.na(num_contigs) || num_contigs <= 0,
      NA_real_, min(1, 1 / num_contigs)
    ),
    circular_ratio,
    ifelse(is.na(frac_multi_copy), NA_real_, 1 - frac_multi_copy),
    frac_canon_spliced
  )
)

radar <- radar %>%
  mutate(score = pmin(pmax(score, 0), 1))

p_radar <- ggplot(radar, aes(x = metric_label, y = score, group = 1)) +
  geom_polygon(fill = NA, linewidth = 0.9) +
  geom_point(size = 2.2) +
  coord_polar() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = percent_format(accuracy = 1)
  ) +
  labs(title = "Mitochondrial assembly quality (4C + splicing)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title        = element_blank(),
    axis.text.x       = element_text(size = 8),
    panel.grid.major  = element_line(linewidth = 0.3),
    panel.grid.minor  = element_blank()
  )

# ------------------ Summary bar panel ------------------
summary_metrics <- m4c %>%
  filter(
    metric %in% c(
      "genes_total",
      "genes_present",
      "high_integrity_gene_prop",
      "num_contigs",
      "N50",
      "total_length",
      "fragmentation_index",
      "frac_multi_copy_genes",
      "frac_spliced_genes",
      "frac_canonical_trans_spliced_detected",
      "circular_contig_ratio"
    )
  ) %>%
  mutate(
    label = paste(category, metric, sep = " / "),
    label = factor(label, levels = unique(label))
  )

p_bar <- ggplot(summary_metrics, aes(x = label, y = value)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Value",
    title = "Selected 4C + splicing metrics"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7)
  )

# ------------------ Contig scatter ------------------
contig_df <- tryCatch(
  read_tsv(contig_file, col_types = cols(.default = col_guess())),
  error = function(e) NULL
)

if (!is.null(contig_df)) {
  names(contig_df) <- tolower(gsub("[^a-z0-9]+", "_", names(contig_df)))
  if (all(c("length", "conserved_gene_number", "gc_content") %in% names(contig_df))) {
    p_contig <- ggplot(
      contig_df,
      aes(
        x = length,
        y = conserved_gene_number,
        colour = gc_content
      )
    ) +
      geom_point(alpha = 0.8, size = 2) +
      scale_x_continuous(labels = comma) +
      scale_colour_viridis_c(option = "C") +
      labs(
        x = "Contig length (bp)",
        y = "# conserved genes",
        colour = "GC",
        title = "Contig-level structure"
      ) +
      theme_bw(base_size = 9)
  } else {
    p_contig <- ggplot() +
      theme_void() +
      labs(title = "Contig scatter: required columns missing")
  }
} else {
  p_contig <- ggplot() +
    theme_void() +
    labs(title = "Contig scatter: contig_table.tsv not available")
}

# ------------------ Splicing panel (robust) ------------------
spliced_df <- tryCatch(
  read_tsv(spliced_file, col_types = cols(.default = col_guess())),
  error = function(e) NULL
)

if (!is.null(spliced_df) && nrow(spliced_df) > 0) {
  spliced_clean <- spliced_df %>%
    mutate(
      splice_type = dplyr::case_when(
        splice_type %in% c("cis", "trans") ~ splice_type,
        TRUE ~ "none"
      ),
      splice_type = factor(splice_type,
        levels = c("cis", "trans", "none")
      ),
      canonical = canonical_trans_spliced == "TRUE"
    )

  # Stacked bar: cis vs trans by canonical/non-canonical
  spliced_counts <- spliced_clean %>%
    filter(splice_type %in% c("cis", "trans")) %>%
    count(canonical, splice_type)

  if (nrow(spliced_counts) > 0) {
    p_splice_bar <- ggplot(
      spliced_counts,
      aes(x = canonical, y = n, fill = splice_type)
    ) +
      geom_col(position = "stack") +
      scale_fill_brewer(palette = "Set2") +
      scale_x_discrete(labels = c("FALSE" = "non-canonical", "TRUE" = "canonical")) +
      labs(
        x = NULL,
        y = "# genes",
        fill = "splice type",
        title = "Splicing classification"
      ) +
      theme_bw(base_size = 9)
  } else {
    p_splice_bar <- ggplot() +
      theme_void() +
      labs(title = "Splicing classification: no cis/trans genes detected")
  }

  # Exon-count proxy histogram
  exon_df <- spliced_clean %>%
    filter(!is.na(min_exon_count_min), min_exon_count_min > 0)

  if (nrow(exon_df) > 0) {
    max_exon <- max(exon_df$min_exon_count_min, na.rm = TRUE)
    max_exon <- ifelse(is.finite(max_exon), max_exon, NA_real_)

    if (!is.na(max_exon) && max_exon >= 1) {
      p_exon <- ggplot(exon_df, aes(x = min_exon_count_min)) +
        geom_bar() +
        scale_x_continuous(breaks = seq(1, max_exon)) +
        labs(
          x = "Minimal exon count (proxy)",
          y = "# genes",
          title = "Exon-count proxy from cumulative integrity"
        ) +
        theme_bw(base_size = 9)
    } else {
      p_exon <- ggplot() +
        theme_void() +
        labs(title = "Exon-count: no multi-exon candidates")
    }
  } else {
    p_exon <- ggplot() +
      theme_void() +
      labs(title = "Exon-count: no genes with min_exon_count_min > 0")
  }

  p_splicing <- p_splice_bar / p_exon
} else {
  p_splicing <- ggplot() +
    theme_void() +
    labs(title = "Splicing: no spliced_genes.tsv or zero rows")
}

# ------------------ Layout & save ------------------
layout <- (p_radar | p_bar) / (p_contig | p_splicing)
ggsave(output_file, layout, width = 11, height = 8)
