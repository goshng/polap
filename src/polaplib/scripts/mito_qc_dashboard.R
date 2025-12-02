#!/usr/bin/env Rscript
# scripts/mito_qc_dashboard.R
# Version: 0.1.0
#
# Make a compact “3C/4C-style” mitochondrial assembly QC dashboard from HiMT TSV outputs.
#
# Inputs (positional):
#   1) contig_metrics_tsv      - output of contig_3C_metrics.py (two columns: metric\tvalue)
#   2) contig_table_tsv        - per-contig info
#                                Index, Contig name, Conserved gene number, GC content, Length
#   3) gene_completeness_tsv   - per-gene integrity table (or "NA" to skip)
#   4) output_file             - output figure (PDF/PNG)
#   5) expected_gene_count     - expected # conserved genes (optional, numeric)
#
# Example:
#   Rscript scripts/mito_qc_dashboard.R \
#     himt_contig_metrics.tsv \
#     himt_contigs.tsv \
#     himt_core_gene_integrity.tsv \
#     himt_mito_qc.pdf \
#     24
#
# where 24 is the # of conserved mitochondrial protein-coding genes used in HiMT
# (nad5, cox2, ccmB, ..., ccmFN).:contentReference[oaicite:3]{index=3}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: mito_qc_dashboard.R <contig_metrics_tsv> <contig_table_tsv> <gene_completeness_tsv|NA> <output_file> [expected_gene_count]\n",
      file = stderr())
  quit(status = 1)
}

contig_metrics_file    <- args[1]
contig_table_file      <- args[2]
gene_completeness_file <- args[3]
output_file            <- args[4]
expected_gene_count    <- if (length(args) >= 5) suppressWarnings(as.numeric(args[5])) else NA_real_

#-----------------------------
# 1. Read contig-level metrics
#-----------------------------
metrics_long <- read_tsv(
  contig_metrics_file,
  col_types = cols(
    metric = col_character(),
    value  = col_double()
  ),
  comment = "#"
)

metrics_wide <- metrics_long %>%
  mutate(metric = tolower(metric)) %>%
  pivot_wider(names_from = metric, values_from = value)

get_metric <- function(df, nm, default = NA_real_) {
  nm <- tolower(nm)
  if (nm %in% names(df)) {
    as.numeric(df[[nm]][1])
  } else {
    default
  }
}

num_contigs         <- get_metric(metrics_wide, "num_contigs")
total_length        <- get_metric(metrics_wide, "total_length")
n50                 <- get_metric(metrics_wide, "n50_from_awk",
                                 get_metric(metrics_wide, "n50"))
min_length          <- get_metric(metrics_wide, "min_length")
max_length          <- get_metric(metrics_wide, "max_length")
mean_gc             <- get_metric(metrics_wide, "mean_gc")
total_conserved_gen <- get_metric(metrics_wide, "total_conserved_genes")

#-----------------------------
# 2. Read per-contig table
#-----------------------------
contig_df <- read_tsv(
  contig_table_file,
  col_types = cols(.default = col_guess())
)

# normalise column names
names(contig_df) <- tolower(gsub("[^a-z0-9]+", "_", names(contig_df)))

if (!all(c("index", "contig_name", "conserved_gene_number", "gc_content", "length") %in% names(contig_df))) {
  warning("Contig table does not have expected columns; got: ",
          paste(names(contig_df), collapse = ", "))
}

#-----------------------------
# 3. Gene-set completeness
#-----------------------------
gene_set_score    <- NA_real_  # mean integrity (0–1)
gene_set_score_95 <- NA_real_  # fraction of genes with >=95% integrity
p_genes <- NULL

if (!is.na(gene_completeness_file) && tolower(gene_completeness_file) != "na") {
  gene_df <- read_tsv(
    gene_completeness_file,
    col_types = cols(.default = col_guess()),
    comment = "#"
  )

  names(gene_df) <- tolower(gsub("[^a-z0-9]+", "_", names(gene_df)))

  # Try to infer column names
  gene_col <- names(gene_df)[grepl("^gene$", names(gene_df))][1]
  if (is.na(gene_col)) {
    gene_col <- names(gene_df)[grepl("gene", names(gene_df))][1]
  }
  integrity_col <- names(gene_df)[grepl("integrity|identity|percent|pct", names(gene_df))][1]

  if (is.na(gene_col) || is.na(integrity_col)) {
    warning("Could not infer 'gene' and 'integrity' columns from gene_completeness_tsv; skipping gene-level metrics.")
  } else {
    gene_by_gene <- gene_df %>%
      transmute(
        gene      = .data[[gene_col]],
        integrity = as.numeric(.data[[integrity_col]])
      ) %>%
      group_by(gene) %>%
      summarise(max_integrity = max(integrity, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_integrity))

    gene_set_score    <- mean(gene_by_gene$max_integrity, na.rm = TRUE) / 100
    gene_set_score_95 <- mean(gene_by_gene$max_integrity >= 95, na.rm = TRUE)

    p_genes <- ggplot(gene_by_gene,
                      aes(x = reorder(gene, max_integrity),
                          y = max_integrity)) +
      geom_col() +
      coord_flip() +
      scale_y_continuous(labels = percent_format(scale = 1),
                         limits = c(0, 100)) +
      labs(x = "Conserved protein coding gene",
           y = "Best integrity per gene (%)",
           title = "Gene-set completeness") +
      theme_bw(base_size = 10)
  }
}

# Fallback: approximate completeness from total conserved genes
if (is.na(gene_set_score) && !is.na(total_conserved_gen) && !is.na(expected_gene_count)) {
  gene_set_score    <- total_conserved_gen / expected_gene_count
  gene_set_score_95 <- gene_set_score
}

#-----------------------------
# 4. Build 0–1 scores for radar
#-----------------------------
radar_metrics <- list()

if (!is.na(gene_set_score)) {
  radar_metrics[["Gene set completeness"]] <- max(min(gene_set_score, 1), 0)
}

if (!is.na(n50) && !is.na(total_length) && total_length > 0) {
  # ideal is N50 == total length (single perfect contig)
  radar_metrics[["N50 / total length"]] <- max(min(n50 / total_length, 1), 0)
}

if (!is.na(mean_gc)) {
  # scale GC into [0,1] assuming typical mito GC ~30–60%
  gc_score <- (mean_gc - 0.3) / (0.6 - 0.3)
  gc_score <- max(min(gc_score, 1), 0)
  radar_metrics[["GC content (scaled)"]] <- gc_score
}

if (!is.na(num_contigs) && num_contigs > 0) {
  # fewer contigs is better; perfect circular assembly => 1
  contig_score <- 1 / num_contigs
  contig_score <- max(min(contig_score, 1), 0)
  radar_metrics[["1 / # contigs"]] <- contig_score
}

if (length(radar_metrics) == 0) {
  stop("No metrics available to plot in radar chart.")
}

radar_df <- tibble::tibble(
  metric = factor(names(radar_metrics), levels = names(radar_metrics)),
  score  = as.numeric(radar_metrics)
)

p_radar <- ggplot(radar_df, aes(x = metric, y = score, group = 1)) +
  geom_polygon(fill = NA, linewidth = 0.8) +
  geom_point(size = 2) +
  coord_polar() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = percent_format(accuracy = 1)
  ) +
  labs(title = "Mitochondrial assembly quality (3C/4C-style)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title        = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(size = 9)
  )

#-----------------------------
# 5. Contig-level scatter
#-----------------------------
if (all(c("length", "conserved_gene_number", "gc_content") %in% names(contig_df))) {
  p_contig <- ggplot(contig_df,
                     aes(x = length,
                         y = conserved_gene_number,
                         colour = gc_content,
                         size = conserved_gene_number)) +
    geom_point(alpha = 0.8) +
    scale_x_continuous(labels = comma) +
    scale_colour_viridis_c(option = "C") +
    guides(size = "none") +
    labs(
      x = "Contig length (bp)",
      y = "# conserved genes",
      colour = "GC content",
      title = "Contig-level metrics"
    ) +
    theme_bw(base_size = 10)
} else {
  p_contig <- ggplot() + theme_void() +
    labs(title = "Contig-level metrics: contig table missing expected columns")
}

#-----------------------------
# 6. Assemble dashboard & save
#-----------------------------
if (!is.null(p_genes)) {
  dashboard <- p_radar + (p_genes / p_contig) + plot_layout(widths = c(1, 1))
} else {
  dashboard <- p_radar + p_contig + plot_layout(ncol = 2)
}

ggsave(
  filename = output_file,
  plot     = dashboard,
  width    = 10,
  height   = if (!is.null(p_genes)) 7 else 5
)

