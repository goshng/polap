#!/usr/bin/env Rscript
# File: scripts/plot_mt_meta_overview.R
# Version: v0.1.0
#
# Generate overview plots and summary tables for plant mt meta-analysis:
#  - Year x sequencing platform (from GenBank COMMENT)
#  - Year x assembler
#  - Data availability (BioProject / SRA presence) over time
#
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  cat("Usage: plot_mt_meta_overview.R --meta meta.tsv --sra sra_runinfo.csv --outdir outdir\n", file = stderr())
  quit(status = 1)
}

meta_file <- NULL
sra_file <- NULL
outdir <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--meta") {
    meta_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--sra") {
    sra_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--outdir") {
    outdir <- args[i + 1]
    i <- i + 2
  } else {
    stop("Unknown argument: ", args[i])
  }
}

if (is.null(meta_file) || is.null(outdir)) {
  stop("Both --meta and --outdir are required (sra can be empty).")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------- Read assembly_meta --------
meta <- read_tsv(meta_file, show_col_types = FALSE)

meta <- meta %>%
  mutate(
    year = suppressWarnings(as.integer(year)),
    has_bioproject = bioproject != "" & !is.na(bioproject),
    has_sra = sra_accessions != "" & !is.na(sra_accessions)
  ) %>%
  filter(!is.na(year), year >= 1990, year <= as.integer(format(Sys.Date(), "%Y")))

# -------- Year x platform --------
meta <- meta %>%
  mutate(
    seq_platform2 = ifelse(
      seq_platform %in% c("Illumina", "PacBio CLR", "PacBio HiFi", "ONT", "MGI"),
      seq_platform,
      "Other/Unknown"
    )
  )

platform_counts <- meta %>%
  group_by(year, seq_platform2) %>%
  summarise(n_genomes = n(), .groups = "drop")

write_tsv(platform_counts, file.path(outdir, "mt_meta_platform_year_counts.tsv"))

p_platform <- ggplot(platform_counts, aes(x = year, y = n_genomes, fill = seq_platform2)) +
  geom_col(position = "stack") +
  scale_x_continuous(breaks = pretty) +
  labs(
    title = "Plant mitochondrial genomes: sequencing platform over time\n(GenBank COMMENT/Assembly-Data)",
    x = "Year",
    y = "Number of mitogenomes",
    fill = "Platform"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(outdir, "mt_meta_platform_year.pdf"), p_platform, width = 9, height = 6)

# -------- Year x assembler --------
common_assemblers <- c("Flye", "Canu", "hifiasm", "NOVOPlasty", "GetOrganelle", "SPAdes", "MitoZ")
meta <- meta %>%
  mutate(
    assembler2 = ifelse(assembler %in% common_assemblers, assembler, "Other/Unknown")
  )

asm_counts <- meta %>%
  group_by(year, assembler2) %>%
  summarise(n_genomes = n(), .groups = "drop")

write_tsv(asm_counts, file.path(outdir, "mt_meta_assembler_year_counts.tsv"))

p_asm <- ggplot(asm_counts, aes(x = year, y = n_genomes, fill = assembler2)) +
  geom_col(position = "stack") +
  scale_x_continuous(breaks = pretty) +
  labs(
    title = "Plant mitochondrial genomes: assembly methods over time\n(GenBank COMMENT/Assembly-Data)",
    x = "Year",
    y = "Number of mitogenomes",
    fill = "Assembler"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(outdir, "mt_meta_assembler_year.pdf"), p_asm, width = 9, height = 6)

# -------- Data availability over time --------
availability <- meta %>%
  group_by(year) %>%
  summarise(
    n_total = n(),
    n_bioproject = sum(has_bioproject),
    n_sra = sum(has_sra),
    .groups = "drop"
  ) %>%
  mutate(
    frac_bioproject = ifelse(n_total > 0, n_bioproject / n_total, NA_real_),
    frac_sra = ifelse(n_total > 0, n_sra / n_total, NA_real_)
  )

write_tsv(availability, file.path(outdir, "mt_meta_data_availability.tsv"))

p_avail <- ggplot(availability, aes(x = year)) +
  geom_line(aes(y = frac_bioproject, color = "BioProject linked")) +
  geom_point(aes(y = frac_bioproject, color = "BioProject linked"), size = 1) +
  geom_line(aes(y = frac_sra, color = "SRA listed")) +
  geom_point(aes(y = frac_sra, color = "SRA listed"), size = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("BioProject linked" = "black", "SRA listed" = "darkred")) +
  labs(
    title = "Data availability for plant mitochondrial genomes over time",
    x = "Year",
    y = "Fraction of mitogenomes",
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(outdir, "mt_meta_availability.pdf"), p_avail, width = 9, height = 6)

# -------- Optional: minimal SRA usage summary (if any) --------
if (!is.null(sra_file) && file.exists(sra_file) && file.info(sra_file)$size > 0) {
  sra <- read_csv(sra_file, show_col_types = FALSE)
  # Common SRA RunInfo columns: Run, LibraryStrategy, Platform, Model, BioProject
  sra_summary <- sra %>%
    group_by(Platform, Model) %>%
    summarise(n_runs = n(), .groups = "drop") %>%
    arrange(desc(n_runs))
  write_tsv(sra_summary, file.path(outdir, "mt_meta_sra_platform_summary.tsv"))
}
