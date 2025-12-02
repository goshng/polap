#!/usr/bin/env Rscript
# File: scripts/plot_mt_assembly_trend.R
# Version: v0.1.0
#
# Read assembly metadata TSV and produce:
#  - plant_mt.assembly_trend_counts.tsv
#  - plant_mt.assembly_trend.pdf
#
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: plot_mt_assembly_trend.R --meta meta.tsv --out-pdf fig.pdf --out-tsv counts.tsv\n", file = stderr())
  quit(status = 1)
}

meta_file <- NULL
out_pdf <- NULL
out_tsv <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--meta") {
    meta_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--out-pdf") {
    out_pdf <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--out-tsv") {
    out_tsv <- args[i + 1]; i <- i + 2
  } else {
    stop("Unknown argument: ", args[i])
  }
}

if (is.null(meta_file) || is.null(out_pdf) || is.null(out_tsv)) {
  stop("Both --meta, --out-pdf, and --out-tsv are required.")
}

dat <- read_tsv(meta_file, show_col_types = FALSE)

# Clean year (only keep plausible years)
dat <- dat %>%
  mutate(
    year = suppressWarnings(as.integer(year))
  ) %>%
  filter(!is.na(year), year >= 1990, year <= as.integer(format(Sys.Date(), "%Y")))

# For nicer plots, collapse rare assemblers into "Other"
common_assemblers <- c("Flye", "Canu", "hifiasm", "NOVOPlasty", "GetOrganelle", "SPAdes", "MitoZ")
dat <- dat %>%
  mutate(
    assembler2 = ifelse(assembler %in% common_assemblers, assembler, "Other"),
    seq_platform2 = ifelse(seq_platform %in% c("Illumina", "PacBio CLR", "PacBio HiFi", "ONT", "MGI"),
                           seq_platform, "Other/Unknown")
  )

# Aggregate counts by year, platform, and assembler
counts <- dat %>%
  group_by(year, seq_platform2, assembler2) %>%
  summarise(n_genomes = n(), .groups = "drop")

write_tsv(counts, out_tsv)

# Main trend plot:
#   x = year, y = n_genomes, fill = seq_platform2, facet by assembler2
p <- ggplot(counts, aes(x = year, y = n_genomes, fill = seq_platform2)) +
  geom_col(position = "stack") +
  facet_wrap(~ assembler2, scales = "free_y") +
  scale_x_continuous(breaks = pretty) +
  scale_fill_brewer(palette = "Set2", name = "Seq platform") +
  labs(
    title = "Trends in sequencing and assembly methods\nfor plant mitochondrial genomes (NCBI)",
    x = "Year",
    y = "Number of mitochondrial genomes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave(out_pdf, p, width = 10, height = 7)


