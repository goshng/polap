#!/usr/bin/env Rscript
###############################################################################
# polap-r-table-s1-histograms.R
#
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose:
#   Draw histograms for Table S1 metrics:
#     - Bases (Gb)
#     - AvgLen (kb)
#     - AvgQ
#
# Usage:
#   Rscript polap-r-table-s1-histograms.R \
#     man/md/v2-0-auto-table-s1.tsv \
#     man/md/v2-0-auto-table-s1-hist.pdf
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript polap-r-table-s1-histograms.R <table.tsv> <out.pdf>")
}

tsv_file <- args[1]
pdf_file <- args[2]

# ---- Libraries ----
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(stringr)
})

# ---- Read data ----
df <- read_tsv(
  tsv_file,
  show_col_types = FALSE,
  progress = FALSE
)

# ---- Clean columns ----
# Remove markdown italics from species
df$Species <- str_replace_all(df$Species, "\\*", "")

# Ensure numeric types
df$`Bases (Gb)` <- as.numeric(df$`Bases (Gb)`)
df$`AvgLen (kb)` <- as.numeric(df$`AvgLen (kb)`)
df$AvgQ <- as.numeric(df$AvgQ)

# ---- Open PDF ----
pdf(pdf_file, width = 7, height = 6)

# ---- Histogram: Bases (Gb) ----
p_bases <- ggplot(df, aes(x = `Bases (Gb)`)) +
  geom_histogram(bins = 15, color = "black") +
  theme_bw(base_size = 14) +
  labs(
    # title = "Distribution of total bases",
    title = NULL,
    x = "Bases (Gb)",
    y = "Number of datasets"
  )

print(p_bases)

# ---- Histogram: AvgLen (kb) ----
p_avglen <- ggplot(df, aes(x = `AvgLen (kb)`)) +
  geom_histogram(bins = 15, color = "black") +
  theme_bw(base_size = 14) +
  labs(
    # title = "Distribution of average read length",
    title = NULL,
    x = "Average read length (kb)",
    y = "Number of datasets"
  )

print(p_avglen)

# ---- Histogram: AvgQ ----
p_avgq <- ggplot(df, aes(x = AvgQ)) +
  geom_histogram(binwidth = 1, color = "black") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    # title = "Distribution of average read quality",
    title = NULL,
    x = "Average quality",
    y = "Number of datasets"
  )

print(p_avgq)

# ---- Close PDF ----
dev.off()

message("[INFO] Wrote histogram PDF: ", pdf_file)
