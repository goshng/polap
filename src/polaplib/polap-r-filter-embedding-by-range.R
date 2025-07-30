#!/usr/bin/env Rscript

# Usage:
# Rscript filter_embedding_by_range.R input.tsv output.tsv [--x-min -5 --x-max 5 --y-min -2 --y-max 10 --invert]

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# === Define CLI options ===
option_list <- list(
  make_option(c("--x-min"), type = "double", default = -Inf),
  make_option(c("--x-max"), type = "double", default = Inf),
  make_option(c("--y-min"), type = "double", default = -Inf),
  make_option(c("--y-max"), type = "double", default = Inf),
  make_option(c("--invert"), action = "store_true", default = FALSE,
              help = "Invert selection (keep rows outside the range)")
)

# === Parse CLI args ===
parser <- OptionParser(usage = "%prog input.tsv output.tsv [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 2)
opt <- args$options
input_file <- args$args[1]
output_file <- args$args[2]

# === Load and rename ===
df <- read_tsv(input_file, col_types = cols())
colnames(df)[1:2] <- c("X", "Y")

# === Apply filter ===
if (opt$invert) {
  filtered <- df %>%
    filter(X < opt$`x-min` | X > opt$`x-max` | Y < opt$`y-min` | Y > opt$`y-max`)
} else {
  filtered <- df %>%
    filter(X >= opt$`x-min`, X <= opt$`x-max`, Y >= opt$`y-min`, Y <= opt$`y-max`)
}

# === Save ===
write_tsv(filtered, output_file)
cat("✅ Filtered", nrow(filtered), ifelse(opt$invert, "(inverted)", ""), "rows →", output_file, "\n")
