#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tools) # for file_path_sans_ext
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || "-o" %in% args == FALSE) {
  stop("Usage: Rscript run.R file1.tsv file2.tsv ... -o output.pdf")
}

# Extract output file
o_idx <- which(args == "-o")
tsv_files <- args[1:(o_idx - 1)]
output_file <- args[o_idx + 1]

# Read and merge data
merged_df <- NULL
for (f in tsv_files) {
  df <- read.table(f, header = TRUE, sep = "\t")
  sample_name <- file_path_sans_ext(basename(f))
  df <- df %>% rename(!!sample_name := alpha)
  if (is.null(merged_df)) {
    merged_df <- df
  } else {
    merged_df <- merge(merged_df, df, by = "index")
  }
}

# Convert to long format for ggplot
long_df <- merged_df %>%
  pivot_longer(-index, names_to = "Sample", values_to = "Alpha")

# Plot
p <- ggplot(long_df, aes(x = index, y = Alpha, color = Sample)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Alpha values over index", x = "Index", y = "Alpha")

ggsave(output_file, plot = p, width = 7, height = 5)
