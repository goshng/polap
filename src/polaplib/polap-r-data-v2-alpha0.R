#!/usr/bin/env Rscript

################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# This script creates one of the two different figures used in polap's
# subsampling-based approach.
#
# Plot 1.
# Line plot of read-coverage thresholds versus subsample size index in Stage 1
# over different initial values of the read-coverage threshold of
# the subsampling-based assemblies for the _Eucalyptus pauciflora_ dataset.
# Plot 2.
# Line plot of read-coverage thresholds versus subsample size index in Stage 1
# over different different increment sizes of the subsampling-based assemblies
# for the _Eucalyptus pauciflora_ dataset.
#
# Example:
# Rscript ./src/polaplib/polap-r-data-v2-alpha0.R \
#   input/0.25.tsv input/0.50.tsv -l alpha0 -o output/alpha0.pdf
# Rscript ./src/polaplib/polap-r-data-v2-alpha0.R \
#   input/?.??.tsv -l delta -o output/alpha0.pdf
#
# See Also:
# man-figure-alpha_genus_species
# man-figure-delta_genus_species
#
# TODO: rename: polap-r-data-v2-alpha0.R
#
# Check: change the command-line processing
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tools)
})

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for required options
if (length(args) < 3 || !("-o" %in% args)) {
  # stop("Usage: Rscript run.R file1.tsv file2.tsv ... [-l legend_title] -o output.pdf")
  tsv_files <- c("0.00.tsv", "1.00.tsv", "2.00.tsv")
  output_file <- "output.pdf"
  legend_title <- "alpha" # default
} else {
  # Find indices of options
  o_idx <- which(args == "-o")
  l_idx <- which(args == "-l")

  # Parse output file
  output_file <- args[o_idx + 1]

  # Parse legend title
  legend_title <- "Sample" # default
  if (length(l_idx) == 1) {
    legend_title <- args[l_idx + 1]
  }

  # Determine tsv files
  option_indices <- sort(c(o_idx, o_idx + 1, l_idx, l_idx + 1))
  tsv_files <- args[-option_indices]
}

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
  labs(title = "", x = "Index", y = "Alpha", color = legend_title)
p
ggsave(output_file, plot = p, width = 7, height = 5)
