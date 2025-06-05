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

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("-g", "--gfa"),
  action = "store",
  help = "GFA sequence part",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$gfa)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "3-gfa.seq.part.tsv")
  output1 <- file.path(input_dir0, "edges_stats.txt")
  args1 <- parse_args(parser, args = c("--gfa", input1, "-o", output1))
}

# x1 <- read_tsv(args1$gfa, col_names = FALSE, show_col_types = FALSE)

# The R script does this:
# The input file, 3-gfa.seq.part.tsv, has columns separated by tabs. I want to create an output file like contigs_stats.txt. The header of the output file is the same. The second column of the input file goes to the first column called #seq_name. The number at the 4th column in the input file goes to the 2nd column length of the output file. The number at the 5th column goes to the 3rd column, coverage, in the output file. These two numbers are extracted by the pattern LN:i:number for the 4th column in the input and dp:i:number for the 5th column in the input file. Most of the rest columns in the output file are just filled nominally or just the same values except for mult and graph_path columns among those of  circular        repeat  mult    telomere        alt_group       graph_path.    mult column values are computed by dividing the coverage value at the same row by the median of the coverage over all the rows. graph_path is the number extracted from the #seq_name column in the same row. The rest of column values are N N (mult value) both * (graph_path value).  For processing the input and creating the output, use R tidyverse package to write an R script.
#
# Key Steps:
# Reading the Input: Reads the input file and extracts the necessary columns.
# Pattern Matching: Extracts the LN:i:number and dp:i:number values for the length and coverage columns.
# Computing mult: Calculates the mult column as the ratio of coverage to the median coverage across all rows.
# Filling Columns: Fills the other columns as required, with fixed values or extracted numbers from seq_name.
# Output: Writes the output to a new file following the structure of contigs_stats.txt.

# Define file paths
input_file <- args1$gfa
output_file <- args1$out

# Read the input file
input_data <- read_tsv(args1$gfa, col_names = FALSE, show_col_types = FALSE)

# Extract columns with required patterns for the output file
output_data <- input_data %>%
  # Create the #seq_name column (from 2nd column of input)
  mutate(
    `#seq_name` = X2,
    # Extract 'LN:i:number' for length (4th column)
    length = as.numeric(str_extract(X4, "(?<=LN:i:)\\d+")),
    # Extract 'dp:i:number' for coverage (5th column)
    coverage = as.numeric(str_extract(X5, "(?<=dp:i:)\\d+"))
  ) %>%
  # Select relevant columns and keep the order for the output file
  select(`#seq_name`, length, coverage)

# Calculate the median of coverage
median_coverage <- median(output_data$coverage, na.rm = TRUE)

# Add the other columns to the output file
output_data <- output_data %>%
  mutate(
    circular = "N",
    "repeat" = "N",
    # Compute 'mult' as coverage/median_coverage
    mult = round(coverage / median_coverage),
    telomere = "both",
    alt_group = "*",
    # Extract graph_path value from the #seq_name
    graph_path = str_extract(`#seq_name`, "\\d+")
  )

# Write the final output data to the file
write_tsv(output_data, output_file, col_names = TRUE)
