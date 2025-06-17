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
# This script creates edges_stats.txt, the edge version of contigs_stats.txt.
# It is used by these two almost the same functions:
# function _run_polap_edges-stats { # create an edge version of
# contigs_stats.txt
# polap_edges-stats() {
# in run-polap-function-annotate.sh.
#
# It follows a similar procedure used one by polap-r-depthfilter-gfa.R.
# In other words, this script uses the sequence part of a gfa assembly graph.
# step 3-1: creating GFA without sequence data
#   convert gfa to another gfa2 using gfatools view command.
# step 3-2: extracting sequence part of GFA
#   extract the sequence part of the gfa2 using grep ^S command.
# Then, this script is called to create the edges_stats.txt.
#
# Example:
# polap edges-stats
# polap edges-stats view
# Rscript ./src/polaplib/run-polap-r-edges-stats.R \
#   --gfa output/3-gfa.seq.part.tsv \
#   --out output/edges_stats.txt
#
# TEST-SCC: done in test/test.sh
#
# TODO: rename: polap-r-edges-stats.R
#
# Check: 2025-06-17
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
  input_dir0 <- file.path("input")
  input_dir1 <- file.path("output")
  input1 <- file.path(input_dir0, "3-gfa.seq.part.tsv")
  output1 <- file.path(input_dir1, "edges_stats.txt")
  args1 <- parse_args(parser, args = c("--gfa", input1, "-o", output1))
}

# x1 <- read_tsv(args1$gfa, col_names = FALSE, show_col_types = FALSE)

# TAG: SCC delete
# The R script does this:
# The input file, 3-gfa.seq.part.tsv, has columns separated by tabs. I want to create an output file like contigs_stats.txt. The header of the output file is the same. The second column of the input file goes to the first column called #seq_name. The number at the 4th column in the input file goes to the 2nd column length of the output file. The number at the 5th column goes to the 3rd column, coverage, in the output file. These two numbers are extracted by the pattern LN:i:number for the 4th column in the input and dp:i:number for the 5th column in the input file. Most of the rest columns in the output file are just filled nominally or just the same values except for mult and graph_path columns among those of  circular        repeat  mult    telomere        alt_group       graph_path.    mult column values are computed by dividing the coverage value at the same row by the median of the coverage over all the rows. graph_path is the number extracted from the #seq_name column in the same row. The rest of column values are N N (mult value) both * (graph_path value).  For processing the input and creating the output, use R tidyverse package to write an R script.
#
# Key Steps:
# Reading the Input: Reads the input file and extracts the necessary columns.
# Pattern Matching: Extracts the LN:i:number and dp:i:number values for the length and coverage columns.
# Computing mult: Calculates the mult column as the ratio of coverage to the median coverage across all rows.
# Filling Columns: Fills the other columns as required, with fixed values or extracted numbers from seq_name.
# Output: Writes the output to a new file following the structure of contigs_stats.txt.

# contigs_stats.txt
# -----------------
# #seq_name       length  coverage        circular        repeat  mult    telomere        alt_group       graph_path
# contig_2        19223   3       N       N       1       both    *       2
# contig_3        70078   6       N       N       1       both    *       1,3,-1
# contig_1        39736   6       N       Y       1       left    *       1

# input: gfa with only sequence part
# The input file, 3-gfa.seq.part.tsv, has columns separated by tabs.
#
# S       edge_1  *       LN:i:39736      dp:i:6
# S       edge_2  *       LN:i:19223      dp:i:3
# S       edge_3  *       LN:i:21885      dp:i:6
#
input_file <- args1$gfa

# output: edges_stats.txt or the edge version of contigs_stats.txt
# #seq_name       length  coverage        circular        repeat  mult    telomere        alt_group       graph_path
# edge_1  39736   6       N       N       1       both    *       1
# edge_2  19223   3       N       N       0       both    *       2
# edge_3  21885   6       N       N       1       both    *       3
output_file <- args1$out

# Read the gfa file of only sequence part.
input_data <- read_tsv(args1$gfa, col_names = FALSE, show_col_types = FALSE)

# Description of what this script does:
# -------------------------------------
# I want to create an output file like contigs_stats.txt.
# The header of the output file is the same.
# The second column of the input file goes to the first column called #seq_name.
# The number at the 4th column in the input file goes to the 2nd column length
# of the output file.
# The number at the 5th column goes to the 3rd column, coverage, in the output
# file.
# These two numbers are extracted by the pattern
# LN:i:number for the 4th column in the input and
# dp:i:number for the 5th column in the input file.
# Most of the rest columns in the output file are just filled nominally or
# just the same values except for mult and graph_path columns among those of
#  circular        repeat  mult    telomere        alt_group       graph_path.
# mult column values are computed by dividing the coverage value
# at the same row by the median of the coverage over all the rows.
# graph_path is the number extracted from the #seq_name column in the same row.
# The rest of column values are N N (mult value) both * (graph_path value).
# For processing the input and creating the output,
# use R tidyverse package to write an R script.

# Extract three columns to save as output_data.
# 1. The X2, the 2nd column, becomes #seq_name or the first column.
# 2. We have length column with the edge length in each row.
# 3. We have coverage column with the depth in each row.
# Then, we select these 3 columns to save as output_data.
#
# References:
# https://stackoverflow.com/a/35804399
# https://stringr.tidyverse.org/articles/regular-expressions.html#look-arounds
# https://stackoverflow.com/a/35804434
# https://stackoverflow.com/a/57438518
# https://stringr.tidyverse.org/reference/str_extract.html
# Tips:
# str_extract("LN:i:39736", "(?<=LN:i:)\\d+") -> "39736"
# str_extract("dp:i:6", "(?<=dp:i:)\\d+") -> "6"
# library(stringr)
# myStrings <- c("MFG: acme", "something else", "MFG: initech")
# str_extract(myStrings, "(?<=MFG:\\s)\\w+")
#
# NOTE: we could have used these simpler one without look-around feature of
# the stringr's regular expression because we have only numbers.
# However, the pattern matching method may be more specific than these simpler
# ones.
# > str_extract("LN:i:39736", "\\d+")
# [1] "39736"
# > str_extract("dp:i:6", "\\d+")
# [1] "6"
output_data <- input_data %>%
  # Have the #seq_name column (from 2nd column of input)
  # Use the backsticks because of the pound sign being a part of the column.
  mutate(
    `#seq_name` = X2,
    # Get the number that follows LN:i
    # like 'LN:i:number' for length (4th column)
    length = as.numeric(str_extract(X4, "(?<=LN:i:)\\d+")),
    # Get the number that follows dp:i
    # like 'dp:i:number' for coverage (5th column)
    coverage = as.numeric(str_extract(X5, "(?<=dp:i:)\\d+"))
  ) %>%
  select(`#seq_name`, length, coverage)

# Calculate the median of coverage
median_coverage <- median(output_data$coverage, na.rm = TRUE)

# Add more columns as follows.
# circular        repeat  mult    telomere        alt_group       graph_path
# 4. circular: N
# 5. repeat: N
# 6. mult: coverage divided by the median coverage
# 7. telomere: both
# 8. alt_group: *
# 9. graph_path: the edge number
#
# Tips:
# str_extract("edge_1", "\\d+") -> "1"
output_data <- output_data %>%
  mutate(
    circular = "N",
    "repeat" = "N",
    # Compute 'mult' as coverage/median_coverage
    mult = round(coverage / median_coverage),
    telomere = "both",
    alt_group = "*",
    # Extract the number from the #seq_name
    graph_path = str_extract(`#seq_name`, "\\d+")
  )

# We need the column names to match it with contigs_stats.txt.
write_tsv(output_data, output_file, col_names = TRUE)
