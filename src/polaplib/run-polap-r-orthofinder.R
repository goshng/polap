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
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# args = commandArgs(trailingOnly=TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "OrthoFinder Orthogroups GeneCount tsv",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "Hylodesmum_podocarpum/taxonomy/04-orthofinder/OrthoFinder/Results_1/Orthogroups/Orthogroups.GeneCount.tsv")
  output1 <- file.path(input_dir0, "minimum.og.txt")
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-o", output1
  ))
}

# Read the TSV file
data <- read_tsv(args1$table, show_col_types = FALSE)

# Calculate the number of taxa (number of columns minus 2)
num_taxa <- ncol(data) - 2

# Filter rows where the 'Total' column is at least 0.9 times the number of taxa
filtered_data <- data %>%
  filter(Total >= 0.9 * num_taxa)

# Print the number of OG groups after filtering
num_og_groups <- nrow(filtered_data)
cat(num_og_groups, "\n")
