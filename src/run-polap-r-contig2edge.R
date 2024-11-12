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
suppressPackageStartupMessages(library("ggplot2"))

parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle annotation table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "contig-annotation-table.txt")
  output1 <- file.path(input_dir0, "contig-annotation-table-edges.txt")

  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
}

data <- read_delim(args1$table, delim = " ", show_col_types = FALSE)


# Extract all unique numbers from the Edge column
unique_edges <- data %>%
  pull(Edge) %>%
  str_split(",") %>%
  unlist() %>%
  as.numeric() %>%
  abs() %>%
  unique() %>%
  sort(na.last = NA) # Remove NA values, if any

# Create labels in the format "edge_<number>"
edge_labels <- paste0("edge_", unique_edges)

# Write to file, one label per line
writeLines(edge_labels, args1$out)
