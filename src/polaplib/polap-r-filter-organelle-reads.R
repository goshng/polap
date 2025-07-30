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
# Copy this R template file to create a new R script that is used in polap.
# The following template file provides argument processing.
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# args = commandArgs(trailingOnly=TRUE)
# Define command-line options
option_list <- list(
  make_option(c("-t", "--table"), type = "character", help = "Path to annotation table (TSV)"),
  make_option(c("-l", "--length"), type = "numeric", help = "Maximum cumulative length (L)"),
  make_option(c("-o", "--output"), type = "character", help = "Output file to write contig IDs")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
args1 <- parse_args(opt_parser)

# Check required args1ions
# if (is.null(args1$table) || is.null(args1$length) || is.null(args1$output)) {
#   print_help(args1_parser)
#   stop("All arguments -t, -l, and -o are required.", call. = FALSE)
# }

if (is_null(args1$table)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"
  print(s)

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.pdf")

  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-l", 1e+5,
    "-o", output1
  ))
}

x0 <- read_tsv(args1$table, show_col_types = FALSE)
# print(args1$table)

# Filter contigs by cumulative length
selected <- x0 %>%
  mutate(cum_length = cumsum(Length)) %>%
  filter(cum_length <= args1$length)

# Write contig IDs to the output file
write_lines(selected$Contig, args1$output)
