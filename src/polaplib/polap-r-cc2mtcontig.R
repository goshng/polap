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

# polaplib/polap-r-cc2mtcontig.R

# TODO: document it later
#
# Subcommand: seeds
#
# Used by:
# function _run_polap_seeds-graph { # select seed contigs

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--mitochondrial"),
  action = "store_true",
  default = TRUE, help = "Mitochondrial genome assembly"
)
parser <- add_option(parser, c("-p", "--plastid"),
  action = "store_false",
  dest = "mitochondrial", help = "Plastid genome assembly"
)
parser <- add_option(parser, c("-s", "--seed"),
  action = "store",
  help = "GFA sequence part",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("--order"),
  action = "store",
  help = "Depth range",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$seed)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "5-gfa.links.seed.txt")
  input2 <- file.path(input_dir0, "4-gfa.links.order.txt")
  output1 <- file.path(input_dir0, "6-gfa.links.mtcontig.txt")
  args1 <- parse_args(parser, args = c("--seed", input1, "--order", input2, "-o", output1))
}

# Define file paths
numbers_file <- args1$seed
pairs_file <- args1$order

# Read the numbers from the first file
numbers <- readLines(numbers_file)
numbers <- as.character(numbers)

# Read the number-string pairs from the second file
pairs <- read.table(pairs_file, header = FALSE, col.names = c("Number", "String"), sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE)

# Filter pairs to keep only those with numbers in the list
matching_pairs <- pairs |>
  filter(Number %in% numbers) |>
  dplyr::select(String) |>
  write_tsv(args1$out, col_names = FALSE)
