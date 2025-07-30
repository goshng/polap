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

# run-polap-r-assemble-bioproject-3-length-match.R

################################################################################
# Read the BLSATN output to sum the length of subjects longer than 2kb and
# 99% percent identity.
# input columns:
# qseqid qstart sseqid sstart sstrand length pident
#
# TODO: rename: remove run- from the name.
#
# Used By:
# function _run_polap_compare-mtdna {
#
# See Also:
# _run_polap_blast-mtdna : polap-cmd-mtdna.sh
#
# Check: 2025-06-18
################################################################################

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("tidyr"))
args <- commandArgs(trailingOnly = TRUE)

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

if (length(args) > 0) {
  input1 <- args[1]
  output1 <- args[2]
  output2 <- args[3]
} else {
  s <- "Spirodela_polyrhiza"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/1/70-bioproject")
  input1 <- paste0(input_dir0, "/3-blastn3.txt")
  output1 <- paste0(input_dir0, "/3-blastn3.length.txt")
}

x1 <- read_tsv(input1, col_names = FALSE)

x1 |>
  filter(X5 == "plus", X6 > 2000, X7 > 99) |>
  summarize(total = sum(X6)) |>
  write_tsv(output1, col_names = FALSE)
