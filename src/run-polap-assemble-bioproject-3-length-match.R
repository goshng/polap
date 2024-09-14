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

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("tidyr"))
args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  input1 <- args[1]
  output1 <- args[2]
  output2 <- args[3]
} else {
  s="bioprojects"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/bioproject")
  input1 <- paste0(input_dir0, "/3-blastn3.txt")
  output1 <- paste0(input_dir0, "/3-blastn3.length.txt")
}

x1 <- read_tsv(input1, col_names = FALSE)

x1 |>
  filter(X5 == "plus", X6 > 2000, X7 > 99) |>
  summarize(total = sum(X6)) |>
  write_tsv(output1, col_names = FALSE)
