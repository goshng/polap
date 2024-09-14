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
  output3 <- args[4]
} else {
  s="bioprojects"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, "o/bioproject")
  input1 <- file.path(input_dir0, "1-runinfo.tsv")
  output1 <- file.path(input_dir0, "1-sra-long-read.txt")
  output2 <- file.path(input_dir0, "1-sra-short-read.txt")
  output3 <- file.path(input_dir0, "1-species.txt")
}

x1 <- read_tsv(input1, col_names = TRUE)

l1 <- x1 |>
  filter(LibrarySource == "GENOMIC", 
         Platform == "OXFORD_NANOPORE", 
         LibraryLayout == "SINGLE") |>
  select(Run, bases, Platform, ScientificName)

l1 |>
  filter(bases == max(bases)) |>
  write_tsv(output1, col_names = FALSE)

s1 <- x1 |>
  filter(LibrarySource == "GENOMIC", 
         Platform == "ILLUMINA", 
         LibraryLayout == "PAIRED") |>
  select(Run, bases, Platform, ScientificName)

s1 |>
  filter(bases == max(bases)) |>
  write_tsv(output2, col_names = FALSE)

l1 |>
  select(ScientificName) |>
  distinct() |>
  write_tsv(output3, col_names = FALSE)
