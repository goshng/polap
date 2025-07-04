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
# polaplib/run-polap-r-get-bioproject-1.R
#
# This script process BioProject information to find long-read data of different
# platforms. We started with collecting long-read sequencing data.
# Use this to extract BioProject information.
#
# Used by:
# function _run_polap_get-bioproject { # get BioProject info from NCBI
#
# Check: 2025-06-18
################################################################################


suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("tidyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  input1 <- args[1]
  output1 <- args[2]
  output2 <- args[3]
  output3 <- args[4]
  output4 <- args[5]
} else {
  s <- "bioprojects"
  input_dir0 <- file.path("/media/h2/goshng/figshare/bioprojects/test/o/0-bioproject")
  input1 <- file.path(input_dir0, "1-runinfo.tsv")
  output1 <- file.path(input_dir0, "1-sra-long-read.txt")
  output2 <- file.path(input_dir0, "1-sra-short-read.txt")
  output3 <- file.path(input_dir0, "1-species.txt")
  output4 <- file.path(input_dir0, "1-runinfo.per.species.tsv")
}

x1 <- read_tsv(input1, col_names = TRUE)

# Load the data
data <- x1

# Filter and select runs for OXFORD_NANOPORE platform
selected_runs_nano <- data |>
  filter(
    LibraryStrategy == "WGS" | LibraryStrategy == "WGA", # we need this.
    LibrarySource == "GENOMIC",
    Platform == "OXFORD_NANOPORE"
  ) |>
  group_by(ScientificName) |>
  filter(bases == max(bases)) |>
  ungroup()

# Filter and select runs for ILLUMINA platform
selected_runs_illumina <- data |>
  filter(
    LibraryStrategy == "WGS" | LibraryStrategy == "WGA", # we need this.
    LibrarySource == "GENOMIC",
    Platform == "ILLUMINA" | Platform == "DNBSEQ"
  ) |>
  group_by(ScientificName) |>
  filter(bases == max(bases)) |>
  ungroup()

# Find ScientificNames present in both selected runs
common_scientific_names <- intersect(
  selected_runs_nano$ScientificName,
  selected_runs_illumina$ScientificName
)

# Filter rows with ScientificName present in both selected runs
selected_runs_common_nano <- selected_runs_nano |>
  filter(ScientificName %in% common_scientific_names)

selected_runs_common_illumina <- selected_runs_illumina |>
  filter(ScientificName %in% common_scientific_names)

# Create a table with three columns.
# ScientificName, Run from OXFORD_NANOPORE, and Run from ILLUMINA
result_table <- selected_runs_common_nano |>
  select(ScientificName, Run_Nano = Run) |>
  left_join(
    selected_runs_common_illumina |>
      select(ScientificName, Run_Illumina = Run),
    by = "ScientificName"
  ) |>
  write_tsv(output4)

l1 <- x1 |>
  filter(
    LibrarySource == "GENOMIC",
    Platform == "OXFORD_NANOPORE",
    LibraryLayout == "SINGLE"
  ) |>
  select(Run, bases, LibraryStrategy, LibrarySource, Platform, ScientificName)

l1 |>
  filter(bases == max(bases)) |>
  write_tsv(output1, col_names = FALSE)

s1 <- x1 |>
  filter(
    LibrarySource == "GENOMIC",
    Platform == "ILLUMINA",
    LibraryLayout == "PAIRED"
  ) |>
  select(Run, bases, LibraryStrategy, LibrarySource, Platform, ScientificName)

s1 |>
  filter(bases == max(bases)) |>
  write_tsv(output2, col_names = FALSE)

l1 |>
  select(ScientificName) |>
  distinct() |>
  write_tsv(output3, col_names = FALSE)
