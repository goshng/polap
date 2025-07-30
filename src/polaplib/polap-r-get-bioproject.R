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
# polaplib/polap-r-get-bioproject.R
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

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# args = commandArgs(trailingOnly=TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"),
  action = "store",
  help = "BioProject runinfo CSV",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output folder"
)
parser <- add_option(parser, c("--pacbio"),
  action = "store_true",
  default = FALSE, help = "Use PacBio"
)

args1 <- parse_args(parser)

if (is_null(args1$input)) {
  input_dir0 <- file.path("input")
  input_dir1 <- file.path("output")
  input1 <- file.path(input_dir0, "PRJDB3943.csv")
  input1 <- file.path(input_dir0, "PRJDB5905.csv")
  input1 <- file.path(input_dir0, "PRJDB4656.csv")
  output_dir1 <- file.path(input_dir1, "PRJDB3943")
  output_dir1 <- file.path(input_dir1, "PRJDB5905")
  output_dir1 <- file.path(input_dir1, "PRJDB4656")
  args1 <- parse_args(parser, args = c(
    "--input", input1,
    "--pacbio",
    "-o", output_dir1
  ))
}

dir.create(args1$out, recursive = TRUE)
output1 <- file.path(args1$out, "1-sra-long-read.txt")
output2 <- file.path(args1$out, "1-sra-short-read.txt")
output3 <- file.path(args1$out, "1-species.txt")
output4 <- file.path(args1$out, "1-runinfo.per.species.tsv")

x1 <- read_csv(args1$input, col_names = TRUE, show_col_types = FALSE)
# x1 <- read_csv("1.csv", col_names = TRUE, show_col_types = FALSE)

if (ncol(x1) == 0) {
  quit(save = "no", status = 1)
}

# x1 <- read_tsv(input1, col_names = TRUE)

# Load the data
data <- x1

if (args1$pacbio) {
  selected_runs_nano <- data |>
    filter(
      LibraryStrategy == "WGS" | LibraryStrategy == "WGA", # we need this.
      LibrarySource == "GENOMIC",
      Platform == "PACBIO_SMRT"
    ) |>
    group_by(ScientificName) |>
    filter(bases == max(bases)) |>
    ungroup()

  # Create a table with three columns.
  # ScientificName, Run from OXFORD_NANOPORE, and Run from ILLUMINA
  result_table <- selected_runs_nano |>
    select(ScientificName, Run_Nano = Run) |>
    write_tsv(output4)

  l1 <- x1 |>
    filter(
      LibrarySource == "GENOMIC",
      Platform == "PACBIO_SMRT"
      # LibraryLayout == "SINGLE"
    ) |>
    select(Run, bases, avgLength, size_MB, LibraryStrategy, LibrarySource, Platform, Model, ScientificName)

  l1 |>
    filter(bases == max(bases)) |>
    write_tsv(output1, col_names = FALSE)
} else {
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
    select(Run, bases, avgLength, size_MB, LibraryStrategy, LibrarySource, Platform, Model, ScientificName)
  # select(Run, bases, LibraryStrategy, LibrarySource, Platform, ScientificName)

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
}


l1 |>
  select(ScientificName) |>
  distinct() |>
  write_tsv(output3, col_names = FALSE)
