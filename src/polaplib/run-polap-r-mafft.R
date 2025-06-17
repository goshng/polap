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
# Parse MAFFT alignment output file to extract alignment statistics.
# Use: Biostrings
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
  help = "MAFFT aligned FASTA",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$input)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "Anthoceros_agrestis/0/mafft/1/out.mafft")
  output1 <- file.path(input_dir0, "Anthoceros_agrestis/0/mafft/1/out.txt")

  args1 <- parse_args(parser, args = c("--input", input1, "-o", output1))
}

# Function to compute alignment statistics
compute_alignment_stats <- function(aligned_fasta_file, output_text_file) {
  # Read the aligned sequences
  aln <- readDNAStringSet(aligned_fasta_file)

  # Extract sequences
  seq1 <- as.character(aln[[1]])
  seq2 <- as.character(aln[[2]])

  # Ensure both sequences have the same length
  if (nchar(seq1) != nchar(seq2)) {
    stop("Error: Sequences have different lengths in alignment!")
  }

  # Compute statistics
  alignment_length <- nchar(seq1)
  matches <- sum(strsplit(seq1, "")[[1]] == strsplit(seq2, "")[[1]] & strsplit(seq1, "")[[1]] != "-")
  mismatches <- sum(strsplit(seq1, "")[[1]] != strsplit(seq2, "")[[1]] & strsplit(seq1, "")[[1]] != "-" & strsplit(seq2, "")[[1]] != "-")
  gaps <- sum(strsplit(seq1, "")[[1]] == "-" | strsplit(seq2, "")[[1]] == "-")
  percent_identity <- (matches / alignment_length) * 100
  # percent_identity <- round((matches / alignment_length) * 100, 5)

  # Print results
  cat("Alignment:", alignment_length, "\n", file = output_text_file)
  cat("Match:", matches, "\n", file = output_text_file, append = TRUE)
  cat("Mismatch:", mismatches, "\n", file = output_text_file, append = TRUE)
  cat("Gap:", gaps, "\n", file = output_text_file, append = TRUE)
  cat(sprintf("Percent Identity: %.5f%%\n", percent_identity), file = output_text_file, append = TRUE)
}

# Run the function with the MAFFT output file
compute_alignment_stats(args1$input, args1$out)
