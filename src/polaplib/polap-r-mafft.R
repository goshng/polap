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
# This script reads a MAFFT aligned FASTA file, computes the alignment length,
# number of matches, mismatches, gaps, and percent identity between the first
# two sequences, and writes the results to a specified output text file.
# We could not use the following code because it does not work for long
# alignments:
# alignment <- pairwiseAlignment(seq1, seq2, type="global")
# We had to use a simpler approach to count matches, mismatches, and gaps.
# Although this is not the most efficient way, it works for long alignments.
#
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
  
  input_dir0 <- file.path("output")
  input1 <- file.path(input_dir0, "mafft1/out.mafft")
  output1 <- file.path(input_dir0, "mafft1/out.txt")
  
  input_dir0 <- file.path("output")
  input1 <- file.path(input_dir0, "mafft2/out.mafft")
  output1 <- file.path(input_dir0, "mafft2/out.txt")
  
  args1 <- parse_args(parser, args = c("--input", input1, "-o", output1))
}

# To test this script, you can run the following two lines in R:
# aligned_fasta_file <- args1$input
# output_text_file <- args1$out

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
  
  ## This is the original code that works for long alignments.
  ## But, it is too complicated.
  # matches <- sum(strsplit(seq1, "")[[1]] == strsplit(seq2, "")[[1]] & strsplit(seq1, "")[[1]] != "-")
  # mismatches <- sum(strsplit(seq1, "")[[1]] != strsplit(seq2, "")[[1]] & strsplit(seq1, "")[[1]] != "-" & strsplit(seq2, "")[[1]] != "-")
  # gaps <- sum(strsplit(seq1, "")[[1]] == "-" | strsplit(seq2, "")[[1]] == "-")
  
  ## This does not work because of too long alignment
  # alignment <- pairwiseAlignment(seq1, seq2, type="global")
  # matches <- nmatch(alignment)
  # mismatches <- nmismatch(alignment)
  # gaps <- nindel(alignment)
  
  chars1 <- strsplit(seq1, "")[[1]]
  chars2 <- strsplit(seq2, "")[[1]]
  matches <- sum(chars1 == chars2 & chars1 != "-")
  mismatches <- sum(chars1 != chars2 & chars1 != "-" & chars2 != "-")
  gaps <- sum(chars1 == "-" | chars2 == "-")
  
  percent_identity <- (matches / alignment_length) * 100
  # percent_identity <- round((matches / alignment_length) * 100, 5)

  # Print alignment length, match, mismatch, gap, and pident value.
  cat("Alignment:", alignment_length, "\n", file = output_text_file)
  cat("Match:", matches, "\n", file = output_text_file, append = TRUE)
  cat("Mismatch:", mismatches, "\n", file = output_text_file, append = TRUE)
  cat("Gap:", gaps, "\n", file = output_text_file, append = TRUE)
  cat(sprintf("Percent Identity: %.5f%%\n", percent_identity), file = output_text_file, append = TRUE)
}

# Run the function with the MAFFT output file
compute_alignment_stats(args1$input, args1$out)
