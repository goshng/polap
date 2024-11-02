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
suppressPackageStartupMessages(library("mixR"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
# args = commandArgs(trailingOnly=TRUE)
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
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.pdf")

  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
}

# Load necessary library
library(ShortRead)

# Define the path to your FASTQ file
fastq_file <- "path/to/your/file.fastq"

# Read the FASTQ file
fastq_data <- readFastq(fastq_file)

# Extract the sequences and their lengths
seq_lengths <- width(sread(fastq_data))

# Define the weights as proportional to the sequence lengths
weights <- seq_lengths / sum(seq_lengths)

# Define the number of sequences you want to sample
num_samples <- 100 # adjust this to your preference

# Sample indices with probabilities proportional to the weights
sampled_indices <- sample(seq_along(seq_lengths), size = num_samples, prob = weights, replace = TRUE)

# Extract the sampled sequences
sampled_fastq <- fastq_data[sampled_indices]

# Write the sampled sequences to a new FASTQ file
writeFastq(sampled_fastq, "path/to/your/sampled_output.fastq")
