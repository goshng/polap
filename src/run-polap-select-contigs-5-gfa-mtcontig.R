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
  input2 <- args[2]
  output1 <- args[3]
} else {
  option1 <- "locus_tag"
  input1 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/30-contigger/graph_final.gfa"
  input2 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/assembly_info_organelle_annotation_count-all.txt"
  output1 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/mt.contig.name-1"
  output2 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/mt.contig.name-1.stats"

  s="Brassica_rapa"
  s="Vigna_radiata"
  input0 <- paste0("/media/h2/goshng/figshare/", s, "/o/0")
  input0 <- paste0(input0, "/mt.contig.name-1")
  output2 <- paste0(input0, ".gfa.links.order.txt")

  input1 <- paste0(input0, ".gfa.links.seed.txt")
  input2 <- paste0(input0, ".gfa.links.order.txt")
  output1 <- paste0(input0, ".gfa.links.mtcontig.txt")
}

# Load necessary libraries
library(dplyr)

# Define file paths
numbers_file <- input1
pairs_file <- input2

# Read the numbers from the first file
numbers <- readLines(numbers_file)
numbers <- as.character(numbers)

# Read the number-string pairs from the second file
pairs <- read.table(pairs_file, header = FALSE, col.names = c("Number", "String"), sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE)

# Filter pairs to keep only those with numbers in the list
matching_pairs <- pairs |>
  filter(Number %in% numbers) |>
  select(String) |>
  write_tsv(output1, col_names = FALSE)
