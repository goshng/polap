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
  s="Vigna_radiata"
  s="Anthoceros_angustus"
  s="Brassica_rapa"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/0/1/mtcontigs")
  input1 <- paste0(input_dir0, "/gfa.seq.filtered.txt")
  input2 <- paste0(input_dir0, "/gfa.link.part.tsv")
  output1 <- paste0(input_dir0, "/gfa.seq.filtered.link.added.txt")
}

# Step 1: Read the two files
# Adjust file paths as needed

# First file with a single column of strings
file1 <- read_tsv(input1, col_names = "string1")

# Second file with two columns of strings
file2 <- read_tsv(input2, col_names = c("string2_col1", "string2_col2"))

# Step 2: Extract the second column values where the first column string exists in the first file
# Use a semi-join to extract matching rows, then pull the second column

extracted_col2 <- file2 |>
  semi_join(file1, by = c("string2_col1" = "string1")) %>%
  pull(string2_col2)

# Step 3: Extract the first column values where the second column string exists in the first file
# Use a semi-join again and pull the first column

extracted_col1 <- file2 |>
  semi_join(file1, by = c("string2_col2" = "string1")) %>%
  pull(string2_col1)

# Step 4: Combine the extracted values with the values from file1, and remove duplicates
combined_strings <- c(file1$string1, extracted_col1, extracted_col2) |>
  unique()

# Step 5: Create a final tibble with the unique values
final_result <- tibble(unique_strings = combined_strings)

# View the result
# print(final_result)

# Show the first few rows of the result
final_result |>
  write_tsv(output1, col_names = FALSE)
