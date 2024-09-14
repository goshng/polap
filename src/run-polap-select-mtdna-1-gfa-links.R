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
  s="Brassica_rapa"
  s="Vigna_radiata"
  s="Anthoceros_angustus"
  s="Spirodela_polyrhiza"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/1/mtdna")
  input1 <- paste0(input_dir0, "/1-gfa.links.tsv")
  output1 <- paste0(input_dir0, "/1-gfa.links.number.txt")
  output2 <- paste0(input_dir0, "/1-gfa.links.order.txt")
}

# x1 <- read_tsv(input0, col_names = c("string", "depth"))

# Step 1: Read the TSV file (assuming it has no header and two columns)
# Adjust the file path as necessary
df <- read_tsv(input1, 
               col_names = c("type", "string1", "strand1", "string2", "strand2",
                             "match", "copy1", "copy2", "copy3"))

# Step 2: Combine both columns to identify all unique strings
unique_strings <- df |>
  select(string1, string2) |>
  unlist() |>
  unique() |>
  sort()

# Step 3: Create a mapping table: each string corresponds to a unique non-negative integer
mapping_table <- tibble(
  int_value = seq(0, length(unique_strings) - 1),
  string = unique_strings
)

# Step 4: Use left_join to map each string to its corresponding integer
df_mapped <- df |>
  left_join(mapping_table, by = c("string1" = "string")) |>
  rename(int1 = int_value) |>
  left_join(mapping_table, by = c("string2" = "string")) |>
  rename(int2 = int_value) |>
  select(int1, int2)

# Step 5: Now `df_mapped` contains integer values corresponding to each string
# and `mapping_table` shows the string-to-integer mapping

# Show the first few rows of the result
df_mapped |>
  write_tsv(output1, col_names = FALSE)

# Show the mapping table
mapping_table |>
  write_tsv(output2, col_names = FALSE)

