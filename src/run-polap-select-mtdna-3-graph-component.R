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
} else {
  s="Vigna_radiata"
  s="Brassica_rapa"
  s="Spirodela_polyrhiza"
  s="Trifolium_pratense"
  n="2"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/", n, "/mtdna")
  input1 <- paste0(input_dir0, "/3-gfa.component.txt")
  output1 <- paste0(input_dir0, "/3-gfa.component.graph.csv")
}

# no need:
# $ cat o/1/mtdna/3-component.edge.txt
# edge_1
# edge_2
# edge_3
#
# input1:
# $ cat o/1/mtdna/3-gfa.component.txt
# L       edge_1  +       edge_2  -       0M      L1:i:1799       L2:i:48402      RC:i:142
# L       edge_1  +       edge_2  +       0M      L1:i:1799       L2:i:48402      RC:i:135
# L       edge_1  -       edge_3  +       0M      L1:i:1799       L2:i:174640     RC:i:127
# L       edge_1  -       edge_3  -       0M      L1:i:1799       L2:i:174640     RC:i:151
#
# output
# $ cat 2.csv
  # node1,node2,distance
  # 1-,3+,1
  # 1-,3-,1
  # 3+,3-,1
  # 1+,1-,1
  # 1+,2+,1
  # 1+,2-,1
  # 2+,2-,1

# x1 <- read_tsv(input1, col_names = FALSE)

# stop()
#
# y <- x2 |> separate(X5, into = c("dp", "i", "depth"), sep = ":") |>
#   mutate(depth = as.integer(depth))


# Step 1: Read the input TSV
input_data <- read_tsv(input1, col_names = FALSE)

# Step 2: Remove 'edge_' and extract unique numbers from columns X2 and X4
unique_numbers <- input_data |>
  mutate(X2 = gsub("edge_", "", X2),
         X4 = gsub("edge_", "", X4)) |>
  select(X2, X4) |>
  pivot_longer(cols = everything(), values_to = "number") |>
  distinct(number) |>
  pull(number)

# Step 3: Create node1 and node2 based on the input data
output_data <- input_data |>
  mutate(node1 = paste(gsub("edge_", "", X2), X3, sep = ""),
         node2 = paste(gsub("edge_", "", X4), X5, sep = ""),
         distance = 1) |>
  select(node1, node2, distance)

# Step 4: Create additional rows from unique numbers
additional_rows <- tibble(
  node1 = paste0(unique_numbers, "+"),
  node2 = paste0(unique_numbers, "-"),
  distance = 1
)

# Step 5: Combine the original rows and additional rows
final_output <- bind_rows(output_data, additional_rows)

# Step 6: Write the output to a CSV file
write_csv(final_output, output1)


