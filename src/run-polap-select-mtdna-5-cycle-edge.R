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
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("stringr"))
args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  input1 <- args[1]
  output1 <- args[2]
} else {
  s="Vigna_radiata"
  s="Brassica_rapa"
  s="Spirodela_polyrhiza"
  s="Macadamia_tetraphyll"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/1/mtdna")
  input1 <- paste0(input_dir0, "/3-gfa.cycle.component.txt")
  output1 <- paste0(input_dir0, "/3-gfa.cycle.edge.txt")
}

# Define the file paths (update with your actual file paths)
input_file <- input1
output_file <- output1

# Read the input file
# Split each line by commas to separate the nodes
# Convert to a tibble for easier manipulation
graph_data <- read_lines(input_file) |>
  map(~ str_split(.x, ",") |> unlist()) |>
  map_dfr(~ tibble(from = .x[1], to = .x[2]))

# Separate the direction indicators from the node identifiers
graph_tidy <- graph_data |>
  mutate(from_node = str_remove(from, "[+-]"),
         from_direction = str_extract(from, "[+-]"),
         to_node = str_remove(to, "[+-]"),
         to_direction = str_extract(to, "[+-]")) |>
  select(from_node, from_direction, to_node, to_direction)

# Filter rows where from_node and to_node are the same
# Create the edge format as per the required output
filtered_edges <- graph_tidy |>
  filter(from_node == to_node) |>
  mutate(edge = paste(paste0("edge_", from_node), from_direction, sep = '\t'))

# Write the result to the output file
write_lines(filtered_edges$edge, output_file)

# Display the resulting edges (optional)
# filtered_edges


