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
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  input1 <- args[1]
  input2 <- args[2]
  input3 <- args[3]
  output1 <- args[4]
} else {
  s="Vigna_radiata"
  s="Brassica_rapa"
  s="Spirodela_polyrhiza"
  s="Trifolium_pratense"
  input0 <- paste0("/media/h2/goshng/figshare/", s, "/o/1")
  input1 <- paste0(input0, "/assembly_info_organelle_annotation_count-all.txt")
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/1/mtdna")
  input2 <- paste0(input_dir0, "/1-gfa.links.order.txt")
  input3 <- paste0(input_dir0, "/2-gfa.links.seed.txt")
  output1 <- paste0(input_dir0, "/3-annotated.edge.txt")
}

# Read the file
data <- read_delim(input1, delim=' ')

# Process the Edge column to extract absolute numbers and assign MT values
processed_data <- data |>
  # Separate the Edge column into individual numbers
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  # Summarize by taking the maximum MT value for each edge
  summarise(MT = max(MT, na.rm = TRUE)) |>
  ungroup()

# Find the edge(s) with the largest MT value
max_mt_edges <- processed_data |>
  filter(MT == max(MT)) |>
  sample_n(1)

# Function to select a line containing a specific number from a comma-delimited file
select_line_by_number <- function(number, file1) {
  # Read the file and split each line by comma
  lines <- readLines(file1)
  selected_line <- NA
  
  # Search for the line that contains the specific number
  for (line in lines) {
    nums <- unlist(strsplit(line, ","))
    if (as.character(number) %in% nums) {
      selected_line <- nums
      break
    }
  }
  
  return(selected_line)
}

# Function to convert numbers from selected line to string based on second file mapping
convert_numbers_to_strings <- function(selected_line, file2) {
  # Read the second file (tab-delimited)
  mapping <- read.delim(file2, header = FALSE, sep = "\t", col.names = c("number", "string"))
  
  # Convert the numbers in the selected line to corresponding strings
  result <- c()
  for (num in selected_line) {
    str <- mapping[mapping$number == as.numeric(num), "string"]
    if (length(str) > 0) {
      result <- c(result, str)
    }
  }
  
  # Collapse result into end-line delimited string
  return(paste(result, collapse = "\n"))
}

find_first_column <- function(tsv_file, edge_number) {
  # Read the TSV file
  df <- read_tsv(tsv_file, col_names = c("number", "edge"))

  # Construct the edge string based on the input number
  edge_string <- paste0("edge_", edge_number)

  # Filter the dataframe to find the matching row
  result <- df %>%
    filter(edge == edge_string) %>%
    select(number)

  # Return the first column number or a message if not found
  if(nrow(result) > 0) {
    return(result$number)
  } else {
    return(NA)  # Return NA if not found
  }
}

# Example usage
# file1 <- "path/to/your/comma_delimited_file.txt"
# file2 <- "path/to/your/tab_delimited_file.txt"
# selected_line <- select_line_by_number(123, file1)  # Replace 123 with the number you're looking for
# result_string <- convert_numbers_to_strings(selected_line, file2)
# cat(result_string)

# Edge -> number using order
#
max_mt_edge_number <- find_first_column(input2, max_mt_edges$Edge)
selected_line <- select_line_by_number(max_mt_edge_number, input3)  # Replace 123 with the number you're looking for
result_string <- convert_numbers_to_strings(selected_line, input2)
cat(result_string, "\n", file = output1, sep = "")

# |>
#   write_tsv(output1, col_names = FALSE)

# Print the result
# print(max_mt_edges)

