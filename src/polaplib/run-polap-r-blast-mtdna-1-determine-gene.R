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

# Load necessary library
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
args = commandArgs(trailingOnly=TRUE)

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

if (length(args) > 0) {
  input1 <- args[1]
  input2 <- args[2]
  input3 <- args[3]
  output1 <- args[4]
  output2 <- args[5]
  output3 <- args[6]
} else {
  s="bioprojects"
  base_dir0 <- "/media/h2/goshng/figshare/"

  input_dir0 <- paste0(base_dir0, s, "/o/1/60-chloroplot/")
  input_dir1 <- paste0(base_dir0, s, "/")
  input1 <- paste0(input_dir0, "2.bed4")
  input2 <- paste0(input_dir1, "src/polap-mt.1.c70.3.faa.name")
  input3 <- paste0(input_dir1, "src/polap-mtgenes.txt")
  output1 <- paste0(input_dir0, "2.gene")
  output2 <- paste0(input_dir0, "2.bed4.description")
  output3 <- paste0(input_dir0, "2.bed4.count")
  output1.check <- paste0(output1, ".check")
}


# Load the first file (1.bed4) which contains 4 columns with the 4th being the ID
bed4 <- read_tsv(input1, col_names = c("chr", "start", "end", "ID"))

# Load the second file (polap-mt.1.c70.3.faa.name) which contains 2 columns: ID and description
polap <- read_tsv(input2, col_names = c("ID", "description"))

# Merge the two datasets by ID
merged <- left_join(bed4, polap, by = "ID")

# Load the gene names and their alternatives from the mtgenes.txt file
genes <- read_tsv(input3, col_names = c("gene", "alternatives"))

# Split the alternatives column into individual gene names, including the main gene name
genes <- genes %>%
  mutate(alternatives = str_split(alternatives, ",")) %>%
  unnest(alternatives) %>%
  mutate(alternatives = str_trim(alternatives))

# Function to count occurrences of each alternative in the description column
count_gene_occurrences <- function(alt_name) {
  sum(str_detect(merged$description, regex(alt_name, ignore_case = TRUE)))
}

# Count occurrences of each gene name and its alternatives in the description field
gene_counts <- genes %>%
  group_by(gene) %>%
  summarize(count = sum(map_int(alternatives, count_gene_occurrences)))

# Find the gene with the highest count
most_frequent_gene <- gene_counts %>% filter(count == max(count))

# Display results
merged |>
  write_tsv(output2)

gene_counts |>
  write_tsv(output3, col_names = FALSE)

most_frequent_gene |>
  select(gene) |>
  filter(nchar(gene) == max(nchar(gene))) |>
  write_tsv(output1, col_names = FALSE)

most_frequent_gene |>
  select(gene) |>
  arrange(desc(nchar(gene)) |>
  write_tsv(output4, col_names = FALSE)


