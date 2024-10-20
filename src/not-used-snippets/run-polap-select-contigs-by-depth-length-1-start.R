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
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  input1 <- args[1]
  output1 <- args[2]
  output2 <- args[3]
} else {
  s <- "Vigna_radiata"
  s <- "Brassica_rapa"
  s <- "Anthoceros_angustus"

  s <- "bioprojects"
  o <- "PRJNA766769"
  jnum <- "1"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input_dir1 <- file.path(input_dir0, jnum, "mtcontigs")
  output1 <- file.path(input_dir1, "1-mtcontig.mt.stats.txt")
  output2 <- file.path(input_dir1, "1-mtcontig.pt.stats.txt")
}

# Load the data from the TSV file
data <- read_delim(input1, delim = " ")

# Sort the data by the 'Copy' column in decreasing order
sorted_data <- data %>%
  select(Contig, Length, V3, Copy, MT, PT, Edge) %>%
  arrange(desc(Copy))

# Calculate the cumulative sum of the 'Length' column
sorted_data <- sorted_data %>%
  mutate(cumulative_length = cumsum(Length))

# Cut off rows at the given cumulative length of 3,000,000
cutoff_data <- sorted_data %>%
  filter(cumulative_length <= 3000000) %>%
  filter(MT > 0 | PT > 0)

x.mt <- cutoff_data |>
  filter(PT < MT)

x.pt <- cutoff_data |>
  filter(PT > MT)

x.mt |>
  summarise(
    depth_lower_bound = mean(V3) - sd(V3) * 5,
    depth_upper_bound = mean(V3) + sd(V3) * 5,
    depth_min = min(V3),
    depth_max = max(V3),
    depth_median = median(V3),
    depth_mean = mean(V3),
    depth_sd = sd(V3),
  ) |>
  write_tsv(output1)

x.pt |>
  summarise(
    depth_lower_bound = mean(V3) - sd(V3) * 5,
    depth_upper_bound = mean(V3) + sd(V3) * 5,
    depth_min = min(V3),
    depth_max = max(V3),
    depth_median = median(V3),
    depth_mean = mean(V3),
    depth_sd = sd(V3),
  ) |>
  write_tsv(output2)

# Find the 'Copy' value at the cutoff
# copy_at_cutoff <- cutoff_data %>%
#   summarise(copy_cutoff = min(Copy))

# Print the result
# print(copy_at_cutoff)
