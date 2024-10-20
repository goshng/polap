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
  output1 <- args[3]
  output2 <- args[4]
} else {
  s="Vigna_radiata"
  s="Brassica_rapa"
  s="Anthoceros_angustus"

  s="bioprojects"
  o="PRJNA766769"
  jnum="1"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input1 <- file.path(input_dir0, "30-contigger/graph_final.gfa")
  input2 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input_dir1 <- file.path(input_dir0, jnum, "mtcontigs")
  output1 <- file.path(input_dir1, "1-mtcontig.annotated.txt")
  output2 <- file.path(input_dir1, "1-mtcontig.stats.txt")
}

x0 <- read_delim(input2, delim=' ')

# PT selection
# gene density cutoff: 1 in 10 kb
x1 <- x0 |>
  filter(Copy > 1) |>
  filter(PT > 1) |>
  filter(PT > MT) |>
  mutate(RT=as.integer(Length/PT)) |>
  filter(RT < 1e+4) |>
  # Separate the Edge column into individual numbers
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  # Summarize by taking the maximum MT value for each edge
  # summarise(MT = max(Length, na.rm = TRUE)) |>
  ungroup() |>
  distinct()

# MT selection
# gene density cutoff: 1 in 100 kb
x2 <- x0 |>
  filter(Copy > 1) |>
  filter(MT > 1) |>
  filter(PT <= MT) |>
  mutate(RT=as.integer(Length/MT)) |>
  # Separate the Edge column into individual numbers
  filter(RT < 1e+5) |>
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  # Summarize by taking the maximum MT value for each edge
  # summarise(MT = max(Length, na.rm = TRUE)) |>
  ungroup() |>
  distinct()

bind_rows(x1, x2) |> mutate(edgename=paste0("edge_",Edge)) |>
  select(edgename, Length, V3, Copy, MT, PT, RT) |>
  arrange(Copy) |>
  write_tsv(output1, col_names = FALSE)
