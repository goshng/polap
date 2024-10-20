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
  input0 <- paste0("/media/h2/goshng/figshare/", s, "/o/0")
  input1 <- paste0(input0, "/30-contigger/graph_final.gfa")
  input2 <- paste0(input0, "/assembly_info_organelle_annotation_count-all.txt")
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/0/1/mtcontigs")
  output1 <- paste0(input_dir0, "/1-mtcontig.annotated.txt")
  output2 <- paste0(input_dir0, "/1-mtcontig.stats.txt")
}

x0 <- read_delim(input2, delim=' ')

# filtering
#
# Case 1. too long without many MT genes
# Case 2. too big or small copy number
  # filter(Length < 1e+6, MT > 1) |>
  # filter(Copy > 1, PT<MT) |>
x1 <- x0 |>
  filter(Copy > 1, PT<MT) |>
  mutate(R1=as.integer(Length/MT)) |>
  filter(R1 < 1e+5) |>
  # Separate the Edge column into individual numbers
  separate_rows(Edge, sep = ",") |>
  distinct() |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge)))

y <- x1 |>
  filter(Copy < 5 * median(x1$Copy), Copy > 0.2 * median(x1$Copy)) |>
  arrange(MT)

y |> mutate(edgename=paste0("edge_",Edge)) |> select(edgename, V3) |>
  write_tsv(output1, col_names = FALSE)
# cat(paste0("edge_",y$Edge), file = output1, sep="\n")
cat(mean(y$V3),sd(y$V3),median(y$V3), file=output2, sep="\n")
