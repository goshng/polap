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
  input3 <- args[3]
  output1 <- args[4]
} else {
  option1 <- "locus_tag"
  input1 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/30-contigger/graph_final.gfa"
  input2 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/assembly_info_organelle_annotation_count-all.txt"
  output1 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/mt.contig.name-1"
  output2 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/mt.contig.name-1.stats"

  s="Brassica_rapa"
  s="Vigna_radiata"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/0/1/mtcontigs")
  input1 <- paste0(input_dir0, "/gfa.seq.part.tsv")
  input2 <- paste0(input_dir0, "/mtcontig.stats.txt")
  input3 <- paste0(input_dir0, "/gfa.links.mtcontig.txt")
  output1 <- paste0(input_dir0, "/gfa.links.mtcontig.depth.txt")
}

x1 <- read_tsv(input1, col_names = FALSE)
x2 <- scan(input2)
m1 <- x2[3] - 5 * x2[2]
if (m1 < 0) {
  m1 <- 10
}
m2 <- x2[3] + 2 * x2[3]

x3 <- read_tsv(input3, col_names = FALSE)

x3 |> left_join(x1, by = c("X1" = "X2")) |>
  separate(X5, into = c("dp", "i", "depth"), sep = ":") |>
  mutate(depth = as.integer(depth)) |>
  filter(m1 < depth, depth < m2) |>
  select(X1) |>
  write_tsv(output1, col_names = FALSE)

