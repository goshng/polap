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
  output2 <- args[4]
} else {
  s="Vigna_radiata"
  s="Brassica_rapa"

  s="bioprojects"
  o="PRJNA766769"
  jnum="1"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir1 <- file.path(input_dir0, jnum, "mtcontigs")
  input1 <- file.path(input_dir1, "2-gfa.seq.part.tsv")
  input2 <- file.path(input_dir1, "1-mtcontig.mt.stats.txt")
  output1 <- file.path(input_dir1, "2-gfa.seq.filtered.txt")
  output2 <- file.path(input_dir1, "2-gfa.seq.filtered.range.txt")
}

x1 <- read_tsv(input1, col_names = FALSE)
x2 <- read_tsv(input2)

y <- x1 |> separate(X5, into = c("dp", "i", "depth"), sep = ":") |>
  mutate(depth = as.integer(depth))

v1 <- x2$depth_lower_bound
v2 <- x2$depth_upper_bound

y |> filter(v1 <= depth, depth <= v2) |>
  select(X2, depth) |>
  write_tsv(output1, col_names = FALSE)

tibble(c(v1,v2)) |>
  write_tsv(output2, col_names = FALSE)

