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
  option1 <- "locus_tag"
  input1 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/30-contigger/graph_final.gfa"
  input2 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/assembly_info_organelle_annotation_count-all.txt"
  output1 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/mt.contig.name-1"
  output2 <- "/media/h2/goshng/figshare/Brassica_rapa/o2/0/mt.contig.name-1.stats"

  s="Vigna_radiata"
  s="Brassica_rapa"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/0/1/mtcontigs")
  input1 <- paste0(input_dir0, "/2-gfa.seq.part.tsv")
  input2 <- paste0(input_dir0, "/1-mtcontig.annotated.txt")
  output1 <- paste0(input_dir0, "/2-gfa.seq.filtered.txt")
  output2 <- paste0(input_dir0, "/2-gfa.seq.filtered.range.txt")
}

x1 <- read_tsv(input1, col_names = FALSE)
x2 <- read_tsv(input2, col_names = FALSE)
m1 <- mean(x2$X2)
sd1 <- sd(x2$X2)
m2 <- median(x2$X2)

# x2 <- scan(input2)

y <- x1 |> separate(X5, into = c("dp", "i", "depth"), sep = ":") |>
  mutate(depth = as.integer(depth))

# x2[1]: mean
# x2[2]: sd
# x2[3]: median
# m1 <- median - 5 * sd
# m2 <- median + 2 * median or 3 * median
#
# most
# m1 <- m2 - 5 * x2[2]
#
# Brassica_rapa
# m1 <- x2[3] - 2 * x2[2]
#
# Case 1. mean >> sd or mean > 5 x sd
# 
#
# Case 2. mean < sd
# delete the largest-depth edge until mean > sd
#
# Case 3. mean > sd, mean < 3 x sd
#

v1 <- 0
v2 <- 0

##################################################
while (m1 < 1 * sd1 && nrow(x2) > 1) {

x2 <- x2 |> filter(X2 < max(X2))
m1 <- mean(x2$X2)
sd1 <- sd(x2$X2)
m2 <- median(x2$X2)

}

##################################################
if (nrow(x2) == 1) {

m1 <- mean(x2$X2)
sd1 <- m1 / 4
m2 <- median(x2$X2)

}

# if (m2/2 < sd1) {
#   v1 <- m2/2
#   v2 <- m2 * 3
# } else {
#   v1 <- m2 - sd1
# }

##################################################
if (m1 > 5 * sd1) {

v1 <- m2 - 2 * sd1 # median - 2 x sd
v2 <- m2 + 2 * m2  # 3 x median for repeats

##################################################
} else if (m1 > 2 * sd1) {

v1 <- m2 - 1 * sd1 # median - 1 x sd
v2 <- m2 + 2 * m2  # 3 x median for repeats

##################################################
} else if (m1 > 1 * sd1) {

v1 <- m2 - 0.5 * sd1 # median - 0.5 x sd
v2 <- m2 + 2 * m2    # 3 x median for repeats

##################################################
} else if (m1 < 1 * sd1) {

  # delete the largest-depth edge in x2
  stop()

##################################################
} else {
  # what case? m1 == sd1
  stop()
}

y |> filter(v1 < depth, depth < v2) |>
  select(X2, depth) |>
  write_tsv(output1, col_names = FALSE)

tibble(c(v1,v2)) |>
  write_tsv(output2, col_names = FALSE)

# m1 <- x2[3] - 1 * x2[2]
# if (m1 < 0) {
#   m1 <- 10
# }
# m2 <- x2[3] + 2 * x2[3]
# # m2 <- x2[1] + 3 * x2[2]
#
# y |> filter(m1 < depth, depth < m2) |>
#   select(X2, depth) |>
#   write_tsv(output1, col_names = FALSE)


# y <- x |> filter(Copy > 1, PT<MT) |> mutate(R1=as.integer(Length/MT)) |> arrange(MT<=PT)
# cat(paste0("edge_",y$Edge), file = output1, sep="\n")
# cat(mean(y$V3),sd(y$V3),median(y$V3), file=output2, sep="\n")
