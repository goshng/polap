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

# polaplib/run-polap-r-depthfilter-gfa.R
# Check: 2025-06-17

################################################################################
# Filter the assembly graphs in a gfa file by depths.
#
# This is used by seeds-graph and archive.
# _polap_seeds_depthfilter-gfa
# step 3-1: creating GFA without sequence data
#   convert gfa to another gfa2 using gfatools view command.
# step 3-2: extracting sequence part of GFA
#   extract the sequence part of the gfa2 using grep ^S command.
# step 3-3: filtering GFA sequence part using depth range
#   use this script to extract sequences with the depth range
#
# Example:
# Rscript ./src/polaplib/run-polap-r-depthfilter-gfa.R \
#   --gfa output/3-gfa.seq.all.tsv \
#   --lower-bound-depth 40 --upper-bound-depth 90 \
#   --out output/3-gfa.seq.depthfiltered.txt.copy
#
# TEST: done in test/test.sh
# input1: a gfa file
# input2: a depth range TSV
# depth_lower_bound	depth_upper_bound
# 40 90
#
# Step 3-1
# output/3-gfa.all.gfa
# S       edge_1  *       LN:i:18273      dp:i:41
# S       edge_2  *       LN:i:88015      dp:i:38
# S       edge_3  *       LN:i:26202      dp:i:89
# L       edge_1  +       edge_3  -       0M      L1:i:18273      L2:i:26202      RC:i:35
# L       edge_1  -       edge_3  -       0M      L1:i:18273      L2:i:26202      RC:i:42
# L       edge_2  +       edge_3  +       0M      L1:i:88015      L2:i:26202      RC:i:36
# L       edge_2  -       edge_3  +       0M      L1:i:88015      L2:i:26202      RC:i:42
# Step 3-2
# output/3-gfa.seq.all.tsv
# S       edge_1  *       LN:i:18273      dp:i:41
# S       edge_2  *       LN:i:88015      dp:i:38
# S       edge_3  *       LN:i:26202      dp:i:89
# Step 3-3
# output/3-gfa.seq.depthfiltered.txt
# edge_1  41
# edge_3  89
#
# See Also:
# 1.
# polap-cmd-seeds.sh
# function _polap_seeds_depthfilter-gfa
# function _run_polap_seeds-graph { # select seed contigs
# 2.
# polap-cmd-archive.sh
# function _polap_archive_gfa-depth-filtered
#
# TODO: rename: polap-r-depthfilter-gfa.R
#
# Check: 2025-06-17
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--mitochondrial"),
  action = "store_true",
  default = TRUE, help = "Mitochondrial genome assembly"
)
parser <- add_option(parser, c("-p", "--plastid"),
  action = "store_false",
  dest = "mitochondrial", help = "Plastid genome assembly"
)
parser <- add_option(parser, c("-g", "--gfa"),
  action = "store",
  help = "GFA sequence part",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("--depth"),
  action = "store",
  help = "Depth range",
  default = NULL,
  metavar = "<FILE>"
)
parser <- add_option(parser, c("--lower-bound-depth"),
  action = "store",
  type = "integer",
  default = NULL,
  help = "A lower-bound depth"
)
parser <- add_option(parser, c("--upper-bound-depth"),
  action = "store",
  type = "integer",
  default = NULL,
  help = "A upper-bound depth"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$gfa)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "3-gfa.seq.part.tsv")
  input2 <- file.path(input_dir0, "1-depth.range.txt")
  output1 <- file.path(input_dir0, "3-gfa.seq.filtered.txt")
  args1 <- parse_args(parser, args = c("--gfa", input1, "--depth", input2, "-o", output1))
}

x1 <- read_tsv(args1$gfa, col_names = FALSE, show_col_types = FALSE)

if (is_null(args1$depth)) {
  stopifnot(!is_null(args1$`lower-bound-depth`))
  stopifnot(!is_null(args1$`upper-bound-depth`))
  v1 <- args1$`lower-bound-depth`
  v2 <- args1$`upper-bound-depth`
} else {
  x2 <- read_tsv(args1$depth, show_col_types = FALSE)
  v1 <- x2$depth_lower_bound
  v2 <- x2$depth_upper_bound
}

y <- x1 |>
  separate(X5, into = c("dp", "i", "depth"), sep = ":") |>
  mutate(depth = as.integer(depth))


y |>
  filter(v1 <= depth, depth <= v2) |>
  dplyr::select(X2, depth) |>
  write_tsv(args1$out, col_names = FALSE)

# tibble(c(v1,v2)) |>
#   write_tsv(output2, col_names = FALSE)
