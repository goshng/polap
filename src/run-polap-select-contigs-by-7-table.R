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

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mixR"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
# args = commandArgs(trailingOnly=TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle annotation table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-m", "--mtcontigseed"),
  action = "store",
  help = "Organelle contig seeds",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  s <- "Vigna_radiata"
  s <- "Brassica_rapa"
  s <- "Anthoceros_angustus"

  s <- "bioprojects"

  # o <- "PRJNA914763" # type 3
  o <- "PRJNA597121" # type 5 + 3
  # o <- "PRJNA766769"    # type 5 + 3
  # o <- "PRJNA635654"    # type 5 + 3
  # o <- "PRJNA646849" # type 5 + 3, type 3
  # o <- "PRJEB26621"     # type 5 + 3, type 3
  # o <- "PRJNA1074514"   # type 5 + 3
  # o <- "PRJEB19787"     # no mtDNA seeds
  # o <- "PRJEB21270" # type 3; type 5 Error: insufficient data for model estimation.
  # o <- "PRJEB79308"

  jnum <- "1"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input2 <- file.path(input_dir0, "mt.contig.name-")
  input2 <- paste0(input2, jnum)
  input_dir1 <- file.path(input_dir0, jnum, "mtcontigs")
  output1 <- file.path(input_dir1, "1-mtcontig")

  # for Plastid
  # args1 <- parse_args(parser, args = c("--plastid", "--table", input1, "--gene-count", 1, "--gene-compare", "--gene-density", 100))
  #
  # for mitochondrial
  #
  # Try
  # type 1
  # type 2
  # type 3
  # type 4 + {1 or 2 or 3}
  # type 5 + {1 or 2 or 3}
  #
  # Recommended
  # type 3
  args1 <- parse_args(parser, args = c("--table", input1, "-m", input2, "-o", output1))
  # type 5 + 3
  # args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare", "--gene-density", 10, "--save-range-copy", "-o", output1))
  # type 3 + depth graph
  # type 5 + 3 + depth graph
  # type 5 + 3 + graph (using depth range from the distribution)
  #
}

output1 <- paste0(args1$out, ".table.tsv")

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)
x1 <- read_delim(args1$mtcontigseed, delim = " ", col_names = c("edgename"), show_col_types = FALSE)

xall <- x0 |>
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  ungroup() |>
  distinct() |>
  dplyr::select(-V4, -V5, -V7, -V8) |>
  mutate(edgename = paste0("edge_", Edge)) |>
  relocate(edgename) |>
  filter(MT > PT, MT > 0) |>
  arrange(Copy)

left_join(x1, xall) |>
  filter(!is.na(MT)) |>
  arrange(desc(MT > PT), desc(MT)) |>
  write_tsv(output1, col_names = FALSE)
