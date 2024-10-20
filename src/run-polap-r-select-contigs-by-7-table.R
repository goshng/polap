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
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))

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
parser <- add_option(parser, c("-l", "--length"),
                     type = "integer",
                     default = -1,
                     help = "Maximum length of seed contigs [default off]",
                     metavar = "number"
)
args1 <- parse_args(parser)



if (is_null(args1$table)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"
  
  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input2 <- file.path(input_dir0, "mt.contig.name-8")
  output1 <- file.path(input_dir0, "1-mtcontig.table.tsv")
  args1 <- parse_args(parser, args = c("--table", input1, "-m", input2, "-o", output1))
  
}

output1 <- args1$out

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
  arrange(Copy)

if (args1$length > 0) {
  xall <- xall |>
    filter(Length < args1$length)
}

# Not sure whether we need this:
# Case: PRJNA838254-Prunus_padus
# The results may need this because of no annotation but linked to 
# those annotated.
# filter(MT > PT, MT > 0) |>
# arrange(Copy)


left_join(x1, xall) |>
  filter(!is.na(MT)) |>
  arrange(desc(MT > PT), desc(MT)) |>
  write_tsv(output1, col_names = FALSE)
