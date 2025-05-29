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

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("--mitochondrial"),
  action = "store_true",
  default = TRUE, help = "Mitochondrial genome assembly"
)
parser <- add_option(parser, c("--plastid"),
  action = "store_false",
  dest = "mitochondrial", help = "Plastid genome assembly"
)
parser <- add_option(parser, c("--copy-number"),
  type = "integer",
  default = -1,
  help = "Minimum contig copy number [default off]",
  metavar = "number"
)
parser <- add_option(parser, c("--gene-count"),
  type = "integer",
  default = -1,
  help = "Minimum gene count [default off]",
  metavar = "number"
)
parser <- add_option(parser, c("-c", "--compare-mt-pt"),
  action = "store_true",
  default = FALSE,
  help = "Compare MT and PT gene counts"
)
# max gene distance: 100 kb for mitochondrial and 10 kb for plastid
parser <- add_option(parser, c("-g", "--gene-density"),
  type = "integer",
  default = -1,
  help = "Minimum gene density or the number of gene per 1 Mb [default off; 10 for mitochondrial and 100 for plastid]",
  metavar = "number"
)
parser <- add_option(parser, c("-d", "--depth-range"),
  type = "character",
  default = "0,0",
  help = "Range of contig depth [default off]",
  metavar = "number"
)
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle annotation table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "1-annotation.txt")
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "--depth-range", "50,100",
    "-c",
    "-g", 10,
    "-o", output1
  ))
}

# filter by the depth range
depth.range <- as.numeric(strsplit(args1$`depth-range`, ",")[[1]])

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# MT selection: x0 -> x1
# gene density cutoff: 1 in 100 kb
x1 <- x0
if (args1$mitochondrial == TRUE) {
  # for mitochondria
  # filter by the depth range
  if (depth.range[1] > 0) {
    x1 <- x1 |> filter(depth.range[1] <= Depth, Depth <= depth.range[2])
  }

  # filter by the minimum copy number
  if (args1$`copy-number` > 0) {
    x1 <- x1 |> filter(Copy >= args1$`copy-number`)
  }

  # filter by the minimum MT gene count
  if (args1$`gene-count` > 0) {
    x1 <- x1 |> filter(MT >= args1$`gene-count`)
  }

  # depth range by mixture model
  # no code

  # filter by the MT/PT gene comparison among genes with MT > 0
  if (args1$`compare-mt-pt` == TRUE) {
    x1 <- x1 |> filter(MT > PT, MT > 0)
  }

  # filter by gene density: 1 MT gene per 1 Mb
  if (args1$`gene-density` > 0) {
    rt1 <- 1e+6 / args1$`gene-density`
    x1 <- x1 |>
      filter(MT > 0) |>
      mutate(RT = as.integer(Length / MT)) |>
      filter(RT < rt1)
  }
  # output: .annotated.txt
  if (nrow(x1) > 0) {
    x1 |>
      arrange(desc(MT > PT), desc(MT)) |>
      write_tsv(args1$out, col_names = FALSE)
  } else {
    tibble() |> write_tsv(args1$out)
  }
} else {
  # for plastid
  # filter by the depth range
  if (depth.range[1] > 0) {
    x1 <- x1 |> filter(depth.range[1] <= Depth, Depth <= depth.range[2])
  }

  # for plastid
  # filter by the minimum copy number
  if (args1$`copy-number` > 0) {
    x1 <- x1 |> filter(Copy >= args1$`copy-number`)
  }

  # for plastid
  # filter by the minimum PT gene count
  if (args1$`gene-count` > 0) {
    x1 <- x1 |> filter(PT >= args1$`gene-count`)
  }

  # depth range by mixture model
  # no code

  # for plastid
  # filter by the MT/PT gene comparison among genes with PT > 0
  if (args1$`compare-mt-pt` == TRUE) {
    x1 <- x1 |> filter(PT > MT, PT > 0)
  }

  # for plastid
  # filter by gene density: 1 MT gene per 100 kb
  if (args1$`gene-density` > 0) {
    # for plastid
    rt1 <- 1e+5 / args1$`gene-density`
    x1 <- x1 |>
      filter(PT > 0) |>
      mutate(RT = as.integer(Length / PT)) |>
      filter(RT < rt1)
  }
  # for plastid
  # output: .annotated.txt
  if (nrow(x1) > 0) {
    x1 |>
      arrange(desc(PT > MT), desc(PT)) |>
      write_tsv(args1$out, col_names = FALSE)
  } else {
    tibble() |> write_tsv(args1$out)
  }
}
