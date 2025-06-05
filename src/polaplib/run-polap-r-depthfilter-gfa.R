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
