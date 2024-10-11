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
parser <- add_option(parser, c("-c", "--copy"),
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
  input2 <- file.path(input_dir0, "1-manual.copy.range.txt")
  output1 <- file.path(input_dir0, "1-manual.depth.range.txt")
  args1 <- parse_args(parser, args = c("--table", input1, "--copy", input2, "-o", output1))
  
}

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)
x1 <- scan(args1$copy)

x0 |>
  filter(x1[1] <= Copy, Copy <= x1[2]) |>
  dplyr::select(V3) |>
  summarise(
    depth_lower_bound = min(V3),
    depth_upper_bound = max(V3)
  ) |>
  write_tsv(args1$out, col_names = TRUE)
