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

parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle annotation table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
parser <- add_option(parser, c("--mixfit"),
  action = "store",
  help = "Output mixfit"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "3-depth.range.mixture.txt")
  output2 <- file.path(input_dir0, "3-mixfit.txt")
  # no output3 for mixfit file
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-o", output1,
    "--mixfit", output2
  ))
}

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# MT selection: x0 -> x1
# gene density cutoff: 1 in 100 kb

# depth range by mixture model
x1 <- x0 |> filter(Copy > 1)
m1 <- mixfit(x1$Copy, ncomp = 3, family = "gamma", max_iter = 30)
if (is_null(m1)) {
  # nothing for x1 if mixfit fails to converge.
  x0 <- x0 |> filter(Copy < 0)
  x1 <- x0
} else {
  # if (args1$`range-type` == 1) {
  #   x1 <- x0 |> filter(m1$mu[1] < Copy, Copy < m1$mu[2] + m1$sd[2] * 3)
  # } else if (args1$`range-type` == 2) {
  #   x1 <- x0 |> filter(m1$mu[1] < Copy, Copy < m1$mu[3] + m1$sd[3])
  # }

  # x1 <- x0 |> filter(m1$mu[1] < Copy, Copy < m1$mu[2] + m1$sd[2] * 3)
  x1 <- x0 |> filter(m1$mu[1] < Copy, Copy < m1$mu[3] + m1$sd[3])
}

if (nrow(x1) > 0) {
  tibble(
    depth_lower_bound = min(x1$V3) - 1,
    depth_upper_bound = max(x1$V3) + 1
  ) |>
    write_tsv(args1$out)
} else {
  tibble() |> write_tsv(args1$out) # create an empty file
}

# .mixfit
if (exists("m1")) {
  if (!is_null(m1)) {
    sink(args1$mixfit)
    print(m1)
    sink()
  }
}
