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
suppressPackageStartupMessages(library("ggplot2"))
# args = commandArgs(trailingOnly=TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("-a", "--dista"),
  action = "store",
  help = "Summary table a: summary2.txt from the disassemble multiple runs",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-b", "--distb"),
                     action = "store",
                     help = "Summary table b: summary2.txt from disassemble single run",
                     metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output file"
)
args1 <- parse_args(parser)

if (is_null(args1$dista)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "jvalidus/disassemble/summary2.txt")
  input2 <- file.path(input_dir0, "jvalidus2/disassemble/summary2.txt")
  output1 <- file.path(input_dir0, "lower-percentage.txt")

  args1 <- parse_args(parser, args = c("--dista", input1, "--distb", input2, "-o", output1))
}

data1 <- read_tsv(args1$dista, show_col_types = FALSE)
data2 <- read_tsv(args1$distb, show_col_types = FALSE)

# Extract the 'length' column from 1.txt and remove zeros
filtered_lengths1 <- data1 |>
  filter(length > 0) |>
  pull(length)

# Extract the single 'length' value from 2.txt
length2 <- data2 |>
  pull(length)

# Compute the lower percentage
lower_percentage <- sum(filtered_lengths1 < length2[1]) / length(filtered_lengths1) * 100

# Display the result
write(lower_percentage, file = args1$out)

