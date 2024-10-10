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
parser <- add_option(parser, c("-m", "--mitochondrial"),
  action = "store_true",
  default = TRUE, help = "Mitochondrial genome assembly"
)
parser <- add_option(parser, c("-p", "--plastid"),
  action = "store_false",
  dest = "mitochondrial", help = "Plastid genome assembly"
)
parser <- add_option(parser, c("-n", "--number-copy"),
  type = "integer",
  default = -1,
  help = "Minimum contig copy number [default off]",
  metavar = "number"
)
parser <- add_option(parser, c("-g", "--gene-count"),
  type = "integer",
  default = -1,
  help = "Minimum gene count [default off]",
  metavar = "number"
)
parser <- add_option(parser, c("-c", "--gene-compare"),
  action = "store_true",
  default = FALSE,
  help = "Compare MT and PT gene counts"
)
# max gene distance: 100 kb for mitochondrial and 10 kb for plastid
parser <- add_option(parser, c("-d", "--gene-density"),
  type = "integer",
  default = -1,
  help = "Minimum gene density or the number of gene per 1 Mb [default off; 10 for mitochondrial and 100 for plastid]",
  metavar = "number"
)
parser <- add_option(parser, c("-r", "--range-copy"),
  action = "store_true",
  default = FALSE,
  help = "Range of contig copy numbers using the copy number distribution"
)
parser <- add_option(parser, c("-s", "--save-range-copy"),
  action = "store_true",
  default = FALSE,
  help = "Save the range of contig copy numbers using the copy number distribution"
)
parser <- add_option(parser, c("-u", "--upper-number-copy"),
  type = "integer",
  default = -1,
  help = "Maximum contig copy number [default off]",
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
  s <- "Vigna_radiata"
  s <- "Brassica_rapa"
  s <- "Anthoceros_angustus"
  s <- "bioprojects"

  o <- "PRJNA597121"
  o <- "PRJEB79308"
  o <- "PRJDB10540a"
  o <- "PRJEB26621a"
  o <- "PRJNA644206-Populus_x_sibirica"
  o <- "PRJEB42431-Ophrys_insectifera_subsp._aymoninii"

  jnum <- "2"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input_dir1 <- file.path(input_dir0, jnum, "mtcontigs")
  output1 <- file.path(input_dir1, "1-mtcontig")

  if (jnum == "1" || jnum == "3") {
    args1 <- parse_args(parser, args = c("--table", input1, "--gene-compare", "--gene-density", 10, "-o", output1))
  } else if (jnum == "2" || jnum == "4") {
    args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare", "--gene-density", 10, "-o", output1))
  } else if (jnum == "5") {
    args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare", "--gene-density", 10, "--save-range-copy", "-o", output1))
  }
}

output2 <- paste0(args1$out, ".depth.stats.txt")

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# Sort the data by the 'Copy' column in decreasing order
cutoff_data <- x0 |>
  dplyr::select(Contig, Length, V3, Copy, MT, PT, Edge) |>
  arrange(desc(Copy)) |>
  # Calculate the cumulative sum of the 'Length' column
  mutate(cumulative_length = cumsum(Length)) |>
  # Cut off rows at the given cumulative length of 3,000,000
  filter(cumulative_length <= 3e+6) |>
  filter(MT > 0 | PT > 0)

# MT selection: x0 -> cutoff_data -> xt
# gene density cutoff: 1 in 100 kb
if (args1$mitochondrial == TRUE) {
  xt <- cutoff_data |>
    filter(PT < MT)
} else {
  xt <- cutoff_data |>
    filter(PT > MT)
}

mean1 <- mean(xt$V3)
sd1 <- sd(xt$V3)
while (mean1 < 1 * sd1 && nrow(xt) > 1) {
  xt <- xt |> filter(V3 < max(V3))
  mean1 <- mean(xt$V3)
  sd1 <- sd(xt$V3)
}

# xt -> .stats
if (nrow(xt) > 1) {
  xt |>
    summarise(
      depth_lower_bound = max(min(V3), mean(V3) - sd(V3) * 3),
      depth_upper_bound = min(max(V3), mean(V3) + sd(V3) * 3),
      depth_min = min(V3),
      depth_max = max(V3),
      depth_median = median(V3),
      depth_mean = mean(V3),
      depth_sd = sd(V3),
    ) |>
    write_tsv(output2)
} else if (nrow(xt) == 1) {
  xt |>
    summarise(
      depth_lower_bound = min(V3),
      depth_upper_bound = min(V3),
      depth_min = min(V3),
      depth_max = max(V3),
      depth_median = median(V3),
      depth_mean = mean(V3),
      depth_sd = sd(V3),
    ) |>
    write_tsv(output2)
} else {
  tibble() |> write_tsv(output2)
}
