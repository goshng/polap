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
parser <- add_option(parser, c("-m", "--mitochondrial"),
  action = "store_true",
  default = TRUE, help = "Mitochondrial genome assembly"
)
parser <- add_option(parser, c("-p", "--plastid"),
  action = "store_false",
  dest = "mitochondrial", help = "Plastid genome assembly"
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
parser <- add_option(parser, c("-n", "--number-copy"),
  type = "integer",
  default = -1,
  help = "Minimum contig copy number [default off]",
  metavar = "number"
)
parser <- add_option(parser, c("-u", "--upper-number-copy"),
  type = "integer",
  default = -1,
  help = "Maximum contig copy number [default off]",
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

  # o <- "PRJNA914763" # type 3
  o <- "PRJNA597121" # type 5 + 3
  # o <- "PRJNA766769"    # type 5 + 3
  # o <- "PRJNA635654"    # type 5 + 3
  # o <- "PRJNA646849" # type 5 + 3, type 3
  # o <- "PRJEB26621"     # type 5 + 3, type 3
  # o <- "PRJNA1074514"   # type 5 + 3
  # o <- "PRJEB19787"     # no mtDNA seeds
  # o <- "PRJEB21270" # type 3; type 5 Error: insufficient data for model estimation.
  o <- "PRJEB79308"
  o <- "PRJDB10540a"
  o <- "PRJEB26621a"
  o <- "PRJNA644206-Populus_x_sibirica"

  jnum <- "2"
  input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
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
  if (jnum == "1" || jnum == "3") {
    args1 <- parse_args(parser, args = c("--table", input1, "--gene-compare", "--gene-density", 10, "-o", output1))
  } else if (jnum == "2" || jnum == "4") {
    args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare", "--gene-density", 10, "-o", output1))
  } else if (jnum == "5") {
    args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare", "--gene-density", 10, "--save-range-copy", "-o", output1))
  }

  # type 3 + depth graph
  # type 5 + 3 + depth graph
  # type 5 + 3 + graph (using depth range from the distribution)
  #
  # Gene count only
  # args1 <- parse_args(parser, args = c("--table", input1, "--gene-count", 1))
  #
  # type 1
  # args1 <- parse_args(parser, args = c("--table", input1, "--gene-compare"))
  # type 2
  # args1 <- parse_args(parser, args = c("--table", input1, "--gene-density", 10))
  # type 3
  # args1 <- parse_args(parser, args = c("--table", input1, "--gene-compare", "--gene-density", 10, "-o", output1))
  #
  # Depth only
  # type 4
  # args1 <- parse_args(parser, args = c("--table", input1, "--number-copy", 2))
  #
  # Depth & Gene count
  # type 4 + 1
  # args1 <- parse_args(parser, args = c("--table", input1, "--number-copy", 2, "--gene-compare"))
  # type 4 + 2
  # args1 <- parse_args(parser, args = c("--table", input1, "--number-copy", 2, "--gene-density", 10))
  # type 4 + 3
  # args1 <- parse_args(parser, args = c("--table", input1, "--number-copy", 2, "--gene-compare", "--gene-density", 10))
  #
  # type 5
  # args1 <- parse_args(parser, args = c("--table", input1, "--range-copy"))
  # type 5 + 4
  # args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--number-copy", 1))
  # type 5 + 1
  # args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare"))
  # type 5 + 2
  # args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-density", 10))
  # type 5 + 3
  # args1 <- parse_args(parser, args = c("--table", input1, "--range-copy", "--gene-compare", "--gene-density", 10, "--save-range-copy", "-o", output1))

  # args1 <- parse_args(parser, args = c(
  #   "--table", input1,
  #   "--gene-count", 1,
  #   "--gene-compare",
  #   "--gene-density", 10,
  #   "--number-copy", 1,
  #   "--range-copy"
  # ))
  #
}

output1 <- paste0(args1$out, ".annotated.txt")
output2 <- paste0(args1$out, ".stats.txt")
output3 <- paste0(args1$out, ".mixfit.txt")

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

xall <- x0 |>
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  ungroup() |>
  distinct() |>
  dplyr::select(-V4, -V5, -V7, -V8)

# Sort the data by the 'Copy' column in decreasing order
cutoff_data <- x0 |>
  dplyr::select(Contig, Length, V3, Copy, MT, PT, Edge) |>
  arrange(desc(Copy)) |>
  # Calculate the cumulative sum of the 'Length' column
  mutate(cumulative_length = cumsum(Length)) |>
  # Cut off rows at the given cumulative length of 3,000,000
  filter(cumulative_length <= 3e+6) |>
  filter(MT > 0 | PT > 0)

# MT selection
# gene density cutoff: 1 in 100 kb
x1 <- x0

if (args1$mitochondrial == TRUE) {
  xt <- cutoff_data |>
    filter(PT < MT)

  if (args1$`number-copy` > 0) {
    x1 <- x1 |> filter(Copy > args1$`number-copy`)
  }
  if (args1$`range-copy` == TRUE) {
    x2 <- x1 |> filter(Copy > 1)
    m1 <- mixfit(x2$Copy, ncomp = 3, family = "gamma")
    # m1 <- NULL
    if (is_null(m1)) {
      # nothing for x1 if mixfit fails to converge.
      x1 <- x1 |> filter(Copy < 0)
    } else {
      # x2 <- x1 |> filter(m1$mu[1] < Copy, Copy < m1$mu[3])
      x2 <- x1 |> filter(m1$mu[1] < Copy, Copy < m1$mu[3] + m1$sd[3])
      x1 <- x2
    }
  }
  if (args1$`gene-count` > 0) {
    x1 <- x1 |> filter(MT >= args1$`gene-count`)
  }
  if (args1$`gene-compare` == TRUE) {
    x1 <- x1 |> filter(MT >= PT, MT > 0)
  }
  if (args1$`gene-density` > 0) {
    rt1 <- 1e+6 / args1$`gene-density`
    x1 <- x1 |>
      filter(MT > 0) |>
      mutate(RT = as.integer(Length / MT)) |>
      filter(RT < rt1)
  }
} else {
  xt <- cutoff_data |>
    filter(PT > MT)

  if (args1$`number-copy` > 0) {
    x1 <- x1 |> filter(Copy > args1$`number-copy`)
  }
  if (args1$`range-copy` == TRUE) {
    x2 <- x1 |> filter(Copy > 1)
    m1 <- mixfit(x2$Copy, ncomp = 3, family = "gamma")
    if (is_null(m1)) {
      x1 <- x1 |> filter(Copy < 0)
    } else {
      x2 <- x1 |> filter(Copy > m1$mu[3])
      x1 <- x2
    }
  }
  if (args1$`gene-count` > 0) {
    x1 <- x1 |> filter(PT >= args1$`gene-count`)
  }
  if (args1$`gene-compare` == TRUE) {
    x1 <- x1 |> filter(PT >= MT, PT > 0)
  }
  if (args1$`gene-density` > 0) {
    rt1 <- 1e+6 / args1$`gene-density`
    x1 <- x1 |>
      filter(PT > 0) |>
      mutate(RT = as.integer(Length / PT)) |>
      filter(RT < rt1)
  }
}

x1 <- x1 |>
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  ungroup() |>
  distinct() |>
  dplyr::select(-V4, -V5, -V7, -V8)

# print(x1)

# if (is_null(args1$out)) {
#   args1$out <- output1
# }

# select(edgename, Length, V3, Copy, MT, PT, RT) |>

if (nrow(x1) > 0) {
  x1 |>
    mutate(edgename = paste0("edge_", Edge)) |>
    relocate(edgename) |>
    arrange(desc(MT > PT), desc(MT)) |>
    distinct() |>
    write_tsv(output1, col_names = FALSE)
} else {
  tibble() |> write_tsv(output1)
}

if (args1$`save-range-copy` == TRUE) {
  if (nrow(x1) > 0) {
    tibble(
      depth_lower_bound = min(x1$V3) - 1,
      depth_upper_bound = max(x1$V3) + 1
    ) |>
      write_tsv(output2)
  } else {
    tibble() |> write_tsv(output2) # create an empty file
  }
} else {
  mean1 <- mean(xt$V3)
  sd1 <- sd(xt$V3)
  while (mean1 < 1 * sd1 && nrow(xt) > 1) {
    xt <- xt |> filter(V3 < max(V3))
    mean1 <- mean(xt$V3)
    sd1 <- sd(xt$V3)
  }

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
}

if (exists("m1")) {
  if (!is_null(m1)) {
    sink(output3)
    print(m1)
    sink()
  }
}
