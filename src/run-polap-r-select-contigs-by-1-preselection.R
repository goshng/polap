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

compare_and_choose_smaller <- function(x, y) {
  if (x < 0 && y < 0) {
    stop("Both values cannot be negative; at least one must be positive.")
  }
  smaller <- min(x, y)
  larger <- max(x, y)

  if (smaller < 0) {
    return(larger)
  } else {
    return(smaller)
  }
}

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
parser <- add_option(parser, c("-r", "--range-type"),
  type = "integer",
  default = 0,
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
parser <- add_option(parser, c("-x", "--depth"),
  action = "store",
  help = "Output depth range"
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
  jnum <- 1
  if (jnum == 1) {
    input_dir0 <- file.path(".")
    input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
    output1 <- file.path(input_dir0, "1-preselection.by.gene.density.txt")
    output2 <- file.path(input_dir0, "1-depth.range.by.gene.density.txt")
    # no output3 for mixfit file
    args1 <- parse_args(parser, args = c("--table", input1, "-c", "-d", 10, "-o", output1, "--depth", output2))
  } else if (jnum == 2) {
    input_dir0 <- file.path(".")
    input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
    output1 <- file.path(input_dir0, "1-preselection.by.depth.mixture.txt")
    output2 <- file.path(input_dir0, "1-depth.range.by.depth.mixture.txt")
    output3 <- file.path(input_dir0, "1-mixfit.txt")
    args1 <- parse_args(parser, args = c("--table", input1, "-c", "-d", 10, "-r", 2, "--out", output1, "--depth", output2, "--mixfit", output3))
  }
}

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# MT selection: x0 -> x1
# gene density cutoff: 1 in 100 kb
x1 <- x0
if (args1$mitochondrial == TRUE) {
  # filter by the minimum copy number
  if (args1$`number-copy` > 0) {
    x1 <- x1 |> filter(Copy >= args1$`number-copy`)
  }

  # filter by the minimum MT gene count
  if (args1$`gene-count` > 0) {
    x1 <- x1 |> filter(MT >= args1$`gene-count`)
  }

  # depth range by mixture model
  if (args1$`range-type` > 0) {
    x2 <- x1 |> filter(Copy > 1)
    m1 <- mixfit(x2$Copy, ncomp = 3, family = "gamma", max_iter = 30)
    if (is_null(m1)) {
      # nothing for x1 if mixfit fails to converge.
      x1 <- x1 |> filter(Copy < 0)
      x2 <- x1
    } else {
      if (args1$`range-type` == 1) {
        x2 <- x1 |> filter(m1$mu[1] < Copy, Copy < m1$mu[2] + m1$sd[2] * 3)
      } else if (args1$`range-type` == 2) {
        x2 <- x1 |> filter(m1$mu[1] < Copy, Copy < m1$mu[3] + m1$sd[3])
      }
      x1 <- x2
    }
  }

  # filter by the MT/PT gene comparison among genes with MT > 0
  if (args1$`gene-compare` == TRUE) {
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
} else {
  x1 <- x1 |> filter(Copy > args1$`number-copy`)
  x1 <- x1 |> filter(PT >= args1$`gene-count`)
  if (args1$`range-copy` == TRUE) {
    x2 <- x1 |> filter(Copy > 1)
    m1 <- mixfit(x2$Copy, ncomp = 3, family = "gamma")
    if (is_null(m1)) {
      x1 <- x1 |> filter(Copy < 0)
      x2 <- x1
    } else {
      x2 <- x1 |> filter(Copy > m1$mu[3])
      x1 <- x2
    }
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

# list by edge
x1 <- x1 |>
  separate_rows(Edge, sep = ",") |>
  # Convert the numbers to absolute values
  mutate(Edge = abs(as.numeric(Edge))) |>
  # Group by the absolute value of Edge
  group_by(Edge) |>
  ungroup() |>
  distinct() |>
  dplyr::select(-V4, -V5, -V7, -V8)

# output: .annotated.txt
if (nrow(x1) > 0) {
  x1 |>
    mutate(edgename = paste0("edge_", Edge)) |>
    relocate(edgename) |>
    arrange(desc(MT > PT), desc(MT)) |>
    distinct() |>
    write_tsv(args1$out, col_names = FALSE)
} else {
  tibble() |> write_tsv(args1$out)
}

# outfile: .stats
if (args1$`range-type` == 0) {
  if (nrow(x1) > 0) {
    tibble(
      depth_lower_bound = max(min(x1$V3) - 1, 1), # CHECK # FIXME
      depth_upper_bound = max(x1$V3) + 1
    ) |>
      write_tsv(args1$depth)
  } else {
    tibble() |> write_tsv(args1$depth) # create an empty file
  }
} else {
  if (nrow(x2) > 0) {
    tibble(
      depth_lower_bound = min(x2$V3) - 1,
      depth_upper_bound = max(x2$V3) + 1
    ) |>
      write_tsv(args1$depth)
  } else {
    tibble() |> write_tsv(args1$depth) # create an empty file
  }
}

# .mixfit
if (exists("m1")) {
  if (!is_null(m1)) {
    sink(args1$mixfit)
    print(m1)
    sink()
  }
}
