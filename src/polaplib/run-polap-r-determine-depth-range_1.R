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

# polaplib/run-polap-r-determine-depth-range-1.R
# Check: 2025-06-17

################################################################################
# This script selects seed contigs using a gene annotation table.
# We have multiple scripts with a similar name except for the last number.
# Each script has a different scheme to select seed contigs.
# Although we have option for plastid flag, this script is for mitochondrial
# only. See polap-r-plastid-determine-depth-range.R for plastid seed selection.
# See the code for detail.
# This is created from the template: polap-r-determine-depth-range.R.
#
# input:
# assembly_info_organelle_annotation_count-all.txt
#
# output:
# 2-depth.range.by.cdf.copy.number.txt
# contig-annotation-cdf-table.txt
#
# TODO: rename: polap-r-determine-depth-range_1.R
#
# See Also:
# polap-function-disassemble-seeds.sh
# polap-cmd-seeds.sh
#
# Check: 2025-06-17
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
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle annotation table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-c", "--cdf"),
  action = "store",
  help = "Organelle annotation table sorted by depths",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename",
  metavar = "<FILE>"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "2-depth.range.by.cdf.copy.number.txt")
  output2 <- file.path(input_dir0, "contig-annotation-cdf-table.txt")
  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1, "-c", output2))
}

# Steps;
# ptDNA contigs are used for the lower bound of ptDNA contig depth.
# filter by the lower bound of ptcontig.
# gene dispersion is a measure of how far two genes are apart.
# filter by pt gene dispersion of 10 kb
# filter by mt gene dispersion of 100 kb
# filter in upto 3 Mb of CDF among the top depth contigs
#
# then
# filter by PT < MT, Depth < ptDNA lower bound, Copy > 0
#
# make the data so that mean < sd;
# remove large copy number contig so that SD is less than MEAN.
#
# if multiple rows:
# lower bound: larger one between mean - 2xsd or min
# upper bound: 3 times the smaller one between mean + 2xsd or max
# otherwise:
# bound is equal to the single value.

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# PT > MT: the top most three -> max depth
# this max depth -> range or lower baound of ptDNA
# FIXME: test this with Anthoroces angustus
# FIXME: Now, this would not work for ptDNA.
pt_lower_bound <- x0 |>
  dplyr::select(Contig, Length, Depth, Copy, MT, PT, Edge) |>
  filter(PT > MT) |>
  arrange(desc(PT)) |>
  slice_head(n = 3) |>
  summarise(max_depth = max(Depth, na.rm = TRUE) * 0.9) |>
  pull(max_depth)

# Sort the data by the 'Copy' column in decreasing order
cutoff_data <- x0 |>
  dplyr::select(Contig, Length, Depth, Copy, MT, PT, Edge) |>
  arrange(desc(Depth)) |>
  mutate(
    pseudo_MT = if_else(MT == 0, MT + 1, MT),
    pseudo_PT = if_else(PT == 0, PT + 1, PT)
  ) |>
  # Calculate the cumulative sum of the 'Length' column
  filter(Depth < pt_lower_bound) |>
  mutate(dispersion_MT = as.integer(Length / pseudo_MT)) |>
  mutate(dispersion_PT = as.integer(Length / pseudo_PT)) |>
  select(-pseudo_MT, -pseudo_PT) |>
  filter(!(MT <= PT & dispersion_PT > 1e+4)) |>
  filter(!(MT >= PT & dispersion_MT > 1e+5)) |>
  # Cut off rows at the given cumulative length of 3,000,000
  filter(MT > 0 | PT > 0) |>
  mutate(cumulative_length = cumsum(Length)) |>
  filter(cumulative_length <= 3e+6)

# MT selection: x0 -> cutoff_data -> xt
# gene density cutoff: 1 in 100 kb
if (args1$mitochondrial == TRUE) {
  xt <- cutoff_data |>
    filter(
      PT < MT,
      Depth < pt_lower_bound,
      Copy > 0
    )
} else {
  xt <- cutoff_data |>
    filter(
      PT > MT,
      Depth > pt_lower_bound,
      Copy > 0
    )
}

# remove large copy number contig so that SD is less than MEAN.
mean1 <- mean(xt$Depth)
sd1 <- sd(xt$Depth)
while (mean1 < 1 * sd1 && nrow(xt) > 1) {
  xt <- xt |> filter(Depth < max(Depth))
  mean1 <- mean(xt$Depth)
  sd1 <- sd(xt$Depth)
}

write_tsv(xt, args1$cdf)

# Carex pseudochinensis
lbound <- function(x) {
  depth_lower_bound <- 0
  if (sd(x) > 1) {
    depth_lower_bound <- round(max(
      min(x), mean(x) - sd(x) * 2
    )) # CHECK # FIXME
  } else {
    # if it is too narrow in the distribution
    # the number 'divided by 3' is critical for C. pseudochinensis case
    depth_lower_bound <- round(mean(x) / 3)
  }
  return(depth_lower_bound)
}

# lower bound: 1/3 of the original value
# upper bound: x3 of the original value
# xt -> .stats
if (nrow(xt) > 1) {
  xt |>
    summarise(
      # 2024-10-19
      # depth_lower_bound = round(max(min(Depth), mean(Depth) - sd(Depth) * 3) / 3), # CHECK # FIXME
      # depth_upper_bound = round(min(max(Depth), mean(Depth) + sd(Depth) * 3) * 3),
      depth_lower_bound = lbound(Depth), # CHECK # FIXME
      depth_upper_bound = round(min(
        max(Depth), mean(Depth) + sd(Depth) * 2
      ) * 3),
      depth_min = min(Depth),
      depth_min2 = (mean(Depth) - sd(Depth) * 2) / 1,
      depth_max = max(Depth) * 3,
      depth_max2 = (mean(Depth) + sd(Depth) * 2) * 3,
      depth_median = median(Depth),
      depth_mean = mean(Depth),
      depth_sd = sd(Depth),
    ) |>
    write_tsv(args1$out)
} else if (nrow(xt) == 1) {
  xt |>
    summarise(
      depth_lower_bound = min(Depth),
      depth_upper_bound = min(Depth),
      depth_min = min(Depth),
      depth_min2 = 0,
      depth_max = max(Depth),
      depth_max2 = 0,
      depth_median = median(Depth),
      depth_mean = mean(Depth),
      depth_sd = sd(Depth),
    ) |>
    write_tsv(args1$out)
} else {
  tibble() |> write_tsv(args1$out)
}
