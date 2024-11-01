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
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle annotation table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-c", "--cdf"),
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
  output1 <- file.path(input_dir0, "2-depth.range.by.cdf.copy.number.txt")
  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
}

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# Sort the data by the 'Copy' column in decreasing order
cutoff_data <- x0 |>
  dplyr::select(Contig, Length, Depth, Copy, MT, PT, Edge) |>
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

write_tsv(xt, args1$cdf)

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
      min(x), (mean(x) - sd(x) * 2) / 1
    )) # CHECK # FIXME
  } else {
    # if it is too narrow in the distribution
    # the number 'divided by 3' is critical for C. pseudochinensis case
    depth_lower_bound <- round(mean(x) / 3)
  }
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
