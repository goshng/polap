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
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.pdf")

  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
}

# output1 <- paste0(args1$out, ".table.tsv")
output1 <- paste0(args1$out)

x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

stats_copy <- x0 |>
  filter(MT > 0, MT > PT, Copy > 1) |>
  summarise(
    average_copy = mean(Copy, na.rm = TRUE),
    sd_copy = sd(Copy, na.rm = TRUE),
    cutoff_copy = pmax(average_copy - 3 * sd_copy, 1) # Ensure the result is at least 1
  )

x1 <- x0 |> filter(Copy > stats_copy$cutoff_copy)

# Create a curved plot of a histogram using density
p1 <- ggplot(x1, aes(x = Copy)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = "Curved Histogram with Density Plot",
    x = paste("Copy cutoffed by", round(stats_copy$cutoff_copy, digits = 1)),
    y = "Density"
  ) +
  theme_minimal()

# Save the plot as a PDF
ggsave(output1, plot = p1, device = "pdf", width = 8, height = 6)
