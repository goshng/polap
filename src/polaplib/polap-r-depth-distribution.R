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

################################################################################
# Create a figure of depth distribution of edges with Copy values greater than
# mean - 3xsd so that we exclude too extreme depths.
#
# Example:
# polap depth-distribution -o Spirodela_polyrhiza/o1
# Spirodela_polyrhiza/o1/0/assembly_info_organelle_annotation_count-all.txt.pdf
#
# Code:
# Rscript "${_POLAPLIB_DIR}"/polap-r-depth-distribution.R \
# 	-t "${_polap_var_ga_annotation_all}" \
# 	-o "${_polap_var_ga_annotation_all}".pdf \
#
# input: assembly_info_organelle_annotation_count-all.txt
# output: assembly_info_organelle_annotation_count-all.txt.pdf
#
# See Also:
# polap-cmd-annotate.sh
# function _run_polap_depth-distribution
#
# References:
# https://how.dev/answers/what-is-the-pmax-function-in-r
#
# TODO: rename: polap-r-depth-distribution.R
#
# Check: 2025-06-17
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

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

# We have a default set of command-line options for debugging purpose.
# Put the input file at the current directory and source this R script
# to execute this R script. Then, it will read the input gene
# annotation table to create a PDF file.
# input: assembly_info_organelle_annotation_count-all.txt
if (is_null(args1$table)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.pdf")

  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
}

# Set the output filename.
# output1 <- paste0(args1$out, ".table.tsv")
output1 <- paste0(args1$out)

# The gene annotation table is delimited by a space.
x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)

# Compute the cut-off for excluding edges with too low depths
# Cut-off is Average - 3 x SD. We use edges with depths greater than the cut-off.
stats_copy <- x0 |>
  filter(MT > 0, MT > PT, Copy > 1) |>
  summarise(
    average_copy = mean(Copy, na.rm = TRUE),
    sd_copy = sd(Copy, na.rm = TRUE),
    cutoff_copy = pmax(average_copy - 3 * sd_copy, 1)
  )
# Tips: https://how.dev/answers/what-is-the-pmax-function-in-r
# a <- c(10, 8, 3, 9, 0, 5)
# b <- c(15, 4, 6, 9, 8, 4)
# pmax(a, b)
# [1] 15  8  6  9  8  5

# Filter out edges with too low depths
x1 <- x0 |> filter(Copy > stats_copy$cutoff_copy)

# A plot of a histogram of edge depths of the gene annotation table
p1 <- ggplot(x1, aes(x = Copy)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = "Histogram of edge depths of the gene annotation table",
    x = paste("Copy or depths cutoffed by", round(stats_copy$cutoff_copy, digits = 1)),
    y = "Density"
  ) +
  theme_minimal()

# Save the plot as a PDF
ggsave(output1, plot = p1, device = "pdf", width = 8, height = 6)
