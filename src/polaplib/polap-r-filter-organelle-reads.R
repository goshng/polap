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
# Copy this R template file to create a new R script that is used in polap.
# The following template file provides argument processing.
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))

# Define command-line options
option_list <- list(
  make_option(c("-t", "--table"), type = "character", help = "Path to annotation table (TSV)"),
  make_option(c("-l", "--length"),
    type = "numeric", default = 3e+7,
    help = "Maximum cumulative length of selected reads [default: %default]"
  ),
  make_option(c("-o", "--output"), type = "character", help = "Output file to write selected contig IDs"),
  make_option(c("--output-all"),
    type = "character", default = NULL,
    help = "Optional output file to write contigs passing dispersion filtering"
  ),
  make_option(c("--diff"),
    type = "numeric", default = 3,
    help = "minimum difference between MT and PT [default: %default]"
  ),
  make_option(c("--min-dispersion-pt"),
    type = "numeric", default = 1e+4,
    help = "Maximum dispersion_PT to retain contig (default: %default)"
  ),
  make_option(c("--min-dispersion-mt"),
    type = "numeric", default = 1e+5,
    help = "Random seed for reproducibility. If <= 0, seed from time."
  ),
  make_option(c("--rng-seed"),
    type = "integer", default = 42,
    help = "Maximum dispersion_MT to retain contig (default: %default)"
  ),
  make_option(c("--plot"),
    type = "character", default = NULL,
    help = "Optional output file (PDF/PNG) to save gene density scatter plot"
  )
)

opt_parser <- OptionParser(option_list = option_list)
args1 <- parse_args(opt_parser)

# Determine RNG seed
if (args1$`rng-seed` <= 0) {
  seed_value <- as.integer(as.numeric(Sys.time()) %% .Machine$integer.max)
  message("Using random RNG seed: ", seed_value)
} else {
  seed_value <- args1$`rng-seed`
  message("Using fixed RNG seed: ", seed_value)
}
set.seed(seed_value)

# Load input table
x0 <- read_tsv(args1$table, show_col_types = FALSE)

# Apply pseudo counts and dispersion filtering
filtered <- x0 %>%
  mutate(
    pseudo_MT = if_else(MT == 0, 1L, MT),
    pseudo_PT = if_else(PT == 0, 1L, PT),
    dispersion_MT = as.integer(Length / pseudo_MT),
    dispersion_PT = as.integer(Length / pseudo_PT)
  ) %>%
  select(-pseudo_MT, -pseudo_PT) %>%
  filter(!(MT <= PT & dispersion_PT > args1$`min-dispersion-pt`)) %>%
  filter(!(MT >= PT & dispersion_MT > args1$`min-dispersion-mt`)) %>%
  filter(MT > 0 | PT > 0) %>%
  filter(abs(MT - PT) > args1$diff)

# Optionally output all contigs after filtering
if (!is.null(args1$`output-all`)) {
  write_lines(filtered$Contig, args1$`output-all`)
}

# Optional plot
if (!is.null(args1$plot)) {
  p <- ggplot(filtered, aes(x = dispersion_MT, y = dispersion_PT)) +
    geom_point(alpha = 0.6, color = "#3366CC") +
    geom_hline(yintercept = args1$`min-dispersion-pt`, linetype = "dashed", color = "red") +
    geom_vline(xintercept = args1$`min-dispersion-mt`, linetype = "dashed", color = "red") +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Gene Density Filtering (Dispersion_MT vs Dispersion_PT)",
      x = "Dispersion_MT (Length / MT)",
      y = "Dispersion_PT (Length / PT)"
    ) +
    theme_minimal()
  ggsave(args1$plot, plot = p, width = 6, height = 5)
}

# Shuffle and select up to cumulative length
shuffled <- filtered[sample(nrow(filtered)), ]
selected <- shuffled %>%
  mutate(cumulative_length = cumsum(Length)) %>%
  filter(cumulative_length <= args1$length)

# Write selected contigs
write_lines(selected$Contig, args1$output)
