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

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"),
  action = "store",
  help = "Input file: x, fragments, bases, depth",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--output"),
  action = "store",
  help = "Output graph"
)
parser <- add_option(parser, c("-s", "--sizes"),
  action = "store",
  type = "character",
  help = "Comma-separated list of sizes, e.g., '5000,7000,9000'"
)
args1 <- parse_args(parser)

if (is_null(args1$input)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "summary.tsv")
  output1 <- file.path(input_dir0, "summary.pdf")
  input2 <- "5000,7000,9000,11000,13000"

  args1 <- parse_args(parser, args = c("--input", input1, "--output", output1, "--sizes", input2))
}

# Split sizes into a vector if provided
if (!is.null(args1$sizes)) {
  selected_sizes <- as.numeric(unlist(strsplit(args1$sizes, ",")))
} else {
  selected_sizes <- NULL
  cat("No sizes provided.\n")
}

# Load the data from the combined TSV file using readr
data <- read_tsv(args1$input, show_col_types = FALSE)

# Vector with selected sizes
# selected_sizes <- c(5000, 7000, 9000, 11000, 13000, 15000, 17000)

# Filter the data to only include rows with the selected sizes
if (is.null(selected_sizes)) {
  filtered_data <- data
} else {
  filtered_data <- data %>% filter(size %in% selected_sizes)
}

# Normalize bases and depth by the maximum fragments value
fragments_max <- max(data$fragments)
data <- filtered_data %>%
  mutate(
    Normalized_Bases = bases / max(bases) * fragments_max,
    Normalized_Depth = depth / max(depth) * fragments_max
  )

# Select all three metrics for plotting
plot_data <- data %>%
  select(size, fragments, Normalized_Bases, Normalized_Depth) %>%
  pivot_longer(
    cols = c("fragments", "Normalized_Bases", "Normalized_Depth"),
    names_to = "Dataset", values_to = "Value"
  )

# Plot with categorical X-axis, using bars for each property and adding text labels to depth bars
p1 <- ggplot(plot_data, aes(x = factor(size), y = Value, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = subset(plot_data, Dataset == "Normalized_Depth"),
    aes(label = round(data$depth, 0)),
    vjust = -0.3,
    hjust = -1.5,
    size = 3.5
  ) +
  labs(x = "Size", y = "Fragments Value") +
  scale_y_continuous(
    name = "Fragments",
    sec.axis = sec_axis(~ . * max(data$bases) / fragments_max, name = "Bases")
  ) +
  scale_fill_manual(values = c("fragments" = "gray80", "Normalized_Bases" = "gray70", "Normalized_Depth" = "gray60")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")
  )

# Save the plot as a PDF
ggsave(args1$output, plot = p1, device = "pdf", width = 8, height = 6)
