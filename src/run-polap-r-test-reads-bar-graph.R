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
args1 <- parse_args(parser)

if (is_null(args1$input)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "output.tsv")
  output1 <- file.path(input_dir0, "output.pdf")

  args1 <- parse_args(parser, args = c("--input", input1, "--output", output1))
}



# Load the data from the combined TSV file using readr
data <- read_tsv(args1$input)


# Normalize the bases data to match the fragments scale
fragments_max <- max(data$fragments)
data <- data %>%
  mutate(Normalized_Bases = bases / max(bases) * fragments_max)

# Select only bases and fragments columns for plotting
plot_data <- data %>%
  select(size, fragments, Normalized_Bases) %>%
  pivot_longer(
    cols = c("fragments", "Normalized_Bases"),
    names_to = "Dataset", values_to = "Value"
  )

# Plot with categorical X-axis and grayscale fill
p1 <- ggplot(plot_data, aes(x = factor(size), y = Value, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Size", y = "Fragments Value") +
  scale_y_continuous(
    name = "Fragments",
    sec.axis = sec_axis(~ . * max(data$bases) / fragments_max, name = "Bases")
  ) +
  scale_fill_manual(values = c("fragments" = "gray70", "Normalized_Bases" = "gray40")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")
  )


# Save the plot as a PDF
ggsave(args1$output, plot = p1, device = "pdf", width = 8, height = 6)
