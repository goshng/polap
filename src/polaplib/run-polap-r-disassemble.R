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
# This script parses summary1.txt in the polap subcommand disassemble.
# Subcommand disassemble consists of three stages; each stage has summary text
# for the stage. The summary shows information about iterations.
# One of the iterations is selected for the next stage or the final choice.
# This script does this step. We use the length distribution to select one
# of the iterations so that we choose one with the length nearest the mode of
# the distribution.
#
# Example:
# Rscript ./src/polaplib/run-polap-r-disassemble.R \
#   --table input/summary1-3.txt \
#   --out output/summary1-3-ordered.txt \
#   --plot input/summary1-3-ordered.pdf
# output: summary1-3-ordered.txt
#
# TEST-SCC: try this out with different summary1.txt.
#
# See Also:
# run-polap-function-disassemble.sh
#
# Check: 2025-06-18
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

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
parser <- add_option(parser, c("-p", "--plot"),
  action = "store",
  help = "Plot file"
)
parser <- add_option(parser, c("-c", "--coverage"),
  action = "store_true",
  default = FALSE, help = "use alignment coverage"
)
parser <- add_option(parser, c("-s", "--scatter"),
  action = "store_true",
  default = FALSE, help = "use alignment coverage"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  input_dir0 <- file.path("Juncus_roemerianus")
  input_dir0 <- file.path("Eucalyptus_pauciflora")
  input1 <- file.path(input_dir0, "0/disassemble/infer-1/3-check/summary1.txt")
  
  input_dir0 <- file.path("input")
  input_dir1 <- file.path("output")
  input1 <- file.path(input_dir0, "summary1-0.txt")
  output1 <- file.path(input_dir1, "summary1-0-ordered.txt")
  output2 <- file.path(input_dir1, "summary1-0-ordered.pdf")
  
  input1 <- file.path(input_dir0, "summary1-1.txt")
  output1 <- file.path(input_dir1, "summary1-1-ordered.txt")
  output2 <- file.path(input_dir1, "summary1-1-ordered.pdf")

  input1 <- file.path(input_dir0, "summary1-2.txt")
  output1 <- file.path(input_dir1, "summary1-2-ordered.txt")
  output2 <- file.path(input_dir1, "summary1-2-ordered.pdf")
  
  input1 <- file.path(input_dir0, "summary1-2-1.txt")
  output1 <- file.path(input_dir1, "summary1-2-1-ordered.txt")
  output2 <- file.path(input_dir1, "summary1-2-1-ordered.pdf")
  
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "--out", output1,
    "--plot", output2,
    "--coverage"
  ))
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "--out", output1,
    "--plot", output2,
    "--scatter"
  ))
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "--out", output1,
    "--plot", output2
  ))
}

df <- read_tsv(args1$table, show_col_types = FALSE)

if (args1$coverage) {
  # Remove rows with length equal to zero
  df <- df %>%
    filter(length > 0) %>%
    mutate(coverage = (coverage_ref + coverage_target) / 2)
} else {
  df <- df %>%
    filter(length > 0)
}

# No iteration or single iteration
if (nrow(df) == 0) {
  print(paste("Number of iterations:", round(nrow(df))))
  quit(save = "no")
} else if (nrow(df) == 1) {
  print(paste("Number of iterations:", round(nrow(df))))
  print(df |> select(index, length))
  write_tsv(df, args1$out)

  message <- paste("#mode: NA")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  message <- paste("#sd: NA")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  # 2025-06-05
  # message <- paste("#index: 0")
  message <- paste("#index: -1")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  # 2025-06-05
  # message <- paste("#n: 0")
  message <- paste("#n: 1")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  quit(save = "no")
}

# We have at least 2 iterations in the data with length > 0.

# Function to remove outliers based on IQR (Interquartile Range)
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  return(x >= (Q1 - 1.5 * IQR) & x <= (Q3 + 1.5 * IQR))
}

# Filter out outliers in assembly length because the mode could be weird
# when we have outliers.
df_no_outliers <- df %>% filter(remove_outliers(length))

if (nrow(df_no_outliers) == 0) {
  print(paste("Number of iterations:", round(nrow(df_no_outliers))))
  quit(save = "no")
} else if (nrow(df_no_outliers) == 1) {
  print(paste("Number of iterations:", round(nrow(df_no_outliers))))
  print(df_no_outliers |> select(index, length))
  write_tsv(df_no_outliers, args1$out)
  
  message <- paste("#mode: NA")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  message <- paste("#sd: NA")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  # 2025-06-05
  # message <- paste("#index: 0")
  message <- paste("#index: -1")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  # 2025-06-05
  # message <- paste("#n: 0")
  message <- paste("#n: 1")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  quit(save = "no")
}

# what if df_no_outliers has no or one element?
# consider cases:
# 0, 1, 2, 3, ...

# Density of the lengths without outliers
density_length <- density(df_no_outliers$length)

# Length nearest the mode of the length distribution
mode_length <- density_length$x[which.max(density_length$y)]

# Date: 2025-06-17
# We comment this out because we get it after sorting by diff from the mode.
# The iteration with length nearest the length mode
#selected_row_nearest_mode <- df_no_outliers %>%
#  slice_min(abs(length - mode_length), n = 1)

# We now have the data to print out as a result.
# We need to consider whether the data has zero, 1, or more points.
# Depending on the count, we need to save differently.

# Column diff for the absolute difference between mode_length and length
df_no_outliers <- df_no_outliers %>%
  mutate(diff = round(abs(length - mode_length)))

# Sort the rows by the absolute difference
df_sorted <- df_no_outliers %>% arrange(diff)

# Check?
# We have a new row selected. This should be the same?
# The first iteration at top of the df_sorted is the selected one.
selected_row_nearest_mode <- df_sorted %>% slice(1)

sd_length <- sd(df_no_outliers$length)
n_length <- length(df_no_outliers$length)

# Print the selected row
# print(selected_row_nearest_mode)
# Print the mode length, SD, and the index selected.
print(paste("Number of iterations:", round(n_length)))
print(paste("Length at the highest point in the density curve (mode):", round(mode_length)))
print(paste("Length standard deviation:", round(sd_length)))
print(paste("Row with the length nearest the highest point in the density curve (mode):", selected_row_nearest_mode$index))

# Print the sorted data
print(df_sorted |> select(index, length))

# Now, we save the result in summary1-ordered.txt.
# Save the sorted data to a new TSV file: summary1-ordered.txt
write_tsv(df_sorted, args1$out)

message <- paste("#mode:", round(mode_length))
cat(message, file = args1$out, append = TRUE, sep = "\n")
message <- paste("#sd:", round(sd_length))
cat(message, file = args1$out, append = TRUE, sep = "\n")
message <- paste("#index:", selected_row_nearest_mode$index)
cat(message, file = args1$out, append = TRUE, sep = "\n")
message <- paste("#n:", n_length)
cat(message, file = args1$out, append = TRUE, sep = "\n")

# We do not use this plotting much.
if (args1$coverage || args1$scatter) {
  # Plot: Scatter plot of length vs. coverage
  p1 <- ggplot(df_no_outliers, aes(x = length, y = coverage)) +
    geom_point(alpha = 0.6) +
    geom_point(data = selected_row_nearest_mode, aes(x = length, y = coverage), color = "red", linewidth = 3) +
    labs(title = "Scatter Plot of Length vs. Coverage", x = "Length", y = "Coverage") +
    theme_minimal()
}

# Plot the trace plot of index vs. length with a horizontal line at mode_length
p2 <- ggplot(df_no_outliers, aes(x = index, y = length)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = selected_row_nearest_mode$length, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste0("(A) Trace Plot: Index vs. Length (Mode: ", round(mode_length), ")"),
    x = "Index",
    y = "Length"
  ) +
  theme_minimal()

# Plot the histogram of length without outliers and indicate Q3 as a vertical line
p3 <- ggplot(df_no_outliers, aes(x = length)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = selected_row_nearest_mode$length, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste0("(B) Histogram of Length (SD: ", round(sd_length), ")"),
    x = "Length",
    y = "Frequency"
  ) +
  theme_minimal()

# Plot the density curve and indicate the mode with a vertical line
p4 <- ggplot(df_no_outliers, aes(x = length)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  geom_vline(xintercept = selected_row_nearest_mode$length, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Density Curve of Length (Mode indicated)",
    x = "Length",
    y = "Density"
  ) +
  theme_minimal()

# Save them in a PDF
if (args1$coverage) {
  combined_plot <- (p1 | p2) / (p3 | p4)
  ggsave(args1$plot, combined_plot, width = 12, height = 10)
} else if (args1$scatter) {
  ggsave(args1$plot, p1, width = 12, height = 10)
} else {
  combined_plot <- (p2) / (p3)
  ggsave(args1$plot, combined_plot, width = 12, height = 10)
}
