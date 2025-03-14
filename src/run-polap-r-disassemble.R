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
suppressPackageStartupMessages(library("patchwork"))
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
  input1 <- file.path(input_dir0, "disassemble/infer-1/1/summary1.txt")
  output1 <- file.path(input_dir0, "out.txt")
  output2 <- file.path(input_dir0, "out.pdf")
  
  args1 <- parse_args(parser, args = c("--table", input1, 
                                       "--out", output1,
                                       "--plot", output2,
                                       "--coverage"))
  args1 <- parse_args(parser, args = c("--table", input1, 
                                       "--out", output1,
                                       "--plot", output2,
                                       "--scatter"))
  args1 <- parse_args(parser, args = c("--table", input1, 
                                       "--out", output1,
                                       "--plot", output2))
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

if (nrow(df) == 0) {
  quit(save = "no")
} else if (nrow(df) == 1) {
  write_tsv(df, args1$out)

  message <- paste("#mode: NA")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  message <- paste("#sd: NA")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  message <- paste("#index: 0")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  message <- paste("#n: 0")
  cat(message, file = args1$out, append = TRUE, sep = "\n")
  quit(save = "no")
}



# Define a function to remove outliers based on IQR (Interquartile Range)
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  return(x >= (Q1 - 1.5 * IQR) & x <= (Q3 + 1.5 * IQR))
}

# Filter out outliers in the length column
df_no_outliers <- df %>% filter(remove_outliers(length))

# Compute the density of the length
density_length <- density(df_no_outliers$length)

# Find the length corresponding to the highest point in the density curve (mode)
mode_length <- density_length$x[which.max(density_length$y)]

# Find the row with length nearest to Q3
selected_row_nearest_mode <- df_no_outliers %>%
  slice_min(abs(length - mode_length), n = 1)

sd_length <- sd(df_no_outliers$length)
n_length <- length(df_no_outliers$length)

# Print the mode length
print(paste("Length at the highest point in the density curve (mode):", round(mode_length)))
print(paste("Length standard deviation:", round(sd_length)))
print(paste("Row with the length nearest the highest point in the density curve (mode):", selected_row_nearest_mode$index))

# Print the selected row
print(selected_row_nearest_mode)

# Add a column for the absolute difference between mode_length and length
df_no_outliers <- df_no_outliers %>%
  mutate(diff = round(abs(length - mode_length)))

# Sort the rows by the absolute difference
df_sorted <- df_no_outliers %>% arrange(diff)

# Print the sorted data
print(df_sorted)

# Save the sorted data to a new TSV file
write_tsv(df_sorted, args1$out)

message <- paste("#mode:", round(mode_length))
cat(message, file = args1$out, append = TRUE, sep = "\n")
message <- paste("#sd:", round(sd_length))
cat(message, file = args1$out, append = TRUE, sep = "\n")
message <- paste("#index:", selected_row_nearest_mode$index)
cat(message, file = args1$out, append = TRUE, sep = "\n")
message <- paste("#n:", n_length)
cat(message, file = args1$out, append = TRUE, sep = "\n")

if (args1$coverage || args1$scatter) {
  # Plot: Scatter plot of length vs. coverage
  p1 <- ggplot(df_no_outliers, aes(x = length, y = coverage)) +
    geom_point(alpha = 0.6) +
    geom_point(data = selected_row_nearest_mode, aes(x = length, y = coverage), color = "red", size = 3) +
    labs(title = "Scatter Plot of Length vs. Coverage", x = "Length", y = "Coverage") +
    theme_minimal()
}

# Plot the trace plot of index vs. length with a horizontal line at mode_length
p2 <- ggplot(df_no_outliers, aes(x = index, y = length)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = selected_row_nearest_mode$length, color = "red", linetype = "dashed", size = 1) +
  labs(title = paste0("(A) Trace Plot: Index vs. Length (Mode: ", round(mode_length), ")"),
       x = "Index",
       y = "Length") +
  theme_minimal()

# Plot the histogram of length without outliers and indicate Q3 as a vertical line
p3 <- ggplot(df_no_outliers, aes(x = length)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = selected_row_nearest_mode$length, color = "red", linetype = "dashed", size = 1) +
  labs(title = paste0("(B) Histogram of Length (SD: ", round(sd_length), ")"),
       x = "Length",
       y = "Frequency") +
  theme_minimal()

# Plot the density curve and indicate the mode with a vertical line
p4 <- ggplot(df_no_outliers, aes(x = length)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  geom_vline(xintercept = selected_row_nearest_mode$length, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Density Curve of Length (Mode indicated)",
       x = "Length",
       y = "Density") +
  theme_minimal()

# Display the combined plot
# print(combined_plot)

if (args1$coverage) {
  # Combine the four plots into a single plot
  combined_plot <- (p1 | p2) / (p3 | p4)
  
  # Save the combined plot as a PDF
  ggsave(args1$plot, combined_plot, width = 12, height = 10)
} else if (args1$scatter) {
  ggsave(args1$plot, p1, width = 12, height = 10)
} else {
  # Combine the four plots into a single plot
  combined_plot <- (p2) / (p3)

  # Save the combined plot as a PDF
  ggsave(args1$plot, combined_plot, width = 12, height = 10)
}

