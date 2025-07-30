#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(imager)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)

if (interactive() && length(args) == 0) {
  args <- c("--dilate-size", "1")  # Set default or test values for RStudio
}

# === Command-line options ===
option_list <- list(
  make_option("--dilate-size", type = "integer", default = 3,
              help = "Size of dilation window (must be odd, e.g., 3, 5, 7) [default %default]")
)
#opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = FALSE)
opt <- parse_args(OptionParser(option_list = option_list), args = args)

# === Select input file interactively ===
#infile <- file.choose()
infile <- "/Volumes/run_thorne_hdd/h3/labshare/goshng/hifi1/Arabidopsis_thaliana-2/kmer/k2.umap.tsv"
output_prefix <- tools::file_path_sans_ext(basename(infile))

# === Load input UMAP/embedding file ===
df <- read_tsv(infile, col_types = cols())

# Standardize coordinate names
coords <- df[, 1:2]
colnames(coords) <- c("col1", "col2")

# === Compute bin2d density ===
binplot <- ggplot(coords, aes(col1, col2)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(count)))
bin_data <- ggplot_build(binplot)$data[[1]]

# Convert bins to a matrix grid
bin_data <- bin_data[bin_data$count > 0, ]
x_breaks <- sort(unique(c(bin_data$xmin, bin_data$xmax)))
y_breaks <- sort(unique(c(bin_data$ymin, bin_data$ymax)))

nx <- length(x_breaks) - 1
ny <- length(y_breaks) - 1
mat <- matrix(0, nrow = nx, ncol = ny)

# Map bins to matrix indices
for (i in seq_len(nrow(bin_data))) {
  xi <- findInterval(bin_data$xmin[i], x_breaks)
  yi <- findInterval(bin_data$ymin[i], y_breaks)
  if (xi >= 1 && xi <= nx && yi >= 1 && yi <= ny) {
    mat[xi, yi] <- bin_data$count[i]
  }
}

# === Convert to binary and dilate ===
binary <- mat > 0
img <- as.cimg(binary)
img_dilated <- imager::dilate_square(img, size = opt$`dilate-size`)
mask <- as.matrix(img_dilated) > 0

# === Label connected components ===
labeled <- as.matrix(label(img_dilated))
label_counts <- as.data.frame(table(labeled[labeled > 0]))
largest_label <- as.integer(label_counts$Var1[which.max(label_counts$Freq)])

# === Find grid centers of the largest component ===
idx <- which(labeled == largest_label, arr.ind = TRUE)
x_centers <- (x_breaks[idx[, 1]] + x_breaks[idx[, 1] + 1]) / 2
y_centers <- (y_breaks[idx[, 2]] + y_breaks[idx[, 2] + 1]) / 2

# === Match points in df to these grid regions ===
selected_flags <- mapply(function(x, y) {
  xi <- findInterval(x, x_breaks)
  yi <- findInterval(y, y_breaks)
  xi >= 1 && xi <= nx && yi >= 1 && yi <= ny && labeled[xi, yi] == largest_label
}, df[[1]], df[[2]])

df$selected <- ifelse(selected_flags, "Inside", "Outside")

# === Save selected/deselected ===
write_tsv(df[df$selected == "Inside", ], paste0(output_prefix, ".selected.tsv"))
write_tsv(df[df$selected == "Outside", ], paste0(output_prefix, ".deselected.tsv"))

# === Plot result ===
p <- ggplot(df, aes(x = df[[1]], y = df[[2]], color = selected)) +
  geom_point(alpha = 0.6, size = 0.7) +
  theme_minimal() +
  labs(title = paste0("Selected Region (Dilation size = ", opt$`dilate-size`, ")")) +
  scale_color_manual(values = c("Inside" = "red", "Outside" = "grey70"))

ggsave(paste0(output_prefix, ".selected_plot.pdf"), plot = p, width = 8, height = 6)

print(p)

cat("✅ Done. Selected points saved with prefix:", output_prefix, "\n")
cat("✅ Dnfile:", infile, "\n")
cat("Number of selected points (Inside):", sum(df$selected == "Inside"), "\n")
cat("Number of selected points (Outside):", sum(df$selected == "Outside"), "\n")
