#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
output_prefix <- args[2]
k <- as.integer(args[3])
label_file <- ifelse(length(args) >= 4, args[4], NA)
ellipse_size <- ifelse(length(args) >= 5, args[5], "medium")
get_ellipse_params <- function(size = "medium") {
  size <- tolower(size)
  switch(size,
    "tight" = list(top_percent = 0.10, ellipse_level = 0.90, scale_factor = 1.0),
    "medium" = list(top_percent = 0.20, ellipse_level = 0.95, scale_factor = 1.1),
    "loose" = list(top_percent = 0.30, ellipse_level = 0.99, scale_factor = 1.2),
    { # fallback/default
      warning("Unknown ellipse size preset. Using 'medium'.")
      list(top_percent = 0.20, ellipse_level = 0.95, scale_factor = 1.1)
    }
  )
}

params <- get_ellipse_params(ellipse_size)
top_percent <- params$top_percent
ellipse_level <- params$ellipse_level
scale_factor <- params$scale_factor

# conda install -y r-pracma r-patchwork r-rmarkdown r-ellipse r-sp r-robustbase r-mass bioconductor-genomeinfodb r-uwot r-ggforce r-r.utils r-biocmanager genomescope2 bbmap pbsim3

library(Biostrings)
library(ggplot2)
library(gridExtra)
# library(Rtsne)
library(uwot)
library(readr)
library(R.utils)
suppressPackageStartupMessages(library(MASS)) # for ellipse fitting
suppressPackageStartupMessages(library(ggforce)) # optional: alternative ellipse display
suppressPackageStartupMessages(library(ellipse))
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(robustbase)) # for `cov.rob`


is_fasta <- grepl(".fa$|.fasta$", infile, ignore.case = TRUE)

cat("📖 Reading reads...\n")
if (grepl(".gz$", infile)) {
  tmp <- tempfile(fileext = if (is_fasta) ".fa" else ".fq")
  gunzip(infile, tmp, overwrite = TRUE, remove = FALSE)
  infile <- tmp
}

reads <- if (is_fasta) {
  readDNAStringSet(infile)
} else {
  readDNAStringSet(infile, format = "fastq")
}
names(reads) <- sub(" .*", "", names(reads))

cat("🔢 Counting ", k, "-mers...\n")
kmers <- oligonucleotideFrequency(reads, width = k)
rownames(kmers) <- names(reads)
abund <- rowSums(kmers)

# Optional: read labels
if (!is.na(label_file) && file.exists(label_file)) {
  labels <- read_tsv(label_file, col_names = c("id", "label"))
  labels <- labels[labels$id %in% rownames(kmers), ]
  kmers <- kmers[labels$id, , drop = FALSE]
  abund <- abund[labels$id]
  group <- labels$label
} else {
  group <- rep("Unlabeled", nrow(kmers))
}

# === PCA
cat("🧬 PCA...\n")
pca <- prcomp(kmers, center = TRUE, scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], label = group, abundance = abund, id = rownames(kmers))
write_tsv(pca_df, paste0(output_prefix, ".pca.tsv"))

# === UMAP
cat("🧬 UMAP...\n")
umap_xy <- umap(kmers, n_neighbors = 15)
umap_df <- data.frame(UMAP1 = umap_xy[, 1], UMAP2 = umap_xy[, 2], label = group, abundance = abund, id = rownames(kmers))
write_tsv(umap_df, paste0(output_prefix, ".umap.tsv"))

# === Plot all
p1 <- ggplot(pca_df, aes(PC1, PC2, color = abundance)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_viridis_c() +
  labs(title = "PCA", color = "k-mer abundance")

p2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = abundance)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_viridis_c() +
  labs(title = "UMAP", color = "k-mer abundance")

p1_label <- ggplot(pca_df, aes(PC1, PC2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (by label)", color = "Label")

p2_label <- ggplot(umap_df, aes(UMAP1, UMAP2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP (by label)", color = "Label")

pdf(paste0(output_prefix, ".plot.pdf"), width = 10, height = 8)
grid.arrange(p1, p2, ncol = 2)
dev.off()

pdf(paste0(output_prefix, ".label_plot.pdf"), width = 10, height = 8)
grid.arrange(p1_label, p2_label, ncol = 2)
dev.off()


# === Heatmap Plots with Contours ===
p1_density0 <- ggplot(pca_df, aes(PC1, PC2)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(count))) +
  geom_density_2d(color = "white", size = 0.3) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() +
  labs(title = "PCA Density + Contours", fill = "Count")

# === Step 1: Get bin2d stats ===
bin_stats <- ggplot_build(
  ggplot(pca_df, aes(PC1, PC2)) +
    stat_bin2d(bins = 100)
)$data[[1]]

# === Step 2: Keep nonzero bins only ===
nonzero_bins <- bin_stats[bin_stats$count > 0, ]

# === Step 3: Top 10% threshold ===
# threshold <- quantile(nonzero_bins$count, 0.75) # size of ellipse
threshold <- quantile(nonzero_bins$count, 1 - top_percent)


# === Step 4: Extract centroids of top bins ===
top_bins <- nonzero_bins[nonzero_bins$count > threshold, ]
top_bin_centers <- data.frame(
  x = (top_bins$xmin + top_bins$xmax) / 2,
  y = (top_bins$ymin + top_bins$ymax) / 2
)

# === Step 5: Fit ellipse (robust covariance) ===
fit <- cov.rob(top_bin_centers, quantile.used = floor(0.9 * nrow(top_bin_centers)))
# ellipse_pts <- ellipse::ellipse(fit$cov, centre = fit$center, level = 0.99, npoints = 200) # size of ellipse
ellipse_pts <- ellipse::ellipse(fit$cov, centre = fit$center, level = ellipse_level, npoints = 200)
ellipse_df <- as.data.frame(ellipse_pts)

# size of ellipse
scale_ellipse <- function(df, center, factor = 1.2) {
  data.frame(
    x = (df$x - center[1]) * factor + center[1],
    y = (df$y - center[2]) * factor + center[2]
  )
}
ellipse_df <- scale_ellipse(ellipse_df, fit$center, factor = scale_factor)


# === Step 6: Select points inside ellipse ===
inside <- sp::point.in.polygon(pca_df$PC1, pca_df$PC2, ellipse_df$x, ellipse_df$y) > 0
pca_df$selected <- ifelse(inside, "Inside", "Outside")

# Create output file paths
selected_file <- paste0(output_prefix, ".pca.ellipse.selected.tsv")
unselected_file <- paste0(output_prefix, ".pca.ellipse.deselected.tsv")

# Subset and write
write_tsv(pca_df[pca_df$selected == "Inside", ], selected_file)
write_tsv(pca_df[pca_df$selected == "Outside", ], unselected_file)

cat("✅ Ellipse-selected points written to:", selected_file, "\n")
cat("✅ Non-selected points written to:", unselected_file, "\n")


# === Step 7: Plot everything ===
p1_density <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = selected), alpha = 0.7, size = 0.8) +
  geom_polygon(data = ellipse_df, aes(x = x, y = y), color = "red", fill = NA, size = 0.8) +
  theme_minimal() +
  labs(title = "PCA Ellipse Around Top 10% Dense Regions", color = "Selection") +
  scale_color_manual(values = c("Inside" = "red", "Outside" = "grey70"))

p2_density <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(count))) +
  geom_density_2d(color = "white", size = 0.3) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() +
  labs(title = "UMAP Density + Contours", fill = "Count")

pdf(paste0(output_prefix, ".density_heatmap.pdf"), width = 10, height = 5)
grid.arrange(p1_density0, p2_density, ncol = 2)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

pdf(paste0(output_prefix, ".density_heatmap2.pdf"), width = 10, height = 5)
grid.arrange(p1_density, p2_density, ncol = 2)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("✅ Done. Plots and TSVs written with prefix:", output_prefix, "\n")
