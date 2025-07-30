#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
output_prefix <- args[2]
k <- as.integer(args[3])
label_file <- ifelse(length(args) >= 4, args[4], NA)

library(Biostrings)
library(ggplot2)
library(gridExtra)
library(Rtsne)
library(uwot)
library(readr)
library(R.utils)

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

# === t-SNE
# cat("🧬 t-SNE...\n")
# tsne <- Rtsne(kmers, dims = 2, perplexity = 30, check_duplicates = FALSE)
# tsne_df <- data.frame(tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2], label = group, abundance = abund, id = rownames(kmers))
# write_tsv(tsne_df, paste0(output_prefix, ".tsne.tsv"))

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

# p3 <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = abundance)) +
#   geom_point(alpha = 0.8) +
#   theme_minimal() +
#   scale_color_viridis_c() +
#   labs(title = "t-SNE", color = "k-mer abundance")

p1_label <- ggplot(pca_df, aes(PC1, PC2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (by label)", color = "Label")

p2_label <- ggplot(umap_df, aes(UMAP1, UMAP2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP (by label)", color = "Label")

# p3_label <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = label)) +
#   geom_point(alpha = 0.8) +
#   theme_minimal() +
#   labs(title = "t-SNE (by label)", color = "Label")

pdf(paste0(output_prefix, ".plot.pdf"), width = 10, height = 8)
grid.arrange(p1, p2, ncol = 2)
dev.off()

pdf(paste0(output_prefix, ".label_plot.pdf"), width = 10, height = 8)
grid.arrange(p1_label, p2_label, ncol = 2)
dev.off()


# === Heatmap Plots with Contours ===
p1_density <- ggplot(pca_df, aes(PC1, PC2)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(count))) +
  geom_density_2d(color = "white", size = 0.3) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() +
  labs(title = "PCA Density + Contours", fill = "Count")

p2_density <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = 100, aes(fill = after_stat(count))) +
  geom_density_2d(color = "white", size = 0.3) +
  scale_fill_viridis_c(trans = "log10") +
  theme_minimal() +
  labs(title = "UMAP Density + Contours", fill = "Count")

pdf(paste0(output_prefix, ".density_heatmap.pdf"), width = 10, height = 5)
grid.arrange(p1_density, p2_density, ncol = 2)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("✅ Done. Plots and TSVs written with prefix:", output_prefix, "\n")
