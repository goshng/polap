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
library(dbscan)
library(dplyr)

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

# === Optional: read labels ===
if (!is.na(label_file) && file.exists(label_file)) {
  labels <- read_tsv(label_file, col_names = c("id", "label"))
  labels <- labels[labels$id %in% rownames(kmers), ]
  kmers <- kmers[labels$id, , drop = FALSE]
  abund <- abund[labels$id]
  group <- labels$label
} else {
  group <- rep("Unlabeled", nrow(kmers))
}

# === PCA ===
cat("🧬 PCA...\n")
pca <- prcomp(kmers, center = TRUE, scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], label = group, abundance = abund, id = rownames(kmers))
write_tsv(pca_df, paste0(output_prefix, ".pca.tsv"))

# === UMAP ===
# n_neighbors in UMAP
# Definition: Number of nearest neighbors used to construct the local graph around each point.
#
# Default: Typically 15.
#
# Effect:
#
# Low value (e.g. 5–10): focuses on local structure — small clusters and fine detail.
#
# High value (e.g. 30–100): incorporates more global structure — broader relationships.
#
# UMAP builds a weighted graph of points based on these neighbors and tries to preserve both local and some global topology.
#
cat("🧬 UMAP...\n")
umap_xy <- umap(kmers, n_neighbors = 50)
umap_df <- data.frame(UMAP1 = umap_xy[, 1], UMAP2 = umap_xy[, 2], label = group, abundance = abund, id = rownames(kmers))

# === t-SNE ===
# perplexity in t-SNE
# Definition: Roughly the number of effective nearest neighbors. It balances attention between local and global aspects.
#
# Default: Often 30.
#
# Effect:
#
# Low perplexity (5–10): more emphasis on fine-grained local structure, but can be noisy.
#
# High perplexity (30–50): better for global cluster relationships, but may miss small details.
#
# t-SNE uses perplexity to determine the size of the Gaussian kernels when modeling distances between points.
#
cat("🧬 t-SNE...\n")
tsne <- Rtsne(kmers, dims = 2, perplexity = 50, check_duplicates = FALSE)
tsne_df <- data.frame(tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2], label = group, abundance = abund, id = rownames(kmers))
write_tsv(tsne_df, paste0(output_prefix, ".tsne.tsv"))

cat("📏 Plotting k-NN distances for DBSCAN eps estimation...\n")
kminPts <- 30 # Number of neighbors for DBSCAN
pdf(paste0(output_prefix, ".dbscan_kNNdistplot.pdf"), width = 6, height = 5)
kNNdistplot(umap_xy, k = kminPts)
abline(h = 1.0, col = "red", lty = 2) # You can adjust 1.0 to your initial eps guess
dev.off()

# === DBSCAN and plastid cluster selection ===
cat("🔍 DBSCAN clustering...\n")
db <- dbscan(umap_xy, eps = 0.7, minPts = kminPts)
umap_df$cluster <- db$cluster

# Find clusters that contain ANY plastid-labeled reads
cluster_has_plastid <- umap_df %>%
  group_by(cluster) %>%
  summarise(has_plastid = any(label == "nuclear")) %>%
  filter(has_plastid) %>%
  pull(cluster)

selected_df <- umap_df %>% filter(cluster %in% cluster_has_plastid)
write_lines(selected_df$id, paste0(output_prefix, ".selected_ids.txt"))
write_tsv(selected_df, paste0(output_prefix, ".selected.tsv"))

cat("✅ Selected", nrow(selected_df), "reads from", length(cluster_has_plastid), "cluster(s) with plastid labels.\n")

# Save UMAP output
write_tsv(umap_df, paste0(output_prefix, ".umap.tsv"))

# === Plots ===
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

# t-SNE
p3 <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = abundance)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "t-SNE", color = "kmer abundance")

p1_label <- ggplot(pca_df, aes(PC1, PC2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (by label)", color = "Label")

p2_label <- ggplot(umap_df, aes(UMAP1, UMAP2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP (by label)", color = "Label")

# t-SNE
p3_label <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = label)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "t-SNE (by label)", color = "Label")

p2_cluster <- ggplot(umap_df, aes(UMAP1, UMAP2, color = factor(cluster))) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP (by DBSCAN cluster)", color = "Cluster")

pdf(paste0(output_prefix, ".plot.pdf"), width = 10, height = 8)
grid.arrange(p1, p2, p3, ncol = 2)
dev.off()

pdf(paste0(output_prefix, ".label_plot.pdf"), width = 10, height = 8)
grid.arrange(p1_label, p2_label, p3_label, p2_cluster, ncol = 2)
dev.off()

cat("✅ Done. Plots and TSVs written with prefix:", output_prefix, "\n")
