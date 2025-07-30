#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(gatepoints))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))  # for p1 + p2 layout

# Input file prompt or hardcoded path
#input_file <- file.choose()  # Or hardcode
input_file <- "/Volumes/run_thorne_hdd/h3/labshare/goshng/hifi1/Dunaliella_tertiolecta-1/kmer/k6.umap.tsv"

input_type <- if (grepl("umap\\.tsv$", input_file)) "umap" else if (grepl("pca\\.tsv$", input_file)) "pca" else stop("Unknown input type")
alt_type <- if (input_type == "umap") "pca" else "umap"

output_file <- sub("\\.tsv$", ".selected.tsv", input_file)
output2_file <- sub("\\.tsv$", ".deselected.tsv", input_file)
plot_file <- sub("\\.tsv$", ".selected.pdf", input_file)

# Record selected k
base_dir <- dirname(input_file)
selected_file <- file.path(base_dir, "selected.txt")
k <- as.integer(sub(".*k([0-9]+)\\..*\\.tsv$", "\\1", input_file))
writeLines(as.character(k), selected_file)

# === Load input file and interactively select ===
df <- read_tsv(input_file, col_types = cols())
coords <- as.data.frame(lapply(df[, 1:2], as.numeric))
colnames(coords) <- c("X", "Y")
plot(coords, pch = 20, col = "grey", main = "")
selected_indices <- gatepoints::fhs(coords)
selected_indices <- as.integer(selected_indices)

selected_points <- df[selected_indices, ]
deselected_points <- df[-selected_indices, ]

# Save selected/deselected
write.table(selected_points, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deselected_points, output2_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Selected", nrow(selected_points), "points saved to:", output_file, "\n")
cat("✅ DeSelected", nrow(deselected_points), "points saved to:", output2_file, "\n")

# === Overlay plot on original projection (selection plot)
df$selected <- ifelse(df$id %in% deselected_points$id, "Deselected", "Other")
colnames(df)[1:2] <- c("X", "Y")

p1 <- ggplot(df, aes(X, Y, color = selected)) +
  geom_point(alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("Deselected" = "red", "Other" = "grey70")) +
  labs(title = paste("Deselected on", toupper(input_type), "(selection space)"), color = "") +
  theme_minimal()

ggsave(sub("\\.tsv$", paste0(".", input_type, ".deselected_selection_plot.pdf"), input_file), p1, width = 7, height = 7)

# === Overlay plot on alternate projection
alt_file <- sub(input_type, alt_type, input_file)
if (!file.exists(alt_file)) {
  stop("⚠️ Corresponding file not found: ", alt_file)
}
alt_df <- read_tsv(alt_file, col_types = cols())
alt_df$selected <- ifelse(alt_df$id %in% deselected_points$id, "Deselected", "Other")
colnames(alt_df)[1:2] <- c("X", "Y")

p2 <- ggplot(alt_df, aes(X, Y, color = selected)) +
  geom_point(alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("Deselected" = "red", "Other" = "grey70")) +
  labs(title = paste("Deselected on", toupper(alt_type), "(alternate space)"), color = "") +
  theme_minimal()

ggsave(sub("\\.tsv$", paste0(".", alt_type, ".deselected_plot.pdf"), input_file), p2, width = 7, height = 7)

cat("📄 Deselected selection plot saved to:", sub("\\.tsv$", paste0(".", input_type, ".deselected_selection_plot.pdf"), input_file), "\n")
cat("📄 Deselected overlay plot saved to:", sub("\\.tsv$", paste0(".", alt_type, ".deselected_plot.pdf"), input_file), "\n")

p3 <- p1 + p2 + plot_layout(ncol = 2)
ggsave(sub("\\.tsv$", ".deselected_combined_plot.pdf", input_file), p3, width = 14, height = 7)

cat("📄 Combined plot saved to:", sub("\\.tsv$", ".deselected_combined_plot.pdf", input_file), "\n")
