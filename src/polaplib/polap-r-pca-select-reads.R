#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(gatepoints))
suppressPackageStartupMessages(library(readr))

# Prompt user for input file
#input_file <- file.choose()
input_file <- "/Volumes/run_thorne_hdd/h3/labshare/goshng/hifi1/Cucumis_sativus_var_hardwickii-1/kmer/k6.umap.tsv"

output_file <- sub("\\.tsv$", ".selected.tsv", input_file)
output2_file <- sub("\\.tsv$", ".deselected.tsv", input_file)
plot_file <- sub("\\.tsv$", ".selected.pdf", input_file)

# Record the selected K in selected.txt
# Extract full path directory
base_dir <- dirname(input_file)
# Create the output file path
selected_file <- file.path(base_dir, "selected.txt")
# Extract the number after 'k' and before '.umap.tsv'
k <- as.integer(sub(".*k([0-9]+)\\.umap\\.tsv$", "\\1", input_file))
# Write k to the file
writeLines(as.character(k), selected_file)

# Load input file
df <- read_tsv(input_file, col_types = cols())

# Use first two columns for coordinates
coords <- as.data.frame(lapply(df[, 1:2], as.numeric))

# === Interactive selection ===
#plot(coords, pch = 20, col = "grey", main = "", xlim=c(-25,25), ylim=c(-50,50))
plot(coords, pch = 20, col = "grey", main = "")
selected_indices <- gatepoints::fhs(coords)
selected_indices <- as.integer(selected_indices)

# Extract selected and deselected
selected_points <- df[selected_indices, ]
deselected_points <- df[-selected_indices, ]

# === Save static figure as PDF (no title) ===
pdf(file = plot_file, width = 7, height = 7)
#plot(coords, pch = 20, col = "grey", main = "", xlim=c(-25,25), ylim=c(-50,50))  # no title
plot(coords, pch = 20, col = "grey", main = "")  # no title
points(coords[selected_indices, ], col = "red", pch = 20)
points(coords[-selected_indices, ], col = "blue", pch = 20)
dev.off()

# Save selected/deselected tables
write.table(selected_points, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deselected_points, output2_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Selected", nrow(selected_points), "points saved to:", output_file, "\n")
cat("✅ DeSelected", nrow(deselected_points), "points saved to:", output2_file, "\n")
cat("📄 PDF figure saved to:", plot_file, "\n")
