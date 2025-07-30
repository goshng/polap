#!/usr/bin/env Rscript

library(zoo)
library(ggplot2)
library(pracma)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
histo_file <- args[1]
output_file <- args[2]
xmin <- ifelse(length(args) >= 3, as.numeric(args[3]), 0)
xmin <- if (is.na(xmin)) 0 else xmin
xmax <- ifelse(length(args) >= 4, as.numeric(args[4]), NA)

message("length(args): ", length(args))
message("xmin: ", xmin)
message("xmax: ", xmax)

# Load histogram
histo <- read.table(histo_file, col.names = c("Multiplicity", "Count"))
histo <- histo[histo$Multiplicity > 1, ]
histo$logCount <- log10(histo$Count + 1)

# Smoothing
histo$smoothed <- rollmean(histo$Count, k = 11, fill = NA, align = "center")

# Detect peaks
peaks <- findpeaks(histo$smoothed, nups = 1, ndowns = 1, sortstr = TRUE)
peak_locs <- if (!is.null(peaks)) peaks[, 2] else integer(0)
peak_df <- histo[peak_locs, ]

# Get top 3 peaks by Count
top_peaks <- head(peak_df[order(-peak_df$Count), ], 3)
top_peaks <- top_peaks[order(top_peaks$Multiplicity), ] # order by x
top_peaks$Label <- c("Nuclear", "Mitochondrial", "Plastid")
top_peaks$Type <- "Peak"
top_peaks$Rank <- 1:3

# === Valley detection between peaks ===
valley_between <- function(start, end, smoothed, depth) {
  range_idx <- which(depth >= start & depth <= end)
  local_min_idx <- range_idx[which.min(smoothed[range_idx])]
  return(local_min_idx)
}

valley1_idx <- valley_between(top_peaks$Multiplicity[1], top_peaks$Multiplicity[2], histo$smoothed, histo$Multiplicity)
valley2_idx <- valley_between(top_peaks$Multiplicity[2], top_peaks$Multiplicity[3], histo$smoothed, histo$Multiplicity)
valleys <- histo[c(valley1_idx, valley2_idx), ]
valleys$Label <- c("Valley1", "Valley2")
valleys$Type <- "Valley"
valleys$Rank <- NA

# === Combine annotations ===
annotations <- bind_rows(top_peaks, valleys)

# === Plot ===
p <- ggplot(histo, aes(x = Multiplicity, y = Count)) +
  geom_line(aes(y = smoothed), color = "blue") +
  scale_y_log10() +
  labs(
    title = "k-mer Abundance Histogram",
    x = "k-mer multiplicity (coverage)",
    y = "Count (log10 scale)"
  ) +
  theme_minimal()

if (!is.na(xmin) & !is.na(xmax)) {
  p <- p + xlim(xmin, xmax)
}

# Add peak/valley annotations
p <- p +
  geom_vline(data = annotations, aes(xintercept = Multiplicity, linetype = Type), color = "red") +
  geom_text(
    data = annotations, aes(x = Multiplicity, y = Count, label = Label),
    vjust = -0.8, size = 3.5, color = "darkred"
  ) +
  scale_linetype_manual(values = c(Peak = "dashed", Valley = "dotted"))

ggsave(output_file, plot = p, width = 8, height = 5)

# Output table to console
annotations_out <- annotations %>%
  select(Rank, Label, Type, Multiplicity, Count, logCount)

# print(annotations_out, row.names = FALSE)

# Assign valleys to sequential ranks after the last peak
next_rank <- max(annotations_out$Rank, na.rm = TRUE)
annotations_out$Rank[is.na(annotations_out$Rank)] <- seq(next_rank + 1, by = 1, length.out = sum(is.na(annotations_out$Rank)))

# Set Rank as rownames and print
rownames(annotations_out) <- annotations_out$Rank
print(annotations_out)
