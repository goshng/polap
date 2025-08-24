#!/usr/bin/env Rscript

library(zoo)

args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 2) {
#   cat("Usage: plot_kmer_histogram_peaks.R <reads.histo> <output.pdf>\n")
#   quit("no")
# }

histo_file <- args[1]
output_file <- args[2]
xmin <- ifelse(length(args) >= 3, as.numeric(args[3]), 0)
xmin <- if (is.na(xmin)) 0 else xmin
xmax <- ifelse(length(args) >= 4, as.numeric(args[4]), NA)

message("length(args): ", length(args))
message("xmin: ", xmin)
message("xmax: ", xmax)

suppressPackageStartupMessages({
  library(ggplot2)
  library(pracma) # for findpeaks
})

# Load histogram
histo <- read.table(histo_file, col.names = c("Multiplicity", "Count"))
histo <- histo[histo$Multiplicity > 1, ] # remove error peak

# Log scale for better visibility
histo$logCount <- log10(histo$Count + 1)

# Find peaks
histo$smoothed <- rollmean(histo$Count, k = 11, fill = NA, align = "center")
peaks <- findpeaks(histo$smoothed, nups = 1, ndowns = 1, sortstr = TRUE)

# counts_vec <- histo$Count
# peaks <- findpeaks(counts_vec, nups = 1, ndowns = 1, sortstr = TRUE)


# Get peak positions and annotate
peak_locs <- if (!is.null(peaks)) peaks[, 2] else integer(0)
peak_df <- histo[peak_locs, ]
top_peaks <- head(peak_df[order(-peak_df$Count), ], 3)

# Label and rank them
top_peaks$Label <- c("Nuclear", "Mitochondrial", "Plastid")
top_peaks$Rank <- 1:nrow(top_peaks)
top_peaks <- top_peaks[, c("Rank", "Multiplicity", "Count", "logCount", "Label")]

# Plot
p <- ggplot(histo, aes(x = Multiplicity, y = Count)) +
  geom_line(aes(y = smoothed), color = "blue") +
  # geom_line(color = "blue") +
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

# Annotate peaks
if (nrow(top_peaks) > 0) {
  p <- p + geom_vline(data = top_peaks, aes(xintercept = Multiplicity), linetype = "dashed", color = "red") +
    geom_text(
      data = top_peaks, aes(x = Multiplicity, y = Count, label = paste0("Peak @ ", Multiplicity)),
      vjust = -0.5, hjust = 0.5, size = 3.5, color = "darkred"
    )
}

ggsave(output_file, plot = p, width = 8, height = 5)

# cat("✅ Plot saved to:", output_file, "\n")
if (nrow(top_peaks) > 0) {
  # cat("📈 Top detected peaks:\n")
  # write.table(top_peaks, file = "peak_thresholds.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  rownames(top_peaks) <- NULL
  print(top_peaks)
}
