#!/usr/bin/env Rscript
# scripts/paf_plot_core.R
# Version: v0.1.0
suppressPackageStartupMessages({
  ok_dt <- requireNamespace("data.table", quietly = TRUE)
  ok_gp <- requireNamespace("ggplot2", quietly = TRUE)
  ok_hb <- requireNamespace("hexbin", quietly = TRUE)
})
if (!ok_dt || !ok_gp || !ok_hb) {
  stop("Requires R packages: data.table, ggplot2, hexbin")
}
library(data.table)
library(ggplot2)
library(hexbin)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: paf_plot_core.R in.core.csv outdir prefix")
}
incsv <- args[[1]]
outdir <- args[[2]]
prefix <- args[[3]]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

DT <- fread(incsv)
# coerce
numcols <- c("qlen", "qstart", "qend", "tlen", "tstart", "tend", "nmatch", "alen", "mapq", "identity", "q_aln_frac")
for (c in numcols) set(DT, j = c, value = as.numeric(DT[[c]]))
DT[, tp := ifelse(tp == "", "NA", tp)]
DT[, strand := factor(strand, levels = c("+", "-"))]

# 1) hist_alen
p1 <- ggplot(DT, aes(x = alen)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  labs(title = "Aligned length (alen)", x = "Aligned length (bp, log10)", y = "Count")
if (length(unique(na.omit(DT$tp))) > 0) p1 <- p1 + facet_grid(rows = vars(tp))
ggsave(file.path(outdir, paste0(prefix, "_hist_alen.png")), p1, width = 7, height = 5, dpi = 180)

# 2) hist_identity
p2 <- ggplot(DT, aes(x = identity)) +
  geom_histogram(binwidth = 0.005, boundary = 0) +
  coord_cartesian(xlim = c(0.5, 1.0)) +
  labs(title = "Per-alignment identity (nmatch/alen)", x = "Identity", y = "Count")
if (length(unique(na.omit(DT$strand))) > 0) p2 <- p2 + facet_grid(cols = vars(strand))
ggsave(file.path(outdir, paste0(prefix, "_hist_identity.png")), p2, width = 7, height = 5, dpi = 180)

# 3) hex_alen_identity
DT[, log10_alen := log10(alen)]
p3 <- ggplot(DT, aes(x = log10_alen, y = identity)) +
  stat_binhex(bins = 50) +
  labs(title = "alen vs identity (hexbin)", x = "log10(alen)", y = "Identity")
ggsave(file.path(outdir, paste0(prefix, "_hex_alen_identity.png")), p3, width = 7, height = 5, dpi = 180)

# 4) hist_mapq
p4 <- ggplot(DT, aes(x = mapq)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title = "MAPQ distribution", x = "MAPQ", y = "Count")
if (length(unique(na.omit(DT$tp))) > 0) p4 <- p4 + facet_grid(rows = vars(tp))
ggsave(file.path(outdir, paste0(prefix, "_hist_mapq.png")), p4, width = 7, height = 5, dpi = 180)

# 5) ECDF of query aligned fraction
p5 <- ggplot(DT, aes(x = q_aln_frac)) +
  stat_ecdf(geom = "step") +
  labs(title = "ECDF of query aligned fraction", x = "(qend - qstart)/qlen", y = "ECDF")
if (length(unique(na.omit(DT$tp))) > 0) p5 <- p5 + facet_grid(rows = vars(tp))
ggsave(file.path(outdir, paste0(prefix, "_ecdf_q_aligned_fraction.png")), p5, width = 7, height = 5, dpi = 180)
