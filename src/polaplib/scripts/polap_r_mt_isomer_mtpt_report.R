#!/usr/bin/env Rscript
# Version: v0.8.0
# Purpose: Compact PDF summary for mt isomer evidence + MTPT load
# Args:
#   1 repeats.tsv
#   2 junctions.tsv
#   3 junction_support.tsv
#   4 mtpt.tsv
#   5 out.pdf

args <- commandArgs(trailingOnly = TRUE)
rep_tsv <- args[1]
junc_tsv <- args[2]
supp_tsv <- args[3]
mtpt_tsv <- args[4]
out_pdf <- args[5]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

rep <- fread(rep_tsv) # pair_id r_name r_s r_e q_name q_s q_e orient len pid
junc <- fread(junc_tsv) # junction_id pair_id type r_name r_pos q_name q_pos flank
supp <- fread(supp_tsv) # junction_id reads_support reads_left reads_right
mtpt <- tryCatch(fread(mtpt_tsv), error = function(e) data.table())

# per-pair total bridge support (sum of all its junctions)
jsum <- merge(junc, supp, by = "junction_id", all.x = TRUE)
jsum[is.na(reads_support), reads_support := 0]
perpair <- jsum[, .(bridges = sum(reads_support)), by = pair_id]
rep2 <- merge(rep, perpair, by = "pair_id", all.x = TRUE)
rep2[is.na(bridges), bridges := 0]

# Plots
p1 <- ggplot(rep2, aes(len, pid, color = orient)) +
  geom_point(alpha = .7) +
  theme_minimal() +
  labs(title = "Repeat landscape (len vs %identity)", x = "Repeat length (bp)", y = "% identity")

p2 <- ggplot(rep2, aes(reorder(pair_id, bridges), bridges, fill = orient)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Total junction-bridging reads per repeat", x = "Repeat pair", y = "Bridging reads")

p3 <- if (nrow(mtpt)) {
  ggplot(mtpt, aes(reorder(mt_contig, length), length / 1e3, fill = class)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title = "MTPT tracts (per contig)", x = "mt contig", y = "MTPT length (kb)")
} else {
  ggplot() +
    theme_void() +
    ggtitle("No MTPT tracts detected")
}

pdf(out_pdf, width = 11, height = 8.5)
grid.arrange(p1, p2, ncol = 2)
print(p3)
dev.off()
