#!/usr/bin/env Rscript
# Version: v0.1.0
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  write("Usage: plot_hist.R <lengths.tsv> <out.pdf>", stderr())
  quit(status = 2)
}
x <- read.table(args[1], sep = "\t", header = FALSE, col.names = c("id", "len"))
pdf(args[2], width = 6, height = 4)
hist(x$len, main = "Lengths", xlab = "bp")
dev.off()
