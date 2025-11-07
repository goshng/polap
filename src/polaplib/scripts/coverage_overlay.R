# FILE: scripts/coverage_overlay.R
# VERSION: 0.1.0
# Overlay mosdepth binned coverages for all ONT vs mt-assigned ONT; shade MTPT regions.
# Usage: Rscript coverage_overlay.R allONT.regions.bed.gz mtAssigned.regions.bed.gz tracts.bed out.pdf
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: coverage_overlay.R all.regions.bed.gz mt.regions.bed.gz tracts.bed out.pdf")
}
allf <- args[1]
mtf <- args[2]
bedf <- args[3]
outf <- args[4]
all <- fread(allf, col.names = c("chr", "start", "end", "all_cov"))
mt <- fread(mtf, col.names = c("chr", "start", "end", "mt_cov"))
x <- merge(all, mt, by = c("chr", "start", "end"))
x[, mid := (start + end) / 2]
tracts <- fread(bedf, col.names = c("chr", "s", "e", "id", "score", "strand"))
pdf(outf, width = 10, height = 4)
gg <- ggplot(x, aes(x = mid)) +
  geom_line(aes(y = all_cov)) +
  geom_line(aes(y = mt_cov)) +
  theme_bw() +
  xlab("Position (bp)") +
  ylab("Coverage") +
  ggtitle("ONT coverage overlay (all vs mt-assigned)")
for (i in 1:nrow(tracts)) {
  gg <- gg + annotate("rect", xmin = tracts$s[i], xmax = tracts$e[i], ymin = -Inf, ymax = Inf, alpha = 0.1)
}
print(gg)
dev.off()
