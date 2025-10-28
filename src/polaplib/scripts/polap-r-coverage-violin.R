#!/usr/bin/env Rscript
# Version: v0.4.1
# Name: polap-r-coverage-violin.R
# Purpose: Plot per-base depth as violin(s) by contig + pooled "ALL".
suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: polap-r-coverage-violin.R depth.tsv out.png 'Title'")
inpath <- args[[1]]
outpng <- args[[2]]
ttl <- args[[3]]
if (!ok) stop("ggplot2 not available")

openg <- function(p) if (grepl("\\.gz$", p)) gzfile(p, "rt") else file(p, "rt")
con <- openg(inpath)
dat <- tryCatch(
  {
    read.table(con,
      header = FALSE, sep = "\t", stringsAsFactors = FALSE,
      col.names = c("contig", "pos", "depth")
    )
  },
  finally = close(con)
)
MAX_POINTS <- 5e6
if (nrow(dat) > MAX_POINTS) {
  set.seed(123)
  dat <- dat[sample.int(nrow(dat), MAX_POINTS), , drop = FALSE]
}

library(ggplot2)
dat_all <- data.frame(contig = "ALL", pos = dat$pos, depth = dat$depth)
d2 <- rbind(dat[, c("contig", "pos", "depth")], dat_all)

p <- ggplot(d2, aes(x = contig, y = depth)) +
  geom_violin(trim = TRUE, scale = "width") +
  stat_summary(fun = median, geom = "point", size = 0.8) +
  labs(title = ttl, x = "Contig", y = "Depth (×)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = outpng, plot = p, width = 12, height = 6, dpi = 180)
