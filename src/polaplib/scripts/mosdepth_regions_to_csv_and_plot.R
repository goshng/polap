#!/usr/bin/env Rscript
# Version: v0.2.0
suppressPackageStartupMessages({ ok <- requireNamespace("ggplot2", quietly=TRUE) })
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) stop("Usage: mosdepth_regions_to_csv_and_plot.R in.regions.bed.gz out.csv out.png 'Title'")
inbed <- args[[1]]; outcsv <- args[[2]]; outpng <- args[[3]]; ttl <- args[[4]]
con <- gzfile(inbed, "rt")
dat <- read.table(con, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                  col.names=c("chrom","start","end","depth"))
close(con)
dat$mid <- (dat$start + dat$end)/2
write.table(dat, file=outcsv, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
if (ok) {
  library(ggplot2)
  p <- ggplot(dat, aes(x=mid, y=depth)) + geom_line() +
    labs(title=ttl, x="Position (bp)", y="Mean depth per window") +
    facet_wrap(~chrom, scales="free_x") + theme_bw(base_size = 12)
  ggsave(outpng, plot=p, width=10, height=6, dpi=150)
}

