#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: make_uniformity_plots.R depth.tsv bins.tsv metrics.tsv outdir binsize")
}
depth.tsv <- args[1]; bins.tsv <- args[2]; metrics.tsv <- args[3]; outdir <- args[4]; binsize <- as.integer(args[5])

suppressPackageStartupMessages(library(ggplot2))

# Read
dep <- read.table(depth.tsv, header=FALSE, sep="\t",
                  col.names=c("contig","pos","depth"))
bins <- read.table(bins.tsv, header=FALSE, sep="\t",
                   col.names=c("contig","start","end","n","sum","mean","sd"))
met <- read.table(metrics.tsv, header=TRUE, sep="\t")

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# Join mean per contig
contig_mean <- met[met$contig != "ALL", c("contig","mean_depth","fold80","gini")]
bins2 <- merge(bins, contig_mean, by="contig", all.x=TRUE)
bins2$mid <- (bins2$start + bins2$end)/2

## (A) Flatness lines along contigs
pA <- ggplot(bins2, aes(x=mid, y=mean)) +
  geom_line() +
  geom_hline(aes(yintercept=mean_depth), data=contig_mean, linetype="dashed") +
  geom_hline(aes(yintercept=mean_depth*1.2), data=contig_mean, linetype="dotted") +
  geom_hline(aes(yintercept=mean_depth*0.8), data=contig_mean, linetype="dotted") +
  labs(x=sprintf("Position (bin=%dbp)", binsize),
       y="Bin mean depth",
       title="Coverage flatness per contig (mean ±20% guides)") +
  facet_wrap(~contig, scales="free_x", ncol=1)
ggsave(file.path(outdir, "A_flatness_lines.pdf"), pA, width=10, height=2+1.6*length(unique(bins2$contig)))

## Prepare normalized coverage per base for ECDF & Lorenz
dep2 <- merge(dep, contig_mean[,c("contig","mean_depth")], by="contig")
dep2$norm <- with(dep2, ifelse(mean_depth>0, depth/mean_depth, NA))

## (B) Uniformity ECDF (normalized coverage)
# ECDF is P(X <= x). Many uniformity plots show P(X >= x); you can read 1-ECDF if desired.
pB <- ggplot(dep2, aes(x=norm)) +
  stat_ecdf(geom="step") +
  coord_cartesian(xlim=c(0,3), ylim=c(0,1)) +
  labs(x="Normalized coverage (depth / contig mean)",
       y="Fraction of bases <= x",
       title="Uniformity ECDF (perfectly uniform would step at x=1)") +
  facet_wrap(~contig, scales="free_y")
ggsave(file.path(outdir, "B_uniformity_ecdf.pdf"), pB, width=8, height=2+1.6*length(unique(dep2$contig)))

## (C) Lorenz curve + diagonal (Gini)
lorenz_list <- lapply(split(dep2, dep2$contig), function(d) {
  x <- sort(d$depth)
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(NULL)
  cb <- seq_len(n)/n
  cc <- cumsum(x)/sum(x)
  data.frame(contig = d$contig[1], cum_bases = cb, cum_coverage = cc)
})
lorenz <- do.call(rbind, lorenz_list)
pC <- ggplot(lorenz, aes(x=cum_bases, y=cum_coverage)) +
  geom_line() + geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x="Cumulative fraction of bases", y="Cumulative fraction of coverage",
       title="Lorenz curves (inequality of coverage; farther from diagonal = less uniform)") +
  facet_wrap(~contig)
ggsave(file.path(outdir, "C_lorenz_gini.pdf"), pC, width=8, height=2+1.6*length(unique(lorenz$contig)))

