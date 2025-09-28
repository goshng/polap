#!/usr/bin/env Rscript
# polap-r-pick-high-recruit.R  v0.0.2
# Input: contig.reads.tsv -> keep.ids (>= median reads)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) quit(status=2)
inp <- args[1]; outp <- args[2]
suppressWarnings(df <- try(read.table(inp, sep="\t", header=FALSE)))
if (inherits(df, "try-error") || nrow(df) == 0) {
  file.create(outp); quit(status=0)
}
colnames(df) <- c("ctg","n")
thr <- median(df$n, na.rm=TRUE)
keep <- df$ctg[df$n >= thr]
write.table(keep, file=outp, quote=FALSE,
            row.names=FALSE, col.names=FALSE)
