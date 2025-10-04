#!/usr/bin/env Rscript
# r_pt_summary.R  v0.1.0
# args: all.ids pt.ids keep.ids
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<3) q(save="no", status=0)
all <- length(readLines(args[1]))
pt  <- length(readLines(args[2]))
keep<- length(readLines(args[3]))
cat(paste0("all=",all,", pt=",pt,", keep=",keep,"\n"))
