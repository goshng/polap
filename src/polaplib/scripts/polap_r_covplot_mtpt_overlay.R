#!/usr/bin/env Rscript
# Version: v0.1.0
# Read two coverage TSVs (allONT vs assignedONT), optional MTPT BED, and plot overlay.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
opt <- list(fasta=NULL, cov_all=NULL, cov_assigned=NULL,
            mask=NULL, bin=200, out_pdf="coverage_overlay.pdf",
            out_tsv="coverage_summary.tsv", title="Coverage overlay")
i <- 1
while (i <= length(args)) {
  key <- args[i]
  if (key=="--fasta") { opt$fasta <- args[i+1]; i <- i+2
  } else if (key=="--cov-all") { opt$cov_all <- args[i+1]; i <- i+2
  } else if (key=="--cov-assigned") { opt$cov_assigned <- args[i+1]; i <- i+2
  } else if (key=="--mask") { opt$mask <- args[i+1]; i <- i+2
  } else if (key=="--bin") { opt$bin <- as.integer(args[i+1]); i <- i+2
  } else if (key=="--out-pdf") { opt$out_pdf <- args[i+1]; i <- i+2
  } else if (key=="--out-tsv") { opt$out_tsv <- args[i+1]; i <- i+2
  } else if (key=="--title") { opt$title <- args[i+1]; i <- i+2
  } else stop("Unknown arg: ", key)
}
stopifnot(!is.null(opt$fasta), !is.null(opt$cov_all), !is.null(opt$cov_assigned))

# lengths
fai <- fread(paste0(opt$fasta, ".fai"), header=FALSE)
setnames(fai, c("contig","len","off","line_bases","line_bytes"))
lens <- fai[, .(contig, len)]

# coverage
A <- fread(opt$cov_all); setnames(A, c("contig","start","end","mean_depth","source"))
B <- fread(opt$cov_assigned); setnames(B, c("contig","start","end","mean_depth","source"))
D <- rbindlist(list(A,B), use.names=TRUE)
D[, kb := (start + end)/2/1000]

# mask (optional)
M <- data.table(); has_mask <- FALSE
if (!is.null(opt$mask) && file.exists(opt$mask) && file.info(opt$mask)$size > 0) {
  M <- fread(opt$mask, header=FALSE)
  setnames(M, c("contig","mstart","mend"))
  has_mask <- TRUE
}

# summary metrics per contig/source
summ <- D[, .(
  n=.N,
  mean=mean(mean_depth),
  sd=sd(mean_depth),
  cv = ifelse(mean(mean_depth)>0, sd(mean_depth)/mean(mean_depth), NA_real_),
  median=median(mean_depth),
  max=max(mean_depth),
  max_over_median = ifelse(median(mean_depth)>0, max(mean_depth)/median(mean_depth), NA_real_)
), by=.(contig, source)]

# If mask exists: on-mask vs off-mask means
if (has_mask) {
  setkey(D, contig, start, end)
  setkey(M, contig, mstart, mend)
  # expand bin endpoints to intervals matching foverlaps expectations
  Dov <- foverlaps(D, M, by.x=c("contig","start","end"), by.y=c("contig","mstart","mend"), nomatch=0L)
  D[, onmask := FALSE]
  if (nrow(Dov)>0) {
    D[ Dov[, .(contig,start,end)], on=.(contig,start,end), onmask := TRUE]
  }
  S2 <- D[, .(mean_onmask = mean(mean_depth[onmask]),
              mean_offmask = mean(mean_depth[!onmask]),
              ratio_on_off = mean(mean_depth[onmask]) / mean(mean_depth[!onmask])),
          by=.(contig, source)]
  summ <- merge(summ, S2, by=c("contig","source"), all.x=TRUE)
}

fwrite(summ[order(contig, source)], file=opt$out_tsv, sep="\t")

# Plot
pdf(opt$out_pdf, width=11, height=8.5)
gg <- ggplot(D, aes(x=kb, y=mean_depth, color=source)) +
  geom_line(size=0.3) +
  facet_wrap(~contig, scales="free_x", ncol=1) +
  labs(x=sprintf("Position (kb, bin=%d bp)", opt$bin),
       y="Mean depth per bin",
       title=opt$title, color="") +
  theme_minimal() +
  theme(strip.text=element_text(face="bold"),
        legend.position="top")
if (has_mask) {
  # overlay mask as light rectangles (per contig facet). Use annotate per contig via geom_rect + data
  M[, y0 := -Inf][, y1 := Inf]
  gg <- gg + geom_rect(data=M, inherit.aes=FALSE,
                       aes(xmin=mstart/1000, xmax=mend/1000, ymin=y0, ymax=y1),
                       alpha=0.15)
}
print(gg)
dev.off()

cat("[covplot] Wrote:", opt$out_pdf, "and", opt$out_tsv, "\n")

