#!/usr/bin/env Rscript
# Version: v0.4.0
# Mask QC: summarize total masked bp / allowed bp and plot a simple stacked bar.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
opt <- list(fasta=NULL, bed=NULL, outpdf="mask_qc.pdf", outtsv="mask_qc.tsv")
i <- 1
while (i <= length(args)) {
  if (args[i]=="--fasta") { opt$fasta <- args[i+1]; i <- i+2
  } else if (args[i]=="--bed") { opt$bed <- args[i+1]; i <- i+2
  } else if (args[i]=="--outpdf") { opt$outpdf <- args[i+1]; i <- i+2
  } else if (args[i]=="--outtsv") { opt$outtsv <- args[i+1]; i <- i+2
  } else { stop("Unknown arg: ", args[i]) }
}
if (is.null(opt$fasta) || is.null(opt$bed)) stop("Required: --fasta, --bed")

# contig lengths from .fai or compute ad hoc
get_lengths <- function(fa) {
  fai <- paste0(fa, ".fai")
  if (file.exists(fai)) {
    dt <- fread(fai, header=FALSE)
    setnames(dt, c("contig","len","off","line_bases","line_bytes"))
    return(dt[, .(contig, len)])
  } else {
    # naive length compute
    con <- file(opt$fasta, "r"); contig <- NULL; L <- 0L
    lens <- list()
    while (length(ln <- readLines(con, n=1, warn=FALSE))>0) {
      if (substr(ln,1,1)==">") {
        if (!is.null(contig)) lens[[length(lens)+1]] <- data.table(contig=contig, len=L)
        contig <- strsplit(sub("^>","",ln), "\\s+")[[1]][1]; L <- 0L
      } else L <- L + nchar(gsub("\\s+","",ln))
    }
    close(con)
    if (!is.null(contig)) lens[[length(lens)+1]] <- data.table(contig=contig, len=L)
    return(rbindlist(lens))
  }
}

lens <- get_lengths(opt$fasta)
mask <- if (file.exists(opt$bed) && file.info(opt$bed)$size>0) fread(opt$bed, header=FALSE) else data.table()
if (nrow(mask)>0) setnames(mask, c("contig","start","end"))

# merge per contig
merge_intervals <- function(dt) {
  if (nrow(dt)==0) return(dt)
  setorder(dt, start, end)
  out <- data.table(start=integer(), end=integer())
  s <- dt$start[1]; e <- dt$end[1]
  for (i in 2:nrow(dt)) {
    if (dt$start[i] <= e) {
      e <- max(e, dt$end[i])
    } else {
      out <- rbind(out, data.table(start=s, end=e))
      s <- dt$start[i]; e <- dt$end[i]
    }
  }
  out <- rbind(out, data.table(start=s, end=e))
  out
}

maskm <- mask[, merge_intervals(.SD), by=contig]
maskbp <- maskm[, .(masked_bp = sum(pmax(0L, end-start))), by=contig]
dt <- merge(lens, maskbp, by="contig", all.x=TRUE)
dt[is.na(masked_bp), masked_bp := 0L]
dt[, allowed_bp := pmax(0L, len - masked_bp)]
fwrite(dt, file=opt$outtsv, sep="\t")

# single stacked bar
sumdt <- dt[, .(masked=sum(masked_bp), allowed=sum(allowed_bp))]
plotdt <- melt(sumdt, measure.vars=c("allowed","masked"), variable.name="class", value.name="bp")
plotdt[, class := factor(class, levels=c("allowed","masked"))]

pdf(opt$outpdf, width=5, height=4)
ggplot(plotdt, aes(x="assembly", y=bp, fill=class)) +
  geom_col(width=0.6) +
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(labels=function(x) sprintf("%.1f Mb", x/1e6)) +
  labs(x=NULL, y="bp", fill=NULL, title="Mask vs Allowed (bp)") +
  theme_minimal() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

cat("[mask_qc] Wrote:", opt$outtsv, "and", opt$outpdf, "\n")

