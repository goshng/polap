args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_merqury_qv.R <out.png> <A.qv> [B.qv...]")
out <- args[1]; qv_files <- args[-1]
labs <- sub("\\.qv$","",basename(qv_files))
qvs <- sapply(qv_files, function(f) as.numeric(readLines(f)[1]))
png(out, 800, 400, res=120)
barplot(qvs, names.arg=labs, ylab="QV (Phred)", main="Merqury QV")
dev.off()
