#!/usr/bin/env Rscript
#
# gene_set_completeness_summary.R
#
# Read PREFIX.gene_max_integrity.tsv and summarise:
#   - counts and proportions of complete / partial / fragment / missing genes
#   - mean and median of MaxGeneIntegrity
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: gene_set_completeness_summary.R <gene_max_integrity.tsv> <out_tsv>", call. = FALSE)
}

gm_file <- args[1]
out_file <- args[2]

gm <- read.table(gm_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

statuses <- c("complete", "partial", "fragment", "missing")
counts <- sapply(statuses, function(s) sum(gm$Status == s, na.rm = TRUE))
total <- sum(counts)

if (total > 0) {
  props <- counts / total
} else {
  props <- rep(NA_real_, length(counts))
}

mean_int <- mean(gm$MaxGeneIntegrity, na.rm = TRUE)
median_int <- median(gm$MaxGeneIntegrity, na.rm = TRUE)

con <- file(out_file, open = "wt")

writeLines("Status\tCount\tProportion", con)
for (i in seq_along(statuses)) {
  line <- paste(statuses[i], counts[i], props[i], sep = "\t")
  writeLines(line, con)
}

writeLines("\nMetric\tValue", con)
writeLines(paste("mean_max_gene_integrity", mean_int, sep = "\t"), con)
writeLines(paste("median_max_gene_integrity", median_int, sep = "\t"), con)

close(con)
