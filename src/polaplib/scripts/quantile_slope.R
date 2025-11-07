# FILE: scripts/quantile_slope.R
# VERSION: 0.1.0
# Collects MTPTs across species; L = (end-start); d = 1 - mean_pid; fits rq(L ~ d) by clade.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(quantreg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: quantile_slope.R <mtpt_calls_dir> <outdir> <clade_map.tsv|''> <taus_csv> <boots>")
}
indir <- args[1]
outdir <- args[2]
clmap <- args[3]
taus <- as.numeric(strsplit(args[4], ",")[[1]])
B <- as.integer(args[5])

# Load all mtpt.tsv
files <- list.files(indir, pattern = "^.*?/mtpt\\.tsv$", full.names = TRUE, recursive = TRUE)
DT <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  sp <- basename(dirname(f))
  x[, species := sp]
  x[, len := end - start]
  x[, d := 1 - as.numeric(mean_pid)]
  x
}), fill = TRUE)

if (nrow(DT) == 0) stop("No mtpt.tsv found under: ", indir)

# Clade mapping (optional)
if (nzchar(clmap) && file.exists(clmap)) {
  M <- fread(clmap) # columns: species clade
  DT <- merge(DT, M, by = "species", all.x = TRUE)
  DT[is.na(clade), clade := species]
} else {
  DT[, clade := species]
}

# Fit rq per clade per tau
res <- list()
for (cl in sort(unique(DT$clade))) {
  dsub <- DT[clade == cl & is.finite(len) & is.finite(d)]
  if (nrow(dsub) < 10) next
  for (tau in taus) {
    fit <- rq(len ~ d, tau = tau, data = dsub)
    # bootstrap CI
    s <- summary(fit, se = "boot", R = B)
    slope <- coef(fit)[["d"]]
    lo <- s$coefficients["d", "Lower bd"]
    hi <- s$coefficients["d", "Upper bd"]
    res[[length(res) + 1]] <- data.table(clade = cl, tau = tau, slope = slope, lo = lo, hi = hi, n = nrow(dsub))
  }
}
RES <- rbindlist(res)
fwrite(RES, file = file.path(outdir, "erosion_slopes.tsv"), sep = "\t")

# Plot
pdf(file.path(outdir, "erosion_plot.pdf"), width = 10, height = 7)
ggplot(RES, aes(x = factor(tau), y = slope)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap(~clade, scales = "free_y") +
  theme_bw() +
  xlab("Quantile (tau)") +
  ylab("Slope of L ~ (1-PID)") +
  ggtitle("MTPT erosion slopes by clade")
dev.off()
cat("Erosion slopes -> ", file.path(outdir, "erosion_slopes.tsv"), "\n")
