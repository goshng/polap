#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-recency-erosion.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Scatter L vs d=1-PID with 50% and 90% quantile regression per clade.
# Input TSV must have: species, CID, mtpt_id, len_bp, pid (0..1 or %), [clade]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(quantreg)
})

args <- commandArgs(trailingOnly = TRUE)
getA <- function(f, d = NULL, L = FALSE) {
  i <- which(args == f)
  if (!length(i)) {
    return(d)
  }
  if (L) TRUE else if (i == length(args)) d else args[i + 1]
}
inf <- getA("--metrics")
outf <- getA("--out")
facet <- getA("--facet-by", "clade")
if (is.null(inf) || is.null(outf)) stop("Usage: --metrics <mtpt_metrics.tsv> --out <pdf> [--facet-by clade|species]")

D <- fread(inf)
need <- c("species", "CID", "mtpt_id", "len_bp", "pid")
if (!all(need %in% names(D))) stop("metrics file must contain: ", paste(need, collapse = ", "))
D[, pid := ifelse(pid > 1, pid / 100, pid)]
D[, d := 1 - pid] # divergence
D <- D[is.finite(len_bp) & is.finite(d)]
if (!facet %in% names(D)) D[, (facet) := "all"]

# Fit quantile regression per facet
levs <- unique(D[[facet]])
curves <- list()
for (f in levs) {
  sub <- D[get(facet) == f]
  if (nrow(sub) < 10) next
  fit50 <- rq(len_bp ~ d, data = sub, tau = 0.5)
  fit90 <- rq(len_bp ~ d, data = sub, tau = 0.9)
  xs <- seq(min(sub$d), max(sub$d), length.out = 100)
  curves[[length(curves) + 1]] <- data.table(!!facet := f,
    d = xs,
    q50 = predict(fit50, data.frame(d = xs)),
    q90 = predict(fit90, data.frame(d = xs))
  )
}
C <- rbindlist(curves, fill = TRUE)

p <- ggplot(D, aes(x = d, y = len_bp)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_line(data = C, aes(x = d, y = q50), linetype = "solid") +
  geom_line(data = C, aes(x = d, y = q90), linetype = "dashed") +
  facet_wrap(reformulate(facet), scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    x = "Divergence to PT (d = 1 − PID)", y = "Fragment length (bp)",
    title = "Recency vs erosion (quantile regression)"
  )

ggsave(outf, p, width = 11, height = 8.5, units = "in")
cat("Wrote:", outf, "\n")
