#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-repeat-enrichment.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Lollipop of endpoint repeat enrichment per species.
# Input: TSV with columns: species, enrichment (ratio), pval

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
getA <- function(f, d = NULL, L = FALSE) {
  i <- which(args == f)
  if (!length(i)) {
    return(d)
  }
  if (L) TRUE else if (i == length(args)) d else args[i + 1]
}
inf <- getA("--enrichment")
outf <- getA("--out")
alpha <- as.numeric(getA("--alpha", "0.05"))
if (is.null(inf) || is.null(outf)) stop("Usage: --enrichment <repeat_enrichment.tsv> --out <pdf> [--alpha 0.05]")

E <- fread(inf)
need <- c("species", "enrichment", "pval")
if (!all(need %in% names(E))) stop("enrichment file must have: species,enrichment,pval")
E[, sig := pval < alpha]

p <- ggplot(E, aes(y = reorder(species, enrichment), x = enrichment)) +
  geom_segment(aes(x = 1, xend = enrichment, y = species, yend = species), linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16), name = paste0("p <", alpha)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    x = "Enrichment ratio (observed / expected)", y = NULL,
    title = "MTPT endpoint proximity to mt repeats"
  )

ggsave(outf, p, width = 8, height = max(4, 0.3 * nrow(E)), units = "in")
cat("Wrote:", outf, "\n")
