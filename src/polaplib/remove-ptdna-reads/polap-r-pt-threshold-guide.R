#!/usr/bin/env Rscript
# polap-r-pt-threshold-guide.R  v0.1.0
# Plot & threshold guide for selecting ptDNA reads from a PAF
# Metrics per row:
#   ident = nmatch/alen
#   qcov  = (qend-qstart)/qlen
#   alen  = PAF col 11
#   mapq  = PAF col 12
#
# Optional labels.tsv: qname \t label   (label in {PT, nonPT})
# Outputs:
#   outdir/form_metrics.csv                  (metrics table)
#   outdir/hist_ident_qcov_alen_mapq.pdf     (histograms/densities)
#   outdir/hex_ident_qcov.pdf                (hexbin/scatter)
#   outdir/threshold_sweep.csv               (grid TPR/FPR table)
#   outdir/threshold_pick.txt                (suggested cutoffs)
#
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(hexbin)
})

opt_list <- list(
  make_option(c("-p", "--paf"),
    type = "character",
    help = "PAF file (formA.paf or similar)"
  ),
  make_option(c("-l", "--labels"),
    type = "character", default = NULL,
    help = "labels.tsv: qname\\tlabel (PT/nonPT) [optional]"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "ptguide_out",
    help = "Output directory [default: %default]"
  ),
  make_option(c("--fpr"),
    type = "double", default = 0.02,
    help = "Target max FPR for nonPT [default: %default]"
  ),
  make_option(c("--tpr"),
    type = "double", default = 0.95,
    help = "Target min TPR for PT [default: %default]"
  ),
  make_option(c("--min_alen"),
    type = "integer", default = 0,
    help = "Minimum aligned length prefilter (optional)"
  ),
  make_option(c("--min_mapq"),
    type = "integer", default = 0,
    help = "Minimum MAPQ prefilter (optional)"
  ),
  make_option(c("--prefix"),
    type = "character", default = "form",
    help = "Filename prefix [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(file.exists(opt$paf))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Read PAF and compute metrics ----
# PAF cols (0..11): qname, qlen, qstart, qend, strand,
#                   tname, tlen, tstart, tend, nmatch, alen, mapq
message("[info] reading PAF: ", opt$paf)
paf <- suppressWarnings(
  read_tsv(opt$paf,
    col_names = FALSE, comment = "#",
    col_types = cols(
      X1 = col_character(), X2 = col_integer(), X3 = col_integer(),
      X4 = col_integer(), X5 = col_character(), X6 = col_character(),
      X7 = col_integer(), X8 = col_integer(), X9 = col_integer(),
      X10 = col_integer(), X11 = col_integer(), X12 = col_integer()
    )
  )
)

if (nrow(paf) == 0) {
  stop("PAF appears empty. Aborting.")
}

df <- paf %>%
  transmute(
    qname = X1,
    qlen = X2,
    qs = X3, qe = X4,
    tname = X6,
    tl = X7,
    ts = X8, te = X9,
    nmatch = X10, alen = X11, mapq = X12,
    ident = if_else(alen > 0, nmatch / alen, 0),
    qcov = if_else(qlen > 0, (qe - qs) / qlen, 0),
    tcov = if_else(tl > 0, (te - ts) / tl, 0)
  ) %>%
  filter(alen >= opt$min_alen, mapq >= opt$min_mapq)

# If multiple alignments per read, keep the best by MAPQ for plotting
df_best <- df %>%
  arrange(desc(mapq), desc(alen)) %>%
  distinct(qname, .keep_all = TRUE)

# Optional labels
if (!is.null(opt$labels) && file.exists(opt$labels)) {
  message("[info] reading labels: ", opt$labels)
  lab <- read_tsv(opt$labels,
    col_names = FALSE,
    col_types = cols(X1 = col_character(), X2 = col_character())
  ) %>%
    rename(qname = X1, label = X2)
  df_best <- df_best %>%
    left_join(lab, by = "qname") %>%
    mutate(label = coalesce(label, "unlabeled"))
} else {
  df_best <- df_best %>% mutate(label = "unlabeled")
}

# Save metrics CSV
metrics_csv <- file.path(opt$outdir, paste0(opt$prefix, "_metrics.csv"))
write_csv(
  df_best %>%
    select(qname, qlen, ident, qcov, alen, mapq, tcov, label),
  metrics_csv
)
message("[info] wrote: ", metrics_csv)

# ---- 2) Plots: histograms/densities ----
p_hist <- function(var, lbl = "label") {
  ggplot(df_best, aes(x = .data[[var]], fill = .data[[lbl]])) +
    geom_histogram(aes(y = after_stat(density)),
      bins = 60, alpha = 0.35, position = "identity"
    ) +
    geom_density(alpha = 0.7, adjust = 1.2) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = var, y = "density", fill = "label") +
    theme_bw(base_size = 11)
}

pdf(file.path(opt$outdir, paste0(opt$prefix, "_hist_ident_qcov_alen_mapq.pdf")),
  width = 9, height = 7
)
print(p_hist("ident"))
print(p_hist("qcov"))
print(p_hist("alen"))
print(p_hist("mapq"))
dev.off()

# ---- 3) Plot: ident vs qcov (hexbin) ----
hex_plot <- function() {
  ggplot(df_best, aes(x = ident, y = qcov, color = label)) +
    stat_binhex(aes(fill = ..count..), bins = 40, alpha = 0.7, show.legend = FALSE) +
    geom_point(
      data = df_best %>% sample_n(min(5000, nrow(df_best))),
      size = 0.3, alpha = 0.2
    ) +
    scale_fill_viridis_c() +
    labs(x = "identity (nmatch/alen)", y = "query coverage ((qe-qs)/qlen)") +
    theme_bw(base_size = 11) +
    facet_wrap(~label, ncol = 2, scales = "free")
}
ggsave(
  filename = file.path(opt$outdir, paste0(opt$prefix, "_hex_ident_qcov.pdf")),
  plot = hex_plot(), width = 9, height = 6
)

# ---- 4) Threshold sweep (grid) to suggest cutoffs ----
# Only meaningful if PT/nonPT labels exist
has_labels <- any(df_best$label %in% c("PT", "nonPT"))
sweep_csv <- file.path(opt$outdir, paste0(opt$prefix, "_threshold_sweep.csv"))
pick_txt <- file.path(opt$outdir, paste0(opt$prefix, "_threshold_pick.txt"))

if (has_labels) {
  # grid over plausible ONT ranges
  ident_grid <- seq(0.88, 0.98, by = 0.005)
  qcov_grid <- seq(0.30, 0.70, by = 0.05)
  alen_grid <- c(600, 800, 1000, 1200, 1500, 2000)
  mapq_grid <- c(10, 15, 20, 25, 30)

  labd <- df_best %>%
    mutate(
      pos = (label == "PT"),
      neg = (label == "nonPT")
    )

  res <- list()
  for (i in ident_grid) {
    for (q in qcov_grid) {
      for (a in alen_grid) {
        for (m in mapq_grid) {
          pass <- (labd$ident >= i) & (labd$qcov >= q) &
            (labd$alen >= a) & (labd$mapq >= m)
          p_pos <- sum(pass & labd$pos)
          p_neg <- sum(pass & labd$neg)
          n_pos <- sum(labd$pos)
          n_neg <- sum(labd$neg)
          TPR <- ifelse(n_pos > 0, p_pos / n_pos, NA_real_)
          FPR <- ifelse(n_neg > 0, p_neg / n_neg, NA_real_)
          res[[length(res) + 1]] <- data.frame(
            ident = i, qcov = q, alen = a, mapq = m,
            TPR = TPR, FPR = FPR, TP = p_pos, FP = p_neg,
            P = n_pos, N = n_neg
          )
        }
      }
    }
  }
  sweep <- bind_rows(res)
  write_csv(sweep, sweep_csv)
  message("[info] wrote: ", sweep_csv)

  # pick minimal cutoffs that satisfy targets (TPR>=, FPR<=)
  pick <- sweep %>%
    filter(!is.na(TPR), !is.na(FPR)) %>%
    filter(TPR >= opt$tpr, FPR <= opt$fpr) %>%
    arrange(alen, mapq, qcov, ident) %>% # prefer lower gates first
    slice_head(n = 1)

  if (nrow(pick) == 0) {
    warning("No grid point meets TPR/FPR targets; writing the best trade-off.")
    pick <- sweep %>%
      filter(!is.na(TPR), !is.na(FPR)) %>%
      arrange(FPR, desc(TPR), alen, mapq, qcov, ident) %>%
      slice_head(n = 1)
  }

  # write human + shell-friendly suggestions
  cat(sprintf(
    "Suggested thresholds (targets: TPR>=%.3f, FPR<=%.3f):\n",
    opt$tpr, opt$fpr
  ), file = pick_txt)

  cat(sprintf(
    "ident>=%.3f qcov>=%.2f alen>=%d mapq>=%d  (TPR=%.3f, FPR=%.3f)\n",
    pick$ident, pick$qcov, pick$alen, pick$mapq, pick$TPR, pick$FPR
  ), file = pick_txt, append = TRUE)

  # also emit shell lines you can source
  cat(sprintf(
    "ident_min=%.3f\nqcov_min=%.2f\nalen_min=%d\nmapq_min=%d\n",
    pick$ident, pick$qcov, pick$alen, pick$mapq
  ), file = pick_txt, append = TRUE)

  message("[info] wrote: ", pick_txt)
} else {
  # no labels → still save the metrics; user can inspect plots
  writeLines("# no labels provided; sweep skipped", pick_txt)
}

message("[done] plots and suggestions are in: ", opt$outdir)
