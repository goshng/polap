#!/usr/bin/env Rscript
# polap-r-scan-autotune-qc.R  v0.1.0
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

opt_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    help = "scan.overlapness.tsv"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "scan_qc",
    help = "output directory [scan_qc]"
  ),
  make_option(c("--target"),
    type = "double", default = 0.80,
    help = "cum-wdegree target [0.80]"
  )
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(file.exists(opt$input))
dir.create(opt$outdir, showWarnings = F)

df <- read_tsv(opt$input,
  comment = "#",
  col_names = c("read_id", "degree", "wdegree"),
  col_types = "cid"
)

stopifnot(nrow(df) > 0)

# Histograms
p1 <- ggplot(df, aes(degree)) +
  geom_histogram(bins = 80) +
  scale_x_continuous(limits = c(0, quantile(df$degree, 0.99))) +
  ggtitle("Degree (scan)") +
  theme_bw(base_size = 11)

p2 <- ggplot(df, aes(wdegree)) +
  geom_histogram(bins = 80) +
  scale_x_continuous(limits = c(0, quantile(df$wdegree, 0.99))) +
  ggtitle("Wdegree (scan)") +
  theme_bw(base_size = 11)

ggsave(file.path(opt$outdir, "scan_degree_hist.pdf"), p1, width = 6, height = 4)
ggsave(file.path(opt$outdir, "scan_wdegree_hist.pdf"), p2, width = 6, height = 4)

# Cum-wdegree curve
d <- df %>%
  arrange(desc(wdegree)) %>%
  mutate(
    cum = cumsum(wdegree) / sum(wdegree),
    frac = row_number() / n()
  )
# target 80%
hit <- d %>%
  filter(cum >= opt$target) %>%
  slice_head(n = 1)
scan_keep_frac <- hit$frac[1]

# knee by max dist to diagonal
d$dist <- abs(d$cum - d$frac)
knee <- d %>% slice_max(order_by = dist, n = 1)

p3 <- ggplot(d, aes(frac, cum)) +
  geom_line() +
  geom_hline(yintercept = opt$target, lty = 2, color = "red") +
  geom_vline(xintercept = scan_keep_frac, lty = 2, color = "red") +
  geom_point(data = knee, aes(frac, cum), color = "blue") +
  labs(
    x = "fraction of reads (sorted by wdegree)",
    y = "cumulative wdegree",
    title = sprintf(
      "Cum-wdegree (80%% at %.3f; knee at %.3f)",
      scan_keep_frac, knee$frac[1]
    )
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(opt$outdir, "scan_cum_wdegree.pdf"), p3, width = 6, height = 4)

# topk_per_q suggestion from p95(degree)
deg_p95 <- quantile(df$degree, 0.95, na.rm = TRUE)
topk <- round(0.75 * as.numeric(deg_p95))
if (topk < 15) topk <- 0

cat(sprintf("scan_keep_frac=%.4f\n", scan_keep_frac))
cat(sprintf("topk_per_q=%d\n", topk))
