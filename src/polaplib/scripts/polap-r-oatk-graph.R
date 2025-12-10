#!/usr/bin/env Rscript
# polap-r-oatk-graph.R
# Version : v2.0.0  (2025-12-09)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Read a unified OATK metrics CSV (from polap-py-oatk-graph-table.py v0.6.0)
# and draw:
#
# 1) total_mem_gb vs cov      for patterns: h-t, t-h, t-x, x-h, x-x, x
# 2) total_time_hours vs cov  for patterns: h-t, t-h, t-x, x-h, x-x, x
# 3) completeness_prop vs cov for patterns: h-t, t-h, t-x, x-h, x-x, x
# 4) num_contigs vs cov       for patterns: h-t, t-h, t-x, x-h, x-x, x
# 5) total_length vs cov      for patterns: h-t, t-h, t-x, x-h, x-x, x
# 6) N50 vs cov               for patterns: h-t, t-h, t-x, x-h, x-x, x
#
# 7) total_mem_gb vs pattern (x-axis pattern)
# 8) total_time_hours vs pattern
# 9) completeness_prop vs pattern
# 10) num_contigs vs pattern
# 11) total_length vs pattern
# 12) N50 vs pattern
#
# 13) Histogram of avg_qual (per species)
# 14) Histogram of mean_length (per species)
#
# Input CSV columns (from v0.6.0 Python):
#   species, code2, avg_qual, mean_length,
#   tool, pattern, cov,
#   mem_gb, time_hours,
#   total_mem_gb, total_time_hours,
#   completeness_prop, num_contigs, total_length, N50,
#   mt_gfa, mt_png, pt_gfa, pt_png
#
# Only rows with tool == "oatk" are used for the pattern/cov graphs.

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

opt <- OptionParser(option_list = list(
  make_option(
    c("--data"),
    type  = "character",
    help  = "CSV with metrics (species, code2, avg_qual, mean_length, tool, pattern, cov, ...)"
  ),
  make_option(
    c("--out"),
    type  = "character",
    help  = "Output PDF file"
  ),
  make_option(
    c("--title"),
    type    = "character",
    default = "OATK metrics by coverage and pattern"
  ),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"), type = "double", default = 0.5),
  make_option(
    c("--ncol-per-row"),
    type = "integer", default = 3,
    help = "Number of panels per row [default: 3]"
  )
)) |> parse_args()

stopifnot(!is.null(opt$data), file.exists(opt$data), !is.null(opt$out))

# -----------------------------------------------------------------------------#
# Read data
# -----------------------------------------------------------------------------#
df_raw <- suppressMessages(read_csv(
  opt$data,
  show_col_types = FALSE
))

# -----------------------------------------------------------------------------#
# Convert types (robust to numeric or character input)
# -----------------------------------------------------------------------------#
to_num <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  readr::parse_number(as.character(x))
}

df <- df_raw %>%
  mutate(
    cov               = to_num(cov),
    avg_qual          = to_num(avg_qual),
    mean_length       = to_num(mean_length),
    mem_gb            = to_num(mem_gb),
    time_hours        = to_num(time_hours),
    total_mem_gb      = to_num(total_mem_gb),
    total_time_hours  = to_num(total_time_hours),
    completeness_prop = to_num(completeness_prop),
    num_contigs       = to_num(num_contigs),
    total_length      = to_num(total_length),
    N50               = to_num(N50),
    pattern           = factor(pattern),
    tool              = as.character(tool)
  )

# OATK-only subset for pattern/cov graphs
df_oatk <- df %>% filter(tool == "oatk")

# Distinct species rows for histograms (avoid over-weighting by patterns/cov)
df_species <- df %>%
  distinct(species, code2, avg_qual, mean_length)

# -----------------------------------------------------------------------------#
# Metric definitions and labels
# -----------------------------------------------------------------------------#
metrics_cov <- c(
  "total_mem_gb",
  "total_time_hours",
  "completeness_prop",
  "num_contigs",
  "total_length",
  "N50"
)

metric_cov_labels <- c(
  total_mem_gb      = "(A) total_mem_gb vs cov",
  total_time_hours  = "(B) total_time_hours vs cov",
  completeness_prop = "(C) completeness_prop vs cov",
  num_contigs       = "(D) number of contigs vs cov",
  total_length      = "(E) total_length (bp) vs cov",
  N50               = "(F) N50 (bp) vs cov"
)

metrics_pat <- metrics_cov

metric_pat_labels <- c(
  total_mem_gb      = "(G) total_mem_gb vs pattern",
  total_time_hours  = "(H) total_time_hours vs pattern",
  completeness_prop = "(I) completeness_prop vs pattern",
  num_contigs       = "(J) number of contigs vs pattern",
  total_length      = "(K) total_length (bp) vs pattern",
  N50               = "(L) N50 (bp) vs pattern"
)

# -----------------------------------------------------------------------------#
# Helper: metric vs coverage (for OATK patterns)
# -----------------------------------------------------------------------------#
make_vs_cov_plot <- function(metric_name) {
  if (!metric_name %in% names(df_oatk)) {
    return(NULL)
  }
  y <- df_oatk[[metric_name]]
  if (all(is.na(y))) {
    return(NULL)
  }

  ggplot(df_oatk, aes(x = cov, y = .data[[metric_name]], color = pattern, group = pattern)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = sort(unique(df_oatk$cov))) +
    labs(
      x     = "Coverage (cov)",
      y     = metric_cov_labels[[metric_name]],
      title = metric_cov_labels[[metric_name]],
      color = "Pattern"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0, face = "bold"),
      axis.title.x    = element_text(size = 9),
      axis.title.y    = element_text(size = 9),
      legend.position = "bottom"
    )
}

# -----------------------------------------------------------------------------#
# Helper: metric vs pattern (pattern as x, values pooled across cov)
# -----------------------------------------------------------------------------#
make_vs_pattern_plot <- function(metric_name) {
  if (!metric_name %in% names(df_oatk)) {
    return(NULL)
  }
  y <- df_oatk[[metric_name]]
  if (all(is.na(y))) {
    return(NULL)
  }

  ggplot(df_oatk, aes(x = pattern, y = .data[[metric_name]], color = pattern)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.8) +
    labs(
      x     = "Pattern",
      y     = metric_pat_labels[[metric_name]],
      title = metric_pat_labels[[metric_name]],
      color = "Pattern"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0, face = "bold"),
      axis.title.x    = element_text(size = 9),
      axis.title.y    = element_text(size = 9),
      legend.position = "none"
    )
}

# -----------------------------------------------------------------------------#
# Histograms for avg_qual and mean_length
# -----------------------------------------------------------------------------#
make_hist_avg_qual <- function() {
  if (!"avg_qual" %in% names(df_species)) {
    return(NULL)
  }
  y <- df_species$avg_qual
  if (all(is.na(y))) {
    return(NULL)
  }

  ggplot(df_species, aes(x = avg_qual)) +
    geom_histogram(bins = 30, fill = "grey70", color = "black") +
    labs(
      x     = "Average quality (avg_qual)",
      y     = "Count (species)",
      title = "(M) Histogram of avg_qual (per species)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

make_hist_mean_length <- function() {
  if (!"mean_length" %in% names(df_species)) {
    return(NULL)
  }
  y <- df_species$mean_length
  if (all(is.na(y))) {
    return(NULL)
  }

  ggplot(df_species, aes(x = mean_length)) +
    geom_histogram(bins = 30, fill = "grey70", color = "black") +
    labs(
      x     = "Mean read length (mean_length)",
      y     = "Count (species)",
      title = "(N) Histogram of mean_length (per species)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

# -----------------------------------------------------------------------------#
# Build all plots
# -----------------------------------------------------------------------------#
plots_vs_cov <- lapply(metrics_cov, make_vs_cov_plot)
plots_vs_cov <- Filter(Negate(is.null), plots_vs_cov)

plots_vs_pat <- lapply(metrics_pat, make_vs_pattern_plot)
plots_vs_pat <- Filter(Negate(is.null), plots_vs_pat)

hist1 <- make_hist_avg_qual()
hist2 <- make_hist_mean_length()

all_plots <- c(plots_vs_cov, plots_vs_pat, list(hist1, hist2))
all_plots <- Filter(Negate(is.null), all_plots)

if (length(all_plots) == 0) {
  stop("No non-NA metrics available to plot.")
}

# -----------------------------------------------------------------------------#
# Combine plots and write PDF
# -----------------------------------------------------------------------------#
ncol_layout <- opt$`ncol-per-row`
combined <- wrap_plots(all_plots, ncol = ncol_layout) +
  plot_annotation(title = opt$title)

pdf(opt$out,
  width  = opt$`page-width-in`,
  height = opt$`page-height-in`
)
grid::grid.draw(combined)
dev.off()

cat(opt$out, "\n")
