#!/usr/bin/env Rscript
# polap-r-oatk-graph.R
# Version : v1.2.1  (2025-12-08)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Read a POLAP OATK coverage-metrics CSV and draw metric-vs-coverage graphs
# for each of the 5 patterns (h-t, t-h, t-x, x-h, x-x).

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
    help  = "CSV with metrics (species, code2, pattern, cov, ...)"
  ),
  make_option(
    c("--out"),
    type  = "character",
    help  = "Output PDF file"
  ),
  make_option(
    c("--title"),
    type    = "character",
    default = "OATK–Tiara–HiMT metrics vs coverage"
  ),
  make_option(c("--page-width-in"),  type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"),      type = "double", default = 0.5),
  make_option(
    c("--ncol-per-row"),
    type    = "integer", default = 3,
    help    = "Number of panels per row [default: 3]"
  )
)) |> parse_args()

stopifnot(!is.null(opt$data), file.exists(opt$data), !is.null(opt$out))

# -----------------------------------------------------------------------------#
# Read data
# -----------------------------------------------------------------------------#
df <- suppressMessages(read_csv(
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

df <- df %>%
  mutate(
    cov                     = to_num(cov),
    geneset_completeness_prop = to_num(geneset_completeness_prop),
    num_contigs             = to_num(num_contigs),
    total_length            = to_num(total_length),
    N50                     = to_num(N50),
    oatk_mem_gb             = to_num(oatk_mem_gb),
    oatk_time_hours         = to_num(oatk_time_hours),
    himt_mem_gb             = to_num(himt_mem_gb),
    himt_time_hours         = to_num(himt_time_hours),
    tiara_mem_gb            = to_num(tiara_mem_gb),
    tiara_time_hours        = to_num(tiara_time_hours),
    pattern                 = factor(pattern)
  )

# -----------------------------------------------------------------------------#
# Metric definitions and labels
# -----------------------------------------------------------------------------#
metrics <- c(
  "geneset_completeness_prop",
  "num_contigs",
  "total_length",
  "N50",
  "oatk_mem_gb",
  "oatk_time_hours",
  "himt_mem_gb",
  "himt_time_hours",
  "tiara_mem_gb",
  "tiara_time_hours"
)

metric_labels <- c(
  geneset_completeness_prop = "(A) Gene-set completeness",
  num_contigs               = "(B) Number of contigs",
  total_length              = "(C) Total length (bp)",
  N50                       = "(D) N50 (bp)",
  oatk_mem_gb               = "(E) OATK peak memory (GB)",
  oatk_time_hours           = "(F) OATK runtime (hours)",
  himt_mem_gb               = "(G) HiMT peak memory (GB)",
  himt_time_hours           = "(H) HiMT runtime (hours)",
  tiara_mem_gb              = "(I) Tiara peak memory (GB)",
  tiara_time_hours          = "(J) Tiara runtime (hours)"
)

# -----------------------------------------------------------------------------#
# Helper to generate one metric plot
# -----------------------------------------------------------------------------#
make_metric_plot <- function(metric_name) {
  if (!metric_name %in% names(df)) {
    return(NULL)
  }
  y <- df[[metric_name]]
  if (all(is.na(y))) {
    return(NULL)
  }

  ggplot(df, aes(x = cov, y = .data[[metric_name]], color = pattern, group = pattern)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = sort(unique(df$cov))) +
    labs(
      x     = "Coverage (cov)",
      y     = metric_labels[[metric_name]],
      title = metric_labels[[metric_name]],
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

plots <- lapply(metrics, make_metric_plot)
plots <- Filter(Negate(is.null), plots)

if (length(plots) == 0) {
  stop("No non-NA metrics available to plot.")
}

# -----------------------------------------------------------------------------#
# Combine plots and write PDF
# -----------------------------------------------------------------------------#
ncol_layout <- opt$`ncol-per-row`
combined <- wrap_plots(plots, ncol = ncol_layout) +
  plot_annotation(title = opt$title)

pdf(opt$out,
    width  = opt$`page-width-in`,
    height = opt$`page-height-in`)
grid::grid.draw(combined)
dev.off()

cat(opt$out, "\n")

