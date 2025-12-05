#!/usr/bin/env Rscript
# polap-r-hifi-graph.R
# Version : v1.1.4  (2025-12-04)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Boxplots of mitogenome assembly QC metrics across 4 pipelines (pmat, tippo, himt, oatk).
# Layout now configurable using --ncol-per-row.
#
# Input CSV (via --data) must have columns:
#   species, pipeline,
#   mem_gb, time_hours, disk_gb,
#   geneset_completeness_prop,
#   num_contigs, N50,
#
# Remove the following:
#   fragmentation_index,
#   contig_length_cv, max_contig_prop, total_length

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

opt <- OptionParser(option_list = list(
  make_option(c("--data"),
    type = "character",
    help = "CSV with metrics (species, pipeline, mem_gb, time_hours, disk_gb, ...)"
  ),
  make_option(c("--out"),
    type = "character",
    help = "Output PDF file"
  ),
  make_option(c("--title"),
    type = "character",
    default = "Mitogenome assembly QC metrics (PMAT, TIPPo, HiMT, Oatk)"
    # default = ""
  ),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"), type = "double", default = 0.5),

  # ⭐ NEW OPTION: number of plots per row
  make_option(c("--ncol-per-row"),
    type = "integer", default = 3,
    help = "Number of boxplots per row [default: 3]"
  )
)) |> parse_args()

stopifnot(!is.null(opt$data), file.exists(opt$data), !is.null(opt$out))

# -----------------------------------------------------------------------------#
# Read data
# -----------------------------------------------------------------------------#
df <- suppressMessages(read_csv(
  opt$data,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
))

# Ensure pipeline factor order
df <- df %>%
  mutate(
    pipeline = factor(pipeline, levels = c("pmat", "tippo", "himt", "oatk"))
  )

# -----------------------------------------------------------------------------#
# Convert metrics to numeric
# -----------------------------------------------------------------------------#
to_num <- function(x) {
  readr::parse_number(x)
}

df <- df %>%
  mutate(
    mem_gb = to_num(mem_gb),
    time_hours = to_num(time_hours),
    disk_gb = to_num(disk_gb),
    geneset_completeness_prop = to_num(geneset_completeness_prop),
    num_contigs = to_num(num_contigs),
    N50 = to_num(N50)
    # fragmentation_index = to_num(fragmentation_index),
    # contig_length_cv = to_num(contig_length_cv),
    # max_contig_prop = to_num(max_contig_prop),
    # total_length = to_num(total_length)
  )

# -----------------------------------------------------------------------------#
# Metric definitions and pretty labels
#
# 2025-12-04
# Remove the following:
# "fragmentation_index",
# "contig_length_cv",
# "max_contig_prop",
# "total_length"
# -----------------------------------------------------------------------------#
metrics <- c(
  "mem_gb",
  "time_hours",
  "disk_gb",
  "geneset_completeness_prop",
  "num_contigs",
  "N50"
)

# 2025-12-04
# Remove the following:
# fragmentation_index = "Fragmentation index",
# contig_length_cv = "Contig length CV",
# max_contig_prop = "Max contig proportion",
# total_length = "Total length (bp)"
metric_labels <- c(
  mem_gb = "(A) Peak memory (GB)",
  time_hours = "(B) Runtime (hours)",
  disk_gb = "(C) Disk (GB)",
  geneset_completeness_prop = "(D) Gene set completeness",
  num_contigs = "(E) Number of contigs",
  N50 = "(F) N50 (bp)"
)

# -----------------------------------------------------------------------------#
# Helper to generate one boxplot
# -----------------------------------------------------------------------------#
make_boxplot <- function(metric_name) {
  y <- df[[metric_name]]
  if (all(is.na(y))) {
    return(NULL)
  }

  ggplot(df, aes(x = pipeline, y = .data[[metric_name]])) +
    geom_boxplot(outlier.shape = 21) +
    scale_x_discrete(labels = c(
      pmat = "PMAT2",
      tippo = "TIPPo",
      himt = "HiMT",
      oatk = "Oatk"
    )) +
    labs(
      x = "Pipeline",
      y = metric_labels[[metric_name]],
      title = metric_labels[[metric_name]]
    ) +
    theme_bw(base_size = 10) +
    theme(
      # hjust = 0 left
      # hjust = 0.5 center
      # hjust = 1 right
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

plots <- lapply(metrics, make_boxplot)
plots <- Filter(Negate(is.null), plots)

if (length(plots) == 0) {
  stop("No non-NA metrics available to plot.")
}

# -----------------------------------------------------------------------------#
# Combine using user-configurable columns per row
# -----------------------------------------------------------------------------#
ncol_layout <- opt$`ncol-per-row`
combined <- wrap_plots(plots, ncol = ncol_layout) +
  plot_annotation(title = opt$title)

# -----------------------------------------------------------------------------#
# Output PDF
# -----------------------------------------------------------------------------#
pdf(opt$out, width = opt$`page-width-in`, height = opt$`page-height-in`)
grid::grid.draw(combined)
dev.off()

cat(opt$out, "\n")
