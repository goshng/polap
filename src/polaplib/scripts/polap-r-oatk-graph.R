#!/usr/bin/env Rscript
# polap-r-oatk-graph.R
# Version : v2.1.2  (2025-12-10)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# OATK / Tiara / HiMT / TIPPo coverage graphs.
#
# Input CSV (via --data) is the long-format table produced by
# polap-py-oatk-graph-table.py (v0.5.x+), with (at least) columns:
#
#   species, code2,
#   avg_qual, mean_length,
#   tool         ∈ { "oatk", "tippo", "himt" }
#   pattern      ∈ { "t-h", "h-t", "t-x", "x-h", "x-x", "x", ... }
#   cov          (numeric / integer, e.g. 2,5,8,11,20,30,40,50; NA for some tools)
#   mem_gb, time_hours,
#   total_mem_gb, total_time_hours, disk_gb,
#   completeness_prop, num_contigs, total_length, N50
#
# For flat OATK (oatk-x 10/20/30), and for HiMT and TIPPo, total_mem_gb and
# total_time_hours are empty, so we fall back to mem_gb and time_hours for
# those cases.

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------#
# CLI
# -----------------------------------------------------------------------------#
opt <- OptionParser(option_list = list(
  make_option(
    c("--data"),
    type  = "character",
    help  = "CSV with metrics (species, code2, tool, pattern, cov, ...)"
  ),
  make_option(
    c("--out"),
    type  = "character",
    help  = "Output PDF file"
  ),
  make_option(c("--out-tsv"),
    type = "character",
    help = "Dummy for polap-r-oatk-graph-summary.R"
  ),
  make_option(
    c("--title"),
    type    = "character",
    default = "OATK / Tiara / HiMT / TIPPo metrics"
  ),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"), type = "double", default = 0.5),
  make_option(
    c("--ncol-per-row"),
    type    = "integer",
    default = 2,
    help    = "Number of panels per row for boxplot grids [default: 2]"
  ),
  make_option(
    c("--complete-threshold-strict"),
    type    = "double",
    default = 1.0,
    help    = "Strict completeness threshold (default: 1.0)"
  ),
  make_option(
    c("--complete-threshold-relaxed"),
    type    = "double",
    default = 0.95,
    help    = "Relaxed completeness threshold (default: 0.95)"
  ),
  make_option(c("--complete-threshold"),
    type = "double",
    default = 0,
    help = "Dummy for polap-r-oatk-graph-summary.R"
  )
)) |>
  parse_args()

stopifnot(!is.null(opt$data), file.exists(opt$data), !is.null(opt$out))

# -----------------------------------------------------------------------------#
# Read data
# -----------------------------------------------------------------------------#
df <- suppressMessages(read_csv(
  opt$data,
  show_col_types = FALSE
))

required_cols <- c(
  "species", "code2", "tool", "pattern", "cov",
  "total_mem_gb", "total_time_hours", "total_disk_gb",
  "completeness_prop", "num_contigs", "total_length", "N50",
  "mem_gb", "time_hours"
)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in --data CSV: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------#
# Type conversion (robust to numeric/character)
# -----------------------------------------------------------------------------#
to_num <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  suppressWarnings(readr::parse_number(as.character(x)))
}

df <- df %>%
  mutate(
    tool              = as.character(tool),
    pattern           = as.character(pattern),
    cov               = to_num(cov),
    mem_gb            = to_num(mem_gb),
    time_hours        = to_num(time_hours),
    total_mem_gb      = to_num(total_mem_gb),
    total_time_hours  = to_num(total_time_hours),
    total_disk_gb     = to_num(total_disk_gb),
    disk_gb           = to_num(disk_gb),
    completeness_prop = to_num(completeness_prop),
    num_contigs       = to_num(num_contigs),
    total_length      = to_num(total_length),
    N50               = to_num(N50)
  )

# -----------------------------------------------------------------------------#
# For oatk-x (10/20/30), himt, and tippo, use mem_gb/time_hours if total_* are NA
# -----------------------------------------------------------------------------#
df <- df %>%
  mutate(
    is_oatk_x_flat = tool == "oatk" &
      pattern %in% c("x", "") &
      !is.na(cov) &
      as.integer(round(cov)) %in% c(10L, 20L, 30L),
    is_himt_tippo = tool %in% c("himt", "tippo")
  ) %>%
  mutate(
    total_mem_gb = ifelse(
      !is.na(total_mem_gb),
      total_mem_gb,
      ifelse(is_oatk_x_flat | is_himt_tippo, mem_gb, total_mem_gb)
    ),
    total_time_hours = ifelse(
      !is.na(total_time_hours),
      total_time_hours,
      ifelse(is_oatk_x_flat | is_himt_tippo, time_hours, total_time_hours)
    ),
    total_disk_gb = ifelse(
      !is.na(total_disk_gb),
      total_disk_gb,
      ifelse(is_oatk_x_flat | is_himt_tippo, disk_gb, total_disk_gb)
    )
  ) %>%
  select(-is_oatk_x_flat, -is_himt_tippo)

# -----------------------------------------------------------------------------#
# Define pattern / coverage / case layout
# -----------------------------------------------------------------------------#
patterns_core <- c("t-h", "h-t", "t-x", "x-h", "x-x")
patterns_main <- c("t-h", "h-t", "t-x", "x-h") # x-x handled separately
cov_all <- c(2, 5, 8, 11, 20, 30, 40, 50)
cov_xx <- c(11, 20, 30, 40, 50) # x-x only for these

# Pattern-first case labels (desired grid/bar order)
cases_pattern_first <- c(
  unlist(lapply(patterns_main, function(p) sprintf("%s-%02d", p, cov_all))),
  sprintf("x-x-%02d", cov_xx),
  "tippo",
  "himt",
  "x-30" # base OATK (pattern 'x', cov=30)
)

# Cov-first case labels
cases_cov_first <- {
  out <- character(0)
  for (cv in cov_all) {
    out <- c(out, sprintf("t-h-%02d", cv))
    out <- c(out, sprintf("h-t-%02d", cv))
    out <- c(out, sprintf("t-x-%02d", cv))
    out <- c(out, sprintf("x-h-%02d", cv))
    if (cv %in% cov_xx) {
      out <- c(out, sprintf("x-x-%02d", cv))
    }
  }
  out <- c(out, "tippo", "himt", "x-30")
  out
}

# Case levels for grid & barplots: pattern-first ordering
cases_grid_levels <- cases_pattern_first

# Aggregated group levels
pattern_group_levels <- c("t-h", "h-t", "t-x", "x-h", "x-x", "tippo", "himt", "x-30")
cov_group_levels <- c(
  "02", "05", "08", "11", "20", "30", "40", "50",
  "tippo", "himt", "x-30"
)

# -----------------------------------------------------------------------------#
# Build case_id and grouping fields
# -----------------------------------------------------------------------------#
df_cases <- df %>%
  mutate(
    cov_num = as.numeric(cov),
    case_id = case_when(
      tool == "oatk" & pattern %in% patterns_core & !is.na(cov_num) ~
        sprintf("%s-%02d", pattern, as.integer(round(cov_num))),
      tool == "oatk" & pattern %in% c("x", "") & !is.na(cov_num) &
        as.integer(round(cov_num)) == 30L ~ "x-30",
      tool == "tippo" ~ "tippo",
      tool == "himt" ~ "himt",
      TRUE ~ NA_character_
    ),
    pattern_group = case_when(
      tool == "oatk" & pattern %in% patterns_core ~ pattern,
      case_id == "x-30" ~ "x-30",
      tool == "tippo" ~ "tippo",
      tool == "himt" ~ "himt",
      TRUE ~ NA_character_
    ),
    cov_group = case_when(
      tool == "oatk" & pattern %in% patterns_core & !is.na(cov_num) ~
        sprintf("%02d", as.integer(round(cov_num))),
      case_id == "x-30" ~ "x-30",
      tool == "tippo" ~ "tippo",
      tool == "himt" ~ "himt",
      TRUE ~ NA_character_
    ),
    pattern_group = factor(pattern_group, levels = pattern_group_levels),
    cov_group = factor(cov_group, levels = cov_group_levels),
    case_pattern_first = factor(case_id, levels = cases_pattern_first),
    case_cov_first = factor(case_id, levels = cases_cov_first)
  ) %>%
  filter(!is.na(case_id)) # only rows relevant to these plots

if (nrow(df_cases) == 0) {
  stop("No rows with recognized tool/pattern/cov combinations (df_cases is empty).")
}

# For coloring we use a display pattern / tool label
df_cases <- df_cases %>%
  mutate(
    pattern_display = case_when(
      pattern %in% patterns_core ~ pattern,
      case_id == "x-30" ~ "x-30",
      tool == "tippo" ~ "tippo",
      tool == "himt" ~ "himt",
      TRUE ~ pattern
    )
  )

# -----------------------------------------------------------------------------#
# Metrics and labels
# -----------------------------------------------------------------------------#
metrics_box <- c(
  "total_mem_gb",
  "total_time_hours",
  "total_disk_gb",
  "completeness_prop"
)

metric_labels <- c(
  total_mem_gb      = "Total peak memory (GB)",
  total_time_hours  = "Total runtime (hours)",
  total_disk_gb     = "Total Disk usage (GB)",
  completeness_prop = "Gene set completeness"
)

# ASCII-only subtitles
style_labels <- list(
  pattern_first = "by pattern-coverage (pattern-first: t-h-02, ..., x-x-50, tippo, himt, x-30)",
  cov_first     = "by coverage-pattern (cov-first: t-h-02, h-t-02, ..., x-x-50, tippo, himt, x-30)",
  pattern_group = "aggregated by pattern (t-h, h-t, t-x, x-h, x-x, tippo, himt, x-30)",
  cov_group     = "aggregated by coverage (02,05,08,11,20,30,40,50, tippo, himt, x-30)"
)

# -----------------------------------------------------------------------------#
# Helpers to make boxplots
# -----------------------------------------------------------------------------#
make_box_pattern_first <- function(metric_name) {
  d <- df_cases %>% filter(!is.na(.data[[metric_name]]))
  if (nrow(d) == 0) {
    return(NULL)
  }

  ggplot(
    d,
    aes(
      x    = case_pattern_first,
      y    = .data[[metric_name]],
      fill = pattern_display
    )
  ) +
    geom_boxplot(outlier.shape = 21, size = 0.2) +
    labs(
      x     = "Pattern-coverage case (pattern-first)",
      y     = metric_labels[[metric_name]],
      title = paste(metric_labels[[metric_name]], style_labels$pattern_first),
      fill  = "Pattern / tool"
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0, face = "bold", size = 9),
      legend.position = "bottom"
    )
}

make_box_cov_first <- function(metric_name) {
  d <- df_cases %>% filter(!is.na(.data[[metric_name]]))
  if (nrow(d) == 0) {
    return(NULL)
  }

  ggplot(
    d,
    aes(
      x    = case_cov_first,
      y    = .data[[metric_name]],
      fill = pattern_display
    )
  ) +
    geom_boxplot(outlier.shape = 21, size = 0.2) +
    labs(
      x     = "Pattern-coverage case (cov-first)",
      y     = metric_labels[[metric_name]],
      title = paste(metric_labels[[metric_name]], style_labels$cov_first),
      fill  = "Pattern / tool"
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0, face = "bold", size = 9),
      legend.position = "bottom"
    )
}

make_box_pattern_group <- function(metric_name) {
  d <- df_cases %>%
    filter(!is.na(.data[[metric_name]]), !is.na(pattern_group))
  if (nrow(d) == 0) {
    return(NULL)
  }

  ggplot(
    d,
    aes(
      x    = pattern_group,
      y    = .data[[metric_name]],
      fill = pattern_group
    )
  ) +
    geom_boxplot(outlier.shape = 21, size = 0.2) +
    labs(
      x     = "Pattern / tool group",
      y     = metric_labels[[metric_name]],
      title = paste(metric_labels[[metric_name]], style_labels$pattern_group),
      fill  = "Pattern / tool"
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0, face = "bold", size = 9),
      legend.position = "none"
    )
}

make_box_cov_group <- function(metric_name) {
  d <- df_cases %>%
    filter(!is.na(.data[[metric_name]]), !is.na(cov_group))
  if (nrow(d) == 0) {
    return(NULL)
  }

  ggplot(
    d,
    aes(
      x    = cov_group,
      y    = .data[[metric_name]],
      fill = cov_group
    )
  ) +
    geom_boxplot(outlier.shape = 21, size = 0.2) +
    labs(
      x     = "Coverage / tool group",
      y     = metric_labels[[metric_name]],
      title = paste(metric_labels[[metric_name]], style_labels$cov_group),
      fill  = "Coverage / tool"
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0, face = "bold", size = 9),
      legend.position = "none"
    )
}

# -----------------------------------------------------------------------------#
# Completeness grid & barplots (strict / relaxed thresholds)
# -----------------------------------------------------------------------------#
# code2_levels <- sort(unique(df_cases$code2))

# -----------------------------------------------------------------------------#
# Y-axis ordering: sort code2 by numeric suffix (e.g. Au01, Al02, Ac03, ...)
# -----------------------------------------------------------------------------#
code2_levels <- df_cases %>%
  distinct(code2) %>%
  mutate(
    code2_chr = as.character(code2),
    code_num  = suppressWarnings(as.integer(sub("^.*?(\\d+)$", "\\1", code2_chr)))
  ) %>%
  arrange(code_num, code2_chr) %>% # numeric primary, alpha secondary
  pull(code2_chr)


# One row per dataset/case with completeness
df_comp <- df_cases %>%
  select(code2, case = case_id, completeness_prop) %>%
  distinct()

# Full grid: all datasets × all cases (pattern-first ordering)
grid_base <- expand_grid(
  code2 = factor(code2_levels, levels = code2_levels),
  case  = factor(cases_grid_levels, levels = cases_grid_levels)
)

grid_df <- grid_base %>%
  left_join(df_comp, by = c("code2", "case")) %>%
  mutate(
    case = factor(as.character(case), levels = cases_grid_levels),
    complete_strict = !is.na(completeness_prop) &
      completeness_prop >= opt$`complete-threshold-strict` - 1e-9,
    complete_relaxed = !is.na(completeness_prop) &
      completeness_prop >= opt$`complete-threshold-relaxed` - 1e-9
  )

# make_grid_plot <- function(flag_col, title_suffix) {
#   d_yes <- grid_df %>% filter(.data[[flag_col]])
#   ggplot(d_yes, aes(x = case, y = code2)) +
#     geom_point(size = 0.8) +
#     scale_x_discrete(drop = FALSE) +
#     scale_y_discrete(drop = FALSE) +
#     labs(
#       x = "Pattern-coverage / tool case",
#       y = "Dataset (code2)",
#       title = paste("Completeness", title_suffix)
#     ) +
#     theme_bw(base_size = 7) +
#     theme(
#       axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       panel.grid   = element_blank(),
#       plot.title   = element_text(hjust = 0, face = "bold", size = 9)
#     )
# }

# working version
# make_grid_plot <- function(flag_col, title_suffix) {
#   d_yes <- grid_df %>% filter(.data[[flag_col]])
#
#   # numeric y positions for horizontal guide lines
#   y_levels <- levels(grid_df$code2)
#   y_pos <- seq_along(y_levels) + 0.5
#
#   ggplot(d_yes, aes(x = case, y = code2)) +
#     # thin horizontal separators
#     geom_hline(
#       yintercept = y_pos,
#       linewidth = 0.25,
#       colour = "grey80"
#     ) +
#     geom_point(size = 0.8) +
#     scale_x_discrete(drop = FALSE) +
#     scale_y_discrete(
#       drop = FALSE,
#       limits = rev(code2_levels) # bottom → top ordering
#     ) +
#     labs(
#       x = "Pattern-coverage / tool case",
#       y = "Dataset (code2)",
#       title = paste("Completeness", title_suffix)
#     ) +
#     theme_bw(base_size = 7) +
#     theme(
#       axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       axis.text.y  = element_text(size = 7),
#       panel.grid   = element_blank(), # keep only our custom lines
#       plot.title   = element_text(hjust = 0, face = "bold", size = 9)
#     )
# }

# version: o x
make_grid_plot <- function(flag_col, title_suffix) {
  # complete cases (●)
  d_yes <- grid_df %>% filter(.data[[flag_col]])

  # incomplete-but-present cases (×)
  d_no <- grid_df %>%
    filter(
      !.data[[flag_col]] &
        !is.na(completeness_prop) # only where data exists
    )

  # numeric y positions for horizontal guide lines
  y_levels <- levels(grid_df$code2)
  y_pos <- seq_along(y_levels) + 0.5

  ggplot() +

    # horizontal separators
    geom_hline(
      yintercept = y_pos,
      linewidth = 0.25,
      colour = "grey80"
    ) +

    # ● complete
    geom_point(
      data = d_yes,
      aes(x = case, y = code2),
      shape = 16,
      size = 0.8
    ) +

    # × incomplete
    geom_point(
      data = d_no,
      aes(x = case, y = code2),
      shape = 4, # X mark
      size = 0.9,
      stroke = 0.6
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(
      drop = FALSE,
      limits = rev(code2_levels)
    ) +
    labs(
      x = "Pattern-coverage / tool case",
      y = "Dataset (code2)",
      title = paste("Completeness", title_suffix)
    ) +
    theme_bw(base_size = 7) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y  = element_text(size = 7),
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0, face = "bold", size = 9)
    )
}



make_bar_count_plot <- function(flag_col, title_suffix) {
  d_counts <- grid_df %>%
    group_by(case) %>%
    summarise(
      n_complete = sum(.data[[flag_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      case = factor(as.character(case), levels = cases_grid_levels)
    )

  ggplot(d_counts, aes(x = case, y = n_complete)) +
    geom_col() +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Pattern-coverage / tool case",
      y = "Number of complete datasets",
      title = paste("Count of complete datasets", title_suffix)
    ) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title   = element_text(hjust = 0, face = "bold", size = 9)
    )
}

grid_strict_plot <- make_grid_plot(
  "complete_strict",
  sprintf("== %.2f", opt$`complete-threshold-strict`)
)
grid_relaxed_plot <- make_grid_plot(
  "complete_relaxed",
  sprintf(">= %.2f", opt$`complete-threshold-relaxed`)
)
bar_strict_plot <- make_bar_count_plot(
  "complete_strict",
  sprintf("for completeness == %.2f", opt$`complete-threshold-strict`)
)
bar_relaxed_plot <- make_bar_count_plot(
  "complete_relaxed",
  sprintf("for completeness >= %.2f", opt$`complete-threshold-relaxed`)
)

# -----------------------------------------------------------------------------#
# Output PDF
# -----------------------------------------------------------------------------#
pdf(
  opt$out,
  width  = opt$`page-width-in`,
  height = opt$`page-height-in`
)

# 1–4) Boxplots (four styles) for each metric
for (metric_name in metrics_box) {
  p1 <- make_box_pattern_first(metric_name)
  p2 <- make_box_cov_first(metric_name)
  p3 <- make_box_pattern_group(metric_name)
  p4 <- make_box_cov_group(metric_name)

  plots_metric <- list(p1, p2, p3, p4)
  plots_metric <- Filter(Negate(is.null), plots_metric)

  if (length(plots_metric) == 0) {
    next
  }

  combined <- wrap_plots(plots_metric, ncol = opt$`ncol-per-row`) +
    plot_annotation(
      title = paste0("Metric: ", metric_labels[[metric_name]])
    ) &
    theme(legend.position = "bottom")

  print(combined)
}

# 5–8) Completeness grids and barplots
print(grid_strict_plot)
print(grid_relaxed_plot)
print(bar_strict_plot)
print(bar_relaxed_plot)

dev.off()

cat(opt$out, "\n")
