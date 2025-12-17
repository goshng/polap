#!/usr/bin/env Rscript
# polap-r-oatk-graph-summary.R
# Version : v1.0.0  (2025-12-10)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# SUMMARY SCRIPT:
#
# Reads the long-format OATK/Tiara/HiMT/TIPPo CSV produced by
# polap-py-oatk-graph-table.py (v0.5.x+).
#
# Instead of drawing graphs (as in polap-r-oatk-graph.R), this script writes
# a TSV listing *which species have missing or zero completeness* for each
# pattern–coverage/tool case.
#
# Output TSV columns:
#   species   code2   case_id
#
# where case_id is one of:
#   t-h-02 ... x-x-50, tippo, himt, x-30
#
# (Exact ordering as in polap-r-oatk-graph.R grid plots.)
#

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
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
  make_option(c("--out-tsv"),
    type = "character",
    help = "Output TSV listing species with 0 or missing completeness"
  ),
  make_option(c("--complete-threshold"),
    type = "double",
    default = 0,
    help = "Report species where completeness_prop <= threshold (default: 0)"
  ),
  make_option(
    c("--out"),
    type  = "character",
    help  = "Output PDF file"
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
  )
)) |> parse_args()

stopifnot(!is.null(opt$data), file.exists(opt$data))
stopifnot(!is.null(opt$`out-tsv`))

# -----------------------------------------------------------------------------#
# Load data
# -----------------------------------------------------------------------------#
df <- suppressMessages(read_csv(opt$data, show_col_types = FALSE))

required <- c("species", "code2", "tool", "pattern", "cov", "completeness_prop")
miss <- setdiff(required, names(df))
if (length(miss) > 0) {
  stop("Missing required columns: ", paste(miss, collapse = ", "))
}

# Numeric conversion
to_num <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  suppressWarnings(readr::parse_number(as.character(x)))
}
df <- df %>%
  mutate(
    cov = to_num(cov),
    completeness_prop = to_num(completeness_prop)
  )

# -----------------------------------------------------------------------------#
# Define pattern/cov cases (same order as plotting script)
# -----------------------------------------------------------------------------#
patterns_main <- c("t-h", "h-t", "t-x", "x-h")
patterns_xx <- c("x-x")
cov_all <- c(2, 5, 8, 11, 20, 30, 40, 50)
# 2025-12-17
# cov_xx <- c(11, 20, 30, 40, 50)
cov_xx <- c(5, 8, 11, 20, 30, 40, 50)

cases_pattern_first <- c(
  unlist(lapply(patterns_main, function(p) sprintf("%s-%02d", p, cov_all))),
  sprintf("x-x-%02d", cov_xx),
  "tippo",
  "himt",
  "x-30"
)

# -----------------------------------------------------------------------------#
# Assign case_id for each row
# -----------------------------------------------------------------------------#
df2 <- df %>%
  mutate(
    case_id = case_when(
      tool == "oatk" & pattern %in% c("t-h", "h-t", "t-x", "x-h", "x-x") & !is.na(cov) ~
        sprintf("%s-%02d", pattern, as.integer(round(cov))),
      tool == "oatk" & pattern %in% c("x", "") &
        !is.na(cov) & as.integer(round(cov)) == 30L ~ "x-30",
      tool == "tippo" ~ "tippo",
      tool == "himt" ~ "himt",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(case_id))

# Ensure ordered factor
df2 <- df2 %>%
  mutate(case_id = factor(case_id, levels = cases_pattern_first))

# -----------------------------------------------------------------------------#
# Identify species whose completeness_prop <= threshold or is NA
# -----------------------------------------------------------------------------#
bad_rows <- df2 %>%
  filter(is.na(completeness_prop) | completeness_prop <= opt$`complete-threshold`) %>%
  arrange(case_id, species, code2)

# -----------------------------------------------------------------------------#
# Write TSV
# -----------------------------------------------------------------------------#
write_tsv <- function(x, file) {
  write.table(
    x,
    file = file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
}

write_tsv(bad_rows %>% select(species, code2, case_id), opt$`out-tsv`)

cat("Wrote summary TSV:", opt$`out-tsv`, "\n")
