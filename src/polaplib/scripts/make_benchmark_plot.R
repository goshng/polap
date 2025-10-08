#!/usr/bin/env Rscript
###############################################################################
# scripts/make_benchmark_plots-v0.1.0.R
#
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Produce time/memory benchmark plots from md/table-<set>-0.tsv.
#
# Usage   :
#   Rscript make_benchmark_plots-v0.1.0.R --table md/table-some-0.tsv --out man/v0.5.4/figures --type time|memory
###############################################################################
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

opt <- OptionParser(option_list = list(
  make_option(c("--table"), type="character", help="input TSV"),
  make_option(c("--out"),   type="character", help="output directory"),
  make_option(c("--type"),  type="character", default="time") # time|memory
)) |> parse_args()

if (is.null(opt$table) || is.null(opt$out)) {
  stop("Missing --table or --out")
}

df <- tryCatch(read_tsv(opt$table, show_col_types=FALSE),
               error=function(e) stop(paste("Failed to read:", e)))

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

if (opt$type == "time") {
  cols <- c(Polap="_total_hours_polap_assemble_ont_pt",
            TIPPo="_total_hours_tippo_nextdenovo_hifi",
            Oatk ="_total_hours_oatk_nextdenovo_30")
  ylab <- "Wall time (hours)"
  ofn  <- file.path(opt$out, "benchmark-time.pdf")
} else {
  cols <- c(Polap="_memory_gb_polap_assemble_ont_pt",
            TIPPo="_memory_gb_tippo_nextdenovo_hifi",
            Oatk ="_memory_gb_oatk_nextdenovo_30")
  ylab <- "Peak memory (GB)"
  ofn  <- file.path(opt$out, "benchmark-memory.pdf")
}

sel <- df |>
  select(any_of(unname(cols))) |>
  rename(Polap=all_of(cols["Polap"]),
         TIPPo=all_of(cols["TIPPo"]),
         Oatk =all_of(cols["Oatk"])) |>
  pivot_longer(cols=c(Polap,TIPPo,Oatk), names_to="Tool", values_to="Value")

p <- ggplot(sel, aes(x=Tool, y=as.numeric(Value))) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.2, outlier.size=0.5) +
  labs(x=NULL, y=ylab, title=paste("Benchmark:", ifelse(opt$type=="time","Performance","Memory"))) +
  theme_bw()

ggsave(ofn, p, width=7, height=5, device=cairo_pdf)
cat(ofn, "\n")
