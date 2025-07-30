#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)
library(IRanges)

# ────────────────────────────────────────────────────────────────────────────────
# Command-line options
option_list <- list(
  make_option(c("-m", "--mt"), type = "character", help = "Path to mt.paf file"),
  make_option(c("-p", "--pt"), type = "character", help = "Path to pt.paf file"),
  make_option(c("-o", "--output"), type = "character", help = "Output folder"),
  make_option(c("--min-mapq"),
    type = "integer", default = 1, dest = "min_mapq",
    help = "Minimum mapping quality [default: %default]"
  ),
  make_option(c("--min-pt"),
    type = "integer", default = 0, dest = "min_pt",
    help = "Minimum PT gene count [default: %default]"
  ),
  make_option(c("--min-identity"),
    type = "double", default = 0.15, dest = "min_identity",
    help = "Minimum alignment identity C10/C11 [default: %default]"
  )
)
# Release
opt <- parse_args(OptionParser(option_list = option_list))

# Test
# args <- commandArgs(trailingOnly = TRUE)
#
# if (interactive() && length(args) == 0) {
#   args <- c("--mt", "mt.paf",
#             "--pt", "pt.paf",
#             "--output", ".",
#             "--min-mapq", "1",
#             "--min-identity", "0.5")
# }
# opt <- parse_args(OptionParser(option_list = option_list), args = args)

output2 <- file.path(opt$output, "assembly_info_organelle_annotation_count-all.txt")
output1 <- file.path(opt$output, "assembly_info_organelle_annotation_count.txt")
output4 <- file.path(opt$output, "contig-annotation-depth-table.txt")
output5 <- file.path(opt$output, "pt-contig-annotation-depth-table.txt")

# ────────────────────────────────────────────────────────────────────────────────
# Read and filter PAF files (truncate to 12 columns)
read_and_filter_paf <- function(file, min_mapq, min_identity) {
  lines <- readLines(file)
  split_lines <- strsplit(lines, "\t")
  trimmed_lines <- lapply(split_lines, function(x) if (length(x) >= 12) x[1:12] else NULL)
  trimmed_lines <- Filter(Negate(is.null), trimmed_lines)

  df <- as_tibble(do.call(rbind, trimmed_lines), .name_repair = "minimal")
  colnames(df) <- c(
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend",
    "nmatch", "alen", "mapq"
  )

  df |>
    mutate(across(c(qlen, qstart, qend, nmatch, alen, mapq), as.numeric)) |>
    filter(!if_any(c(qlen, qstart, qend, nmatch, alen, mapq), is.na)) |>
    filter(mapq >= min_mapq, nmatch / alen >= min_identity)
}

# ────────────────────────────────────────────────────────────────────────────────
# Count merged (non-overlapping) intervals per contig
count_merged_regions <- function(df) {
  if (nrow(df) == 0) {
    return(tibble(Contig = character(), Length = numeric(), Count = integer()))
  }
  df |>
    group_by(qname, qlen) |>
    summarize(
      Contig = dplyr::first(qname),
      Length = dplyr::first(qlen),
      Count = {
        ir <- IRanges(start = qstart + 1, end = qend)
        length(reduce(ir))
      },
      .groups = "drop"
    )
}

# ────────────────────────────────────────────────────────────────────────────────
# Process mt and pt data
mt_data <- read_and_filter_paf(opt$mt, opt$min_mapq, opt$min_identity) |>
  count_merged_regions() |>
  dplyr::rename(MT = Count, Length_mt = Length)

pt_data <- read_and_filter_paf(opt$pt, opt$min_mapq, opt$min_identity) |>
  count_merged_regions() |>
  dplyr::rename(PT = Count, Length_pt = Length)

# ────────────────────────────────────────────────────────────────────────────────
# Combine and finalize
all_contigs <- union(mt_data$Contig, pt_data$Contig)

df <- tibble(Contig = all_contigs) |>
  left_join(mt_data, by = "Contig") |>
  left_join(pt_data, by = "Contig") |>
  mutate(
    Length = coalesce(Length_mt, Length_pt),
    MT = replace_na(MT, 0),
    PT = replace_na(PT, 0),
    Depth = 1,
    Copy = 1,
    Edge = str_extract(Contig, "(?<=\\.)\\d+$")
  ) |>
  select(Contig, Length, Depth, Copy, MT, PT, Edge) |>
  filter(PT >= opt$min_pt)

# ────────────────────────────────────────────────────────────────────────────────
# Write output
df |>
  arrange(desc(MT)) |>
  write.table(file = output2, sep = "\t", row.names = FALSE, quote = FALSE)

df |>
  arrange(desc(MT)) |>
  arrange(MT <= PT) |>
  write.table(file = output1, sep = "\t", row.names = FALSE, quote = FALSE)

df |>
  arrange(desc(PT)) |>
  arrange(MT <= PT) |>
  filter(MT > PT) |>
  write.table(file = output4, sep = "\t", row.names = FALSE, quote = FALSE)

df |>
  arrange(desc(MT)) |>
  arrange(PT <= MT) |>
  filter(PT > MT) |>
  write.table(file = output5, sep = "\t", row.names = FALSE, quote = FALSE)
