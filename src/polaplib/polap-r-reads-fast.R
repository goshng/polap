#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

opt_list <- list(
  make_option(c("-m", "--mt"), type = "character", help = "Path to mt.paf"),
  make_option(c("-p", "--pt"), type = "character", help = "Path to pt.paf"),
  make_option(c("-o", "--output"), type = "character", help = "Output folder"),
  make_option(c("--min-mapq"), type = "integer", default = 1, dest = "min_mapq"),
  make_option(c("--min-pt"), type = "integer", default = 0, dest = "min_pt"),
  make_option(c("--min-identity"), type = "double", default = 0.15, dest = "min_identity")
)
opt <- parse_args(OptionParser(option_list = opt_list))

out2 <- file.path(opt$output, "assembly_info_organelle_annotation_count-all.txt")
out1 <- file.path(opt$output, "assembly_info_organelle_annotation_count.txt")
out4 <- file.path(opt$output, "contig-annotation-depth-table.txt")
out5 <- file.path(opt$output, "pt-contig-annotation-depth-table.txt")

# Read just the columns we need: 1,2,3,4,10,11,12 (1-based)
# qname, qlen, qstart, qend, nmatch, alen, mapq
read_paf_fast <- function(path, min_mapq, min_id) {
  if (is.null(path) || path == "") {
    return(data.table())
  }
  cols <- c("qname", "qlen", "qstart", "qend", "nmatch", "alen", "mapq")
  dt <- suppressWarnings(fread(
    path,
    sep = "\t", header = FALSE, select = c(1, 2, 3, 4, 10, 11, 12),
    col.names = cols, showProgress = FALSE, data.table = TRUE
  ))
  if (nrow(dt) == 0L) {
    return(dt)
  }
  dt <- dt[mapq >= min_mapq & (as.numeric(nmatch) / as.numeric(alen)) >= min_id]
  # PAF is 0-based, half-open in qstart,qend; we want 1-based closed
  dt[, `:=`(qstart = as.integer(qstart) + 1L, qend = as.integer(qend))]
  dt[, qlen := as.integer(qlen)]
  dt[]
}

# Count merged non-overlapping intervals per qname
count_merged_per_qname <- function(dt) {
  if (nrow(dt) == 0L) {
    return(data.table(Contig = character(), Length = integer(), Count = integer()))
  }
  setorder(dt, qname, qstart, qend)
  # Run-length encode groups of qname
  r <- rleid(dt$qname)
  ans <- dt[,
    {
      # merge intervals for one qname
      s <- qstart
      e <- qend
      merged <- 0L
      cur_s <- s[1L]
      cur_e <- e[1L]
      if (length(s) > 1L) {
        for (i in 2L:length(s)) {
          if (s[i] <= cur_e) {
            if (e[i] > cur_e) cur_e <- e[i]
          } else {
            merged <- merged + 1L
            cur_s <- s[i]
            cur_e <- e[i]
          }
        }
      }
      merged <- merged + 1L
      .(Contig = qname[1L], Length = qlen[1L], Count = merged)
    },
    by = r
  ]
  ans[, r := NULL][]
}

mt_dt <- read_paf_fast(opt$mt, opt$min_mapq, opt$min_identity)
pt_dt <- read_paf_fast(opt$pt, opt$min_mapq, opt$min_identity)

mt_tab <- count_merged_per_qname(mt_dt)
setnames(mt_tab, c("Contig", "Length", "MT"))
pt_tab <- count_merged_per_qname(pt_dt)
setnames(pt_tab, c("Contig", "Length", "PT"))

# Combine
allc <- merge(mt_tab, pt_tab, by = "Contig", all = TRUE, suffixes = c("_mt", "_pt"))
allc[, `:=`(
  Length = fifelse(!is.na(Length_mt), Length_mt, Length_pt),
  MT = fifelse(!is.na(MT), MT, 0L),
  PT = fifelse(!is.na(PT), PT, 0L)
)]
allc[, `:=`(Depth = 1L, Copy = 1L)]
allc[, Edge := {
  # numeric tail after last '.'
  p <- as.integer(regexpr("\\.[0-9]+$", Contig))
  ifelse(p > 0, sub("^.*\\.", "", Contig), NA_character_)
}]
res <- allc[, .(Contig, Length, Depth, Copy, MT, PT, Edge)]
res <- res[PT >= opt$min_pt]

# Output 2: sort by desc(MT)
fwrite(res[order(-MT)], out2, sep = "\t", quote = FALSE)

# Output 1: arrange(desc(MT)); then stable arrange(MT <= PT)
tmp <- res[order(-MT)]
tmp[, key2 := (MT <= PT)]
setorder(tmp, key2) # stable in data.table
tmp[, key2 := NULL]
fwrite(tmp, out1, sep = "\t", quote = FALSE)

# Output 4: MT > PT by (MT-PT) desc
fwrite(res[MT > PT][order(-(MT - PT))], out4, sep = "\t", quote = FALSE)

# Output 5: PT > MT (mimic two-step arrange)
tmp2 <- res[order(-MT)]
tmp2[, key3 := (PT <= MT)]
setorder(tmp2, key3)
tmp2 <- tmp2[PT > MT]
tmp2[, key3 := NULL]
fwrite(tmp2, out5, sep = "\t", quote = FALSE)
