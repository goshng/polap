#!/usr/bin/env Rscript
# polap-r-anno-span.R
# Join Length/Depth with PT/MT merged bed-bp, compute Copy, write 3 TSVs
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: polap-r-anno-span.R OUTDIR id_len_depth.tsv pt.bp.tsv mt.bp.tsv MASTER.tsv MT_OUT.tsv PT_OUT.tsv")
}
outdir <- args[1]
lenfile <- args[2]
ptbp <- args[3]
mtbp <- args[4]
master <- args[5]
mt_out <- args[6]
pt_out <- args[7]

# Length/Depth table; Depth integer; Copy integer
ld <- read_tsv(lenfile, show_col_types = FALSE)
stopifnot(all(c("id", "Length", "Depth") %in% names(ld)))
ld <- ld %>%
  mutate(
    Depth = as.integer(round(Depth))
  )
med <- median(ld$Depth, na.rm = TRUE)
if (is.na(med) || med <= 0) med <- 1L
ld <- ld %>%
  mutate(Copy = as.integer(round(Depth / med)))

PT <- if (file.exists(ptbp) && file.info(ptbp)$size > 0) {
  readr::read_tsv(ptbp, col_names = c("id", "PT"), show_col_types = FALSE)
} else {
  tibble::tibble(id = character(), PT = integer())
}

MT <- if (file.exists(mtbp) && file.info(mtbp)$size > 0) {
  readr::read_tsv(mtbp, col_names = c("id", "MT"), show_col_types = FALSE)
} else {
  tibble::tibble(id = character(), MT = integer())
}

tab <- ld %>%
  left_join(MT, by = "id") %>%
  left_join(PT, by = "id") %>%
  mutate(
    MT = replace_na(MT, 0L),
    PT = replace_na(PT, 0L),
    Edge = if_else(str_detect(id, "^edge_[0-9]+$"),
      parse_number(id) %>% as.integer(),
      0L
    )
  ) %>%
  select(Contig = id, Length, Depth, Copy, MT, PT, Edge)

# Master table
write_tsv(tab, file = master)

# MT-first: MT>0 & MT>PT
tab %>%
  filter(MT > 0L, MT > PT) %>%
  arrange(desc(MT), desc(PT), desc(Length)) %>%
  bind_rows(tab[0, ]) %>%
  write_tsv(file = mt_out)

# PT-first: PT>0 & PT>MT
tab %>%
  filter(PT > 0L, PT > MT) %>%
  arrange(desc(PT), desc(MT), desc(Length)) %>%
  bind_rows(tab[0, ]) %>%
  write_tsv(file = pt_out)
