#!/usr/bin/env Rscript
# polap-r-paf-diagnose.R  v0.0.2
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: polap-r-paf-diagnose.R <mt.paf> [out_prefix] [fixed_qcov=0.50] [fixed_id=0.80] [min_alnlen=0] [primary_only=false]\n")
  quit(status = 2)
}
paf_file <- args[[1]]
out_prefix <- ifelse(length(args) >= 2, args[[2]], "paf_diag")
fixed_qcov <- ifelse(length(args) >= 3, as.numeric(args[[3]]), 0.50)
fixed_id <- ifelse(length(args) >= 4, as.numeric(args[[4]]), 0.80)
min_alnlen <- ifelse(length(args) >= 5, as.numeric(args[[5]]), 0)
primary_only <- ifelse(length(args) >= 6, tolower(args[[6]]) %in% c("1", "t", "true", "yes"), FALSE)

cols <- cols(
  X1 = col_character(), X2 = col_double(), X3 = col_double(), X4 = col_double(),
  X5 = col_character(), X6 = col_character(), X7 = col_double(), X8 = col_double(),
  X9 = col_double(), X10 = col_double(), X11 = col_double(), X12 = col_double(),
  .default = col_character()
)
paf <- suppressWarnings(read_tsv(paf_file, col_names = FALSE, col_types = cols, progress = FALSE))
if (nrow(paf) == 0) {
  cat("Empty PAF: ", paf_file, "\n", sep = "")
  quit(status = 0)
}

extract_tp <- function(row) {
  if (ncol(row) < 13) {
    return(NA_character_)
  }
  for (i in 13:ncol(row)) {
    v <- row[[i]]
    if (!is.na(v) && str_starts(v, "tp:A:")) {
      return(substr(v, 6, 6))
    }
  }
  NA_character_
}

# --- BEGIN PATCH: v0.0.3 primary_only handling ---
if (primary_only) {
  # Work only on tag columns (13..n); some rows may be shorter (NA)
  if (ncol(paf) < 13) {
    # No tags at all → nothing to filter; keep all as primary
    paf$tp <- "P"
  } else {
    tag_mat <- as.matrix(paf[, 13:ncol(paf), drop = FALSE])
    # Ensure it's character for grepl
    storage.mode(tag_mat) <- "character"

    # For each row, pick the first tag matching ^tp:A:
    tp <- apply(tag_mat, 1, function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      hit <- v[grepl("^tp:A:", v, perl = TRUE)]
      if (length(hit)) substring(hit[1], 6, 6) else NA_character_
    })

    paf$tp <- tp
  }

  # Keep only primary hits when tp is available; if all NA, this keeps none
  # If you want to **keep all** when tp is missing, change `tp == "P"` to `is.na(tp) | tp == "P"`.
  paf <- paf %>% dplyr::filter(.data$tp == "P")
}
# --- END PATCH ---



df <- paf %>%
  transmute(
    qname = X1, qlen = X2, qstart = X3, qend = X4, strand = X5,
    tname = X6, tlen = X7, tstart = X8, tend = X9,
    nmatch = X10, alnlen = X11, mapq = X12
  ) %>%
  mutate(
    id   = if_else(alnlen > 0, nmatch / alnlen, NA_real_),
    qcov = if_else(qlen > 0, (qend - qstart) / qlen, NA_real_)
  ) %>%
  filter(!is.na(id), !is.na(qcov))

if (!is.na(min_alnlen) && min_alnlen > 0) df <- df %>% filter(alnlen >= min_alnlen)
total_n <- nrow(df)
if (total_n == 0) {
  cat("No usable PAF rows after prefilter.\n")
  quit(status = 0)
}

sweep_vec <- function(from = 0.30, to = 0.90, by = 0.01) round(seq(from, to, by), 3)

# 1) Sweep identity; hold qcov fixed
id_grid <- sweep_vec(0.30, 0.90, 0.01)
res_id <- tibble(thresh_id = id_grid) %>%
  mutate(
    fixed_qcov = fixed_qcov,
    pass = map_int(thresh_id, ~ sum(df$id >= .x & df$qcov >= fixed_qcov)),
    pct = 100 * pass / total_n,
    drop = c(NA, diff(pct))
  )
elbow_id <- res_id %>%
  filter(!is.na(drop)) %>%
  arrange(drop) %>%
  slice(1)

# 2) Sweep qcov; hold identity fixed
qc_grid <- sweep_vec(0.30, 0.90, 0.01)
res_qc <- tibble(thresh_qcov = qc_grid) %>%
  mutate(
    fixed_id = fixed_id,
    pass = map_int(thresh_qcov, ~ sum(df$qcov >= .x & df$id >= fixed_id)),
    pct = 100 * pass / total_n,
    drop = c(NA, diff(pct))
  )
elbow_qc <- res_qc %>%
  filter(!is.na(drop)) %>%
  arrange(drop) %>%
  slice(1)

# 3) Combined
res_long <- bind_rows(
  res_id %>% transmute(metric = "identity", thresh = thresh_id, fixed = fixed_qcov, pass, pct, total = total_n, drop),
  res_qc %>% transmute(metric = "qcov", thresh = thresh_qcov, fixed = fixed_id, pass, pct, total = total_n, drop)
)

# Write TSVs
write_tsv(res_id, paste0(out_prefix, ".sweep_identity.tsv"))
write_tsv(res_qc, paste0(out_prefix, ".sweep_qcov.tsv"))
write_tsv(res_long, paste0(out_prefix, ".sweep_combined.tsv"))

# Elbow messages
cat(sprintf(
  "[identity v0.0.2] total=%d fixed_qcov=%.2f  biggest_drop @ id=%.2f  drop=%.2f%%  pct_after=%.2f%% (pass=%d)\n",
  total_n, fixed_qcov, elbow_id$thresh_id, elbow_id$drop, elbow_id$pct, elbow_id$pass
))
cat(sprintf(
  "[qcov     v0.0.2] total=%d fixed_id=%.2f    biggest_drop @ qcov=%.2f drop=%.2f%%  pct_after=%.2f%% (pass=%d)\n",
  total_n, fixed_id, elbow_qc$thresh_qcov, elbow_qc$drop, elbow_qc$pct, elbow_qc$pass
))

# Plots with elbow annotations
p1 <- ggplot(res_id, aes(thresh_id, pct)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_vline(xintercept = elbow_id$thresh_id, linetype = "dashed", color = "red") +
  annotate("text",
    x = elbow_id$thresh_id, y = max(res_id$pct, na.rm = TRUE),
    label = sprintf("elbow id=%.2f\n%%=%.1f", elbow_id$thresh_id, elbow_id$pct),
    hjust = -0.05, vjust = 1.1, size = 3, color = "red"
  ) +
  labs(title = "Percent passing vs identity (with elbow)", x = "identity threshold", y = "% passing") +
  theme_bw(base_size = 12)
ggsave(paste0(out_prefix, ".sweep_identity.png"), p1, width = 6.5, height = 4.0, dpi = 150)

p2 <- ggplot(res_qc, aes(thresh_qcov, pct)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_vline(xintercept = elbow_qc$thresh_qcov, linetype = "dashed", color = "red") +
  annotate("text",
    x = elbow_qc$thresh_qcov, y = max(res_qc$pct, na.rm = TRUE),
    label = sprintf("elbow qcov=%.2f\n%%=%.1f", elbow_qc$thresh_qcov, elbow_qc$pct),
    hjust = -0.05, vjust = 1.1, size = 3, color = "red"
  ) +
  labs(title = "Percent passing vs qcov (with elbow)", x = "qcov threshold", y = "% passing") +
  theme_bw(base_size = 12)
ggsave(paste0(out_prefix, ".sweep_qcov.png"), p2, width = 6.5, height = 4.0, dpi = 150)
