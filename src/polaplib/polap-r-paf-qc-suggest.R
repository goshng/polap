#!/usr/bin/env Rscript
# polap-r-paf-qc-suggest.R  v0.1.1
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

opt_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "allvsall.paf[.gz]"),
  make_option(c("-o", "--outdir"),
    type = "character", default = "paf_qc",
    help = "output directory [paf_qc]"
  ),
  make_option(c("--ident_shift"),
    type = "double", default = 0.05,
    help = "subtract this from center to set min_ident [0.05]"
  ),
  make_option(c("--min_olen_floor"),
    type = "integer", default = 800,
    help = "minimum floor for min_olen [800]"
  )
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input)) stop("Usage: polap-r-paf-qc-suggest.R -i allvsall.paf[.gz] [-o outdir]")

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- Robust read of exactly the first 12 PAF columns (pad missing with NA) ---
# Use fixed names V1..V12 at read-time to avoid name-length issues.
ct <- cols(
  V1 = col_character(), # qname
  V2 = col_double(), # qlen
  V3 = col_double(), # qstart
  V4 = col_double(), # qend
  V5 = col_character(), # strand
  V6 = col_character(), # tname
  V7 = col_double(), # tlen
  V8 = col_double(), # tstart
  V9 = col_double(), # tend
  V10 = col_double(), # nmatch
  V11 = col_double(), # alen
  V12 = col_double() # mapq
)

df <- suppressMessages(
  read_tsv(
    file = opt$input,
    col_types = ct,
    col_names = paste0("V", 1:12), # force 12 names
    col_select = 1:12, # and only those 12 cols
    comment = "#",
    progress = FALSE
  )
)

# If nothing read, emit defaults and exit gracefully
if (nrow(df) == 0) {
  cat("min_olen=1000\nmin_ident=0.84\nw_floor=0.12\n", sep = "")
  quit(save = "no", status = 0)
}

# Now safely rename (we know we have exactly 12 columns)
names(df) <- c(
  "qname", "qlen", "qstart", "qend", "strand",
  "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"
)

# Drop rows that cannot be used (NA numerics or zero/negative lengths)
df <- df %>%
  filter(!is.na(qlen), !is.na(tlen), !is.na(alen), !is.na(nmatch)) %>%
  filter(qlen > 0, tlen > 0, alen > 0)

# If still empty after cleaning, emit defaults
if (nrow(df) == 0) {
  cat("min_olen=1000\nmin_ident=0.84\nw_floor=0.12\n", sep = "")
  quit(save = "no", status = 0)
}

# Metrics
df <- df %>%
  mutate(
    identity = pmax(0, pmin(1, nmatch / alen)),
    frac     = pmax(0, pmin(1, alen / pmin(qlen, tlen))),
    weight   = identity * frac
  )

# Summary table
summ <- df %>%
  summarise(
    n_align         = n(),
    mean_identity   = mean(identity, na.rm = TRUE),
    median_identity = median(identity, na.rm = TRUE),
    mean_alen       = mean(alen, na.rm = TRUE),
    median_alen     = median(alen, na.rm = TRUE),
    mean_weight     = mean(weight, na.rm = TRUE),
    median_weight   = median(weight, na.rm = TRUE)
  )

write_tsv(summ, file.path(opt$outdir, "paf_summary.tsv"))

# Plots
p_alen <- ggplot(df, aes(alen)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(trans = "log10") +
  theme_bw(base_size = 11) +
  labs(title = "Alignment length", x = "aligned length (bp, log10)", y = "count")

p_ident <- ggplot(df, aes(identity)) +
  geom_histogram(bins = 100) +
  theme_bw(base_size = 11) +
  labs(title = "Identity", x = "nmatch/alen", y = "count")

p_weight <- ggplot(df, aes(weight)) +
  geom_histogram(bins = 100) +
  theme_bw(base_size = 11) +
  labs(title = "Edge weight = identity * (alen / min(qlen,tlen))", x = "weight", y = "count")

p_hex <- ggplot(df, aes(alen, identity)) +
  geom_hex(bins = 80) +
  scale_x_continuous(trans = "log10") +
  theme_bw(base_size = 11) +
  labs(title = "Identity vs aligned length", x = "aligned length (bp, log10)", y = "identity")

ggsave(file.path(opt$outdir, "alen_hist.pdf"), p_alen, width = 6, height = 4)
ggsave(file.path(opt$outdir, "identity_hist.pdf"), p_ident, width = 6, height = 4)
ggsave(file.path(opt$outdir, "weight_hist.pdf"), p_weight, width = 6, height = 4)
ggsave(file.path(opt$outdir, "identity_vs_alen.pdf"), p_hex, width = 7, height = 5)

# Suggestions (data-driven)
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))

min_olen_sug <- max(opt$min_olen_floor, floor(q(df$alen, 0.50)))
center_ident <- q(df$identity, 0.65)
min_ident_sug <- max(0.78, min(0.95, center_ident - opt$ident_shift))
w25 <- q(df$weight, 0.25)
w_floor_sug <- max(0.05, min(0.20, w25))

# Emit shell-friendly suggestions to stdout
cat(sprintf("min_olen=%d\n", min_olen_sug))
cat(sprintf("min_ident=%.3f\n", min_ident_sug))
cat(sprintf("w_floor=%.3f\n", w_floor_sug))

# Also save suggestions TSV
tibble(
  param = c("min_olen", "min_ident", "w_floor"),
  value = c(min_olen_sug, min_ident_sug, w_floor_sug)
) %>% write_tsv(file.path(opt$outdir, "suggestions.tsv"))
