#!/usr/bin/env Rscript
# Version: v0.2.0
# mtpt_turnover_fig1.R
# Ingest all mtpt.tsv files under a root, summarize per species, and generate Fig. 1:
#   (A) MTPT load (kb) per species
#   (B) Recency mix (recent/intermediate/ancient fractions) per species
#
# Usage:
#   Rscript mtpt_turnover_fig1.R \
#     --root Anthoceros_agrestis-0  --pattern mtpt.tsv \
#     --out out_mtpt_summary  \
#     [--mtlen-tsv mt_lengths.tsv] \
#     [--species-from dirname] \
#     [--strip-suffix "-0$"]
#
# Notes:
#  - Each mtpt.tsv must have at least columns: mt_contig, mt_start, mt_end, length, median_pid, class
#  - Species name is derived from the parent folder name by default (dirname of mtpt.tsv).
#  - If you pass --mtlen-tsv (tab-separated: species, mt_len), script will add %mtDNA.
#
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  root = ".",
  pattern = "mtpt.tsv",
  out = "mtpt_summary_out",
  mtlen = NULL,
  species_from = "dirname", # or "column" if you pre-annotate
  strip_suffix = NULL
)
# naive arg parse
i <- 1
while (i <= length(args)) {
  key <- args[i]
  if (key == "--root") {
    opt$root <- args[i + 1]
    i <- i + 2
  } else if (key == "--pattern") {
    opt$pattern <- args[i + 1]
    i <- i + 2
  } else if (key == "--out") {
    opt$out <- args[i + 1]
    i <- i + 2
  } else if (key == "--mtlen-tsv") {
    opt$mtlen <- args[i + 1]
    i <- i + 2
  } else if (key == "--species-from") {
    opt$species_from <- args[i + 1]
    i <- i + 2
  } else if (key == "--strip-suffix") {
    opt$strip_suffix <- args[i + 1]
    i <- i + 2
  } else {
    stop("Unknown arg: ", key)
  }
}

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# 1) find all mtpt.tsv files
all_files <- list.files(opt$root, pattern = opt$pattern, recursive = TRUE, full.names = TRUE)
if (length(all_files) == 0) stop("No files matching pattern under --root: ", opt$root)

# 2) read & stamp species (species = top-level folder directly under --root)
normalize_slash <- function(x) gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = TRUE))
root_abs <- normalize_slash(opt$root)
# drop trailing slash to make substring math stable
if (substr(root_abs, nchar(root_abs), nchar(root_abs)) == "/") {
  root_abs <- substr(root_abs, 1, nchar(root_abs) - 1)
}

mtpt_list <- list()
for (f in all_files) {
  dt <- fread(f)
  req <- c("mt_contig", "mt_start", "mt_end", "length", "median_pid", "class")
  if (!all(req %in% names(dt))) stop("File missing required columns: ", f)

  f_abs <- normalize_slash(f)
  # literal (non-regex) removal of the --root prefix + one slash
  rel <- sub(paste0(root_abs, "/"), "", f_abs, fixed = TRUE)
  parts <- strsplit(rel, "/")[[1]]
  if (length(parts) < 1L) stop("Cannot derive species from path (no components): ", rel)

  # species = first component under --root (e.g., Brassica_rapa in Brassica_rapa/.../mtpt.tsv)
  sp <- parts[1]
  if (!is.null(opt$strip_suffix)) sp <- sub(opt$strip_suffix, "", sp)

  dt[, species := sp]
  mtpt_list[[length(mtpt_list) + 1]] <- dt
}
mtpt <- rbindlist(mtpt_list, use.names = TRUE, fill = TRUE)

# 3) per-species metrics
# normalize class labels
mtpt[, class := tolower(class)]
mtpt[, class := fifelse(class %chin% c("recent", "intermediate", "ancient"), class, "other")]

sum_by_class <- mtpt[, .(
  n_tracts = .N,
  mtpt_bp = sum(length, na.rm = TRUE),
  recent_bp = sum(length[class == "recent"], na.rm = TRUE),
  interm_bp = sum(length[class == "intermediate"], na.rm = TRUE),
  ancient_bp = sum(length[class == "ancient"], na.rm = TRUE),
  max_tract_len = max(length, na.rm = TRUE),
  max_tract_pid = max(median_pid, na.rm = TRUE)
), by = species]

sum_by_class[, mtpt_kb := mtpt_bp / 1000]
sum_by_class[, total_class_bp := recent_bp + interm_bp + ancient_bp]
sum_by_class[total_class_bp > 0, `:=`(
  recent_frac = recent_bp / total_class_bp,
  interm_frac = interm_bp / total_class_bp,
  ancient_frac = ancient_bp / total_class_bp
)]
sum_by_class[total_class_bp == 0, `:=`(
  recent_frac = 0, interm_frac = 0, ancient_frac = 0
)]

# Optional: add mt genome length to compute %mtDNA
if (!is.null(opt$mtlen)) {
  len_tab <- fread(opt$mtlen)
  if (!all(c("species", "mt_len") %in% names(len_tab))) {
    stop("--mtlen-tsv must have columns: species, mt_len")
  }
  sum_by_class <- merge(sum_by_class, len_tab[, .(species, mt_len = as.numeric(mt_len))],
    by = "species", all.x = TRUE
  )
  sum_by_class[, mtpt_pct_mt := ifelse(is.finite(mt_len) & mt_len > 0, 100 * mtpt_bp / mt_len, NA_real_)]
} else {
  sum_by_class[, mtpt_pct_mt := NA_real_]
}

# write summary
fwrite(sum_by_class[order(-mtpt_kb)], file = file.path(opt$out, "mtpt_summary.tsv"), sep = "\t")

# 4) Fig. 1 — (A) load per species; (B) recency mix
# order species by load (descending)
sum_by_class <- sum_by_class[order(-mtpt_kb)]
sum_by_class[, species_f := factor(species, levels = species)]

# Panel A
pA <- ggplot(sum_by_class, aes(x = species_f, y = mtpt_kb)) +
  geom_col() +
  labs(x = NULL, y = "MTPT load (kb)", title = "Fig. 1A - MTPT load per species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Panel B
mix_long <- melt(sum_by_class[, .(species_f, recent_frac, interm_frac, ancient_frac)],
  id.vars = "species_f",
  variable.name = "class", value.name = "frac"
)
mix_long[, class := factor(class,
  levels = c("recent_frac", "interm_frac", "ancient_frac"),
  labels = c("recent (>=97%)", "intermediate (90-97%)", "ancient (<90%)")
)]

pB <- ggplot(mix_long, aes(x = species_f, y = frac, fill = class)) +
  geom_col(width = 0.95) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL, y = "Recency composition", fill = "Class",
    title = "Fig. 1B - Recency mix per species"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# 5) save a single PDF with both panels
pdf(file.path(opt$out, "fig1_mtpt_turnover.pdf"), width = 12, height = 7)
grid.arrange(pA, pB, ncol = 1, heights = c(1, 1))
dev.off()

# also save individual panels if desired
ggsave(file.path(opt$out, "fig1A_mtpt_load.pdf"), pA, width = 12, height = 4)
ggsave(file.path(opt$out, "fig1B_recency_mix.pdf"), pB, width = 12, height = 4)

cat("[OK] Wrote:\n",
  "  ", file.path(opt$out, "mtpt_summary.tsv"), "\n",
  "  ", file.path(opt$out, "fig1_mtpt_turnover.pdf"), "\n",
  sep = ""
)
