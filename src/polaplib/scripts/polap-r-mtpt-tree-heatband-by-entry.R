#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-tree-heatband-by-entry.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Heatmap (rows=CIDs, cols=species in tree order), grouped by entry edge.
# Optional metadata sidebars:
#   --meta <tsv> with columns: CID, len_bp, pid (0..1 or %; we convert to 0..1)

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(gridExtra)
})
# suppressWarnings(source("scripts/polap-r-mtpt-common.R"))

# ------------------------------------------------------------------------------
# BEGIN: Locate and source the shared helpers (polap-r-mtpt-common.R)
# ------------------------------------------------------------------------------
polap_scripts_dir <- local({
  # 1) If caller exported _POLAPLIB_DIR (root), prefer its scripts/ subdir.
  env_root <- Sys.getenv("_POLAPLIB_DIR", unset = NA)
  if (!is.na(env_root) && nzchar(env_root)) {
    cand <- file.path(normalizePath(env_root), "scripts")
    if (file.exists(cand)) {
      return(cand)
    }
  }
  # 2) Resolve the directory of this running script (Rscript --file=...).
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cmdArgs)
  if (length(m) > 0) {
    this_file <- normalizePath(sub("^--file=", "", cmdArgs[m[length(m)]]))
    return(dirname(this_file))
  }
  # 3) Fallback for sourced/interactive contexts that set 'ofile'.
  of <- tryCatch(normalizePath(attr(sys.frames()[[1]], "ofile")), error = function(e) NA)
  if (!is.na(of)) {
    return(dirname(of))
  }
  # 4) Last resort: current working directory.
  normalizePath(".")
})

common_path <- file.path(polap_scripts_dir, "polap-r-mtpt-common.R")
if (!file.exists(common_path)) {
  stop(
    "Cannot locate polap-r-mtpt-common.R at: ", common_path,
    "\nHint: export _POLAPLIB_DIR=/path/to/polaplib or run via Rscript <polap-lib>/scripts/<file>.R"
  )
}
suppressWarnings(source(common_path))
# ------------------------------------------------------------------------------
# END: Locate and source the shared helpers (polap-r-mtpt-common.R)
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
getA <- function(f, d = NULL, L = FALSE) {
  i <- which(args == f)
  if (!length(i)) {
    return(d)
  }
  if (L) TRUE else if (i == length(args)) d else args[i + 1]
}
treef <- getA("--tree")
presf <- getA("--presence")
outf <- getA("--out")
model <- tolower(getA("--model", "fitch"))
strip_ts <- getA("--strip-tree-suffix", "")
strip_ps <- getA("--strip-presence-suffix", "")
meta <- getA("--meta", "")
chronos_on <- isTRUE(getA("--chronos", L = TRUE))
lambda <- as.numeric(getA("--lambda", "1.0"))
max_rows <- as.integer(getA("--max-rows", "300"))
tip_size <- as.numeric(getA("--tip-size", "2.2"))

if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--model fitch|dollo|mk_ER] [--meta meta.tsv] [--chronos] [--lambda 1.0] [--max-rows 300]")
}

al <- align_tree_presence(treef, presf, strip_ts, strip_ps)
tr <- al$tr
if (chronos_on) {
  tr <- tryCatch(chronos(tr, lambda = lambda, quiet = TRUE), error = function(e) {
    tr
  })
}
sc <- tree_scaffold(tr)
tips <- tr$tip.label
Pres <- copy(al$Pres)

# Y <- as.matrix(Pres[, ..al$species_norm_order])
# mode(Y) <- "integer"
# Y[is.na(Y)] <- 0L
# Y[Y != 0L] <- 1L
# rownames(Y) <- Pres$CID
Y <- presence_matrix_01(al$Pres, al$species_norm_order)
tips <- tr$tip.label

# entry edge per CID
entries <- vector("list", nrow(Pres))
for (i in seq_len(nrow(Pres))) {
  cid <- Pres$CID[i]
  y <- Y[cid, ]
  names(y) <- tips
  e <- entry_edge_for_cluster(y, tr, sc, model = model)
  entries[[i]] <- data.table(CID = cid, entry_edge = if (is.null(e)) NA_integer_ else e$edge)
}
ENT <- rbindlist(entries)
Pres <- merge(Pres, ENT, by = "CID", all.x = TRUE)
Pres[is.na(entry_edge), entry_edge := -1L] # 'unknown' group

# Optional metadata
META <- NULL
if (nzchar(meta) && file.exists(meta)) {
  META <- fread(meta)
  if (!all(c("CID", "len_bp", "pid") %in% names(META))) {
    warning("META missing columns (CID,len_bp,pid); ignoring")
    META <- NULL
  } else {
    META[, pid := ifelse(pid > 1, pid / 100, pid)]
    Pres <- merge(Pres, META, by = "CID", all.x = TRUE)
  }
}

# limit rows for plotting (top by entry edge then CID)
if (nrow(Pres) > max_rows) {
  setorder(Pres, -entry_edge, CID)
  Pres <- Pres[1:max_rows]
}

# Long format for tiles
long <- melt(Pres,
  id.vars = c("CID", "entry_edge", if (!is.null(META)) c("len_bp", "pid") else NULL),
  variable.name = "species", value.name = "presence"
)
long[, presence := as.integer(presence)]

# Heatmap
heat <- ggplot(long, aes(x = species, y = reorder(CID, entry_edge), fill = presence)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", guide = "none") +
  facet_grid(rows = vars(entry_edge), scales = "free_y", space = "free_y", labeller = label_both) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = "CIDs grouped by entry_edge", title = sprintf("CID presence by species (%s)", toupper(model)))

if (!is.null(META)) {
  # Sidebars as small heat tiles at left (two columns: length and pid)
  side <- unique(long[, .(CID, entry_edge, len_bp, pid)])
  side_long <- melt(side, id.vars = c("CID", "entry_edge"), variable.name = "metric", value.name = "value")
  side_long[metric == "pid", value := 1 - value] # encode erosion (1-PID) darker
  side_plot <- ggplot(side_long, aes(x = metric, y = reorder(CID, entry_edge), fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(na.value = "grey90") +
    facet_grid(rows = vars(entry_edge), scales = "free_y", space = "free_y") +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      legend.position = "right"
    )
  g <- grid.arrange(side_plot, heat, ncol = 2, widths = c(1.2, 6))
} else {
  g <- heat
}

ggsave(outf, g, width = 12, height = min(10, 2 + 0.02 * nrow(Pres)), units = "in")
cat("Wrote:", outf, "\n")
