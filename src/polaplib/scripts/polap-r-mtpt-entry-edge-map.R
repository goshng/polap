#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-entry-edge-map.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Plot one symbol per CID at its inferred entry edge (first 0->1), under selected model.

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)
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
model <- tolower(getA("--model", "fitch")) # fitch|dollo|mk_ER
chronos_on <- isTRUE(getA("--chronos", L = TRUE))
lambda <- as.numeric(getA("--lambda", "1.0"))
strip_ts <- getA("--strip-tree-suffix", "")
strip_ps <- getA("--strip-presence-suffix", "")
label_cids <- strsplit(getA("--label-cids", ""), ",")[[1]]
show_legend <- isTRUE(getA("--legend", L = TRUE))
tip_size <- as.numeric(getA("--tip-size", "2.5"))

if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--model fitch|dollo|mk_ER] [--chronos] [--lambda 1.0] [--strip-tree-suffix -0] [--strip-presence-suffix -0] [--label-cids cid1,cid2] [--legend]")
}

al <- align_tree_presence(treef, presf, strip_ts, strip_ps)
tr <- al$tr
Pres <- al$Pres
if (chronos_on) {
  tr <- tryCatch(chronos(tr, lambda = lambda, quiet = TRUE), error = function(e) {
    message("[chronos] failed; using original.")
    tr
  })
}
sc <- tree_scaffold(tr)
tips <- tr$tip.label

# ------------------------------------------------------------------------------
# Fixed
# ------------------------------------------------------------------------------
# Y <- as.matrix(Pres[, ..al$species_norm_order])
# mode(Y) <- "integer"
# Y[is.na(Y)] <- 0L
# Y[Y != 0L] <- 1L
# rownames(Y) <- Pres$CID
Y <- presence_matrix_01(al$Pres, al$species_norm_order)
tips <- tr$tip.label

entries <- vector("list", nrow(Pres))
for (i in seq_len(nrow(Pres))) {
  cid <- Pres$CID[i]
  y <- Y[cid, ]
  names(y) <- tips
  e <- entry_edge_for_cluster(y, tr, sc, model = model)
  entries[[i]] <- data.table(CID = cid, edge = if (is.null(e)) NA_integer_ else e$edge)
}
ENT <- rbindlist(entries)
ENT <- ENT[!is.na(edge)]

mid <- edge_midpoints_df(tr)
ENT <- merge(ENT, mid[, c("edge", "x_mid", "y_mid")], by = "edge", all.x = TRUE)

p <- ggtree(tr, size = 0.4) + theme_tree2() + geom_tiplab(size = tip_size)
p <- p + geom_point(data = ENT, aes(x = x_mid, y = y_mid, color = CID), shape = 21, stroke = 0.8, fill = NA)
if (!show_legend) p <- p + theme(legend.position = "none")
p <- p + ggtitle(sprintf(
  "Entry edges per MTPT cluster (%s)%s",
  toupper(model), if (chronos_on) " [chronos]" else ""
))

ggsave(outf, p, width = 11, height = 8.5, units = "in")
cat("Wrote:", outf, "\n")
