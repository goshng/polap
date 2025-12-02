#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-cluster-panels.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Multi-page PDF, one page per chosen CID with:
#  - tree, arrow at entry edge, and tip presence strip
#  - optional text box: genes, len_bp, pid (from --meta)

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
cids <- getA("--cids", "")
cids_file <- getA("--cid-file", "")
model <- tolower(getA("--model", "fitch"))
strip_ts <- getA("--strip-tree-suffix", "")
strip_ps <- getA("--strip-presence-suffix", "")
meta <- getA("--meta", "")
if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--cids CID1,CID2|--cid-file cids.txt] [--model fitch|dollo|mk_ER] [--meta meta.tsv]")
}

al <- align_tree_presence(treef, presf, strip_ts, strip_ps)
tr <- al$tr
sc <- tree_scaffold(tr)
tips <- tr$tip.label
Pres <- copy(al$Pres)

# Y <- as.matrix(Pres[, ..al$species_norm_order])
# mode(Y) <- "integer"
# Y[Y != 0L] <- 1L
# rownames(Y) <- Pres$CID
Y <- presence_matrix_01(al$Pres, al$species_norm_order)

want <- character(0)
if (nzchar(cids)) want <- strsplit(cids, ",")[[1]]
if (nzchar(cids_file) && file.exists(cids_file)) want <- c(want, scan(cids_file, what = "character", quiet = TRUE))
want <- unique(want)
if (!length(want)) want <- Pres$CID

META <- NULL
if (nzchar(meta) && file.exists(meta)) {
  META <- fread(meta)
}

pdf(outf, width = 11, height = 8.5, onefile = TRUE)
for (cid in want) {
  if (!cid %in% Pres$CID) next
  y <- Y[cid, ]
  names(y) <- tips
  e <- entry_edge_for_cluster(y, tr, sc, model = model)
  p <- ggtree(tr, size = 0.5) + theme_tree2()
  p <- p + geom_tiplab(size = 2.6)
  # presence strip at tip
  D <- data.frame(sp = tips, pres = as.integer(y))
  p <- p + geom_point(data = D, aes(x = max(p$data$x) + 0.02, y = match(sp, tips), shape = factor(pres)), size = 2)
  if (!is.null(e)) {
    mid <- edge_midpoints_df(tr)
    pt <- mid[mid$edge == e$edge, ]
    p <- p + annotate("segment",
      x = pt$x_parent, xend = pt$x_child, y = pt$y_parent, yend = pt$y_child,
      colour = "red", size = 1.1, alpha = 0.8
    )
  }
  title <- cid
  if (!is.null(META) && cid %in% META$CID) {
    mr <- META[META$CID == cid][1]
    title <- sprintf(
      "%s  | len=%s, pid=%s, genes=%s",
      cid,
      if ("len_bp" %in% names(mr)) mr$len_bp else "NA",
      if ("pid" %in% names(mr)) mr$pid else "NA",
      if ("genes" %in% names(mr)) mr$genes else "NA"
    )
  }
  print(p + ggtitle(title))
}
dev.off()
cat("Wrote:", outf, "\n")
