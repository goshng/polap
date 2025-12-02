#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-concordance-map.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Show, for each edge, which models (Fitch, Dollo-like, Mk ER) support at least one entry.

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
strip_ts <- getA("--strip-tree-suffix", "")
strip_ps <- getA("--strip-presence-suffix", "")
chronos_on <- isTRUE(getA("--chronos", L = TRUE))
lambda <- as.numeric(getA("--lambda", "1.0"))
if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--chronos] [--lambda 1.0]")
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

# ------------------------------------------------------------------------------
# Fixed
# ------------------------------------------------------------------------------
# Y <- as.matrix(al$Pres[, ..al$species_norm_order])
# mode(Y) <- "integer"
# Y[Y != 0L] <- 1L
# rownames(Y) <- al$Pres$CID
Y <- presence_matrix_01(al$Pres, al$species_norm_order)
tips <- tr$tip.label

edge_has <- function(model) {
  g <- logical(nrow(tr$edge))
  for (i in seq_len(nrow(al$Pres))) {
    y <- Y[i, ]
    names(y) <- tips
    e <- entry_edge_for_cluster(y, tr, sc, model = model)
    if (!is.null(e)) g[e$edge] <- TRUE
  }
  g
}
g_fitch <- edge_has("fitch")
g_dollo <- edge_has("dollo")
g_mk <- edge_has("mk_ER")

mid <- edge_midpoints_df(tr)
DF <- mid
DF$fitch <- as.integer(g_fitch)
DF$dollo <- as.integer(g_dollo)
DF$mk <- as.integer(g_mk)

p <- ggtree(tr, size = 0.4) + theme_tree2() + geom_tiplab(size = 2.5)
# place shapes with offsets
if (any(DF$fitch == 1)) p <- p + geom_point(data = DF[DF$fitch == 1, ], aes(x = x_mid, y = y_mid + 0.02), shape = 24)
if (any(DF$dollo == 1)) p <- p + geom_point(data = DF[DF$dollo == 1, ], aes(x = x_mid, y = y_mid), shape = 22)
if (any(DF$mk == 1)) p <- p + geom_point(data = DF[DF$mk == 1, ], aes(x = x_mid, y = y_mid - 0.02), shape = 21)
p <- p + guides(shape = "none") +
  ggtitle(sprintf(
    "Model concordance of entry edges (▲ Fitch, ■ Dollo-like, ● Mk)%s",
    if (chronos_on) " [chronos]" else ""
  ))

ggsave(outf, p, width = 11, height = 8.5, units = "in")
cat("Wrote:", outf, "\n")
