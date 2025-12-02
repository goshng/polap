#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-timeline-density.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Density + rug of MTPT entry times (edge midpoints) on a (optionally) time-scaled tree.

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
})
# suppressWarnings(source("scripts/polap-r-mtpt-common.R"))
# suppressWarnings(source("polap/src/polaplib/scripts/polap-r-mtpt-common.R"))

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
chronos_on <- isTRUE(getA("--chronos", L = TRUE))
lambda <- as.numeric(getA("--lambda", "1.0"))

if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--model fitch|dollo|mk_ER] [--chronos] [--lambda 1.0]")
}

# --- Error handlers (quiet by default; stack prints short call chain) ---------
quiet_handler <- function() {
  quit(save = "no", status = as.integer(ecode))
}

stack_handler <- local({
  print_stack <- function() {
    calls <- sys.calls()
    nc <- length(calls)
    if (!nc) {
      cat("Stack (empty)\n")
      return()
    }
    cat("Stack (oldest → newest):\n")
    for (i in seq_len(nc)) {
      cl <- calls[[i]]
      sr <- attr(cl, "srcref")
      txt <- paste(deparse(cl), collapse = " ")
      if (!is.null(sr)) {
        sf <- attr(sr, "srcfile")
        fn <- if (!is.null(sf) && !is.null(sf$filename)) basename(sf$filename) else "<src>"
        cat(sprintf(
          "#%d %s:%d:%d  %s\n",
          i, fn, as.integer(sr[[1]]), as.integer(sr[[2]]), txt
        ))
      } else {
        cat(sprintf("#%d %s\n", i, txt))
      }
    }
  }
  function() {
    try(print_stack(), silent = TRUE)
    quit(save = "no", status = as.integer(ecode))
  }
})

ecode <- 1
wantStack <- 1
options(error = if (wantStack) stack_handler else quiet_handler)

# --- Context & optional pre-crash side effect ---------------------------------
cat(sprintf(
  "R test_crash: message\n"
))

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

entries <- vector("list", nrow(al$Pres))
for (i in seq_len(nrow(al$Pres))) {
  cid <- al$Pres$CID[i]
  y <- Y[cid, ]
  names(y) <- tips
  e <- entry_edge_for_cluster(y, tr, sc, model = model)
  if (is.null(e)) next
  entries[[i]] <- data.table(CID = cid, edge = e$edge, parent = e$parent, child = e$child)
}
ENT <- rbindlist(entries, fill = TRUE)
times <- entry_times_for_all(tr, sc, ENT)

D <- data.table(time = times)
p <- ggplot(D, aes(x = time)) +
  geom_density(fill = NA, color = "black") +
  geom_rug(alpha = 0.4) +
  theme_bw(base_size = 12) +
  labs(
    x = "Time from root (branch length units)", y = "Density",
    title = sprintf(
      "Timeline of MTPT entries (%s)%s",
      toupper(model), if (chronos_on) " [chronos]" else ""
    )
  )

ggsave(outf, p, width = 8, height = 4, units = "in")
cat("Wrote:", outf, "\n")
