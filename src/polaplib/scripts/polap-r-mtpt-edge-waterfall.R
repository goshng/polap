#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-edge-waterfall.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Barplot per edge: total gains (up) and losses (down) across all CIDs.

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
model <- tolower(getA("--model", "fitch"))
strip_ts <- getA("--strip-tree-suffix", "")
strip_ps <- getA("--strip-presence-suffix", "")
chronos_on <- isTRUE(getA("--chronos", L = TRUE))
lambda <- as.numeric(getA("--lambda", "1.0"))

if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--model fitch|dollo|mk_ER] [--chronos] [--lambda 1.0] [--strip-tree-suffix -0] [--strip-presence-suffix -0]")
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
# Y[is.na(Y)] <- 0L
# Y[Y != 0L] <- 1L
# rownames(Y) <- al$Pres$CID
Y <- presence_matrix_01(al$Pres, al$species_norm_order)

gains <- integer(nrow(tr$edge))
losses <- integer(nrow(tr$edge))
for (i in seq_len(nrow(al$Pres))) {
  y <- Y[i, ]
  names(y) <- tips
  assign <- switch(model,
    fitch = fitch_assign(y, tr, sc, FALSE, tips),
    dollo = fitch_assign(y, tr, sc, TRUE, tips),
    mk_ER = mk_assign(y, tr, sc, tips),
    fitch_assign(y, tr, sc, FALSE, tips)
  )
  if (is.null(assign)) next
  ev <- edge_events_from_assign(assign, tr)
  gains <- gains + ev$gain
  losses <- losses + ev$loss
}

df <- data.frame(edge = seq_len(nrow(tr$edge)), gain = gains, loss = losses)
df$gainloss <- ifelse(df$gain > 0, df$gain, 0) - ifelse(df$loss > 0, df$loss, 0)

p <- ggplot(df, aes(x = factor(edge), y = gainloss)) +
  geom_col(fill = "grey30") +
  geom_hline(yintercept = 0, color = "black") +
  labs(
    x = "Edge index", y = "Gains (up) / Losses (down)",
    title = sprintf(
      "Edge waterfall of MTPT events (%s)%s",
      toupper(model), if (chronos_on) " [chronos]" else ""
    )
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(outf, p, width = 11, height = 4.5, units = "in")
cat("Wrote:", outf, "\n")
