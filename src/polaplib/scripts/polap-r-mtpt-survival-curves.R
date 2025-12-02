#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-survival-curves.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Kaplan–Meier-like survival of MTPT presence after entry.
# Derives event/censor times from parsimony/Mk node states:
#   - For each descendant lineage of the entry node, a 1->0 edge gives a loss time (event=1).
#   - If tip still 1, right-censor at tip time (event=0).
# Group by --by (e.g., functional class) if provided.

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(survival)
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
groupf <- getA("--group", "") # optional TSV: CID, group

if (is.null(treef) || is.null(presf) || is.null(outf)) {
  stop("Usage: --tree <treefile> --presence <presence.tsv> --out <pdf> [--model fitch|dollo|mk_ER] [--chronos] [--lambda 1.0] [--group cid_group.tsv]")
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
depths <- sc$depths
edge_len <- tr$edge.length

# time of node = distance from root
node_time <- depths

# presence matrix
Pres <- al$Pres
# Y <- as.matrix(Pres[, ..al$species_norm_order])
# mode(Y) <- "integer"
# Y[Y != 0L] <- 1L
# rownames(Y) <- Pres$CID
Y <- presence_matrix_01(al$Pres, al$species_norm_order)
tips <- tr$tip.label

# optional grouping
grpmap <- NULL
if (nzchar(groupf) && file.exists(groupf)) {
  G <- fread(groupf)
  if (!all(c("CID", "group") %in% names(G))) {
    warning("group file missing CID,group; ignored")
  } else {
    grpmap <- setNames(G$group, G$CID)
  }
}

collect <- list()
for (i in seq_len(nrow(Pres))) {
  cid <- Pres$CID[i]
  y <- Y[cid, ]
  names(y) <- tips
  assign <- switch(model,
    fitch = fitch_assign(y, tr, sc, FALSE, tips),
    dollo = fitch_assign(y, tr, sc, TRUE, tips),
    mk_ER = mk_assign(y, tr, sc, tips),
    fitch_assign(y, tr, sc, FALSE, tips)
  )
  if (is.null(assign)) next
  # entry edge
  ent <- entry_edge_for_cluster(y, tr, sc, model = model)
  if (is.null(ent)) next
  entry_time <- node_time[ent$parent] + edge_len[ent$edge] / 2
  # For each descendant tip, determine first loss along path
  # Build parent map and paths by backtracking
  for (t in seq_len(sc$Ntip)) {
    # check if tip is descendant of entry child
    # climb from tip to root, check if we pass through ent$child
    v <- t
    seen_child <- FALSE
    last_state <- assign[v]
    loss_time <- NA_real_
    event <- 0L
    prev <- v
    while (!is.na(sc$parent[v])) {
      p <- sc$parent[v]
      if (p == ent$child) seen_child <- TRUE
      # check 1->0 at edge (p->v)
      if (!is.na(assign[p]) && !is.na(assign[v]) && assign[p] == 1L && assign[v] == 0L && seen_child) {
        # loss time at node v (approx: parent time + half edge)
        ekey <- paste(p, v, sep = ">")
        eidx <- sc$edge_index[ekey]
        loss_time <- node_time[p] + (edge_len[eidx] / 2)
        event <- 1L
        break
      }
      v <- p
    }
    if (!seen_child) next # tip outside the gained clade
    # if not lost, censor at tip time
    if (is.na(loss_time)) {
      loss_time <- node_time[prev] # tip time
      event <- 0L
    }
    collect[[length(collect) + 1]] <- data.table(
      CID = cid, time = loss_time - entry_time, event = event,
      group = if (!is.null(grpmap) && cid %in% names(grpmap)) grpmap[[cid]] else "all"
    )
  }
}
DF <- rbindlist(collect, fill = TRUE)
if (!nrow(DF)) stop("No events/censors derived; check inputs.")

# Fit survival per group
fits <- by(DF, DF$group, function(d) {
  survfit(Surv(time, event) ~ 1, data = d)
})
# Convert to df and plot
to_df <- function(sf, grp) {
  s <- summary(sf)
  data.table(group = grp, time = s$time, surv = s$surv, lower = s$lower, upper = s$upper)
}
PL <- rbindlist(lapply(names(fits), function(g) to_df(fits[[g]], g)), fill = TRUE)

p <- ggplot(PL, aes(x = time, y = surv, color = group)) +
  geom_step() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.15, color = NA) +
  theme_bw(base_size = 12) +
  labs(
    x = "Time since entry (branch units or Myr)",
    y = "Fraction of lineages still carrying MTPT",
    title = "Survival of MTPT presence after entry"
  )

ggsave(outf, p, width = 7.5, height = 5, units = "in")
cat("Wrote:", outf, "\n")
