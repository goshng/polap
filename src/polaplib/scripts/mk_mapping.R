#!/usr/bin/env Rscript
# FILE: scripts/mk_mapping.R
# VERSION: 0.3.4
# SPDX-License-Identifier: GPL-3.0-or-later
# Usage: Rscript mk_mapping.R <treefile> <presence.tsv> <outdir> <model>
# Models: mk | fitch | dollo | all

suppressPackageStartupMessages({
  library(data.table)
  library(ape)  # read.tree, ace
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: mk_mapping.R <treefile> <presence.tsv> <outdir> <model>")
}
treef <- args[1]
presf <- args[2]
outd  <- args[3]
model <- tolower(args[4])
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

# -------- I/O --------
tr0 <- read.tree(treef)
if (is.null(tr0) || length(tr0$tip.label) == 0)
  stop("Empty/invalid tree: ", treef)

P0 <- fread(presf)
if (!"CID" %in% names(P0)) stop("presence.tsv must have a 'CID' column.")
species_cols0 <- setdiff(colnames(P0), "CID")
if (!length(species_cols0)) stop("presence.tsv has no species columns.")

# -------- Normalization helpers --------
normalize <- function(x) {
  # Trim and strip a common suffix of the form '-digits' (e.g., "-0")
  x <- trimws(as.character(x))
  sub("-[0-9]+$", "", x, perl = TRUE)
}

# Original and normalized labels
tips_raw   <- tr0$tip.label
tips_norm  <- normalize(tips_raw)
pres_raw   <- species_cols0
pres_norm  <- normalize(pres_raw)

# Intersect in normalized space
in_both_norm <- intersect(tips_norm, pres_norm)
if (length(in_both_norm) == 0L) {
  # Helpful diagnostics
  stop(
    "No overlapping species between tree tips and presence.tsv columns.\n",
    "  tree tips (n=", length(tips_raw), "): ",
    paste(head(tips_raw, 10), collapse=", "), " ...\n",
    "  presence cols (n=", length(pres_raw), "): ",
    paste(head(pres_raw, 10), collapse=", "), " ...\n",
    "  tree tips (normalized, n=", length(tips_norm), "): ",
    paste(head(tips_norm, 10), collapse=", "), " ...\n",
    "  presence cols (normalized, n=", length(pres_norm), "): ",
    paste(head(pres_norm, 10), collapse=", "), " ..."
  )
}

# Keep tree tips that are present after normalization; preserve tree order
keep_idx <- which(tips_norm %in% in_both_norm)
tr <- if (length(keep_idx) < length(tips_norm)) keep.tip(tr0, tips_raw[keep_idx]) else tr0

# Recompute normalized vectors after pruning
tips_raw  <- tr$tip.label
tips_norm <- normalize(tips_raw)

# Build a mapping from normalized -> presence original columns
# If normalization created duplicates on presence side, we'll collapse by OR later
norm_to_pres <- split(pres_raw, f = pres_norm)  # list: norm -> vector of original colnames
# Ordered species (normalized) following the tree
species_norm_order <- tips_norm

# -------- Build presence matrix aligned to tree order --------
# Collapse duplicates (same normalized name) by OR across columns
collapse_by_norm <- function(DT, norm_map, order_norm) {
  out <- vector("list", length(order_norm))
  names(out) <- order_norm
  for (nm in order_norm) {
    cols <- norm_map[[nm]]
    if (is.null(cols)) {
      # species absent in presence table -> all zeros
      out[[nm]] <- integer(nrow(DT))
    } else {
      M <- as.matrix(DT[, ..cols])
      mode(M) <- "integer"
      M[is.na(M)] <- 0L
      M[M != 0L] <- 1L
      # OR across columns
      out[[nm]] <- as.integer(rowSums(M) > 0L)
    }
  }
  as.data.table(out)
}

Pres_aligned <- collapse_by_norm(P0, norm_to_pres, species_norm_order)
Pres_aligned[, CID := P0$CID]
setcolorder(Pres_aligned, c("CID", species_norm_order))

# -------- Tree scaffolding for parsimony --------
Ntip  <- length(tr$tip.label)
Nnode <- tr$Nnode
Ntot  <- Ntip + Nnode

children <- vector("list", Ntot); parent <- integer(Ntot); parent[] <- NA_integer_
for (i in seq_len(nrow(tr$edge))) {
  a <- tr$edge[i,1]; b <- tr$edge[i,2]
  parent[b] <- a
  children[[a]] <- c(children[[a]], b)
}
for (v in seq_len(Ntot)) if (is.null(children[[v]])) children[[v]] <- integer(0)
root <- setdiff(tr$edge[,1], tr$edge[,2])[1]

postorder <- function(root, children) {
  ord <- integer(0)
  stack_v <- c(root); stack_vis <- c(FALSE)
  while (length(stack_v)) {
    v <- tail(stack_v, 1L); vis <- tail(stack_vis, 1L)
    stack_v <- head(stack_v, -1L); stack_vis <- head(stack_vis, -1L)
    if (!vis) {
      stack_v   <- c(stack_v, v, rev(children[[v]]))
      stack_vis <- c(stack_vis, TRUE, rep(FALSE, length(children[[v]])))
    } else ord <- c(ord, v)
  }
  ord
}
PO <- postorder(root, children)

# Fitch set encoding: 1={0}, 2={1}, 3={0,1}
enc_tip <- function(x01) {
  x01 <- ifelse(is.na(x01), 0L, x01)
  ifelse(x01==0L, 1L, 2L)
}

fitch_assign <- function(y01, dollo_bias=FALSE) {
  S <- integer(Ntot); S[] <- 0L
  # y01 is named by tip labels (in tree order species_norm_order)
  tip_map <- setNames(seq_len(Ntip), tr$tip.label)
  S[seq_len(Ntip)] <- enc_tip(y01[ names(tip_map) ])

  # Postorder: compute sets at internal nodes
  for (v in PO[PO > Ntip]) {
    ch <- children[[v]]; if (!length(ch)) next
    inter_val <- S[ch[1L]]
    if (length(ch) > 1L) for (u in ch[-1L]) inter_val <- bitwAnd(inter_val, S[u])
    if (inter_val != 0L) {
      S[v] <- inter_val
    } else {
      uni_val <- S[ch[1L]]
      if (length(ch) > 1L) for (u in ch[-1L]) uni_val <- bitwOr(uni_val, S[u])
      S[v] <- uni_val
    }
  }

  # Preorder: select states
  assign <- integer(Ntot); assign[] <- NA_integer_
  rset <- S[root]
  root_state <- if (rset==3L) { if (dollo_bias) 0L else 1L } else if (rset==2L) 1L else 0L
  assign[root] <- root_state
  queue <- children[[root]]
  while (length(queue)) {
    v <- queue[1]; queue <- queue[-1]
    pstate <- assign[parent[v]]
    vset <- S[v]
    if (pstate==0L) assign[v] <- if (vset %in% c(1L,3L)) 0L else 1L
    else            assign[v] <- if (vset %in% c(2L,3L)) 1L else 0L
    if (length(children[[v]])) queue <- c(queue, children[[v]])
  }

  gains <- 0L; losses <- 0L
  for (iE in seq_len(nrow(tr$edge))) {
    a <- tr$edge[iE,1]; b <- tr$edge[iE,2]
    if (assign[a]==0L && assign[b]==1L) gains  <- gains + 1L
    if (assign[a]==1L && assign[b]==0L) losses <- losses + 1L
  }
  list(gains=as.integer(gains), losses=as.integer(losses))
}

# -------- ASR across clusters --------
species <- tr$tip.label            # original tip order; columns of Pres_aligned are normalized but in same order
Y <- as.matrix(Pres_aligned[, ..species_norm_order])
mode(Y) <- "integer"; Y[is.na(Y)] <- 0L; Y[Y!=0L] <- 1L
rownames(Y) <- Pres_aligned$CID

res_list <- vector("list", nrow(Pres_aligned))
for (i in seq_len(nrow(Pres_aligned))) {
  cid <- Pres_aligned$CID[i]
  y <- Y[cid, ]
  names(y) <- species  # map to actual tip labels; order consistent with species_norm_order

  rows <- list()

  if (model %in% c("fitch","all")) {
    gl <- fitch_assign(y, dollo_bias = FALSE)
    rows[[length(rows)+1]] <- data.table(CID=cid, model="fitch", gains=gl$gains, losses=gl$losses)
  }
  if (model %in% c("dollo","all")) {
    gl <- fitch_assign(y, dollo_bias = TRUE)
    rows[[length(rows)+1]] <- data.table(CID=cid, model="dollo_like", gains=gl$gains, losses=gl$losses)
  }
  if (model %in% c("mk","all")) {
    # Use factor states to enforce columns "0" and "1"
    vec <- factor(as.integer(y), levels=c(0L,1L))
    names(vec) <- species
    ok <- TRUE; mk_g <- NA_integer_; mk_l <- NA_integer_
    suppressWarnings({
      fit <- try(ace(vec, tr, type="discrete", model="ER"), silent=TRUE)
      if (inherits(fit, "try-error")) ok <- FALSE
    })
    if (ok && !is.null(fit$lik.anc)) {
      # Columns correspond to levels(vec): "0","1"
      colnames_la <- colnames(fit$lik.anc)
      if (!all(c("0","1") %in% colnames_la)) {
        # Fallback: assume second column is "1" if unnamed
        colnames_la <- if (ncol(fit$lik.anc)==2) c("0","1") else colnames_la
      }
      assign <- integer(Ntot); assign[] <- NA_integer_
      assign[seq_len(Ntip)] <- as.integer(as.character(vec))  # 0/1 at tips
      for (k in seq_len(Nnode)) {
        nd <- Ntip + k
        p1 <- fit$lik.anc[k, "1"]
        assign[nd] <- as.integer(p1 >= 0.5)
      }
      g <- 0L; l <- 0L
      for (iE in seq_len(nrow(tr$edge))) {
        a <- tr$edge[iE,1]; b <- tr$edge[iE,2]
        if (assign[a]==0L && assign[b]==1L) g <- g + 1L
        if (assign[a]==1L && assign[b]==0L) l <- l + 1L
      }
      mk_g <- as.integer(g); mk_l <- as.integer(l)
    } else {
      warning("Mk (ace) failed for CID=", cid, "; skipping.")
    }
    rows[[length(rows)+1]] <- data.table(CID=cid, model="mk_ER", gains=mk_g, losses=mk_l)
  }

  res_list[[i]] <- data.table::rbindlist(rows, fill=TRUE)
}

EVENTS <- data.table::rbindlist(res_list, fill=TRUE)

# -------- Write outputs --------
# Write a normalized, tree-ordered presence matrix so downstream is deterministic
P_out <- copy(Pres_aligned)
setcolorder(P_out, c("CID", species_norm_order))
fwrite(P_out, file=file.path(outd, "presence.tsv"), sep="\t")
fwrite(EVENTS, file=file.path(outd, "events.tsv"), sep="\t")
# timeline placeholder (dating optional)
fwrite(data.table(CID=P_out$CID, mid_time=NA_real_), file=file.path(outd, "timeline.tsv"), sep="\t")
cat("ASR done: presence.tsv, events.tsv, timeline.tsv\n")
