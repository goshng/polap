#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-common.R
# VERSION: 0.3.1
# SPDX-License-Identifier: GPL-3.0-or-later
# Helpers to (i) normalize labels, (ii) align presence to tree, (iii) build tree
# scaffolds, (iv) parsimony/Mk assignments, (v) entry-edge calls, and (vi) edge midpoints.

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)  # plotting and convenient node coordinates
})

normalize_suffix <- function(x, suffix="") {
  x <- trimws(as.character(x))
  if (nzchar(suffix)) {
    x <- ifelse(endsWith(x, suffix), substr(x, 1, nchar(x) - nchar(suffix)), x)
  } else {
    x <- sub("-[0-9]+$", "", x, perl=TRUE)  # default robustness: strip "-digits"
  }
  x
}

align_tree_presence <- function(treef, presf, strip_tree_suffix="", strip_presence_suffix="") {
  tr0 <- read.tree(treef)
  if (is.null(tr0) || length(tr0$tip.label)==0) stop("Invalid tree: ", treef)

  P0 <- fread(presf)
  if (!"CID" %in% names(P0)) stop("presence.tsv must contain a 'CID' column.")
  species_raw <- setdiff(names(P0), "CID"); if (!length(species_raw)) stop("No species columns in presence.tsv")

  tips_raw <- tr0$tip.label
  tips_norm <- normalize_suffix(tips_raw, strip_tree_suffix)
  pres_norm <- normalize_suffix(species_raw, if (nzchar(strip_presence_suffix)) strip_presence_suffix else strip_tree_suffix)

  in_both <- intersect(tips_norm, pres_norm)
  if (!length(in_both)) {
    stop("No overlapping species after normalization.\n",
         "tree tips (raw)   : ", paste(head(tips_raw, 10), collapse=", "), "\n",
         "presence (raw)    : ", paste(head(species_raw, 10), collapse=", "), "\n",
         "tree tips (norm)  : ", paste(head(tips_norm, 10), collapse=", "), "\n",
         "presence (norm)   : ", paste(head(pres_norm, 10), collapse=", "))
  }
  keep_idx <- which(tips_norm %in% in_both)
  tr <- if (length(keep_idx) < length(tips_norm)) keep.tip(tr0, tips_raw[keep_idx]) else tr0

  # recompute after pruning
  tips_raw <- tr$tip.label
  tips_norm <- normalize_suffix(tips_raw, strip_tree_suffix)

  # collapse duplicates on presence side by OR
  pres_raw <- species_raw
  norm_to_pres <- split(pres_raw, f=pres_norm)
  order_norm <- tips_norm
  collapse_by_norm <- function(DT, norm_map, order_norm) {
    out <- vector("list", length(order_norm)); names(out) <- order_norm
    for (nm in order_norm) {
      cs <- norm_map[[nm]]
      if (is.null(cs)) out[[nm]] <- integer(nrow(DT))
      else {
        M <- as.matrix(DT[, ..cs]); mode(M) <- "integer"; M[is.na(M)] <- 0L; M[M!=0L] <- 1L
        out[[nm]] <- as.integer(rowSums(M) > 0L)
      }
    }
    as.data.table(out)
  }
  Pres <- collapse_by_norm(P0, norm_to_pres, order_norm)
  Pres[, CID := P0$CID]
  setcolorder(Pres, c("CID", order_norm))

  list(tr=tr, Pres=Pres, tips=tr$tip.label, species_norm_order=order_norm)
}

tree_scaffold <- function(tr) {
  Ntip <- length(tr$tip.label); Nnode <- tr$Nnode; Ntot <- Ntip + Nnode
  children <- vector("list", Ntot); parent <- integer(Ntot); parent[] <- NA_integer_
  for (i in seq_len(nrow(tr$edge))) {
    a <- tr$edge[i,1]; b <- tr$edge[i,2]
    parent[b] <- a
    children[[a]] <- c(children[[a]], b)
  }
  for (v in seq_len(Ntot)) if (is.null(children[[v]])) children[[v]] <- integer(0)
  root <- setdiff(tr$edge[,1], tr$edge[,2])[1]
  # iterative postorder
  postorder <- (function(root, children) {
    ord <- integer(0); stack_v <- c(root); stack_vis <- c(FALSE)
    while (length(stack_v)) {
      v <- tail(stack_v, 1L); vis <- tail(stack_vis, 1L)
      stack_v <- head(stack_v, -1L); stack_vis <- head(stack_vis, -1L)
      if (!vis) { stack_v <- c(stack_v, v, rev(children[[v]])); stack_vis <- c(stack_vis, TRUE, rep(FALSE, length(children[[v]]))) }
      else ord <- c(ord, v)
    }
    ord
  })(root, children)
  depths <- node.depth.edgelength(tr) # distance from root
  # map (parent,child) -> edge index
  edge_df <- data.frame(tr$edge); names(edge_df) <- c("parent","child")
  edge_index <- setNames(seq_len(nrow(tr$edge)), paste(edge_df$parent, edge_df$child, sep=">"))
  list(Ntip=Ntip, Nnode=Nnode, Ntot=Ntot, children=children, parent=parent, root=root,
       postorder=postorder, depths=depths, edge_index=edge_index)
}

# encodings: 1={0}, 2={1}, 3={0,1}
.enc_tip <- function(x01) { x01 <- ifelse(is.na(x01), 0L, x01); ifelse(x01==0L, 1L, 2L) }

fitch_assign <- function(y01, tr, sc, dollo_bias=FALSE, tip_order=NULL) {
  if (is.null(tip_order)) tip_order <- tr$tip.label
  Ntip <- sc$Ntip; Ntot <- sc$Ntot; children <- sc$children; root <- sc$root; PO <- sc$postorder
  S <- integer(Ntot); S[] <- 0L
  y01 <- y01[tip_order]
  S[seq_len(Ntip)] <- .enc_tip(y01)
  for (v in PO[PO > Ntip]) {
    ch <- children[[v]]; if (!length(ch)) next
    inter_val <- S[ch[1L]]
    if (length(ch) > 1L) for (u in ch[-1L]) inter_val <- bitwAnd(inter_val, S[u])
    if (inter_val != 0L) S[v] <- inter_val
    else {
      uni_val <- S[ch[1L]]
      if (length(ch) > 1L) for (u in ch[-1L]) uni_val <- bitwOr(uni_val, S[u])
      S[v] <- uni_val
    }
  }
  # preorder
  assign <- integer(Ntot); assign[] <- NA_integer_
  rset <- S[root]
  root_state <- if (rset==3L) { if (dollo_bias) 0L else 1L } else if (rset==2L) 1L else 0L
  assign[root] <- root_state
  queue <- sc$children[[root]]
  while (length(queue)) {
    v <- queue[1]; queue <- queue[-1]
    pstate <- assign[ sc$parent[v] ]
    vset <- S[v]
    if (pstate==0L) assign[v] <- if (vset %in% c(1L,3L)) 0L else 1L
    else            assign[v] <- if (vset %in% c(2L,3L)) 1L else 0L
    if (length(sc$children[[v]])) queue <- c(queue, sc$children[[v]])
  }
  assign
}

mk_assign <- function(y01, tr, sc, tip_order=NULL) {
  if (is.null(tip_order)) tip_order <- tr$tip.label
  vec <- factor(as.integer(y01[tip_order]), levels=c(0L,1L))
  names(vec) <- tip_order
  fit <- suppressWarnings(try(ace(vec, tr, type="discrete", model="ER"), silent=TRUE))
  if (inherits(fit, "try-error") || is.null(fit$lik.anc)) return(NULL)
  assign <- integer(sc$Ntot); assign[] <- NA_integer_
  assign[seq_len(sc$Ntip)] <- as.integer(as.character(vec))
  for (k in seq_len(sc$Nnode)) {
    nd <- sc$Ntip + k
    p1 <- if (!is.null(colnames(fit$lik.anc)) && "1" %in% colnames(fit$lik.anc))
      fit$lik.anc[k,"1"] else fit$lik.anc[k, min(2, ncol(fit$lik.anc))]
    assign[nd] <- as.integer(p1 >= 0.5)
  }
  assign
}

edge_events_from_assign <- function(assign, tr) {
  gains <- integer(nrow(tr$edge)); losses <- integer(nrow(tr$edge))
  for (e in seq_len(nrow(tr$edge))) {
    a <- tr$edge[e,1]; b <- tr$edge[e,2]
    if (assign[a]==0L && assign[b]==1L) gains[e]  <- gains[e]  + 1L
    if (assign[a]==1L && assign[b]==0L) losses[e] <- losses[e] + 1L
  }
  data.frame(edge=seq_len(nrow(tr$edge)),
             parent=tr$edge[,1], child=tr$edge[,2],
             gain=gains, loss=losses, stringsAsFactors=FALSE)
}

# choose the "entry" edge for one cluster as the gain edge with smallest parent depth
entry_edge_for_cluster <- function(y01, tr, sc, model=c("fitch","dollo","mk_ER")) {
  model <- match.arg(model)
  tip_order <- tr$tip.label
  assign <- switch(model,
                   fitch = fitch_assign(y01, tr, sc, dollo_bias=FALSE, tip_order=tip_order),
                   dollo = fitch_assign(y01, tr, sc, dollo_bias=TRUE,  tip_order=tip_order),
                   mk_ER = mk_assign(y01, tr, sc, tip_order=tip_order))
  if (is.null(assign)) return(NULL)
  df <- edge_events_from_assign(assign, tr)
  gain_edges <- which(df$gain > 0)
  if (!length(gain_edges)) return(NULL)
  # pick earliest (closest to root): minimal parent depth
  depths <- sc$depths
  pe <- df$parent[gain_edges]
  idx <- which.min(depths[pe])
  list(edge=gain_edges[idx], parent=df$parent[gain_edges[idx]], child=df$child[gain_edges[idx]])
}

edge_midpoints_df <- function(tr) {
  p <- ggtree(tr, size=0.2)
  dd <- p$data
  e <- data.frame(tr$edge); names(e) <- c("parent","node")
  d1 <- merge(e, dd[, c("node","x","y")], by="node", all.x=TRUE)
  names(d1)[names(d1)=="x"] <- "x_child"; names(d1)[names(d1)=="y"] <- "y_child"
  d2 <- merge(d1, dd[, c("node","x","y")], by.x="parent", by.y="node", all.x=TRUE)
  names(d2)[names(d2)=="x"] <- "x_parent"; names(d2)[names(d2)=="y"] <- "y_parent"
  d2$x_mid <- (d2$x_child + d2$x_parent)/2
  d2$y_mid <- (d2$y_child + d2$y_parent)/2
  d2$edge <- seq_len(nrow(tr$edge))
  d2
}

entry_times_for_all <- function(tr, sc, entries) {
  # entry time = midpoint distance-from-root; requires branch lengths
  edgel <- tr$edge.length
  dpar  <- sc$depths[entries$parent]
  dpar + edgel[entries$edge] * 0.5
}


# ... keep the existing suppressPackageStartupMessages(...) and helpers ...

# Convert presence DT -> strict 0/1 integer matrix in the exact species order.
presence_matrix_01 <- function(Pres, species_order) {
  if (!all(species_order %in% names(Pres))) {
    missing_cols <- setdiff(species_order, names(Pres))
    stop("[presence_matrix_01] species_order contains columns not in Pres: ",
         paste(missing_cols, collapse=", "))
  }
  DT <- as.data.table(Pres[, ..species_order])

  # Robust per-column coercion (handles "0"/"1", factors, NA/?, etc.)
  for (j in seq_along(species_order)) {
    v <- DT[[j]]
    if (is.factor(v)) v <- as.character(v)
    if (is.character(v)) {
      v[v %in% c("", "NA", "NaN", "nan", "?", "null")] <- NA_character_
      v <- suppressWarnings(as.integer(v))
    } else if (!inherits(v, "integer")) {
      v <- suppressWarnings(as.integer(v))
    }
    v[is.na(v)] <- 0L
    v[v != 0L] <- 1L
    DT[[j]] <- v
  }

  Y <- as.matrix(DT)
  storage.mode(Y) <- "integer"

  # Sanity: keep rows = nrow(Pres) and set rownames safely.
  if (nrow(Y) != nrow(Pres)) {
    stop(sprintf("[presence_matrix_01] matrix rows (%d) != Pres rows (%d).",
                 nrow(Y), nrow(Pres)))
  }
  rownames(Y) <- as.character(Pres$CID)
  Y
}
