#!/usr/bin/env Rscript
# FILE: scripts/polap-r-mtpt-plot-tree.R
# VERSION: 0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Draw a (time-scaled, optional) plastid (PT) tree with MTPT gain/loss counts per edge.
# - Reads PT tree (Newick) and presence matrix (CID x species columns, 0/1).
# - Normalizes labels on both sides (strip trailing "-digits", e.g., "-0") unless overridden.
# - Aligns presence species to the tree (and collapses duplicates post-normalization by OR).
# - Aggregates per-edge event counts across all clusters under:
#     * Fitch parsimony (unordered)           --model fitch
#     * Dollo-like parsimony (gain-averse)    --model dollo
#     * Likelihood Mk(ER) via ape::ace        --model mk_ER
# - Optionally time-scales the tree with ape::chronos (--chronos --lambda 1.0).
#
# Dependencies: data.table, ape, ggplot2, ggtree (Bioconductor)
#
# Example:
#   Rscript scripts/polap-r-mtpt-plot-tree.R \
#     --tree man/analysis/pt_tree/concat/iqtree.treefile \
#     --presence man/analysis/gain_loss/presence.tsv \
#     --out man/analysis/results/mtpt_events_fitch.pdf \
#     --model fitch --chronos --lambda 1.0 --strip-tree-suffix "-0"
#
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)   # bioconductor-ggtree
})

# ---------------- CLI ----------------
args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, def=NULL, logical=FALSE) {
  i <- which(args==flag)
  if (!length(i)) return(def)
  if (logical) return(TRUE)
  if (i==length(args)) return(def)
  args[i+1]
}
treef    <- getArg("--tree")
presf    <- getArg("--presence")
outpdf   <- getArg("--out")
model    <- tolower(getArg("--model","fitch"))   # fitch|dollo|mk_ER
chronos_on <- isTRUE(getArg("--chronos", logical=TRUE))
lambda   <- as.numeric(getArg("--lambda","1.0"))
tip_size <- as.numeric(getArg("--tip-size","2.5"))
width_in <- as.numeric(getArg("--width","11"))
height_in<- as.numeric(getArg("--height","8.5"))
strip_tree_suffix    <- getArg("--strip-tree-suffix", "")      # e.g., "-0"
strip_presence_suffix<- getArg("--strip-presence-suffix", "")  # blank => auto-use strip_tree_suffix

if (is.null(treef) || is.null(presf) || is.null(outpdf)) {
  stop("Usage: --tree <file> --presence <file> --out <pdf> ",
       "[--model fitch|dollo|mk_ER] [--chronos] [--lambda 1.0] ",
       "[--tip-size 2.5] [--width 11] [--height 8.5] ",
       "[--strip-tree-suffix \"-0\"] [--strip-presence-suffix \"\"]")
}

if (strip_tree_suffix != "" && strip_presence_suffix == "") {
  # Harmonize stripping unless user explicitly provided presence suffix
  strip_presence_suffix <- strip_tree_suffix
}

# --------------- IO ---------------
tr0 <- read.tree(treef)
if (is.null(tr0) || length(tr0$tip.label)==0)
  stop("Invalid tree: ", treef)

P0 <- fread(presf)
stopifnot("CID" %in% colnames(P0))
species_cols0 <- setdiff(colnames(P0), "CID")
stopifnot(length(species_cols0) > 0)

# --------------- Normalization & Alignment ---------------
normalize_suffix <- function(x, suffix) {
  x <- trimws(as.character(x))
  if (nzchar(suffix)) {
    # strip specific suffix if present
    x <- ifelse(endsWith(x, suffix), substr(x, 1, nchar(x) - nchar(suffix)), x)
  } else {
    # default: strip trailing "-digits" (e.g., "-0", "-1") to be robust
    x <- sub("-[0-9]+$", "", x, perl = TRUE)
  }
  x
}

tips_raw  <- tr0$tip.label
tips_norm <- normalize_suffix(tips_raw, strip_tree_suffix)
pres_raw  <- species_cols0
pres_norm <- normalize_suffix(pres_raw, strip_presence_suffix)

in_both_norm <- intersect(tips_norm, pres_norm)
if (!length(in_both_norm)) {
  stop("No overlapping species after normalization.\n",
       "  tree tips (raw, n=", length(tips_raw), "): ",
       paste(head(tips_raw, 10), collapse=", "), " ...\n",
       "  presence cols (raw, n=", length(pres_raw), "): ",
       paste(head(pres_raw, 10), collapse=", "), " ...\n",
       "  tree tips (norm, n=", length(tips_norm), "): ",
       paste(head(tips_norm, 10), collapse=", "), " ...\n",
       "  presence cols (norm, n=", length(pres_norm), "): ",
       paste(head(pres_norm, 10), collapse=", "), " ...\n",
       "Consider adding --strip-tree-suffix and/or --strip-presence-suffix.")
}

# Keep/prune tree to those that exist in presence after normalization
keep_idx <- which(tips_norm %in% in_both_norm)
tr <- if (length(keep_idx) < length(tips_norm)) keep.tip(tr0, tips_raw[keep_idx]) else tr0

# Recompute after pruning
tips_raw  <- tr$tip.label
tips_norm <- normalize_suffix(tips_raw, strip_tree_suffix)

# Map normalized -> original presence columns (collapse duplicates by OR)
norm_to_pres <- split(pres_raw, f = pres_norm)
species_norm_order <- tips_norm  # follow tree order

collapse_by_norm <- function(DT, norm_map, order_norm) {
  out <- vector("list", length(order_norm))
  names(out) <- order_norm
  for (nm in order_norm) {
    cols <- norm_map[[nm]]
    if (is.null(cols)) {
      out[[nm]] <- integer(nrow(DT))
    } else {
      M <- as.matrix(DT[, ..cols])
      mode(M) <- "integer"; M[is.na(M)] <- 0L; M[M!=0L] <- 1L
      out[[nm]] <- as.integer(rowSums(M) > 0L)
    }
  }
  as.data.table(out)
}

Pres <- collapse_by_norm(P0, norm_to_pres, species_norm_order)
Pres[, CID := P0$CID]
setcolorder(Pres, c("CID", species_norm_order))

# --------------- Optional time scaling ---------------
if (chronos_on) {
  tr <- tryCatch(chronos(tr, lambda=lambda, quiet=TRUE),
                 error=function(e) { message("[chronos] failed; using original branch lengths."); tr })
}

# --------------- Tree scaffolding for parsimony/Mk ---------------
Ntip  <- length(tr$tip.label)
Nnode <- tr$Nnode
Ntot  <- Ntip + Nnode

children <- vector("list", Ntot)
parent   <- integer(Ntot); parent[] <- NA_integer_
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

# Fitch encodings: 1={0}, 2={1}, 3={0,1}
enc_tip <- function(x01) { x01 <- ifelse(is.na(x01), 0L, x01); ifelse(x01==0L, 1L, 2L) }

fitch_assign <- function(y01, dollo_bias=FALSE, tip_order) {
  S <- integer(Ntot); S[] <- 0L
  # y01 is named by tip labels; ensure in tree tip order
  y01 <- y01[tip_order]
  S[seq_len(Ntip)] <- enc_tip(y01)
  # postorder
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
  # preorder picks
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
  assign
}

mk_assign <- function(y01, tip_order) {
  vec <- factor(as.integer(y01[tip_order]), levels=c(0L,1L))
  names(vec) <- tip_order
  fit <- suppressWarnings(try(ace(vec, tr, type="discrete", model="ER"), silent=TRUE))
  if (inherits(fit, "try-error") || is.null(fit$lik.anc)) return(NULL)
  assign <- integer(Ntot); assign[] <- NA_integer_
  assign[seq_len(Ntip)] <- as.integer(as.character(vec))
  # fit$lik.anc rows correspond to internal nodes in order 1..Nnode
  for (k in seq_len(Nnode)) {
    nd <- Ntip + k
    # enforce column "1" (prob of state=1); if unnamed, assume 2nd
    p1 <- if (!is.null(colnames(fit$lik.anc)) && "1" %in% colnames(fit$lik.anc))
      fit$lik.anc[k,"1"] else fit$lik.anc[k, min(2, ncol(fit$lik.anc))]
    assign[nd] <- as.integer(p1 >= 0.5)
  }
  assign
}

# Presence matrix aligned to tree order (normalized columns)
species_norm_order <- colnames(Pres)[colnames(Pres)!="CID"]
tip_labels <- tr$tip.label
stopifnot(length(species_norm_order) == length(tip_labels))
Y <- as.matrix(Pres[, ..species_norm_order])
mode(Y) <- "integer"; Y[is.na(Y)] <- 0L; Y[Y!=0L] <- 1L
rownames(Y) <- Pres$CID

# --------------- Aggregate per-edge events over clusters ---------------
edge_gains <- integer(nrow(tr$edge))
edge_losses<- integer(nrow(tr$edge))

for (i in seq_len(nrow(Pres))) {
  cid <- Pres$CID[i]
  y <- Y[cid, ]
  names(y) <- tip_labels

  assign <- switch(model,
    "fitch" = fitch_assign(y, dollo_bias=FALSE, tip_order=tip_labels),
    "dollo" = fitch_assign(y, dollo_bias=TRUE,  tip_order=tip_labels),
    "mk_er" = mk_assign(y, tip_order=tip_labels),
    fitch_assign(y, dollo_bias=FALSE, tip_order=tip_labels)
  )
  if (is.null(assign)) next

  for (e in seq_len(nrow(tr$edge))) {
    a <- tr$edge[e,1]; b <- tr$edge[e,2]
    if (assign[a]==0L && assign[b]==1L) edge_gains[e]  <- edge_gains[e]  + 1L
    if (assign[a]==1L && assign[b]==0L) edge_losses[e] <- edge_losses[e] + 1L
  }
}

# --------------- Plot with ggtree ---------------
p <- ggtree(tr, size=0.4) + theme_tree2()
dd <- p$data

edge_df <- data.frame(tr$edge)
colnames(edge_df) <- c("parent","node")
edge_df$gain <- edge_gains
edge_df$loss <- edge_losses

# Node coordinates for midpoints
dd2 <- merge(edge_df, dd[, c("node","x","y")], by="node", all.x=TRUE)
colnames(dd2)[colnames(dd2)=="x"] <- "x_child"
colnames(dd2)[colnames(dd2)=="y"] <- "y_child"
dd2 <- merge(dd2, dd[, c("node","x","y")], by.x="parent", by.y="node", all.x=TRUE)
colnames(dd2)[colnames(dd2)=="x"] <- "x_parent"
colnames(dd2)[colnames(dd2)=="y"] <- "y_parent"
dd2$x_mid <- (dd2$x_child + dd2$x_parent)/2
dd2$y_mid <- (dd2$y_child + dd2$y_parent)/2

# Keep only edges with any events
pts_gain <- dd2[dd2$gain > 0, , drop=FALSE]
pts_loss <- dd2[dd2$loss > 0, , drop=FALSE]

# Compose
p <- p + geom_tiplab(size=tip_size, align=FALSE)

# Gains ▲ and losses ▼; sizes reflect counts
if (nrow(pts_gain)) {
  p <- p + geom_point(data=pts_gain, aes(x=x_mid, y=y_mid, size=gain),
                      shape=24, fill=NA, stroke=0.8)
}
if (nrow(pts_loss)) {
  p <- p + geom_point(data=pts_loss, aes(x=x_mid, y=y_mid, size=loss),
                      shape=25, fill=NA, stroke=0.8)
}

p <- p + scale_size_continuous(name="Events", range=c(1.5,5)) +
         ggtitle(sprintf("PT phylogeny with MTPT gain/loss (%s)%s",
                 toupper(model),
                 if (chronos_on) " [chronos time-scale]" else ""))

# Output
ggsave(filename=outpdf, plot=p, width=width_in, height=height_in, units="in")
cat("Wrote:", outpdf, "\n")

# ---------------- References ----------------
# ape / chronos: Paradis & Schliep (2019) Bioinformatics 35:526–528.
# Fitch parsimony: Fitch (1971) Syst Biol 20:406–416.
# Mk likelihood: Lewis (2001) Syst Biol 50:913–925; ape::ace implementation.
# ggtree: Yu et al. (2017) Methods Ecol Evol 8:28–36.
