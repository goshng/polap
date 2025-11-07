# FILE: scripts/plot_mtpt_timescaled.R
#!/usr/bin/env Rscript
# Version: 0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Draw a (time-scaled) PT phylogeny with MTPT gain/loss events per edge.
# - Reads PT tree (Newick) and presence.tsv (CID x species).
# - Computes per-edge gains/losses across all clusters with:
#     * Fitch parsimony (unordered)           --model fitch
#     * Dollo-like parsimony (gain-averse)    --model dollo
#     * Likelihood Mk(ER) via ape::ace        --model mk_ER
# - Optionally time-scales the tree with ape::chronos (—chronos).
#
# Dependencies: data.table, ape, ggtree, ggplot2
#
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)   # bioconductor-ggtree
})

# ---------- CLI ----------
args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, def=NULL, logical=FALSE) {
  i <- which(args==flag)
  if (!length(i)) return(def)
  if (logical) return(TRUE)
  if (i==length(args)) return(def)
  args[i+1]
}
treef   <- getArg("--tree")
presf   <- getArg("--presence")
outpdf  <- getArg("--out")
model   <- tolower(getArg("--model","fitch"))
chronos_on <- isTRUE(getArg("--chronos", logical=TRUE))
lambda  <- as.numeric(getArg("--lambda","1.0"))
tip_size<- as.numeric(getArg("--tip-size","2.5"))

if (is.null(treef) || is.null(presf) || is.null(outpdf))
  stop("Usage: --tree <file> --presence <file> --out <pdf> [--model fitch|dollo|mk_ER] [--chronos] [--lambda 1.0] [--tip-size 2.5]")

# ---------- IO ----------
tr0 <- read.tree(treef)
if (is.null(tr0) || length(tr0$tip.label)==0) stop("Invalid tree: ", treef)

P <- fread(presf) # CID + species cols
stopifnot("CID" %in% colnames(P))
species_all <- setdiff(colnames(P), "CID")
stopifnot(length(species_all)>0)

# Match species
tips_in <- intersect(tr0$tip.label, species_all)
if (!length(tips_in)) {
  stop("No overlapping species between tree tips and presence.tsv columns.\n",
       "tips (n=", length(tr0$tip.label), "): ", paste(head(tr0$tip.label,10),collapse=", "), "\n",
       "presence cols: ", paste(head(species_all,10),collapse=", "))
}
tr <- if (length(tips_in)<length(tr0$tip.label)) keep.tip(tr0, tips_in) else tr0
species <- tr$tip.label

# Optional time scaling
if (chronos_on) {
  # quiet chronos; if it fails, keep the original branch lengths
  tr <- tryCatch(chronos(tr, lambda=lambda, quiet=TRUE), error=function(e) { message("[chronos] failed; using original lengths."); tr })
}

# ---------- Tree scaffolding ----------
Ntip <- length(tr$tip.label); Nnode <- tr$Nnode; Ntot <- Ntip + Nnode
children <- vector("list", Ntot)
parent   <- integer(Ntot); parent[] <- NA_integer_
for (i in seq_len(nrow(tr$edge))) {
  a <- tr$edge[i,1]; b <- tr$edge[i,2]
  parent[b] <- a
  children[[a]] <- c(children[[a]], b)
}
for (v in seq_len(Ntot)) if (is.null(children[[v]])) children[[v]] <- integer(0)
root <- setdiff(tr$edge[,1], tr$edge[,2])[1]

# Robust iterative postorder
postorder <- function(root, children) {
  ord <- integer(0)
  stack_v <- c(root); stack_vis <- c(FALSE)
  while (length(stack_v)) {
    v <- stack_v[[length(stack_v)]]
    vis <- stack_vis[[length(stack_vis)]]
    stack_v   <- stack_v[-length(stack_v)]
    stack_vis <- stack_vis[-length(stack_vis)]
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

fitch_assign <- function(y01, dollo_bias=FALSE) {
  S <- integer(Ntot); S[] <- 0L
  S[seq_len(Ntip)] <- enc_tip(y01[species])
  # postorder sets
  for (v in PO[PO > Ntip]) {
    ch <- children[[v]]; if (!length(ch)) next
    inter_val <- S[ch[1L]]
    for (u in ch[-1L]) inter_val <- bitwAnd(inter_val, S[u])
    if (inter_val != 0L) {
      S[v] <- inter_val
    } else {
      uni_val <- S[ch[1L]]
      for (u in ch[-1L]) uni_val <- bitwOr(uni_val, S[u])
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

mk_assign <- function(y01) {
  # Likelihood Mk(ER) via ape::ace
  vec <- as.integer(y01); names(vec) <- species
  fit <- suppressWarnings(try(ace(vec, tr, type="discrete", model="ER"), silent=TRUE))
  if (inherits(fit, "try-error") || is.null(fit$lik.anc)) return(NULL)
  assign <- integer(Ntot); assign[] <- NA_integer_
  assign[seq_len(Ntip)] <- vec
  for (k in seq_len(Nnode)) {
    nd <- Ntip + k
    p1 <- fit$lik.anc[k, "1"]
    assign[nd] <- as.integer(p1 >= 0.5)
  }
  assign
}

# Presence matrix (matched species order)
Pres <- as.matrix(P[, ..species]); mode(Pres) <- "integer"
Pres[is.na(Pres)] <- 0L; Pres[Pres != 0L] <- 1L
rownames(Pres) <- P$CID

# ---------- Per-edge event aggregation ----------
edge_gains <- integer(nrow(tr$edge)); edge_losses <- integer(nrow(tr$edge))

for (i in seq_len(nrow(P))) {
  cid <- P$CID[i]
  y <- Pres[cid, ]; names(y) <- species

  assign <- switch(model,
    "fitch" = fitch_assign(y, dollo_bias=FALSE),
    "dollo" = fitch_assign(y, dollo_bias=TRUE),
    "mk_er" = mk_assign(y),
    fitch_assign(y, dollo_bias=FALSE)
  )
  if (is.null(assign)) next

  for (e in seq_len(nrow(tr$edge))) {
    a <- tr$edge[e,1]; b <- tr$edge[e,2]
    if (assign[a]==0L && assign[b]==1L) edge_gains[e]  <- edge_gains[e]  + 1L
    if (assign[a]==1L && assign[b]==0L) edge_losses[e] <- edge_losses[e] + 1L
  }
}

# ---------- Build plotting data (edge midpoints) ----------
# ggtree gives node coordinates; place bubbles at edge midpoints
p <- ggtree(tr, size=0.4) + theme_tree2()
dd <- p$data  # contains node positions (x,y), parent

# Map edges to parent/child node rows
edge_df <- data.frame(tr$edge)
colnames(edge_df) <- c("parent","node")
edge_df$gain  <- edge_gains
edge_df$loss  <- edge_losses

# merge coordinates
dd2 <- merge(edge_df, dd[, c("node","x","y")], by="node", all.x=TRUE)
colnames(dd2)[colnames(dd2)=="x"] <- "x_child"
colnames(dd2)[colnames(dd2)=="y"] <- "y_child"
dd2 <- merge(dd2, dd[, c("node","x","y")], by.x="parent", by.y="node", all.x=TRUE)
colnames(dd2)[colnames(dd2)=="x"] <- "x_parent"
colnames(dd2)[colnames(dd2)=="y"] <- "y_parent"
dd2$x_mid <- (dd2$x_child + dd2$x_parent)/2
dd2$y_mid <- (dd2$y_child + dd2$y_parent)/2

# Keep only edges with any events
pts_gain <- dd2[dd2$gain  > 0, , drop=FALSE]
pts_loss <- dd2[dd2$loss  > 0, , drop=FALSE]

# ---------- Compose plot ----------
# Tip labels
p <- p + geom_tiplab(size=tip_size, align=FALSE)

# Gains (▲) and losses (▼) sized by counts
if (nrow(pts_gain)) {
  p <- p + geom_point(data=pts_gain, aes(x=x_mid, y=y_mid, size=gain), shape=24, fill=NA, stroke=0.8)
}
if (nrow(pts_loss)) {
  p <- p + geom_point(data=pts_loss, aes(x=x_mid, y=y_mid, size=loss), shape=25, fill=NA, stroke=0.8)
}

# A minimalist legend title
p <- p + scale_size_continuous(name="Events", range=c(1.5,5)) +
         ggtitle(sprintf("PT phylogeny with MTPT gain/loss (%s)%s",
                         toupper(model),
                         if (chronos_on) " [chronos time-scale]" else ""))

ggsave(filename=outpdf, plot=p, width=11, height=8.5, units="in")
cat("Wrote:", outpdf, "\n")
