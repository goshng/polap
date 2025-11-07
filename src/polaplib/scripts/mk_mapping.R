# FILE: scripts/mk_mapping.R
#!/usr/bin/env Rscript
# VERSION: 0.3.2
# SPDX-License-Identifier: GPL-3.0-or-later
# Usage: Rscript mk_mapping.R <treefile> <presence.tsv> <outdir> <model>

suppressPackageStartupMessages({
  library(data.table)
  library(ape)        # read.tree, ace
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) stop("Usage: mk_mapping.R <treefile> <presence.tsv> <outdir> <model>")
treef <- args[1]; presf <- args[2]; outd <- args[3]; model <- tolower(args[4])
dir.create(outd, showWarnings=FALSE, recursive=TRUE)

tr0 <- read.tree(treef)
if (is.null(tr0) || length(tr0$tip.label) == 0) stop("Empty/invalid tree: ", treef)

P <- fread(presf)
if (!"CID" %in% names(P)) stop("presence.tsv must have a 'CID' column.")
species_all <- setdiff(colnames(P), "CID")
if (!length(species_all)) stop("presence.tsv has no species columns.")

# Match tips
tips_in <- intersect(tr0$tip.label, species_all)
if (!length(tips_in)) {
  stop("No overlapping species between tree tips and presence.tsv columns.\n",
       "  tree tips (n=", length(tr0$tip.label), "): ",
       paste(head(tr0$tip.label, 10), collapse=", "), " ...\n",
       "  presence cols (n=", length(species_all), "): ",
       paste(head(species_all, 10), collapse=", "), " ...")
}
tr <- if (length(tips_in) < length(tr0$tip.label)) keep.tip(tr0, tips_in) else tr0
species <- tr$tip.label

# Build tree scaffolding
Ntip <- length(tr$tip.label); Nnode <- tr$Nnode; Ntot <- Ntip + Nnode
children <- vector("list", Ntot); parent <- integer(Ntot); parent[] <- NA_integer_
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
    stack_v <- stack_v[-length(stack_v)]
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
  gains <- 0L; losses <- 0L
  for (iE in seq_len(nrow(tr$edge))) {
    a <- tr$edge[iE,1]; b <- tr$edge[iE,2]
    if (assign[a]==0L && assign[b]==1L) gains  <- gains + 1L
    if (assign[a]==1L && assign[b]==0L) losses <- losses + 1L
  }
  list(gains=as.integer(gains), losses=as.integer(losses))
}

# Presence matrix for matched species only (0/1)
Pres <- as.matrix(P[, ..species]); mode(Pres) <- "integer"
Pres[is.na(Pres)] <- 0L; Pres[Pres != 0L] <- 1L; rownames(Pres) <- P$CID

events <- vector("list", nrow(P))
for (i in seq_len(nrow(P))) {
  cid <- P$CID[i]
  y <- Pres[cid, ]; names(y) <- species
  rows <- list()

  if (model %in% c("fitch","all")) {
    gl <- fitch_assign(y, dollo_bias=FALSE)
    rows[[length(rows)+1]] <- data.table(CID=cid, model="fitch", gains=gl$gains, losses=gl$losses)
  }
  if (model %in% c("dollo","all")) {
    gl <- fitch_assign(y, dollo_bias=TRUE)
    rows[[length(rows)+1]] <- data.table(CID=cid, model="dollo_like", gains=gl$gains, losses=gl$losses)
  }
  if (model %in% c("mk","all")) {
    # APE's ace for discrete ER Mk
    # ace() expects named vector of states for tips
    lab <- species
    vec <- as.integer(y); names(vec) <- lab
    ok <- TRUE; mk_g <- NA_integer_; mk_l <- NA_integer_
    suppressWarnings({
      fit <- try(ace(vec, tr, type="discrete", model="ER"), silent=TRUE)
      if (inherits(fit, "try-error")) ok <- FALSE
    })
    if (ok && !is.null(fit$lik.anc)) {
      # choose state 1 if P(state=1) >= 0.5
      assign <- integer(Ntot); assign[] <- NA_integer_
      assign[seq_len(Ntip)] <- vec
      for (k in seq_len(Nnode)) {
        nd <- Ntip + k
        # fit$lik.anc rows correspond to internal nodes in node order:
        # internal node numbers are (Ntip+1):(Ntip+Nnode)
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

  events[[i]] <- data.table::rbindlist(rows, fill=TRUE)
}

EVENTS <- data.table::rbindlist(events, fill=TRUE)
fwrite(P, file=file.path(outd, "presence.tsv"), sep="\t")
fwrite(EVENTS, file=file.path(outd, "events.tsv"), sep="\t")
fwrite(data.table(CID=P$CID, mid_time=NA_real_), file=file.path(outd, "timeline.tsv"), sep="\t")
cat("ASR done: presence.tsv, events.tsv, timeline.tsv\n")

