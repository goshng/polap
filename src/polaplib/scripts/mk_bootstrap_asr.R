# FILE: scripts/mk_bootstrap_asr.R
# VERSION: 0.1.0
# For B bootstrap trees, run stoch-mapping (nsim per tree), aggregate expected gains/losses.
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: mk_bootstrap_asr.R <presence.tsv> <ufboot.trees> <out.tsv> <B> <nsim> <model: auto|ER|ARD> <model_compare.tsv>")
}
presf <- args[1]
ufboot <- args[2]
outf <- args[3]
B <- as.integer(args[4])
nsim <- as.integer(args[5])
mode <- args[6]
mcfile <- args[7]

P <- fread(presf)
species <- setdiff(colnames(P), "CID")

trees <- read.tree(ufboot) # multiPhylo
if (length(trees) < 1) stop("No bootstrap trees read from ufboot.")
set.seed(1)
idx <- sample.int(length(trees), size = min(B, length(trees)), replace = FALSE)
trees <- trees[idx]

MC <- NULL
if (mode == "auto" && file.exists(mcfile)) {
  MC <- fread(mcfile)[, .(CID, best_model)]
}

count_events <- function(mm) {
  sims <- if (inherits(mm, "multiSimmap")) mm else list(mm)
  gains <- 0
  losses <- 0
  for (s in sims) {
    for (edge_maps in s$maps) {
      st <- names(edge_maps)
      if (length(st) >= 2) {
        for (i in 1:(length(st) - 1)) {
          if (st[i] == "0" && st[i + 1] == "1") gains <- gains + 1
          if (st[i] == "1" && st[i + 1] == "0") losses <- losses + 1
        }
      }
    }
  }
  c(gains = gains / length(sims), losses = losses / length(sims))
}

rows <- list()
for (cid in P$CID) {
  y <- as.integer(P[CID == cid, ..species])
  names(y) <- species
  st_ref <- factor(y, levels = c(0, 1), labels = c("0", "1"))

  mdl <- if (!is.null(MC)) MC[CID == cid, best_model] else NA_character_
  if (mode == "ER") mdl <- "ER"
  if (mode == "ARD") mdl <- "ARD"
  if (is.na(mdl)) mdl <- "ER"

  Eg <- c()
  El <- c()
  for (t in trees) {
    # Drop tips not in matrix
    keep <- t$tip.label %in% species
    if (!all(keep)) t <- drop.tip(t, t$tip.label[!keep])
    st <- st_ref[match(t$tip.label, names(st_ref))]
    st <- factor(st, levels = c("0", "1"))

    fit <- tryCatch(fitMk(t, st, model = mdl), error = function(e) NULL)
    if (is.null(fit)) next
    if (mdl == "ER") {
      q <- if (fit$rates[1] <= 0) 1e-6 else fit$rates[1]
      Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
      rownames(Q) <- colnames(Q) <- c("0", "1")
    } else {
      q01 <- if (fit$rates[1] <= 0) 1e-6 else fit$rates[1]
      q10 <- if (fit$rates[2] <= 0) 1e-6 else fit$rates[2]
      Q <- matrix(c(-q01, q01, q10, -q10), 2, 2, byrow = TRUE)
      rownames(Q) <- colnames(Q) <- c("0", "1")
    }
    sim <- make.simmap(t, st, Q = Q, nsim = nsim, pi = "estimated", model = "ARD")
    ev <- count_events(sim)
    Eg <- c(Eg, ev["gains"])
    El <- c(El, ev["losses"])
  }
  if (length(Eg) == 0) {
    rows[[length(rows) + 1]] <- data.table(
      CID = cid, Eg_med = NA_real_, Eg_lo = NA_real_, Eg_hi = NA_real_,
      El_med = NA_real_, El_lo = NA_real_, El_hi = NA_real_
    )
  } else {
    rows[[length(rows) + 1]] <- data.table(
      CID = cid,
      Eg_med = median(Eg), Eg_lo = quantile(Eg, 0.025), Eg_hi = quantile(Eg, 0.975),
      El_med = median(El), El_lo = quantile(El, 0.025), El_hi = quantile(El, 0.975)
    )
  }
}
out <- rbindlist(rows)
fwrite(out, file = outf, sep = "\t")
cat("Bootstrap ASR summary -> ", outf, "\n")
