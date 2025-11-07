# FILE: scripts/mk_model_compare.R
# VERSION: 0.1.0
# Fits Mk models (ER & ARD) using phytools::fitMk, compares by AICc,
# stochastically maps to estimate expected gains/losses, and flags "Dollo-like".
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: mk_model_compare.R <treefile> <presence.tsv> <out.tsv> <nsim> <qratio_thr>")
}
treef <- args[1]
presf <- args[2]
outf <- args[3]
nsim <- as.integer(args[4])
qthr <- as.numeric(args[5])

tr <- read.tree(treef)
P <- fread(presf)
stopifnot("CID" %in% names(P))
species <- setdiff(colnames(P), "CID")

# Drop tips not in matrix
keep <- tr$tip.label %in% species
if (!all(keep)) tr <- drop.tip(tr, tr$tip.label[!keep])
species <- tr$tip.label
n <- length(species)

aicc <- function(logL, k, n) {
  2 * k - 2 * logL + (2 * k * (k + 1)) / (n - k - 1)
}

count_events <- function(mm) {
  # mm: multiSimmap or single simmap; return expected gains (0->1) and losses (1->0)
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

res <- list()
for (i in seq_len(nrow(P))) {
  cid <- P$CID[i]
  y <- as.integer(P[i, ..species])
  names(y) <- species
  # states as factor labels "0"/"1"
  st <- factor(y, levels = c(0, 1), labels = c("0", "1"))

  # Fit ER & ARD
  fit_er <- tryCatch(fitMk(tr, st, model = "ER"), error = function(e) NULL)
  fit_ard <- tryCatch(fitMk(tr, st, model = "ARD"), error = function(e) NULL)

  if (is.null(fit_er) && is.null(fit_ard)) {
    res[[length(res) + 1]] <- data.table(
      CID = cid, ll_ER = NA, AICc_ER = NA, q01_ER = NA, q10_ER = NA,
      ll_ARD = NA, AICc_ARD = NA, q01_ARD = NA, q10_ARD = NA,
      best_model = NA, gain_rate_ratio = NA, E_gains = NA, E_losses = NA,
      dollo_like = NA
    )
    next
  }

  # Extract params
  ll_ER <- if (!is.null(fit_er)) fit_er$logLik else NA_real_
  q_ER <- if (!is.null(fit_er)) fit_er$rates[1] else NA_real_
  AICc_ER <- if (!is.null(fit_er)) aicc(ll_ER, k = 1, n = n) else NA_real_

  ll_ARD <- if (!is.null(fit_ard)) fit_ard$logLik else NA_real_
  q01_ARD <- if (!is.null(fit_ard)) fit_ard$rates[1] else NA_real_
  q10_ARD <- if (!is.null(fit_ard)) fit_ard$rates[2] else NA_real_
  AICc_ARD <- if (!is.null(fit_ard)) aicc(ll_ARD, k = 2, n = n) else NA_real_

  # Choose best by AICc (lower better)
  best <- if (is.na(AICc_ER)) "ARD" else if (is.na(AICc_ARD)) "ER" else if (AICc_ARD < AICc_ER) "ARD" else "ER"

  # Build Q for best model
  if (best == "ER") {
    q <- if (is.na(q_ER) || q_ER <= 0) 1e-6 else q_ER
    Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
    rownames(Q) <- colnames(Q) <- c("0", "1")
  } else {
    q01 <- if (is.na(q01_ARD) || q01_ARD <= 0) 1e-6 else q01_ARD
    q10 <- if (is.na(q10_ARD) || q10_ARD <= 0) 1e-6 else q10_ARD
    Q <- matrix(c(-q01, q01, q10, -q10), 2, 2, byrow = TRUE)
    rownames(Q) <- colnames(Q) <- c("0", "1")
  }

  # Stochastic mapping (expected gains/losses)
  sim <- make.simmap(tr, st, Q = Q, nsim = nsim, pi = "estimated", model = "ARD")
  ev <- count_events(sim)
  gr <- if (!is.null(fit_ard)) q01_ARD / q10_ARD else NA_real_
  dollo_like <- if (!is.na(gr) && gr < qthr) "yes" else "no"

  res[[length(res) + 1]] <- data.table(
    CID = cid,
    ll_ER = ll_ER, AICc_ER = AICc_ER, q01_ER = q_ER, q10_ER = q_ER,
    ll_ARD = ll_ARD, AICc_ARD = AICc_ARD, q01_ARD = q01_ARD, q10_ARD = q10_ARD,
    best_model = best, gain_rate_ratio = gr,
    E_gains = ev["gains"], E_losses = ev["losses"],
    dollo_like = dollo_like
  )
}
out <- rbindlist(res, fill = TRUE)
fwrite(out, file = outf, sep = "\t")
cat("Model comparison table written to:", outf, "\n")
