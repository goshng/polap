#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-select-mt-by-x-and-kmer.py  v0.2.0

Select mitochondrial reads (assembly-free) by combining:
  (A) Coverage proxy X = log10(median syncmer multiplicity) bands learned from mt (and optional pt) anchors
  (B) Supervised k-mer composition classifier (Fisher LDA) trained from anchors

Inputs
  --tsv       syncfilter quickview TSV (columns include: read_id, median, ...)
  --reads     FASTQ(.gz) with matching read IDs (first token of header)
  --mt        mt anchor IDs (txt; one id per line)  [required]
  --pt        pt anchor IDs (txt; optional)
  --prefix    output prefix                         [required]

Key options
  --x-tail FLOAT            one-sided tail for band width (default 0.10)
  --sigma-floor FLOAT       min sigma for band (default 0.03)
  --kmer-order {4,6}        k-mer length for composition (default 4)
  --use-hpc                 homopolymer-compress reads before k-mers (for ONT)
  --min-anchor INT          min anchors per class to train LDA (default 25)
  --mix-mode {any,both,weighted}  combine X & LDA decisions (default any)
  --lda-th FLOAT            LDA score threshold after min-max rescale (0..1, default 0.5)
  --lda-strong FLOAT        strong-confidence cutoff for weighted mode (default 0.75)
  --version                 print version and exit

Outputs
  <prefix>.mt.ids
  <prefix>.pt.ids
  <prefix>.nuclear.ids
  <prefix>.bandstats.tsv     (μ/σ & band ranges learned)
  <prefix>.kmer_pca.png      (2D PCA of k-mers; anchors highlighted)
  <prefix>.x_hist.png        (X histogram + bands)

Notes
- LDA is trained: positive = mt anchors; negative = nuclear-like by X (outside mt/pt bands).
- If anchors are scarce, the script falls back to X-only selection (warns).
"""

import sys, os, argparse, csv, gzip, math, collections
from typing import Dict, List, Tuple, Optional

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

VERSION = "0.2.0"


# --------------------------- utils ------------------------------------------
def open_guess(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "rt")


def norm_id(s: str) -> str:
    return (s or "").split()[0]


def read_ids(path: Optional[str]) -> set:
    S = set()
    if not path:
        return S
    with open(path) as fh:
        for line in fh:
            t = line.strip().split()
            if t:
                S.add(t[0])
    return S


def load_X_from_quickview(tsv: str) -> Dict[str, float]:
    X = {}
    with open(tsv, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rid = norm_id(row["read_id"])
            try:
                med = float(row["median"])
                if med > 0:
                    X[rid] = math.log10(med)
            except:
                pass
    if not X:
        print("[error] no usable rows/medians in TSV", file=sys.stderr)
        sys.exit(2)
    return X


def z_for_tail(t):
    # one-sided z
    if t <= 0.01:
        return 2.3263478740408408
    if t <= 0.02:
        return 2.0537489106318225
    if t <= 0.05:
        return 1.6448536269514722
    if t <= 0.10:
        return 1.2815515655446004
    return 1.6448536269514722


def sigma_mad(x: np.ndarray) -> Optional[float]:
    if x.size == 0:
        return None
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return float(np.std(x)) if x.size > 1 else None
    return float(mad / 0.6744897501960817)


def learn_band(Xmap: Dict[str, float], anchors: set, tail=0.10, sigma_floor=0.03):
    xs = np.array([Xmap[a] for a in anchors if a in Xmap], dtype=float)
    if xs.size < 5:
        return None
    mu = float(np.median(xs))
    sg = float(np.std(xs)) if xs.size > 1 else None
    if sg is None or not np.isfinite(sg) or sg <= 0:
        sg = sigma_mad(xs)
    if sg is None or not np.isfinite(sg) or sg <= 0:
        sg = sigma_floor
    else:
        sg = max(sg, sigma_floor)
    z = z_for_tail(tail)
    return dict(mu=mu, sigma=sg, range=(mu - z * sg, mu + z * sg), n=int(xs.size))


# --------------------------- k-mers ------------------------------------------
def hpc_compress(seq: str) -> str:
    out = []
    prev = ""
    for c in seq.upper():
        if c not in "ACGT":
            prev = ""
            continue
        if c != prev:
            out.append(c)
            prev = c
    return "".join(out)


def kmer_index(order: int) -> Dict[str, int]:
    m = {}

    def rec(prefix, depth, idx):
        if depth == order:
            m[prefix] = idx[0]
            idx[0] += 1
            return
        for ch in "ACGT":
            rec(prefix + ch, depth + 1, idx)

    rec("", 0, [0])
    return m


def kmer_vec(seq: str, order: int, table: Dict[str, int]) -> np.ndarray:
    v = np.zeros(len(table), dtype=np.float32)
    n = len(seq)
    if n < order:
        return v
    for i in range(n - order + 1):
        k = seq[i : i + order]
        if k in table:
            v[table[k]] += 1.0
    s = v.sum()
    if s > 0:
        v /= s
    return v


def fastq_iter(path: str):
    with open_guess(path) as fh:
        while True:
            n = fh.readline()
            if not n:
                break
            s = fh.readline()
            _ = fh.readline()
            q = fh.readline()
            if not q:
                break
            rid = norm_id(n[1:].strip())
            yield rid, s.strip()


def build_kmer_matrix(reads: str, idset: set, order=4, use_hpc=False, max_reads=None):
    tbl = kmer_index(order)
    D = len(tbl)
    X = []
    names = []
    for rid, seq in fastq_iter(reads):
        if rid in idset:
            if use_hpc:
                seq = hpc_compress(seq)
            X.append(kmer_vec(seq, order, tbl))
            names.append(rid)
            if max_reads and len(names) >= max_reads:
                break
    if not X:
        return None, None
    return np.vstack(X), names


# --------------------------- LDA ---------------------------------------------
def lda_train(X_pos: np.ndarray, X_neg: np.ndarray, shrink=1e-6):
    """
    Fisher LDA (two-class). Returns (w,b) for score = w·x + b, plus class means.
    """
    mu1 = X_pos.mean(0)
    mu0 = X_neg.mean(0)
    C1 = np.cov(X_pos.T, bias=False) if X_pos.shape[0] > 1 else np.eye(X_pos.shape[1])
    C0 = np.cov(X_neg.T, bias=False) if X_neg.shape[0] > 1 else np.eye(X_neg.shape[1])
    C = (C1 + C0) / 2.0 + np.eye(X_pos.shape[1]) * shrink
    w = np.linalg.solve(C, (mu1 - mu0))
    b = -0.5 * (mu1 + mu0) @ w
    return w, b, mu1, mu0


def lda_score(X: np.ndarray, w: np.ndarray, b: float) -> np.ndarray:
    return X @ w + b


def minmax_rescale(y: np.ndarray) -> np.ndarray:
    lo, hi = np.min(y), np.max(y)
    if hi <= lo:
        return np.zeros_like(y)
    return (y - lo) / (hi - lo)


# --------------------------- plotting ----------------------------------------
def pca_2d(X: np.ndarray, k=2) -> np.ndarray:
    Xc = X - X.mean(0, keepdims=True)
    U, S, VT = np.linalg.svd(Xc, full_matrices=False)
    return Xc @ VT[:k].T


def plot_kmer_pca(X_all, names_all, mt_anchors, pt_anchors, y_pred, path_png):
    Z = pca_2d(X_all, 2)
    plt.figure(figsize=(6.4, 5.0))
    cmap = {
        "mt": "tab:orange",
        "pt": "tab:green",
        "nuclear": "tab:blue",
        "ambig": "gray",
    }
    # background by predicted class
    for cls in ("nuclear", "pt", "mt", "ambig"):
        idx = [i for i, c in enumerate(y_pred) if c == cls]
        if idx:
            plt.scatter(
                Z[idx, 0], Z[idx, 1], s=10, alpha=0.25, label=f"pred {cls}", c=cmap[cls]
            )
    # anchors on top
    idx_mt = [i for i, n in enumerate(names_all) if n in mt_anchors]
    if idx_mt:
        plt.scatter(
            Z[idx_mt, 0], Z[idx_mt, 1], s=20, c="black", marker="x", label="mt anchors"
        )
    idx_pt = [i for i, n in enumerate(names_all) if n in pt_anchors]
    if idx_pt:
        plt.scatter(
            Z[idx_pt, 0], Z[idx_pt, 1], s=20, c="black", marker="+", label="pt anchors"
        )
    plt.legend(loc="best", fontsize=8)
    plt.title("k-mer PCA (2D)")
    plt.tight_layout()
    plt.savefig(path_png, dpi=140)
    plt.close()


def plot_x_hist(Xvals, band_mt, band_pt, out_png):
    plt.figure(figsize=(7.2, 4.0))
    plt.hist(Xvals, bins=140, alpha=0.35, label="all reads (X)")
    if band_mt:
        a, b = band_mt["range"]
        mu = band_mt["mu"]
        plt.axvline(a, color="tab:orange", ls=":", lw=1.6)
        plt.axvline(b, color="tab:orange", ls=":", lw=1.6)
        plt.axvline(mu, color="tab:orange", ls="--", lw=1.6, label="mt band & μ")
    if band_pt:
        a, b = band_pt["range"]
        mu = band_pt["mu"]
        plt.axvline(a, color="tab:green", ls=":", lw=1.6)
        plt.axvline(b, color="tab:green", ls=":", lw=1.6)
        plt.axvline(mu, color="tab:green", ls="--", lw=1.6, label="pt band & μ")
    plt.xlabel("X = log10(median syncmer multiplicity)")
    plt.ylabel("reads")
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=140)
    plt.close()


# --------------------------- main -------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Select mt reads by coverage (X) + LDA on k-mer composition."
    )
    ap.add_argument("--tsv", required=True, help="syncfilter quickview TSV")
    ap.add_argument("--reads", required=True, help="FASTQ(.gz)")
    ap.add_argument("--mt", required=True, help="mt anchor ID list (txt)")
    ap.add_argument("--pt", help="pt anchor ID list (txt)")
    ap.add_argument("--prefix", required=True, help="output prefix")

    ap.add_argument("--x-tail", type=float, default=0.10)
    ap.add_argument("--sigma-floor", type=float, default=0.03)
    ap.add_argument("--kmer-order", type=int, choices=[4, 6], default=4)
    ap.add_argument("--use-hpc", action="store_true")
    ap.add_argument("--min-anchor", type=int, default=25)
    ap.add_argument("--mix-mode", choices=["any", "both", "weighted"], default="any")
    ap.add_argument("--lda-th", type=float, default=0.50)
    ap.add_argument("--lda-strong", type=float, default=0.75)
    ap.add_argument("--version", action="store_true")
    args = ap.parse_args()

    if args.version:
        print(VERSION)
        sys.exit(0)

    os.makedirs(os.path.dirname(args.prefix) or ".", exist_ok=True)

    # 1) Load X and anchors
    Xmap = load_X_from_quickview(args.tsv)
    mt_anchor = read_ids(args.mt)
    pt_anchor = read_ids(args.pt) if args.pt else set()

    # 2) Learn bands from anchors
    band_mt = learn_band(
        Xmap, mt_anchor, tail=args.x_tail, sigma_floor=args.sigma_floor
    )
    if not band_mt:
        print("[error] not enough mt anchors with X; need >=5", file=sys.stderr)
        sys.exit(2)
    band_pt = (
        learn_band(Xmap, pt_anchor, tail=args.x_tail, sigma_floor=args.sigma_floor)
        if pt_anchor
        else None
    )

    # Save band stats
    with open(f"{args.prefix}.bandstats.tsv", "w") as w:
        w.write("type\tn\tmu_log10\tsigma_log10\tf1_log10\tf2_log10\n")
        w.write(
            f"mt\t{band_mt['n']}\t{band_mt['mu']:.5f}\t{band_mt['sigma']:.5f}\t{band_mt['range'][0]:.5f}\t{band_mt['range'][1]:.5f}\n"
        )
        if band_pt:
            w.write(
                f"pt\t{band_pt['n']}\t{band_pt['mu']:.5f}\t{band_pt['sigma']:.5f}\t{band_pt['range'][0]:.5f}\t{band_pt['range'][1]:.5f}\n"
            )

    # 3) X-only provisional label
    def in_band(x, band):
        a, b = band["range"]
        return (x >= a) and (x <= b)

    cls_X = {}
    for rid, x in Xmap.items():
        if band_mt and in_band(x, band_mt):
            cls = "mt"
        elif band_pt and in_band(x, band_pt):
            cls = "pt"
        else:
            cls = "nuclear"
        cls_X[rid] = cls

    # 4) Build k-mer features for all reads having X
    all_ids = set(Xmap.keys())
    Kmat, Knames = build_kmer_matrix(
        args.reads, all_ids, order=args.kmer_order, use_hpc=args.use_hpc
    )
    if Kmat is None:
        print("[warn] no k-mer features built; falling back to X-only", file=sys.stderr)
        use_kmer = False
    else:
        use_kmer = True

    # 5) Train LDA (mt vs nuclear-like) using anchors
    y_pred_kmer = None
    kscore_map = {}
    if use_kmer:
        name_to_idx = {n: i for i, n in enumerate(Knames)}
        pos = [name_to_idx[n] for n in mt_anchor if n in name_to_idx]
        neg = [
            name_to_idx[n]
            for n in Knames
            if (n not in mt_anchor) and (cls_X.get(n, "nuclear") == "nuclear")
        ]
        if len(pos) >= args.min_anchor and len(neg) >= args.min_anchor:
            Xpos = Kmat[pos, :]
            Xneg = Kmat[neg, :]
            w, b, mu1, mu0 = lda_train(Xpos, Xneg, shrink=1e-6)
            raw = lda_score(Kmat, w, b)
            prob = minmax_rescale(raw)  # 0..1
            y_pred_kmer = np.where(prob >= args.lda_th, "mt", "nuclear")
            kscore_map = {Knames[i]: float(prob[i]) for i in range(len(Knames))}
            # PCA plot
            try:
                plot_kmer_pca(
                    Kmat,
                    Knames,
                    mt_anchor,
                    pt_anchor,
                    y_pred_kmer.tolist(),
                    f"{args.prefix}.kmer_pca.png",
                )
            except Exception as e:
                print(f"[warn] PCA plot failed: {e}", file=sys.stderr)
        else:
            print(
                f"[warn] insufficient anchors for LDA (mt={len(pos)}, neg={len(neg)}); k-mer disabled",
                file=sys.stderr,
            )
            use_kmer = False

    # X histogram plot
    try:
        plot_x_hist(
            np.array(list(Xmap.values())), band_mt, band_pt, f"{args.prefix}.x_hist.png"
        )
    except Exception as e:
        print(f"[warn] X histogram plot failed: {e}", file=sys.stderr)

    # 6) Combine decisions
    mt_ids, pt_ids, nuc_ids = [], [], []
    for rid, x in Xmap.items():
        # X decision
        x_cls = cls_X[rid]
        # k-mer decision
        km_cls = None
        if use_kmer and rid in kscore_map:
            km_cls = "mt" if kscore_map[rid] >= args.lda_th else "nuclear"

        # mixing
        if args.mix_mode == "any":
            is_mt = (x_cls == "mt") or (km_cls == "mt")
        elif args.mix_mode == "both":
            is_mt = (x_cls == "mt") and (km_cls == "mt")
        else:  # weighted
            sc = kscore_map.get(rid, None)
            if sc is not None and sc >= args.lda_strong:
                is_mt = True
            elif sc is not None and sc <= (1.0 - args.lda_strong):
                is_mt = False
            else:
                is_mt = x_cls == "mt"

        if is_mt:
            mt_ids.append(rid)
        elif band_pt and in_band(x, band_pt):
            pt_ids.append(rid)
        else:
            nuc_ids.append(rid)

    # 7) Write outputs
    with open(f"{args.prefix}.mt.ids", "w") as w:
        for rid in mt_ids:
            w.write(rid + "\n")
    with open(f"{args.prefix}.pt.ids", "w") as w:
        for rid in pt_ids:
            w.write(rid + "\n")
    with open(f"{args.prefix}.nuclear.ids", "w") as w:
        for rid in nuc_ids:
            w.write(rid + "\n")

    print(f"[done] mt={len(mt_ids)}  pt={len(pt_ids)}  nuclear={len(nuc_ids)}")
    print(
        f"[out] {args.prefix}.mt.ids  {args.prefix}.pt.ids  {args.prefix}.nuclear.ids"
    )
    if use_kmer:
        print(
            f"[info] k-mer LDA used (order={args.kmer_order}, th={args.lda_th}, mix={args.mix_mode})"
        )
    else:
        print("[info] k-mer LDA disabled; X-only selection")


if __name__ == "__main__":
    main()
