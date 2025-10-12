#!/usr/bin/env python3
# Version: v0.7.0
"""
Compute per-repeat recombinant fraction f_hat = R_norm / (P_norm + R_norm),
with optional normalization by local depth around breakpoints and Wilson 95% CI from RAW counts.

Inputs
------
--repeats   repeats.tsv  (pair_id r_name r_s r_e q_name q_s q_e orient len pid)
--parent    junction_support_parent.tsv  (junction_id reads_support reads_left reads_right)
--recomb    junction_support_recomb.tsv  (same columns)
--depth     samtools depth -aa (contig pos depth) [optional but recommended]
--flank     int (bp window for averaging depth at four breakpoints)
--label     sample label
--out-tsv   repeat_fractions.tsv

Output columns
--------------
label pair_id orient P_raw R_raw P_norm R_norm f_hat CI_lo CI_hi flank_mean_depth
"""
from __future__ import annotations
import argparse, logging
from typing import Dict, Tuple
from collections import defaultdict
import math

__version__ = "0.7.0"


def load_support(path: str) -> Dict[str, int]:
    D = {}
    with open(path) as fh:
        hdr = fh.readline()
        for ln in fh:
            if not ln.strip():
                continue
            jid, both, _, _ = ln.rstrip("\n").split("\t")
            D[jid] = int(both)
    return D


def load_meta(meta_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Map junction_id -> (pair_id, type)
    We infer meta path by replacing 'support' with '' is fragile; better pass explicit meta if needed.
    Here assume .tsv sits next to support with known name 'junc_parent.tsv' or 'junc_recomb.tsv'.
    """
    M = {}
    with open(meta_path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            M[f[cols["junction_id"]]] = (f[cols["pair_id"]], f[cols["type"]])
    return M


def mean_flank_depth(
    depth_map: Dict[Tuple[str, int], int], contig: str, pos: int, flank: int
) -> float:
    if not depth_map:
        return 0.0
    s = max(1, pos - flank)
    e = pos + flank
    sm = n = 0
    for x in range(s, e + 1):
        sm += depth_map.get((contig, x), 0)
        n += 1
    return (sm / n) if n > 0 else 0.0


def wilson_ci(k: int, n: int, z: float = 1.96) -> Tuple[float, float, float]:
    if n <= 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) / n) + (z * z / (4 * n * n))) / denom
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return (p, lo, hi)


def main():
    ap = argparse.ArgumentParser(
        description="Compute per-repeat recombinant fraction with Wilson CI."
    )
    ap.add_argument("--repeats", required=True)
    ap.add_argument("--parent", required=True)
    ap.add_argument("--recomb", required=True)
    ap.add_argument(
        "--parent-meta", help="Optional explicit meta for parent (junc_parent.tsv)"
    )
    ap.add_argument(
        "--recomb-meta", help="Optional explicit meta for recomb (junc_recomb.tsv)"
    )
    ap.add_argument("--depth", help="samtools depth -aa output (optional)", default="")
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument(
        "-v",
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    ap.add_argument("--version", action="store_true")
    a = ap.parse_args()
    if a.version:
        print(__version__)
        return
    logging.basicConfig(
        level=getattr(logging, a.log_level), format="%(levelname)s: %(message)s"
    )

    # Load repeats (for orient + breakpoint coordinates)
    rep_info = {}
    with open(a.repeats) as fh:
        fh.readline()
        for ln in fh:
            if not ln.strip():
                continue
            pair, r, rs, re, q, qs, qe, orient, _, _ = ln.rstrip("\n").split("\t")
            rep_info[pair] = (r, int(rs), int(re), q, int(qs), int(qe), orient)

    P = load_support(a.parent)
    R = load_support(a.recomb)

    # Meta mapping
    if a.parent_meta:
        Mp = load_meta(a.parent_meta)
    else:
        pmeta = a.parent.replace("support", "").rstrip(".tsv") + "tsv"
        Mp = load_meta(pmeta)
    if a.recomb_meta:
        Mr = load_meta(a.recomb_meta)
    else:
        rmeta = a.recomb.replace("support", "").rstrip(".tsv") + "tsv"
        Mr = load_meta(rmeta)

    # Group counts by repeat pair
    Praw = defaultdict(int)
    Rraw = defaultdict(int)
    for jid, c in P.items():
        pair, _ = Mp.get(jid, (None, None))
        if pair:
            Praw[pair] += c
    for jid, c in R.items():
        pair, _ = Mr.get(jid, (None, None))
        if pair:
            Rraw[pair] += c

    # Depth map
    depth_map = {}
    if a.depth:
        with open(a.depth) as fh:
            for ln in fh:
                c, p, d = ln.rstrip("\n").split("\t")
                depth_map[(c, int(p))] = int(d)

    with open(a.out_tsv, "w") as oh:
        oh.write(
            "label\tpair_id\torient\tP_raw\tR_raw\tP_norm\tR_norm\tf_hat\tCI_lo\tCI_hi\tflank_mean_depth\n"
        )
        for pair, (r, rs, re, q, qs, qe, orient) in rep_info.items():
            pcount = Praw.get(pair, 0)
            rcount = Rraw.get(pair, 0)
            # normalization depth (mean of four breakpoint windows)
            if depth_map:
                d_rl = mean_flank_depth(depth_map, r, min(rs, re), a.flank)
                d_rr = mean_flank_depth(depth_map, r, max(rs, re), a.flank)
                d_ql = mean_flank_depth(depth_map, q, min(qs, qe), a.flank)
                d_qr = mean_flank_depth(depth_map, q, max(qs, qe), a.flank)
                d_mean = max(1e-6, (d_rl + d_rr + d_ql + d_qr) / 4.0)
            else:
                d_mean = 0.0
            Pn = (pcount / d_mean) if d_mean > 0 else float(pcount)
            Rn = (rcount / d_mean) if d_mean > 0 else float(rcount)
            denom = Pn + Rn
            fhat = (Rn / denom) if denom > 0 else 0.0
            p, lo, hi = wilson_ci(rcount, (pcount + rcount))
            oh.write(
                f"{a.label}\t{pair}\t{orient}\t{pcount}\t{rcount}\t{Pn:.6f}\t{Rn:.6f}\t{fhat:.6f}\t{lo:.6f}\t{hi:.6f}\t{d_mean:.3f}\n"
            )


if __name__ == "__main__":
    main()
