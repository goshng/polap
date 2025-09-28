#!/usr/bin/env python3
# polap-py-overlapness-from-paf.py  v0.0.2
# From a PAF (file or stdin), compute per-read degree and weighted-degree (wdeg).
# Edge weight: ident * (alen / min(qlen, tlen)).
import sys, gzip, argparse
from collections import defaultdict


def openg(p, m="rt"):
    if p == "-":
        return sys.stdin
    return gzip.open(p, m) if p.endswith(".gz") else open(p, m)


ap = argparse.ArgumentParser(description="Compute degree/wdegree from a PAF")
ap.add_argument("paf", help="PAF file (.gz) or '-' for stdin")
ap.add_argument("--min_olen", type=int, default=1200, help="Min aligned length [1200]")
ap.add_argument("--min_ident", type=float, default=0.84, help="Min identity [0.84]")
ap.add_argument("--w_floor", type=float, default=0.12, help="Min edge weight [0.12]")
args = ap.parse_args()

deg, wdeg = defaultdict(int), defaultdict(float)

with openg(args.paf, "rt") as f:
    for ln in f:
        if not ln.strip() or ln[0] == "#":
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 12:
            continue
        q, ql, qs, qe = p[0], int(p[1]), int(p[2]), int(p[3])
        t, tl, ts, te = p[5], int(p[6]), int(p[7]), int(p[8])
        nm, al, mq = int(p[9]), int(p[10]), int(p[11])
        if q == t:
            continue
        ident = (nm / float(al)) if al > 0 else 0.0
        if al < args.min_olen or ident < args.min_ident:
            continue
        den = float(min(ql, tl)) if min(ql, tl) > 0 else 1.0
        w = ident * (al / den)
        if w < args.w_floor:
            continue
        deg[q] += 1
        wdeg[q] += w
        deg[t] += 1
        wdeg[t] += w

print("#read_id\tdegree\twdegree")
for r in wdeg:
    print(f"{r}\t{deg[r]}\t{wdeg[r]:.6f}")
