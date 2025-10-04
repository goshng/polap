#!/usr/bin/env python3
# polap-py-refilter-edges-to-overlapness.py  v0.1.0
# Re-filter a loose edge list by a single edge-weight threshold (eweight)
# and recompute overlapness.tsv (degree, wdegree) WITHOUT re-running minimap2.
#
# INPUT edge file (TSV or .gz) with columns (header optional):
#   u   v   alen   ident   weight
#
# Filtering rule:
#   keep edge (u,v) iff weight >= --w_floor (aka eweight)
#   (min_olen / min_ident are intentionally ignored by design)
#
# OUTPUT (TSV):
#   #read_id  degree  wdegree
#   <node>    <int>   <float>
#
# Usage:
#   python polap-py-refilter-edges-to-overlapness.py \
#     --edges edges_loose.tsv.gz \
#     --w_floor 0.12 \
#     --out overlapness_strict.tsv
#
import sys
import os
import argparse
import gzip
from collections import defaultdict


def openg(path, mode="rt"):
    if path == "-":
        # stdin only allowed for reading text
        return sys.stdin
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def openw(path):
    if path in (None, "-", ""):
        return sys.stdout
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Refilter edge list (u v alen ident weight) by a single "
            "edge-weight threshold and emit overlapness.tsv."
        )
    )
    ap.add_argument(
        "--edges",
        required=True,
        help="edge list TSV(.gz): columns u v alen ident weight (header OK)",
    )
    ap.add_argument(
        "--w_floor",
        type=float,
        default=0.0,
        help="edge-weight floor (keep edges with weight >= this) [0.0]",
    )
    ap.add_argument(
        "--out",
        default="-",
        help="output overlapness.tsv ('.gz' ok, '-' for stdout) [stdout]",
    )
    return ap.parse_args()


def main():
    args = parse_args()

    deg = defaultdict(int)
    wdeg = defaultdict(float)

    with openg(args.edges, "rt") as f:
        first = True
        for ln in f:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            # Skip header if present
            if first and (p[0].lower() in ("u", "source", "node_u")):
                first = False
                continue
            first = False
            if len(p) < 5:
                continue
            u, v = p[0], p[1]
            # ignore self loops just in case
            if u == v:
                continue
            # parse weight; other columns are not used for filtering here
            try:
                w = float(p[4])
            except ValueError:
                continue
            if w < args.w_floor:
                continue
            # keep edge: update both endpoints
            deg[u] += 1
            wdeg[u] += w
            deg[v] += 1
            wdeg[v] += w

    out = openw(args.out)
    close_out = out is not sys.stdout
    try:
        out.write("#read_id\tdegree\twdegree\n")
        # union of keys (covers nodes that appeared only as one endpoint)
        nodes = set(deg.keys()) | set(wdeg.keys())
        for r in nodes:
            out.write(f"{r}\t{deg.get(r,0)}\t{wdeg.get(r,0.0):.6f}\n")
    finally:
        if close_out:
            out.close()


if __name__ == "__main__":
    main()
