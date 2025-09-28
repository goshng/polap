#!/usr/bin/env python3
# polap-py-topfrac-by-wdeg.py  v0.0.1
# Select top fraction of reads by weighted-degree from overlapness TSV.
import sys, argparse
ap = argparse.ArgumentParser()
ap.add_argument('overlapness_tsv')
ap.add_argument('--top_frac', type=float, default=0.2)
args = ap.parse_args()
rows = []
with open(args.overlapness_tsv) as f:
    for ln in f:
        if not ln.strip() or ln.startswith('#'): continue
        rid, deg, wdeg = ln.rstrip('\n').split('\t')
        rows.append((rid, float(wdeg)))
rows.sort(key=lambda x: x[1], reverse=True)
k = max(1, int(len(rows) * args.top_frac))
for i, (rid, _) in enumerate(rows):
    if i < k: print(rid)
