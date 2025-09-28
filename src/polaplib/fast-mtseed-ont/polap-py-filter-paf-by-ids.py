#!/usr/bin/env python3
# polap-py-filter-paf-by-ids.py  v0.0.1
# Keep PAF rows whose q and t are both in keep.ids. Output may be .gz.
import sys, gzip, argparse
def openg(p, m='rt'): return gzip.open(p, m) if p.endswith('.gz') else open(p, m)
ap = argparse.ArgumentParser()
ap.add_argument('paf'); ap.add_argument('ids'); ap.add_argument('-o','--out',required=True)
args = ap.parse_args()
keep = set(x.strip() for x in open(args.ids) if x.strip())
g = gzip.open(args.out, 'wt') if args.out.endswith('.gz') else open(args.out, 'w')
with openg(args.paf, 'rt') as f, g:
    for ln in f:
        if not ln.strip() or ln[0] == '#': continue
        p = ln.split('\t')
        if len(p) < 6: continue
        if p[0] in keep and p[5] in keep: g.write(ln)
