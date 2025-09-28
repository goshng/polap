#!/usr/bin/env python3
# polap-py-cds-coverage-from-paf.py  v0.0.1
# From miniprot PAF and contig .fai, emit contigs with CDS frac <= max_cds_frac.
import sys, argparse
from collections import defaultdict
ap = argparse.ArgumentParser()
ap.add_argument('paf'); ap.add_argument('contig_fai')
ap.add_argument('--max_cds_frac', type=float, default=0.20)
args = ap.parse_args()
length = {}
with open(args.contig_fai) as f:
    for ln in f:
        if not ln.strip(): continue
        p = ln.split('\t'); length[p[0]] = int(p[1])
cov = defaultdict(int)
with open(args.paf) as f:
    for ln in f:
        if not ln.strip() or ln[0] == '#': continue
        p = ln.split('\t')
        if len(p) < 12: continue
        tname = p[5]; tstart = int(p[7]); tend = int(p[8])
        cov[tname] += max(0, tend - tstart)
for ctg, L in length.items():
    frac = (cov.get(ctg, 0) / float(L)) if L > 0 else 0.0
    if frac <= args.max_cds_frac: print(ctg)
