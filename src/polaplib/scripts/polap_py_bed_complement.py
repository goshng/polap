#!/usr/bin/env python3
# Version: v0.4.0
"""
Compute the complement BED: all unmasked regions per contig.

USAGE
  python3 polap_py_bed_complement.py \
    --fasta target.fa \
    --bed mask.bed \
    --out-bed allow.bed
"""
import sys, argparse, os
from collections import defaultdict

def read_fai_lengths(fasta):
    fai = fasta + ".fai"
    if not os.path.exists(fai):
        lens = {}
        name=None; L=0
        with open(fasta) as fh:
            for ln in fh:
                if ln.startswith(">"):
                    if name is not None: lens[name]=L
                    name = ln[1:].strip().split()[0]; L=0
                else: L += len(ln.strip())
            if name is not None: lens[name]=L
        return lens
    lens={}
    with open(fai) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            lens[f[0]] = int(f[1])
    return lens

def merge_intervals(it):
    if not it: return []
    it.sort()
    out=[list(it[0])]
    for s,e in it[1:]:
        if s <= out[-1][1]:
            if e > out[-1][1]: out[-1][1]=e
        else: out.append([s,e])
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True, help="mask BED (0-based half-open)")
    ap.add_argument("--out-bed", required=True)
    a = ap.parse_args()

    lens = read_fai_lengths(a.fasta)
    mask = defaultdict(list)
    if os.path.exists(a.bed) and os.path.getsize(a.bed)>0:
        with open(a.bed) as fh:
            for ln in fh:
                if not ln.strip(): continue
                t,s,e = ln.rstrip("\n").split("\t")[:3]
                s=int(s); e=int(e)
                if e<=s: continue
                if t not in lens: continue
                # clamp within contig
                s=max(0, min(s, lens[t])); e=max(0, min(e, lens[t]))
                if e<=s: continue
                mask[t].append((s,e))
    out=[]
    for t,L in lens.items():
        if t not in mask or not mask[t]:
            out.append((t,0,L)); continue
        iv = merge_intervals(mask[t])
        start=0
        for s,e in iv:
            if s>start: out.append((t,start,s))
            start = max(start, e)
        if start < L: out.append((t,start,L))
    with open(a.out_bed,"w") as oh:
        for t,s,e in sorted(out, key=lambda x:(x[0],x[1],x[2])):
            if s<e: oh.write(f"{t}\t{s}\t{e}\n")
    sys.stderr.write(f"[bed_complement] wrote {len(out)} intervals to {a.out_bed}\n")

if __name__ == "__main__":
    main()

