#!/usr/bin/env python3
# Version: v0.4.0
"""
Create a merged, padded BED mask on the TARGET assembly from a PAF of assembly-vs-assembly hits.

- Reads PAF (OTHER -> TARGET)
- (Optionally) filters by identity/length (nmatch/alnLen, alnLen >=)
- Extracts [tname, tstart, tend], merges overlapping/adjacent intervals,
- Pads each merged interval by --pad, clamps to contig lengths from FASTA (or .fai),
- Re-merges after padding, writes BED (0-based, half-open).

USAGE
  python3 polap_py_mask_from_paf.py \
    --paf a2a.paf --fasta target.fa \
    --min-ident 0.85 --min-len 150 --pad 1000 \
    --out-bed mask.bed
"""
import sys, argparse, os
from collections import defaultdict

def read_fai_lengths(fasta):
    fai = fasta + ".fai"
    if not os.path.exists(fai):
        # build in-memory lengths from FASTA
        lens = {}
        name=None; L=0
        with open(fasta) as fh:
            for ln in fh:
                if ln.startswith(">"):
                    if name is not None: lens[name]=L
                    name = ln[1:].strip().split()[0]; L=0
                else:
                    L += len(ln.strip())
            if name is not None: lens[name]=L
        return lens
    lens={}
    with open(fai) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            lens[f[0]] = int(f[1])
    return lens

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--fasta", required=True, help="TARGET assembly (for contig lengths)")
    ap.add_argument("--min-ident", type=float, default=0.0, help="filter: nmatch/alnLen")
    ap.add_argument("--min-len",   type=int,   default=0,   help="filter: alnLen")
    ap.add_argument("--pad",       type=int,   default=0,   help="pad merged intervals by ±PAD bp")
    ap.add_argument("--out-bed",   required=True)
    return ap.parse_args()

def merge_intervals(it):
    if not it: return []
    it.sort()
    out=[list(it[0])]
    for s,e in it[1:]:
        if s <= out[-1][1]:  # overlap/adjacent
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s,e])
    return out

def main():
    a = parse_args()
    lens = read_fai_lengths(a.fasta)
    raw = defaultdict(list)  # tname -> list of (s,e)
    n_in=n_keep=0
    with open(a.paf) as fh:
        for ln in fh:
            if not ln.strip() or ln[0]=='#': continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12: continue
            tname = f[5]
            try:
                ts = int(f[7]); te = int(f[8])
                nmatch = int(f[9]); alen = int(f[10])
            except Exception:
                continue
            if ts > te: ts,te = te,ts
            if tname not in lens: continue
            n_in += 1
            keep = True
            if a.min_len > 0 and alen < a.min_len: keep = False
            if a.min_ident > 0.0:
                ident = nmatch/alen if alen>0 else 0.0
                if ident < a.min_ident: keep = False
            if keep:
                n_keep += 1
                raw[tname].append((ts,te))
    # merge, pad, clamp, re-merge
    out = []
    for tname, lst in raw.items():
        if not lst: continue
        merged = merge_intervals(lst)
        padded = []
        L = lens[tname]
        for s,e in merged:
            s2 = max(0, s - a.pad)
            e2 = min(L, e + a.pad)
            padded.append((s2,e2))
        padded.sort()
        merged2 = merge_intervals(padded)
        for s,e in merged2:
            if s<e:
                out.append((tname, s, e))
    with open(a.out_bed, "w") as oh:
        for t,s,e in sorted(out, key=lambda x:(x[0],x[1],x[2])):
            oh.write(f"{t}\t{s}\t{e}\n")
    sys.stderr.write(f"[mask_from_paf] kept {n_keep}/{n_in} PAF, wrote {len(out)} BED intervals to {a.out_bed}\n")

if __name__ == "__main__":
    main()

