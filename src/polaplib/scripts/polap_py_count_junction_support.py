#!/usr/bin/env python3
# Version: v0.8.0
"""
Count reads that bridge junction boundaries in a junction-template FASTA.

Inputs:
  --bam     BAM of ONT reads aligned to junctions.fasta
  --fasta   junctions.fasta with headers containing 'flank=N'
  --meta    junctions.tsv (provenance; not parsed for counting)
  --min-mapq  default 20
  --min-span  default 200
Output:
  junction_id  reads_support  reads_left  reads_right
"""
from __future__ import annotations
import argparse, re, pysam


def parse_flank_headers(path):
    d = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                jid = ln[1:].split("|")[0]
                m = re.search(r"flank=(\d+)", ln)
                if not m:
                    raise SystemExit("Missing flank= in FASTA header")
                d[jid] = int(m.group(1))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-span", type=int, default=200)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    fl = parse_flank_headers(a.fasta)
    bam = pysam.AlignmentFile(a.bam, "rb")
    counts = {k: [0, 0, 0] for k in fl}  # both,left,right

    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.mapping_quality < a.min_mapq:
            continue
        ref = bam.get_reference_name(aln.reference_id)
        if ref not in counts:
            continue
        boundary = fl[ref] - 1  # 0-based
        # build merged ref coverage
        pos = aln.reference_start
        cov = []
        for op, l in aln.cigartuples or []:
            if op in (0, 7, 8):
                cov.append((pos, pos + l))
                pos += l
            elif op == 3:
                pos += l
        cov.sort()
        merged = []
        for s, e in cov:
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        left = right = False
        for s, e in merged:
            if s <= boundary - a.min_span and e >= boundary + a.min_span:
                left = right = True
                break
            if s <= boundary - a.min_span and e > boundary:
                left = True
            if s <= boundary and e >= boundary + a.min_span:
                right = True
        if left and right:
            counts[ref][0] += 1
        if left:
            counts[ref][1] += 1
        if right:
            counts[ref][2] += 1

    with open(a.out_tsv, "w") as oh:
        oh.write("junction_id\treads_support\treads_left\treads_right\n")
        for jid, (b, l, r) in counts.items():
            oh.write(f"{jid}\t{b}\t{l}\t{r}\n")


if __name__ == "__main__":
    main()
