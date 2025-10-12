#!/usr/bin/env python3
# Version: v0.8.0
"""
Collapse cp→mt BLAST hits into MTPT tracts; bin by recency.

Input:
  --blast6    cp->mt BLAST (outfmt 6): qseqid sseqid pident length qstart qend sstart send qcovs
  --mt-fasta  mitochondrial FASTA (for total length if needed)
  --min-len   minimum tract length [100]
  --recent    PID >= recent => 'recent' (default 97)
  --intermediate PID >= inter => 'intermediate' else 'ancient'
Output:
  --out-tsv   mtpt.tsv  (mt_contig mt_start mt_end length median_pid class cp_fragments)
  --out-bed   mtpt.bed  (for browser)
"""
from __future__ import annotations
import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast6", required=True)
    ap.add_argument("--mt-fasta", required=True)
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--recent", type=int, default=97)
    ap.add_argument("--intermediate", type=int, default=90)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-bed", required=True)
    a = ap.parse_args()

    hits = []
    with open(a.blast6) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            q, s, pid, L, qs, qe, ss, se, qcov = ln.rstrip("\n").split("\t")
            pid = float(pid)
            L = int(L)
            qs = int(qs)
            qe = int(qe)
            ss = int(ss)
            se = int(se)
            if L < a.min_len:
                continue
            if ss > se:
                ss, se = se, ss
            hits.append((s, ss, se, pid, L, q, qs, qe))
    hits.sort(key=lambda x: (x[0], x[1], x[2]))

    # collapse per contig
    collapsed = []
    cur = None
    for s, ss, se, pid, L, q, qs, qe in hits:
        if cur is None or s != cur[0] or ss > cur[2] + 10:
            cur = [s, ss, se, [pid], [L], 1]
            collapsed.append(cur)
        else:
            cur[2] = max(cur[2], se)
            cur[3].append(pid)
            cur[4].append(L)
            cur[5] += 1

    def recency(p):
        if p >= a.recent:
            return "recent"
        if p >= a.intermediate:
            return "intermediate"
        return "ancient"

    with open(a.out_tsv, "w") as oh, open(a.out_bed, "w") as bed:
        oh.write(
            "mt_contig\tmt_start\tmt_end\tlength\tmedian_pid\tclass\tcp_fragments\n"
        )
        for c in collapsed:
            contig, start, end, PIDs, Ls, nfrag = c[0], c[1], c[2], c[3], c[4], c[5]
            PIDs.sort()
            med = PIDs[len(PIDs) // 2]
            cls = recency(med)
            length = end - start + 1
            oh.write(f"{contig}\t{start}\t{end}\t{length}\t{med:.2f}\t{cls}\t{nfrag}\n")
            bed.write(f"{contig}\t{start-1}\t{end}\tMTPT|{cls}|{med:.1f}\n")


if __name__ == "__main__":
    main()
