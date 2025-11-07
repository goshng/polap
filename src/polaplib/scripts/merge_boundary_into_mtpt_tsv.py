# FILE: scripts/merge_boundary_into_mtpt_tsv.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Merge boundary support counts into MTPT TSV.
"""
import sys, csv

if len(sys.argv) < 3:
    sys.stderr.write(
        "Usage: merge_boundary_into_mtpt_tsv.py <mtpt.tsv> <boundary.tsv>\n"
    )
    sys.exit(1)
mtpt, bound = sys.argv[1:3]
bmap = {}
with open(bound) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        bmap[r["tract_id"]] = (r["left_span_reads"], r["right_span_reads"])

with open(mtpt) as fh:
    header = fh.readline().rstrip("\n")
    print(header + "\tleft_span_reads\tright_span_reads")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        tid = f[0]
        l, r = bmap.get(tid, ("0", "0"))
        print("\t".join(f + [l, r]))
