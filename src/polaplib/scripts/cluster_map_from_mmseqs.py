#!/usr/bin/env python3
# cluster_map_from_mmseqs.py
# Version: v0.1.1
# Convert "rep<TAB>member" MMseqs TSV to "MTPT_id<TAB>CID" map with CID = C0001, C0002, ...
import sys, csv
if len(sys.argv) != 3:
    sys.stderr.write("Usage: cluster_map_from_mmseqs.py <raw.tsv> <out.tsv>\n")
    sys.exit(2)
raw, out = sys.argv[1], sys.argv[2]
clmap = {}
with open(raw, 'r', encoding='utf-8') as fh:
    for line in fh:
        line=line.strip()
        if not line: continue
        rep, mem = line.split("\t", 1)
        clmap.setdefault(rep, []).append(mem)
with open(out, 'w', encoding='utf-8', newline='') as o:
    w = csv.writer(o, delimiter="\t")
    w.writerow(["MTPT_id", "CID"])
    for idx, (rep, mems) in enumerate(clmap.items(), start=1):
        cid = f"C{idx:04d}"
        for m in mems:
            # MMseqs createtsv yields sequence IDs as they appeared; keep first token, strip leading >
            mtpt = m.split()[0]
            if mtpt.startswith(">"): mtpt = mtpt[1:]
            if mtpt.endswith(".fa"): mtpt = mtpt[:-3]
            w.writerow([mtpt, cid])
