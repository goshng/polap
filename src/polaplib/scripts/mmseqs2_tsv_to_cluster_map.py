# FILE: scripts/mmseqs2_tsv_to_cluster_map.py
#!/usr/bin/env python3
# VERSION: 0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Convert mmseqs2 createtsv (cluster mapping) into cluster_map.tsv with columns:
MTPT_id, CID
Assumes TSV rows: <rep_id>\t<member_id>
"""
import sys, csv

if len(sys.argv) != 3:
    sys.stderr.write(
        "Usage: mmseqs2_tsv_to_cluster_map.py <linclust.tsv> <cluster_map.tsv>\n"
    )
    sys.exit(1)
lin, out = sys.argv[1], sys.argv[2]
rep2cid = {}
rows = []
with open(lin, newline="", encoding="utf-8") as fh:
    for rep, mem, *rest in csv.reader(fh, delimiter="\t"):
        cid = rep2cid.setdefault(rep, f"CID_{len(rep2cid)+1:05d}")
        rows.append((mem, cid))
# ensure reps themselves included
for rep, cid in rep2cid.items():
    rows.append((rep, cid))
# unique
rows = sorted(set(rows))
with open(out, "w", encoding="utf-8") as o:
    o.write("MTPT_id\tCID\n")
    for m, c in rows:
        o.write(f"{m}\t{c}\n")
