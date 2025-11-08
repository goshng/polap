# FILE: scripts/align_clusters_mafft.py
#!/usr/bin/env python3
# VERSION: 0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Extract sequences per cluster and align them with MAFFT L-INS-i.
Inputs:
  all_mtpts.fa, cluster_map.tsv, outdir, threads
Writes:
  clusters/<CID>/align.fa  and clusters/<CID>/align.stats.tsv (simple stats).
"""
import sys, os, csv, subprocess, tempfile
from collections import defaultdict

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: align_clusters_mafft.py <all_mtpts.fa> <cluster_map.tsv> <outdir> <threads>\n"
    )
    sys.exit(1)
allfa, cmap, outdir, threads = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

# load sequences
seqs = {}
with open(allfa, "r", encoding="utf-8") as fh:
    cur = None
    buf = []
    for ln in fh:
        if ln.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf).replace("\n", "")
            cur = ln[1:].strip().split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if cur is not None:
        seqs[cur] = "".join(buf).replace("\n", "")

cid2ids = defaultdict(list)
with open(cmap, newline="", encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        cid2ids[r["CID"]].append(r["MTPT_id"])

for cid, ids in cid2ids.items():
    d = os.path.join(outdir, cid)
    os.makedirs(d, exist_ok=True)
    infa = os.path.join(d, "seqs.fa")
    with open(infa, "w") as o:
        for sid in ids:
            if sid in seqs:
                o.write(f">{sid}\n{seqs[sid]}\n")
    # align
    aln = os.path.join(d, "align.fa")
    subprocess.run(
        [
            "mafft",
            "--thread",
            str(threads),
            "--localpair",
            "--maxiterate",
            "1000",
            infa,
        ],
        check=True,
        stdout=open(aln, "w"),
        stderr=subprocess.DEVNULL,
    )
    # simple stats
    nseq = len(ids)
    totlen = sum(len(seqs.get(sid, "")) for sid in ids)
    with open(os.path.join(d, "align.stats.tsv"), "w") as o:
        o.write("CID\tnseq\ttot_bases\n")
        o.write(f"{cid}\t{nseq}\t{totlen}\n")
