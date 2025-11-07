# FILE: scripts/emit_cluster_alignments.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Create per-cluster FASTAs, align with MAFFT L-INS-i, compute simple stats.
Usage: emit_cluster_alignments.py <pooled.fa> <cluster_map.tsv> <out_cluster_dir> <threads>
Outputs: <out>/<CID>/align.fa and <out>/<CID>/align.stats.tsv
"""
import sys, os, subprocess, tempfile, statistics, math
from Bio import SeqIO
import pandas as pd

if len(sys.argv) < 5:
    sys.stderr.write(
        "Usage: emit_cluster_alignments.py <pooled.fa> <cluster_map.tsv> <outdir> <threads>\n"
    )
    sys.exit(1)
pooled, cmap, outd, threads = sys.argv[1:5]
threads = int(threads)

# load pooled
seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(pooled, "fasta")}

# cluster map
import csv, collections

cl = collections.defaultdict(list)
with open(cmap) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        cl[r["CID"]].append(r["MTPT_id"])

for cid, ids in cl.items():
    cdir = os.path.join(outd, cid)
    os.makedirs(cdir, exist_ok=True)
    fasta = os.path.join(cdir, "seqs.fa")
    with open(fasta, "w") as o:
        for tid in ids:
            if tid in seqs:
                o.write(f">{tid}\n{seqs[tid]}\n")
    aln = os.path.join(cdir, "align.fa")
    subprocess.run(
        [
            "mafft",
            "--thread",
            str(threads),
            "--localpair",
            "--maxiterate",
            "1000",
            fasta,
        ],
        check=True,
        stdout=open(aln, "w"),
    )
    # stats
    # pairwise PID naive
    alns = [str(rec.seq) for rec in SeqIO.parse(aln, "fasta")]
    if len(alns) >= 2:
        import itertools

        pids = []
        for a, b in itertools.combinations(alns, 2):
            matches = sum(1 for x, y in zip(a, b) if x == y and x != "-" and y != "-")
            denom = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
            if denom > 0:
                pids.append(matches / denom)
        occ = sum(
            1 for i in range(len(alns[0])) if any(s[i] != "-" for s in alns)
        ) / len(alns[0])
        with open(os.path.join(cdir, "align.stats.tsv"), "w") as o:
            o.write("CID\tn_seqs\tmean_pairwise_pid\talignment_occupancy\n")
            o.write(
                f"{cid}\t{len(alns)}\t{(sum(pids)/len(pids)) if pids else 1.0:.4f}\t{occ:.4f}\n"
            )
