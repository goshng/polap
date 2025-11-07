# FILE: scripts/partition_builder.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Concatenate masked codon alignments, output concat.fna and IQ-TREE partitions.nex (per gene).
"""
import sys, os
from Bio import SeqIO

if len(sys.argv) < 4:
    sys.stderr.write(
        "Usage: partition_builder.py <aln_dir> <out_concat.fna> <out_partitions.nex>\n"
    )
    sys.exit(1)

adir, outfa, outnex = sys.argv[1:4]
genes = sorted(
    [
        g[: -len(".codon.aln.masked.fna")]
        for g in os.listdir(adir)
        if g.endswith(".codon.aln.masked.fna")
    ]
)
# read all taxa
taxa = set()
for g in genes:
    for rec in SeqIO.parse(os.path.join(adir, g + ".codon.aln.masked.fna"), "fasta"):
        taxa.add(rec.id)
taxa = sorted(taxa)

# build concatenation
seqs = {t: "" for t in taxa}
parts = []
start = 1
for g in genes:
    aln = list(SeqIO.parse(os.path.join(adir, g + ".codon.aln.masked.fna"), "fasta"))
    L = len(aln[0].seq) if aln else 0
    aln_map = {r.id: str(r.seq) for r in aln}
    for t in taxa:
        seqs[t] += aln_map.get(t, "N" * L)
    parts.append((g, start, start + L - 1))
    start += L
with open(outfa, "w") as o:
    for t in taxa:
        o.write(f">{t}\n{seqs[t]}\n")
with open(outnex, "w") as o:
    o.write("#nexus\nbegin sets;\n")
    for g, s, e in parts:
        o.write(f"charset {g} = {s}-{e};\n")
    o.write("end;\n")
