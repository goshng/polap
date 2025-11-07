# FILE: scripts/mask_alignment_by_gap.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Mask (remove) alignment columns with gap fraction > threshold for each codon alignment in a dir.
Writes in-place: .masked.fna
"""
import sys, os
from Bio import AlignIO

if len(sys.argv) < 3:
    sys.stderr.write("Usage: mask_alignment_by_gap.py <aln_dir> <gap_frac>\n")
    sys.exit(1)
adir, thr = sys.argv[1:3]
thr = float(thr)

for f in sorted(os.listdir(adir)):
    if not f.endswith(".codon.aln.fna"):
        continue
    path = os.path.join(adir, f)
    aln = AlignIO.read(path, "fasta")
    n = len(aln)
    L = aln.get_alignment_length()
    keep = []
    for i in range(L):
        gap = sum(1 for r in aln if r[i] == "-")
        if gap / float(n) <= thr:
            keep.append(i)
    # write masked
    with open(path.replace(".codon.aln.fna", ".codon.aln.masked.fna"), "w") as o:
        for r in aln:
            seq = "".join(r.seq[i] for i in keep)
            o.write(f">{r.id}\n{seq}\n")
