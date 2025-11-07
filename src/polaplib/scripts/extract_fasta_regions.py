# FILE: scripts/extract_fasta_regions.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Extract each BED interval as a FASTA record file under outdir, and also write a pooled FASTA.
BED name field is used as record ID if present; otherwise tract_id is taken from BED col4.

Usage: extract_fasta_regions.py <mt.fa> <tracts.bed> <outdir> <pooled.fa>
"""
import sys
from Bio import SeqIO

if len(sys.argv) < 5:
    sys.stderr.write(
        "Usage: extract_fasta_regions.py <mt.fa> <tracts.bed> <outdir> <pooled.fa>\n"
    )
    sys.exit(1)
fa, bed, outdir, pooled = sys.argv[1:5]

# index sequences
seqs = {r.id: r.seq for r in SeqIO.parse(fa, "fasta")}

records = []
with open(bed) as fh:
    for line in fh:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        chrom, s, e, name, score, strand = f[:6]
        s, e = int(s), int(e)
        seq = seqs[chrom][s:e]
        if strand == "-":
            seq = seq.reverse_complement()
        rec_id = f[3]
        # Write per-tract file
        outpath = f"{outdir}/{rec_id}.fa"
        with open(outpath, "w") as o:
            o.write(f">{rec_id}\n")
            for i in range(0, len(seq), 60):
                o.write(str(seq[i : i + 60]) + "\n")
        records.append((rec_id, seq))

with open(pooled, "w") as o:
    for rec_id, seq in records:
        o.write(f">{rec_id}\n")
        for i in range(0, len(seq), 60):
            o.write(str(seq[i : i + 60]) + "\n")
