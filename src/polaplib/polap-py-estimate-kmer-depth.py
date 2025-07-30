#!/usr/bin/env python3
import sys
import gzip
from Bio import SeqIO

if len(sys.argv) < 5:
    print(f"Usage: {sys.argv[0]} <kmer_counts.txt> <reads.fq.gz/fa.gz> <k> <output_tsv>")
    sys.exit(1)

kmer_file, reads_file, k, report_file = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]

print("[1/3] Loading k-mer counts...")
kmer_counts = {}
with open(kmer_file) as f:
    for line in f:
        kmer, count = line.strip().split()
        kmer_counts[kmer] = int(count)

fmt = "fastq" if reads_file.endswith("q") else "fasta"
handle = gzip.open(reads_file, "rt") if reads_file.endswith(".gz") else open(reads_file)

print("[2/3] Estimating mean k-mer depth per read...")
with open(report_file, "w") as out:
    out.write("read_id\tmean_kmer_depth\n")
    for record in SeqIO.parse(handle, fmt):
        seq = str(record.seq)
        kmer_vals = []
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k].upper()
            rev_kmer = kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            canonical = min(kmer, rev_kmer)
            kmer_vals.append(kmer_counts.get(canonical, 0))
        mean = sum(kmer_vals) / len(kmer_vals) if kmer_vals else 0
        out.write(f"{record.id}\t{mean:.2f}\n")

handle.close()
print(f"[3/3] Output written to {report_file}")
