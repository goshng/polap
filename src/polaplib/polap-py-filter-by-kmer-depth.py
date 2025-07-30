#!/usr/bin/env python3
import sys
import gzip
from Bio import SeqIO

if len(sys.argv) < 5:
    print(f"Usage: {sys.argv[0]} <reads.fq.gz/fa.gz> <report.tsv> <min_depth> <prefix>")
    sys.exit(1)

reads_file, report_file, min_depth, prefix = (
    sys.argv[1],
    sys.argv[2],
    float(sys.argv[3]),
    sys.argv[4],
)

print("[1/4] Reading k-mer depth report...")
depths = {}
with open(report_file) as f:
    next(f)
    for line in f:
        read_id, val = line.strip().split("\t")
        depths[read_id] = float(val)

fmt = "fastq" if reads_file.endswith("q") else "fasta"
read_handle = (
    gzip.open(reads_file, "rt") if reads_file.endswith(".gz") else open(reads_file)
)
filtered_file = f"{prefix}_filtered.{fmt}.gz"
discarded_file = f"{prefix}_discarded.{fmt}.gz"
ids_file = f"{prefix}_discarded.ids.txt"

out_filt = gzip.open(filtered_file, "wt")
out_disc = gzip.open(discarded_file, "wt")
id_disc = open(ids_file, "w")

print("[2/4] Filtering reads...")
keep, discard = 0, 0
for record in SeqIO.parse(read_handle, fmt):
    avg = depths.get(record.id, 0)
    if avg > min_depth:
        SeqIO.write(record, out_filt, fmt)
        keep += 1
    else:
        SeqIO.write(record, out_disc, fmt)
        id_disc.write(record.id + "\n")
        discard += 1

print(f"[3/4] Kept: {keep}, Discarded: {discard}")
read_handle.close()
out_filt.close()
out_disc.close()
id_disc.close()
print(f"[4/4] Files: {filtered_file}, {discarded_file}, {ids_file}")
