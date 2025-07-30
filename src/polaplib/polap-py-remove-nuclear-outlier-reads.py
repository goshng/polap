#!/usr/bin/env python3
import argparse
import gzip
from Bio import SeqIO


def canonical_kmer(kmer):
    return min(kmer, kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1])


def load_kmer_counts(kmer_file):
    kmer_counts = {}
    with open(kmer_file) as f:
        for line in f:
            k, c = line.strip().split()
            kmer_counts[k] = int(c)
    return kmer_counts


def filter_reads(reads_file, kmer_counts_file, k, mkc, minlen, output_prefix):
    # fmt = "fastq" if reads_file.endswith("q.gz") else "fasta"
    fmt = "fastq"
    open_func = gzip.open if reads_file.endswith(".gz") else open

    lkc = 0.3 * mkc
    hkc = 5 * mkc

    keep_fp = gzip.open(f"{output_prefix}.fastq.gz", "wt")
    discard_fp = gzip.open(f"{output_prefix}.discarded.fastq.gz", "wt")

    kmer_counts = load_kmer_counts(kmer_counts_file)
    total = kept = filtered = 0

    for record in SeqIO.parse(open_func(reads_file, "rt"), fmt):
        seq = str(record.seq).upper()
        if len(seq) < minlen:
            SeqIO.write(record, discard_fp, fmt)
            filtered += 1
            continue

        kmers = [canonical_kmer(seq[i : i + k]) for i in range(len(seq) - k + 1)]
        counts = [kmer_counts.get(km, 0) for km in kmers]
        too_low = sum(1 for c in counts if c < lkc)
        too_high = sum(1 for c in counts if c > hkc)

        if (too_low + too_high) > len(counts) / 5:
            SeqIO.write(record, discard_fp, fmt)
            filtered += 1
        else:
            SeqIO.write(record, keep_fp, fmt)
            kept += 1

        total += 1

    keep_fp.close()
    discard_fp.close()

    print(f"✅ Filtering summary:")
    print(f" - Total reads:  {total}")
    print(f" - Kept:         {kept}")
    print(f" - Removed:      {filtered}")
    print(f" - Output:       {output_prefix}.fastq.gz")
    print(f" - Discarded:    {output_prefix}.discarded.fastq.gz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True)
    parser.add_argument("--kmer-counts", required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--mkc", type=float, required=True)
    parser.add_argument("--minlen", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    filter_reads(
        args.reads, args.kmer_counts, args.k, args.mkc, args.minlen, args.output
    )
