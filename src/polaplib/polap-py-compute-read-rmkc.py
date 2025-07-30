#!/usr/bin/env python3

import argparse
import gzip
import statistics
from Bio import SeqIO


def canonical_kmer(kmer):
    return min(kmer, kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1])


def load_kmer_counts(kmer_file):
    kmer_counts = {}
    with open(kmer_file) as f:
        for line in f:
            kmer, count = line.strip().split()
            kmer_counts[kmer] = int(count)
    return kmer_counts


def main():
    parser = argparse.ArgumentParser(
        description="Compute per-read median k-mer count (rmkc)"
    )
    parser.add_argument("--reads", required=True, help="Input gzipped FASTQ/FASTA")
    parser.add_argument("--kmer-counts", required=True, help="Jellyfish dump -c output")
    parser.add_argument("--k", type=int, required=True, help="K-mer size")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    kmer_counts = load_kmer_counts(args.kmer_counts)

    # fmt = "fastq" if args.reads.endswith("q.gz") else "fasta"
    fmt = "fastq"
    open_func = gzip.open if args.reads.endswith(".gz") else open
    reads = SeqIO.parse(open_func(args.reads, "rt"), fmt)

    with open(args.output, "w") as out:
        out.write("read_id\trmkc\n")
        for record in reads:
            seq = str(record.seq).upper()
            kmers = [
                canonical_kmer(seq[i : i + args.k])
                for i in range(len(seq) - args.k + 1)
            ]
            values = [kmer_counts.get(k, 0) for k in kmers]
            rmkc = statistics.median(values) if values else 0
            out.write(f"{record.id}\t{rmkc}\n")


if __name__ == "__main__":
    main()
