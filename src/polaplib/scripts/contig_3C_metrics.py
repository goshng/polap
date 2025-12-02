#!/usr/bin/env python3
# scripts/contig_3C_metrics.py
# Version: 0.2.0
#
# Compute simple contig-level 3C-style metrics from a contig table:
#   - num contigs
#   - total length
#   - N50
#   - min / max length
#   - mean GC
#   - total conserved genes
#
# Expected input: TSV with header line:
#   Index\tContig name\tConserved gene number\tGC content\tLength
#
# Usage:
#   contig_3C_metrics.py contig_table.tsv [metrics_out.tsv]
#
# If metrics_out.tsv is omitted, results are written to stdout.

import csv
import sys
from typing import Dict, List


def compute_metrics(path: str) -> Dict[str, float]:
    n_ctg = 0
    total_len = 0
    total_genes = 0
    sum_gc = 0.0
    lengths: List[int] = []

    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)  # skip header
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 5:
                continue

            try:
                genes = float(row[2])
                gc = float(row[3])
                length = int(row[4])
            except ValueError:
                # skip malformed lines
                continue

            n_ctg += 1
            total_len += length
            total_genes += genes
            sum_gc += gc
            lengths.append(length)

    if n_ctg == 0:
        return {
            "num_contigs": 0,
            "total_length": 0,
            "N50_from_awk": 0,
            "min_length": 0,
            "max_length": 0,
            "mean_gc": 0.0,
            "total_conserved_genes": 0,
        }

    lengths_sorted = sorted(lengths)
    min_len = lengths_sorted[0]
    max_len = lengths_sorted[-1]

    half = total_len / 2.0
    csum = 0
    # N50: iterate from largest to smallest
    for L in reversed(lengths_sorted):
        csum += L
        if csum >= half:
            N50 = L
            break

    mean_gc = sum_gc / n_ctg

    return {
        "num_contigs": n_ctg,
        "total_length": total_len,
        "N50_from_awk": N50,  # name preserved for compatibility
        "min_length": min_len,
        "max_length": max_len,
        "mean_gc": mean_gc,
        "total_conserved_genes": total_genes,
    }


def main(argv=None) -> int:
    argv = sys.argv if argv is None else argv

    if len(argv) not in (2, 3):
        print(
            f"Usage: {argv[0]} contig_table.tsv [metrics_out.tsv]",
            file=sys.stderr,
        )
        return 1

    contig_path = argv[1]
    out_path = argv[2] if len(argv) == 3 else None

    try:
        metrics = compute_metrics(contig_path)
    except OSError as e:
        print(f"ERROR: Could not read contig table: {e}", file=sys.stderr)
        return 1

    out_f = sys.stdout if out_path is None else open(out_path, "w", newline="")
    try:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key in [
            "num_contigs",
            "total_length",
            "N50_from_awk",
            "min_length",
            "max_length",
            "mean_gc",
            "total_conserved_genes",
        ]:
            writer.writerow([key, metrics.get(key, "NA")])
    finally:
        if out_f is not sys.stdout:
            out_f.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
