#!/usr/bin/env python3
# scripts/mt_mito_4c_metrics.py
# Version: v0.3.0
#
# Computes 4C metrics:
#   completeness
#   contiguity
#   correctness (proxies)
#   consistency/complexity (proxies)
#
# NEVER exits silently. ALWAYS writes PREFIX.4c_metrics.tsv or errors clearly.

import csv
import math
import os
import sys
import statistics
from typing import Dict, List, Tuple


def fail(msg):
    print(f"[4C ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def read_gene_integrity(path: str) -> Dict[str, List[float]]:
    if not os.path.exists(path):
        fail(f"Missing gene_integrity file: {path}")
    genes = {}
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gene = row["gene"]
            val = float(row["integrity"])
            genes.setdefault(gene, []).append(val)
    return genes


def read_splice_info(path: str) -> Tuple[int, int, int, int]:
    if not os.path.exists(path):
        return (0, 0, 0, 0)

    n_spliced = 0
    n_trans = 0
    n_canon = 0
    n_canon_spliced = 0

    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if "splice_type" not in r.fieldnames:
            fail(f"spliced_genes.tsv missing splice_type column: {path}")
        if "canonical_trans_spliced" not in r.fieldnames:
            fail(f"spliced_genes.tsv missing canonical_trans_spliced column: {path}")

        for row in r:
            stype = row["splice_type"]
            canon = row["canonical_trans_spliced"].upper() == "TRUE"

            is_spliced = stype in ("cis", "trans")
            if is_spliced:
                n_spliced += 1
                if stype == "trans":
                    n_trans += 1

            if canon:
                n_canon += 1
                if is_spliced:
                    n_canon_spliced += 1

    return n_spliced, n_trans, n_canon, n_canon_spliced


def read_basic_info(path: str) -> Dict[str, str]:
    if not os.path.exists(path):
        fail(f"Missing mito_basic_info file: {path}")
    info = {}
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            info[row["metric"]] = str(row["value"])
    return info


def read_contigs(path: str):
    if not os.path.exists(path):
        fail(f"Missing contig_table file: {path}")

    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        header_l = [h.lower() for h in header]

        def find(sub):
            for i, h in enumerate(header_l):
                if sub in h:
                    return i
            fail(f"Column containing '{sub}' not found in: {header}")

        idx_len = find("length")
        idx_gc = find("gc")
        idx_genes = find("conserved")

        lengths, gcs, gene_counts = [], [], []
        for row in reader:
            if not row:
                continue
            try:
                lengths.append(int(float(row[idx_len])))
            except:
                lengths.append(math.nan)
            try:
                gcs.append(float(row[idx_gc]))
            except:
                gcs.append(math.nan)
            try:
                gene_counts.append(float(row[idx_genes]))
            except:
                gene_counts.append(math.nan)

    return lengths, gcs, gene_counts


def gini(values):
    vals = [v for v in values if v >= 0]
    if not vals:
        return math.nan
    vals = sorted(vals)
    n = len(vals)
    cum = 0
    for i, v in enumerate(vals, start=1):
        cum += i * v
    tot = sum(vals)
    if tot == 0:
        return 0.0
    return (2 * cum) / (n * tot) - (n + 1) / n


def compute_4c(prefix: str):
    # 1. Inputs
    gi = f"{prefix}.gene_integrity.tsv"
    sp = f"{prefix}.spliced_genes.tsv"
    bi = f"{prefix}.mito_basic_info.tsv"
    ct = f"{prefix}.contig_table.tsv"

    genes = read_gene_integrity(gi)
    n_genes = len(genes)
    if n_genes == 0:
        fail("No genes found in gene_integrity.tsv")

    best_vals = []
    for g, vals in genes.items():
        best_vals.append(max(float(v) for v in vals))

    # splicing info
    n_sp, n_tr, n_can, n_can_sp = read_splice_info(sp)

    # contigs
    lengths, gcs, gene_counts = read_contigs(ct)
    lengths_clean = [L for L in lengths if isinstance(L, (int, float))]
    tot_len = sum(L for L in lengths_clean if not math.isnan(L))

    # basic info
    basic = read_basic_info(bi)

    # compute metrics
    metrics = []

    def add(cat, name, val, desc):
        metrics.append((cat, name, val, desc))

    genes_present = sum(1 for v in best_vals if v > 0)
    completeness = genes_present / n_genes

    add("completeness", "genes_total", n_genes, "Total genes")
    add("completeness", "genes_present", genes_present, "Present genes")
    add("completeness", "geneset_completeness_prop", completeness, "Proportion present")
    add(
        "completeness",
        "mean_best_integrity",
        statistics.mean(best_vals) / 100,
        "Mean best integrity",
    )
    add(
        "completeness",
        "median_best_integrity",
        statistics.median(best_vals) / 100,
        "Median best integrity",
    )

    # contiguity
    if lengths_clean:
        max_len = max(lengths_clean)
        N50 = math.nan
        half = tot_len / 2 if tot_len else math.nan
        csum = 0
        for L in sorted(lengths_clean, reverse=True):
            csum += L
            if csum >= half:
                N50 = L
                break
        add("contiguity", "num_contigs", len(lengths_clean), "# contigs")
        add("contiguity", "total_length", tot_len, "total length")
        add("contiguity", "N50", N50, "N50")

    # consistency/splicing
    frac_sp = n_sp / n_genes if n_genes else math.nan
    add("consistency", "frac_spliced_genes", frac_sp, "Fraction spliced (cis+trans)")

    if n_sp > 0:
        add(
            "consistency",
            "frac_trans_spliced_among_spliced",
            n_tr / n_sp,
            "Among spliced, fraction trans",
        )
    if n_can > 0:
        add(
            "consistency",
            "frac_canonical_trans_spliced_detected",
            n_can_sp / n_can,
            "Detected canonical trans-spliced",
        )

    return metrics


def main():
    if len(sys.argv) != 2:
        fail("Usage: mt_mito_4c_metrics.py PREFIX")

    prefix = sys.argv[1]
    out_path = f"{prefix}.4c_metrics.tsv"

    metrics = compute_4c(prefix)

    with open(out_path, "w", encoding="utf-8", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["category", "metric", "value", "description"])
        for row in metrics:
            w.writerow(row)

    print(f"[4C] Wrote: {out_path}")


if __name__ == "__main__":
    main()
