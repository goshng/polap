#!/usr/bin/env python3
"""
Parse a HiMT mitochondrial HTML report and extract:
- gene integrity matrix (conserved protein coding genes vs copies)
- mitochondrial basic information table (total length, total contig number, etc.)
- contig information table
- derived gene-set completeness and simple 3C/4C-style metrics.

Usage:
    himt_html_to_metrics.py --input himt_mitochondrial.html --out-prefix sample_prefix \
        [--complete-threshold 95] [--partial-threshold 50]
"""

import argparse
import json
import math
import re
import statistics
from pathlib import Path
from typing import Dict, List, Tuple


def extract_json_block(pattern: str, html: str, label: str) -> Dict:
    """Extract a JSON-like block from HTML using a regex capturing group and parse it."""
    m = re.search(pattern, html, re.S)
    if not m:
        raise RuntimeError(
            f"Could not find JSON block for {label} using pattern: {pattern}"
        )
    block = m.group(1)
    try:
        obj = json.loads(block)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse JSON for {label}: {e}")
    return obj


def write_gene_integrity(
    out_prefix: Path, copies: List[str], genes: List[str], z: List[List[float]]
) -> Path:
    gi_path = out_prefix.with_suffix("")  # strip any extension
    gi_path = Path(str(gi_path) + ".gene_integrity.tsv")
    with gi_path.open("w", encoding="utf-8") as out:
        out.write("Gene\tCopy\tGeneIntegrity\n")
        for gene, row in zip(genes, z):
            for copy, val in zip(copies, row):
                out.write(f"{gene}\t{copy}\t{val}\n")
    return gi_path


def write_gene_max_integrity(
    out_prefix: Path,
    genes: List[str],
    z: List[List[float]],
    complete_thr: float,
    partial_thr: float,
) -> Tuple[Path, Dict[str, float]]:
    gm_path = out_prefix.with_suffix("")
    gm_path = Path(str(gm_path) + ".gene_max_integrity.tsv")

    max_vals = {}
    status_counts = {"complete": 0, "partial": 0, "fragment": 0, "missing": 0}
    dup_genes = 0
    total_copy_nonzero = 0

    with gm_path.open("w", encoding="utf-8") as out:
        out.write("Gene\tMaxGeneIntegrity\tNumCopiesNonZero\tStatus\n")
        for gene, row in zip(genes, z):
            num_nonzero = sum(1 for v in row if v > 0)
            max_val = max(row) if row else 0.0
            if max_val == 0:
                status = "missing"
            elif max_val >= complete_thr:
                status = "complete"
            elif max_val >= partial_thr:
                status = "partial"
            else:
                status = "fragment"

            if num_nonzero > 1:
                dup_genes += 1
            total_copy_nonzero += num_nonzero

            status_counts[status] += 1
            max_vals[gene] = max_val
            out.write(f"{gene}\t{max_val}\t{num_nonzero}\t{status}\n")

    # gene-set stats summary
    total_genes = len(genes)
    completeness_prop = (
        status_counts["complete"] / total_genes if total_genes else float("nan")
    )
    mean_integrity = statistics.mean(max_vals.values()) if max_vals else float("nan")
    mean_integrity_prop = (
        mean_integrity / 100.0 if not math.isnan(mean_integrity) else float("nan")
    )
    avg_copies_per_gene = (
        total_copy_nonzero / total_genes if total_genes else float("nan")
    )

    summary = {
        "total_genes": total_genes,
        "complete_genes": status_counts["complete"],
        "partial_genes": status_counts["partial"],
        "fragment_genes": status_counts["fragment"],
        "missing_genes": status_counts["missing"],
        "gene_set_completeness_prop_complete": completeness_prop,
        "gene_set_mean_max_integrity": mean_integrity,
        "gene_set_mean_max_integrity_prop": mean_integrity_prop,
        "genes_with_multiple_copies": dup_genes,
        "avg_copies_per_gene_with_coverage": avg_copies_per_gene,
    }
    return gm_path, summary


def write_basic_info(out_prefix: Path, labels: List[str], values: List) -> Path:
    bi_path = out_prefix.with_suffix("")
    bi_path = Path(str(bi_path) + ".mito_basic_info.tsv")
    with bi_path.open("w", encoding="utf-8") as out:
        out.write("Metric\tValue\n")
        for label, value in zip(labels, values):
            out.write(f"{label}\t{value}\n")
    return bi_path


def write_contigs(out_prefix: Path, headers: List[str], cols: List[List]) -> Path:
    ctg_path = out_prefix.with_suffix("")
    ctg_path = Path(str(ctg_path) + ".contigs.tsv")
    with ctg_path.open("w", encoding="utf-8") as out:
        out.write("\t".join(headers) + "\n")
        for row in zip(*cols):
            out.write("\t".join(str(x) for x in row) + "\n")
    return ctg_path


def write_summary(
    out_prefix: Path,
    gene_summary: Dict[str, float],
    basic_info: Dict[str, float],
    contigs_cols: Dict[str, List],
) -> Path:
    """Write a simple 3C/4C-style summary TSV."""
    summ_path = out_prefix.with_suffix("")
    summ_path = Path(str(summ_path) + ".summary_metrics.tsv")

    # contig-related metrics
    lengths = [int(x) for x in contigs_cols.get("Length", [])]
    conserved_gene_nums = [
        int(x) for x in contigs_cols.get("Conserved gene number", [])
    ]
    gc_contents = [float(x) for x in contigs_cols.get("GC content", [])]

    num_contigs = len(lengths)
    tot_length_from_contigs = sum(lengths)
    mean_gc = (sum(gc_contents) / len(gc_contents)) if gc_contents else float("nan")
    genes_total = sum(conserved_gene_nums) if conserved_gene_nums else 0
    genes_per_contig_mean = (genes_total / num_contigs) if num_contigs else float("nan")
    genes_per_contig_max = max(conserved_gene_nums) if conserved_gene_nums else 0

    # N50 calculation from lengths (contiguity metric)
    n50_val = float("nan")
    if lengths:
        lengths_sorted = sorted(lengths, reverse=True)
        half = tot_length_from_contigs / 2.0
        csum = 0
        for L in lengths_sorted:
            csum += L
            if csum >= half:
                n50_val = L
                break

    basic_info_map = basic_info

    with summ_path.open("w", encoding="utf-8") as out:

        def w(key: str, value):
            out.write(f"{key}\t{value}\n")

        w("# Completeness-related metrics", "")
        for k, v in gene_summary.items():
            w(k, v)

        w("\n# Contiguity-related metrics", "")
        w("num_contigs_from_table2", num_contigs)
        w("total_length_from_contigs", tot_length_from_contigs)
        w("N50_from_contigs", n50_val)
        w("mean_gc_from_contigs", mean_gc)
        w("total_conserved_genes_in_contigs", genes_total)
        w("mean_conserved_genes_per_contig", genes_per_contig_mean)
        w("max_conserved_genes_in_single_contig", genes_per_contig_max)

        w("\n# Basic information from table1", "")
        for k, v in basic_info_map.items():
            w(k, v)

    return summ_path


def main():
    ap = argparse.ArgumentParser(
        description="Parse HiMT mitochondrial HTML report into tabular metrics."
    )
    ap.add_argument(
        "--input", "-i", required=True, help="HiMT mitochondrial HTML file."
    )
    ap.add_argument(
        "--out-prefix",
        "-o",
        required=True,
        help="Output prefix for generated TSV/metrics files.",
    )
    ap.add_argument(
        "--complete-threshold",
        type=float,
        default=95.0,
        help="Gene integrity threshold (%%) to call a gene 'complete' (default: 95).",
    )
    ap.add_argument(
        "--partial-threshold",
        type=float,
        default=50.0,
        help=(
            "Lower gene integrity threshold (%%) to call a gene 'partial'. "
            "Anything between partial and complete is 'partial', >0 and <partial is 'fragment'. "
            "Default: 50."
        ),
    )
    args = ap.parse_args()

    html_path = Path(args.input)
    out_prefix = Path(args.out_prefix)

    html = html_path.read_text(encoding="utf-8")

    # layout1: gene integrity heatmap
    layout1 = extract_json_block(
        r"var layout1 = (\{.*?\});\s*Plotly\.newPlot\('figure1'",
        html,
        "layout1 (gene integrity heatmap)",
    )
    data0 = layout1["data"][0]
    copies = data0["x"]
    genes = data0["y"]
    z = data0["z"]

    # gene integrity outputs
    gi_path = write_gene_integrity(out_prefix, copies, genes, z)
    gm_path, gene_summary = write_gene_max_integrity(
        out_prefix, genes, z, args.complete_threshold, args.partial_threshold
    )

    # table1: mitochondrial basic information
    table1 = extract_json_block(
        r"//table1\s*var tableLayout = (\{.*?\});\s*Plotly\.newPlot\('table1-chart'",
        html,
        "table1 (mitochondrial basic information)",
    )
    t1data = table1["data"][0]
    labels = t1data["cells"]["values"][0]
    values = t1data["cells"]["values"][1]
    bi_path = write_basic_info(out_prefix, labels, values)
    basic_info_map = {str(k): v for k, v in zip(labels, values)}

    # table2: contig information
    table2 = extract_json_block(
        r"//table2\s*var tableLayout = (\{.*?\});\s*Plotly\.newPlot\('table2-chart'",
        html,
        "table2 (contig information)",
    )
    t2data = table2["data"][0]
    headers = t2data["header"]["values"]
    cols = t2data["cells"]["values"]
    ctg_path = write_contigs(out_prefix, headers, cols)

    contig_cols_map = {h: col for h, col in zip(headers, cols)}

    # summary metrics
    summary_path = write_summary(
        out_prefix, gene_summary, basic_info_map, contig_cols_map
    )

    print("Wrote:")
    for p in [gi_path, gm_path, bi_path, ctg_path, summary_path]:
        print("  ", p)


if __name__ == "__main__":
    main()
