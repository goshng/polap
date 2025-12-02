#!/usr/bin/env python3
# scripts/himt_parse_mito_html.py
# Version: v0.6.0
#
# Parse a HiMT mitochondrial HTML report and extract:
#   - <prefix>.gene_integrity.tsv
#   - <prefix>.gene_integrity_summary.tsv
#   - <prefix>.spliced_genes.tsv       (extended splicing info)
#   - <prefix>.mito_basic_info.tsv
#   - <prefix>.contig_table.tsv
#
# New in v0.6.0:
#   - Extended splice-gene detection:
#       * exon-count proxy (min_exon_count_min = ceil(sum_integrity / 100))
#       * cis- vs trans-splicing classification
#       * canonical trans-spliced gene flag (nad1/2/5/7/9, ccmFc)
#       * optional gene model length input to estimate ORF bp lengths
#
# Usage:
#   himt_parse_mito_html.py [--gene-model-lengths gene_lengths.tsv] report.html prefix
#
# gene_lengths.tsv format:
#   gene<TAB>length_bp
#

import argparse
import csv
import json
import math
import sys
import statistics
from typing import Dict, List, Any, Tuple


CANONICAL_TRANS_SPLICE_GENES = {"nad1", "nad2", "nad5", "nad7", "nad9", "ccmFc"}


def extract_json_object(html: str, var_name: str, plot_id: str) -> Dict[str, Any]:
    """
    Extract JSON object from:
        var <var_name> = {...};
        Plotly.newPlot('<plot_id>', ...)
    """
    anchor = f"Plotly.newPlot('{plot_id}'"
    anchor_idx = html.find(anchor)
    if anchor_idx == -1:
        raise RuntimeError(f"Cannot find Plotly.newPlot('{plot_id}')")

    assign_pat = f"var {var_name}"
    assign_idx = html.rfind(assign_pat, 0, anchor_idx)
    if assign_idx == -1:
        raise RuntimeError(f"Cannot find '{assign_pat}' before {anchor}")

    brace_start = html.find("{", assign_idx)
    if brace_start == -1:
        raise RuntimeError(f"No '{{' after {assign_pat}")

    brace_end = html.find("};", brace_start)
    if brace_end == -1:
        raise RuntimeError("No matching '};' found for JSON object")

    json_str = html[brace_start : brace_end + 1]
    return json.loads(json_str)


def write_tsv(path: str, header: List[str], rows: List[List[Any]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        if header:
            w.writerow(header)
        for row in rows:
            w.writerow(row)


def load_gene_lengths(path: str | None) -> Dict[str, int]:
    if path is None:
        return {}
    gene_lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            gene, length_str = fields[0], fields[1]
            try:
                length = int(length_str)
            except ValueError:
                continue
            gene_lengths[gene] = length
    return gene_lengths


def parse_gene_integrity(html: str, prefix: str, gene_lengths: Dict[str, int]) -> None:
    layout1 = extract_json_object(html, "layout1", "figure1")
    trace = layout1["data"][0]

    genes = trace["y"]
    copies = trace["x"]
    zmat = trace["z"]

    # -------- gene_integrity.tsv (long format) --------
    long_rows: List[List[Any]] = []
    for gene, row in zip(genes, zmat):
        for copy_label, val in zip(copies, row):
            long_rows.append([gene, copy_label, float(val)])
    write_tsv(f"{prefix}.gene_integrity.tsv", ["gene", "copy", "integrity"], long_rows)

    # -------- summary & extended splicing info --------
    best_vals: List[float] = []
    sum_vals: List[float] = []

    spliced_rows: List[List[Any]] = []

    for gene, row in zip(genes, zmat):
        vals = [float(v) for v in row]
        best = max(vals)
        total = sum(vals)
        nonzero = sum(1 for v in vals if v > 0)
        best_vals.append(best)
        sum_vals.append(total)

        # exon-count proxy (minimum number of fragments/exons)
        min_exon_count = math.ceil(total / 100.0) if total > 0 else 0

        # classify splicing type
        splice_type = "none"
        if total >= 90.0 and best < 95.0:
            if nonzero == 1:
                splice_type = "cis"
            elif nonzero > 1:
                splice_type = "trans"

        canonical = gene in CANONICAL_TRANS_SPLICE_GENES

        model_len = gene_lengths.get(gene)
        if model_len is None:
            model_len_bp = ""
            best_orf_bp = ""
            cum_orf_bp = ""
        else:
            model_len_bp = model_len
            best_orf_bp = int(round(best / 100.0 * model_len))
            cum_orf_bp = int(round(total / 100.0 * model_len))

        # Record all genes for transparency; downstream can filter by splice_type
        spliced_rows.append(
            [
                gene,
                best,
                total,
                nonzero,
                min_exon_count,
                splice_type,
                "TRUE" if canonical else "FALSE",
                model_len_bp,
                best_orf_bp,
                cum_orf_bp,
                ",".join(str(v) for v in vals),
            ]
        )

    # summary metrics
    genes_total = len(best_vals)
    genes_present = sum(1 for v in best_vals if v > 0)
    genes_missing = genes_total - genes_present
    genes_high95 = sum(1 for v in best_vals if v >= 95.0)

    mean_best = statistics.mean(best_vals) if best_vals else 0.0
    median_best = statistics.median(best_vals) if best_vals else 0.0
    completeness_prop = genes_present / genes_total if genes_total else 0.0
    high95_prop = genes_high95 / genes_total if genes_total else 0.0

    summary_rows = [
        ["genes_total", genes_total],
        ["genes_present", genes_present],
        ["genes_missing", genes_missing],
        ["genes_high_integrity_ge95", genes_high95],
        ["geneset_completeness_prop", completeness_prop],
        ["high_integrity_gene_prop", high95_prop],
        ["mean_best_integrity", mean_best],
        ["median_best_integrity", median_best],
    ]
    write_tsv(f"{prefix}.gene_integrity_summary.tsv", ["metric", "value"], summary_rows)

    # extended spliced genes TSV
    write_tsv(
        f"{prefix}.spliced_genes.tsv",
        [
            "gene",
            "best_integrity",
            "sum_integrity",
            "nonzero_copies",
            "min_exon_count_min",
            "splice_type",
            "canonical_trans_spliced",
            "model_length_bp",
            "best_orf_bp",
            "cumulative_orf_bp",
            "copy_integrities",
        ],
        spliced_rows,
    )


def parse_basic_info(html: str, prefix: str) -> None:
    table1 = extract_json_object(html, "tableLayout", "table1-chart")
    data = table1["data"][0]
    names = data["cells"]["values"][0]
    vals = data["cells"]["values"][1]

    rows = list(zip(names, vals))
    write_tsv(f"{prefix}.mito_basic_info.tsv", ["metric", "value"], rows)


def parse_contig_table(html: str, prefix: str) -> None:
    table2 = extract_json_object(html, "tableLayout", "table2-chart")
    data = table2["data"][0]
    header = data["header"]["values"]
    columns = data["cells"]["values"]

    n = len(columns[0])
    rows = []
    for i in range(n):
        rows.append([col[i] for col in columns])

    write_tsv(f"{prefix}.contig_table.tsv", header, rows)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Parse HiMT mito HTML into TSVs, with extended splice-gene detection."
    )
    ap.add_argument(
        "--gene-model-lengths",
        help="TSV with columns: gene<TAB>length_bp",
        default=None,
    )
    ap.add_argument("html")
    ap.add_argument("prefix")
    args = ap.parse_args()

    with open(args.html, encoding="utf-8") as f:
        html = f.read()

    gene_lengths = load_gene_lengths(args.gene_model_lengths)

    parse_gene_integrity(html, args.prefix, gene_lengths)
    parse_basic_info(html, args.prefix)
    parse_contig_table(html, args.prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
