#!/usr/bin/env python3
# scripts/report_readassemble_mt0.py
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a JSON-only KPI bundle for MT-0 (protein-guided read selection).
# Inputs are discovered under --base-dir (polap-readassemble/annotate-read-mtseed).
#
# Output:
#   --out <base-dir>/report/mt0-report.json
#
# KPIs:
#   - read_label_counts: PT, MT (from pt.id.all.txt / mt.id.all.txt), AMB/NA left blank unless sources are added
#   - mt_candidates: count & total Length from contig-annotation-depth-table.txt (base)
#   - pt_candidates: count & total Length from pt-contig-annotation-depth-table.txt (base)
#   - gene_tallies (optional): sum MT/PT cols from assembly_info_organelle_annotation_count-all.txt
#   - qc: scatter PDF paths for mt and pt
#
import os
import sys
import csv
import json
import argparse


def safe_line_count(path: str) -> str:
    """Return number of non-empty lines as string, or empty string if missing."""
    if not path or not os.path.isfile(path):
        return ""
    try:
        n = 0
        with open(path, encoding="utf-8") as f:
            for ln in f:
                if ln.strip():
                    n += 1
        return str(n)
    except Exception:
        return ""


def sum_length_and_count(path: str) -> dict:
    """
    From a TSV with a 'Length' column, compute:
      {"count": "<rows>", "sum_length": "<sum Length>"}
    Empty strings if missing or unparsable.
    """
    out = {"count": "", "sum_length": ""}
    if not path or not os.path.isfile(path):
        return out
    try:
        with open(path, newline="", encoding="utf-8") as f:
            rd = csv.reader(f, delimiter="\t")
            hdr = next(rd, None)
            if not hdr:
                return out
            # case-insensitive hdr map
            pos = {name.strip().lower(): i for i, name in enumerate(hdr)}
            idx = pos.get("length")
            if idx is None:
                return out
            c = 0
            s = 0
            for row in rd:
                if not row or idx >= len(row):
                    continue
                try:
                    L = float(row[idx])
                    if L > 0:
                        c += 1
                        s += int(L)
                except Exception:
                    continue
            out["count"] = str(c)
            out["sum_length"] = str(s)
            return out
    except Exception:
        return out


def parse_gene_tallies(path: str) -> dict:
    """
    Sum columns named 'PT' and 'MT' (case-insensitive) if present.
    Returns {"pt_genes_total":"", "mt_genes_total":""} if not available.
    """
    res = {"pt_genes_total": "", "mt_genes_total": ""}
    if not path or not os.path.isfile(path):
        return res
    try:
        with open(path, newline="", encoding="utf-8") as f:
            rd = csv.reader(f, delimiter="\t")
            hdr = next(rd, None)
            if not hdr:
                return res
            pos = {name.strip().lower(): i for i, name in enumerate(hdr)}
            i_pt = pos.get("pt")
            i_mt = pos.get("mt")
            if i_pt is None and i_mt is None:
                return res
            sum_pt = 0
            sum_mt = 0
            for row in rd:
                if i_pt is not None and i_pt < len(row):
                    try:
                        sum_pt += int(float(row[i_pt]))
                    except Exception:
                        pass
                if i_mt is not None and i_mt < len(row):
                    try:
                        sum_mt += int(float(row[i_mt]))
                    except Exception:
                        pass
            if i_pt is not None:
                res["pt_genes_total"] = str(sum_pt)
            if i_mt is not None:
                res["mt_genes_total"] = str(sum_mt)
            return res
    except Exception:
        return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--base-dir",
        required=True,
        help="Path to polap-readassemble/annotate-read-mtseed",
    )
    ap.add_argument("--out", required=True, help="Output JSON path")
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)

    # Files (relative to base)
    assembly_info_all = os.path.join(
        base, "assembly_info_organelle_annotation_count-all.txt"
    )
    mt_candidates_tsv = os.path.join(base, "contig-annotation-depth-table.txt")
    pt_candidates_tsv = os.path.join(base, "pt-contig-annotation-depth-table.txt")

    pt_ids = os.path.join(base, "pt.id.all.txt")
    pt_scatter_pdf = os.path.join(
        base, "pt", "pt-contig-annotation-depth-table.txt.scatter.pdf"
    )

    mt_ids = os.path.join(base, "mt.id.all.txt")
    mt_scatter_pdf = os.path.join(
        base, "mt", "contig-annotation-depth-table.txt.scatter.pdf"
    )

    # Read label counts
    pt_count = safe_line_count(pt_ids)
    mt_count = safe_line_count(mt_ids)
    amb_count = ""  # add a source later if you record AMB ids
    na_count = ""  # add a source later if you record NA ids

    # Candidates (MT table at base/, PT table at base/)
    mt_cand = sum_length_and_count(mt_candidates_tsv)
    pt_cand = sum_length_and_count(pt_candidates_tsv)

    # Gene tallies (optional)
    gene_tally = parse_gene_tallies(assembly_info_all)

    data = {
        "base_dir": base,
        "kpis": {
            "read_label_counts": {
                "PT": pt_count,
                "MT": mt_count,
                "AMB": amb_count,
                "NA": na_count,
            },
            "mt_candidates": {
                "count": mt_cand.get("count", ""),
                "sum_length": mt_cand.get("sum_length", ""),
            },
            "pt_candidates": {
                "count": pt_cand.get("count", ""),
                "sum_length": pt_cand.get("sum_length", ""),
            },
            "gene_tallies": {
                "pt_total": gene_tally.get("pt_genes_total", ""),
                "mt_total": gene_tally.get("mt_genes_total", ""),
            },
        },
        "qc": {
            "mt_scatter_pdf": mt_scatter_pdf if os.path.isfile(mt_scatter_pdf) else "",
            "pt_scatter_pdf": pt_scatter_pdf if os.path.isfile(pt_scatter_pdf) else "",
        },
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
