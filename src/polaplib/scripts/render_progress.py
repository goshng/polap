#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_progress.py
Render a simple HTML table from the aggregate progress JSON.

Usage:
  # Inputs table
  python3 render_progress.py \
    --json progress-report.json \
    --type inputs \
    --html-out progress-inputs.html

  # MT table
  python3 render_progress.py \
    --json progress-report.json \
    --type mt \
    --html-out progress-mt.html

  # PT table
  python3 render_progress.py \
    --json progress-report.json \
    --type pt \
    --html-out progress-pt.html
"""

from __future__ import annotations
import argparse
import html
import json
import os
import sys
from typing import Any, Dict, List

VERSION = "render_progress.py v0.4.0"


def safe_get(d: Dict[str, Any], path: List[str], default="NA") -> Any:
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    if cur in (None, "", []):
        return default
    return cur


def fmt(val: Any) -> str:
    return "NA" if val in (None, "", "NA") else str(val)


def build_rows_inputs(model: Dict) -> (List[str], List[List[str]]):
    header = [
        "species-code",
        "long-read total bases (sum_len)",
        "long-read read count (num_seqs)",
        "long-read AvgQual",
    ]
    rows: List[List[str]] = []
    for rec in model.get("species", []):
        code = fmt(rec.get("code_full"))
        sum_len = fmt(safe_get(rec, ["reads_summary", "long", "sum_len"]))
        num_seqs = fmt(safe_get(rec, ["reads_summary", "long", "num_seqs"]))
        avgq = fmt(safe_get(rec, ["reads_summary", "long", "AvgQual"]))
        rows.append([code, sum_len, num_seqs, avgq])
    return header, rows


def build_rows_mt(model: Dict) -> (List[str], List[List[str]]):
    header = [
        "species-code",
        "long-read SRA accession",
        "long-read total bases (sum_len)",
        "actually used long-read bases (sum_len)",
        "mtDNA NCBI accession",
        "mtDNA NCBI sequence length",
        "polap mt (seq count)",
        "polap mt (total bases)",
    ]
    rows: List[List[str]] = []
    for rec in model.get("species", []):
        code = fmt(rec.get("code_full"))
        long_sra = fmt(safe_get(rec, ["sra", "long"]))
        total_long = fmt(safe_get(rec, ["reads_summary", "long", "sum_len"]))
        used_long = fmt(safe_get(rec, ["reads_used", "long", "sum_len"]))
        mt_acc = fmt(safe_get(rec, ["ncbi", "mt", "accession"]))
        mt_len = fmt(safe_get(rec, ["ncbi", "mt", "length"]))
        mt_count = fmt(safe_get(rec, ["polap", "mt", "seq_count"]))
        mt_total = fmt(safe_get(rec, ["polap", "mt", "total_bases"]))
        rows.append(
            [code, long_sra, total_long, used_long, mt_acc, mt_len, mt_count, mt_total]
        )
    return header, rows


def build_rows_pt(model: Dict) -> (List[str], List[List[str]]):
    header = [
        "species-code",
        "long-read SRA accession",
        "long-read total bases (sum_len)",
        "actually used long-read bases (sum_len)",
        "ptDNA NCBI accession",
        "ptDNA NCBI sequence length",
        "polap pt (seq count)",
        "polap pt (total bases)",
    ]
    rows: List[List[str]] = []
    for rec in model.get("species", []):
        code = fmt(rec.get("code_full"))
        long_sra = fmt(safe_get(rec, ["sra", "long"]))
        total_long = fmt(safe_get(rec, ["reads_summary", "long", "sum_len"]))
        used_long = fmt(safe_get(rec, ["reads_used", "long", "sum_len"]))
        pt_acc = fmt(safe_get(rec, ["ncbi", "pt", "accession"]))
        pt_len = fmt(safe_get(rec, ["ncbi", "pt", "length"]))
        pt_count = fmt(safe_get(rec, ["polap", "pt", "seq_count"]))
        pt_total = fmt(safe_get(rec, ["polap", "pt", "total_bases"]))
        rows.append(
            [code, long_sra, total_long, used_long, pt_acc, pt_len, pt_count, pt_total]
        )
    return header, rows


def html_page(
    title: str, header: List[str], rows: List[List[str]], meta: Dict[str, str]
) -> str:
    # Minimal, dependency-free HTML with simple table styling
    meta_info = f"Generated with {html.escape(VERSION)} from {html.escape(meta.get('version', ''))}"
    s = []
    s.append("<!doctype html>")
    s.append('<html lang="en">')
    s.append("<head>")
    s.append('<meta charset="utf-8">')
    s.append(f"<title>{html.escape(title)}</title>")
    s.append(
        """
<style>
  body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, Noto Sans, Helvetica, Arial, 'Apple Color Emoji', 'Segoe UI Emoji', sans-serif;
         margin: 24px; }
  h1 { margin-bottom: 0.2rem; }
  .meta { color: #555; font-size: 0.9rem; margin-bottom: 1.0rem; }
  table { border-collapse: collapse; width: 100%; font-size: 0.95rem; }
  th, td { border: 1px solid #ddd; padding: 6px 8px; text-align: left; }
  thead th { background: #f5f5f7; position: sticky; top: 0; z-index: 1; }
  tr:nth-child(even) { background: #fafafa; }
  .small { font-size: 0.85rem; color: #666; }
</style>
    """
    )
    s.append("</head>")
    s.append("<body>")
    s.append(f"<h1>{html.escape(title)}</h1>")
    s.append(f"<div class='meta'>{meta_info}</div>")
    s.append("<table>")
    s.append("<thead><tr>")
    for h in header:
        s.append(f"<th>{html.escape(h)}</th>")
    s.append("</tr></thead>")
    s.append("<tbody>")
    for r in rows:
        s.append("<tr>")
        for cell in r:
            s.append(f"<td>{html.escape(str(cell))}</td>")
        s.append("</tr>")
    s.append("</tbody></table>")
    s.append("</body></html>")
    return "\n".join(s)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Render one of three progress tables (inputs|mt|pt) to HTML."
    )
    ap.add_argument(
        "--json", required=True, help="Aggregate JSON from report_progress.py"
    )
    ap.add_argument(
        "--type",
        required=True,
        choices=["inputs", "mt", "pt"],
        help="Which table to render",
    )
    ap.add_argument("--html-out", required=True, help="Output HTML path")
    ap.add_argument("--title", default=None, help="Optional page title")
    args = ap.parse_args()

    with open(args.json, "r", encoding="utf-8") as fh:
        model = json.load(fh)

    if args.type == "inputs":
        header, rows = build_rows_inputs(model)
        title = args.title or "Organelle Assembly Progress – Inputs"
    elif args.type == "mt":
        header, rows = build_rows_mt(model)
        title = args.title or "Organelle Assembly Progress – mtDNA"
    else:
        header, rows = build_rows_pt(model)
        title = args.title or "Organelle Assembly Progress – ptDNA"

    page = html_page(title, header, rows, meta=model.get("meta", {}))
    os.makedirs(os.path.dirname(os.path.abspath(args.html_out)), exist_ok=True)
    with open(args.html_out, "w", encoding="utf-8") as out:
        out.write(page)
    print(f"Wrote {args.html_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
