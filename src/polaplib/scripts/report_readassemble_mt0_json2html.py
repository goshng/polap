#!/usr/bin/env python3
# scripts/report_readassemble_mt0_json2html.py
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Read mt0-report.json (from polap-bash-report-readassemble-mt0.sh)
# and write a compact HTML summary with species-rooted asset paths.
#
# Sections:
#   1) Read label counts (PT/MT)
#   2) Candidates (MT>PT and PT>MT) — count + total bases in Mb and Gb
#   3) QC: embed scatter PDFs for mt and pt if present
#
# Usage:
#   python3 report_readassemble_mt0_json2html.py --json mt0-report.json --html mt0-report.html
#
import os
import re
import sys
import json
import html
import argparse


def esc(x):
    return html.escape("" if x is None else str(x))


def fmt_gb(x):
    try:
        return f"{float(x)/1e9:.2f} Gb"
    except:
        return ""


def fmt_mb(x):
    try:
        return f"{float(x)/1e6:.1f} Mb"
    except:
        return ""


def find_species_name(json_path: str):
    # Guess species folder name from .../<Species>/v6/0/...
    m = re.search(r"/([^/]+)/v\d+/", json_path)
    return m.group(1) if m else None


def rel_from_species(path: str, species: str | None) -> str:
    if not path or not species:
        return path
    i = path.find(species)
    return path[i:] if i != -1 else path


def file_exists_from_json_dir(json_path: str, relpath_from_species: str) -> bool:
    """Check existence of file addressed by path starting with <species>/..."""
    if not relpath_from_species:
        return False
    species = find_species_name(json_path)
    if not species:  # fallback
        return os.path.exists(relpath_from_species)
    # Walk up from the JSON dir until we hit the species dir
    cur = os.path.dirname(os.path.abspath(json_path))
    found = None
    for _ in range(12):
        if os.path.basename(cur) == species:
            found = cur
            break
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    if not found:
        return os.path.exists(relpath_from_species)
    # If relpath already starts with species/, strip the prefix to join to found
    tail = (
        relpath_from_species.split(species + os.sep, 1)[-1]
        if relpath_from_species.startswith(species + os.sep)
        else relpath_from_species
    )
    abs_path = os.path.join(found, tail)
    return os.path.exists(abs_path)


def build_html(data: dict, json_path: str) -> str:
    species = find_species_name(json_path)
    title = f"MT-0 Report ({species})" if species else "MT-0 Report"

    kpis = data.get("kpis", {})
    counts = kpis.get("read_label_counts", {})
    mt_cand = kpis.get("mt_candidates", {})
    pt_cand = kpis.get("pt_candidates", {})
    qc = data.get("qc", {})

    # Read label counts
    pt_count = counts.get("PT", "")
    mt_count = counts.get("MT", "")

    # Candidates
    mt_count_cand = mt_cand.get("count", "")
    mt_sum_len = mt_cand.get("sum_length", "")
    mt_sum_mb = fmt_mb(mt_sum_len)
    mt_sum_gb = fmt_gb(mt_sum_len)

    pt_count_cand = pt_cand.get("count", "")
    pt_sum_len = pt_cand.get("sum_length", "")
    pt_sum_mb = fmt_mb(pt_sum_len)
    pt_sum_gb = fmt_gb(pt_sum_len)

    # QC PDFs (species-rooted)
    mt_pdf = rel_from_species(qc.get("mt_scatter_pdf", ""), species)
    pt_pdf = rel_from_species(qc.get("pt_scatter_pdf", ""), species)

    css = """
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,Arial,sans-serif;margin:20px;line-height:1.45}
h1{font-size:1.5rem;margin-bottom:0.5rem}
h2{font-size:1.2rem;margin-top:1.4rem}
.card{border:1px solid #ddd;border-radius:8px;padding:12px;margin:10px 0;background:#fafafa}
.kv{display:grid;grid-template-columns:260px 1fr;gap:6px 14px;max-width:820px}
.kv div.k{font-weight:600}
.block{margin:8px 0}
object.pdf{width:100%; height:520px; border:1px solid #eee; border-radius:6px; background:#fff}
small{color:#666}
"""
    lines = []
    lines.append("<!DOCTYPE html>")
    lines.append("<html lang='en'><head><meta charset='utf-8'>")
    lines.append(f"<title>{esc(title)}</title>")
    lines.append(f"<style>{css}</style></head><body>")
    lines.append(f"<h1>{esc(title)}</h1>")

    # 1) Read label counts (for next stage)
    lines.append("<h2>1. Read label counts</h2>")
    lines.append("<div class='card kv'>")
    lines.append(
        "<div class='k'>PT-labeled reads</div><div>{}</div>".format(esc(pt_count))
    )
    lines.append(
        "<div class='k'>MT-labeled reads</div><div>{}</div>".format(esc(mt_count))
    )
    lines.append("</div>")

    # 2) Candidate sets
    lines.append("<h2>2. Candidate sets</h2>")
    lines.append("<div class='card kv'>")
    lines.append(
        "<div class='k'>MT candidates (MT > PT)</div><div>{} (total: {} / {})</div>".format(
            esc(mt_count_cand), esc(mt_sum_mb), esc(mt_sum_gb)
        )
    )
    lines.append(
        "<div class='k'>PT candidates (PT > MT)</div><div>{} (total: {} / {})</div>".format(
            esc(pt_count_cand), esc(pt_sum_mb), esc(pt_sum_gb)
        )
    )
    lines.append("</div>")

    # 3) QC PDFs
    lines.append("<h2>3. QC</h2>")
    lines.append("<div class='card'>")

    def embed_pdf(label: str, relp: str):
        if not relp:
            return
        if file_exists_from_json_dir(json_path, relp):
            lines.append(f"<div class='block'><b>{esc(label)}</b></div>")
            # Try to embed via <object>; browsers may render inline
            lines.append(
                f"<object class='pdf' data='{esc(relp)}' type='application/pdf'>"
            )
            lines.append(f"  <p>PDF: <a href='{esc(relp)}'>{esc(relp)}</a></p>")
            lines.append("</object>")

    embed_pdf("MT scatter (contig annotation vs depth)", mt_pdf)
    embed_pdf("PT scatter (contig annotation vs depth)", pt_pdf)
    if not any((mt_pdf, pt_pdf)):
        lines.append("<small>No QC PDFs available.</small>")

    lines.append("</div>")  # card
    lines.append("</body></html>")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--html", required=True)
    args = ap.parse_args()

    with open(args.json, encoding="utf-8") as f:
        data = json.load(f)
    html_doc = build_html(data, os.path.abspath(args.json))
    os.makedirs(os.path.dirname(os.path.abspath(args.html)), exist_ok=True)
    with open(args.html, "w", encoding="utf-8") as f:
        f.write(html_doc)


if __name__ == "__main__":
    sys.exit(main())
