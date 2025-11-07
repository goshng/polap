#!/usr/bin/env python3
# report_assemble.py
# Version: v1.2.1
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a BLOCK JSON for rendering a polap-assemble report.
# The JSON contains sections with "blocks" (figure/table/text/kpi) so a Jinja2 template
# can render a rich report (figures: png/pdf; tables: TSV or space-delimited; text).
#
# Usage:
#   python3 report_assemble.py \
#     --base-dir Brassica_rapa \
#     --out Brassica_rapa/v6/0/polap-assemble/report/assemble-report.json
#
# Notes:
# - This script reads files from the conventional run layout under <Species>/v6/0/polap-assemble
#   and produces a single JSON with the exact blocks you requested.
# - Figure hrefs are species-rooted (e.g., "<Species>/v6/0/...") so the renderer can convert them
#   to links relative to the final HTML location.

from __future__ import annotations
import argparse
import csv
import datetime
import json
import os
import re
from typing import Any, Dict, List, Optional, Tuple

VERSION = "report_assemble.py v1.2.0"


# ------------------------------ basic helpers ------------------------------


def _iso_now() -> str:
    try:
        return datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    except Exception:
        return ""


def _exists(path: str) -> bool:
    return bool(path) and os.path.exists(path)


def _species_and_root(base_dir: str) -> Tuple[str, str]:
    """Return (species_name, species_abs_root)."""
    sp = os.path.basename(os.path.abspath(base_dir.rstrip("/")))
    return sp, os.path.abspath(base_dir)


def _join_species_root(sp_root: str, rel: str) -> str:
    """Join species root and a species-relative path."""
    return os.path.join(sp_root, rel)


def _species_href(abs_path: str, species_root: str, species: str) -> str:
    """
    Produce a href like '<Species>/...' for any absolute file path under the species root.
    Fallback to basename if trimming fails.
    """
    ap = os.path.abspath(abs_path).replace("\\", "/")
    root = os.path.abspath(species_root).replace("\\", "/")
    if ap.startswith(root + "/"):
        tail = ap[len(root) + 1 :]
        return f"{species}/{tail}"
    i = ap.find("/" + species + "/")
    if i != -1:
        return ap[i + 1 :]
    return os.path.basename(ap)


# ------------------------------ table/figure/text blocks ------------------------------


def _read_tsv_select_columns(
    path: str, keep_cols: List[str]
) -> Optional[Dict[str, Any]]:
    """
    TSV reader selecting columns by header (case-sensitive).
    Returns table dict: {title, headers, rows} or None.
    """
    if not _exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd, None)
        if not hdr:
            return None
        hpos = {h: i for i, h in enumerate(hdr)}
        rows: List[List[str]] = []
        for r in rd:
            row = []
            for col in keep_cols:
                i = hpos.get(col)
                row.append(r[i] if i is not None and i < len(r) else "")
            rows.append(row)
    if not rows:
        rows = [[]]
    return {"title": os.path.basename(path), "headers": keep_cols, "rows": rows}


def _read_tsv_generic(path: str, max_rows: int = 2000) -> Optional[Dict[str, Any]]:
    """Read generic TSV into headers, rows."""
    if not _exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd, None)
        if hdr is None:
            # Fallback: split on whitespace, build generic headers
            f.seek(0)
            lines = [ln.strip() for ln in f if ln.strip()]
            if not lines:
                return None
            rows = [ln.split() for ln in lines[:max_rows]]
            width = max(len(r) for r in rows)
            headers = [f"col{i+1}" for i in range(width)]
            return {"title": os.path.basename(path), "headers": headers, "rows": rows}
        rows: List[List[str]] = []
        for r in rd:
            rows.append(r)
            if len(rows) >= max_rows:
                break
    return {"title": os.path.basename(path), "headers": hdr, "rows": rows}


def _read_space_table(path: str, max_rows: int = 2000) -> Optional[Dict[str, Any]]:
    """
    Read a whitespace-delimited table (space/tabs), infer headers from first row if it looks textual.
    """
    if not _exists(path):
        return None
    rows: List[List[str]] = []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split()
            rows.append(parts)
            if len(rows) >= max_rows:
                break
    if not rows:
        return None
    if any(any(ch.isalpha() for ch in cell) for cell in rows[0]):
        headers = rows[0]
        data = rows[1:]
    else:
        width = max(len(r) for r in rows)
        headers = [f"col{i+1}" for i in range(width)]
        data = rows
    return {"title": os.path.basename(path), "headers": headers, "rows": data}


def _text_block_from_file(
    path: str, title: str, limit_kb: int = 256
) -> Optional[Dict[str, Any]]:
    """Read small text file into a text block (truncated message if too large)."""
    if not _exists(path):
        return None
    size = os.path.getsize(path)
    if size > limit_kb * 1024:
        return {
            "type": "text",
            "text": {
                "title": title,
                "content": f"(file too large: {size} bytes)",
                "as_markdown": False,
            },
        }
    with open(path, encoding="utf-8", errors="replace") as f:
        content = f.read()
    return {
        "type": "text",
        "text": {"title": title, "content": content, "as_markdown": False},
    }


def _figure_block_from_path(
    path: str, label: str, species_root: str, species: str
) -> Optional[Dict[str, Any]]:
    if not _exists(path):
        return None
    href = _species_href(path, species_root, species)
    return {"type": "figure", "label": label, "href": href}


def _table_block(table: Dict[str, Any], title: Optional[str] = None) -> Dict[str, Any]:
    t = dict(table)
    if title:
        t["title"] = title
    # NOTE: keep nested .table form for compatibility with existing templates
    return {"type": "table", "table": t}


# ------------------------------ performance parsers ------------------------------


def _parse_summary_polap(path: str) -> Dict[str, str]:
    """
    Parse summary-polap-assemble.txt for:
      - elapsed_hms      "08:30:29"
      - net_increase_kb  "34743552"
      - disk_used_gb     "122"
    """
    out = {
        "elapsed_hms": "NA",
        "net_increase_kb": "NA",
        "disk_used_gb": "NA",
        "source": path if os.path.exists(path) else "NA",
    }
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                ln = line.strip()
                m = re.match(r"^Elapsed time:\s*([0-9]{2}:[0-9]{2}:[0-9]{2})", ln)
                if m:
                    out["elapsed_hms"] = m.group(1)
                    continue
                m = re.match(r"^Net increase:\s*([0-9]+)\s*KB", ln)
                if m:
                    out["net_increase_kb"] = m.group(1)
                    continue
                m = re.match(r"^Disk used:\s*([0-9]+)\s*GB", ln)
                if m:
                    out["disk_used_gb"] = m.group(1)
                    continue
    except FileNotFoundError:
        pass
    return out


def _parse_timing_polap(path: str) -> Dict[str, str]:
    """
    Parse timing-polap-assemble.txt for:
      - hostname
      - cpu_model
      - cpu_cores
    """
    out = {
        "hostname": "NA",
        "cpu_model": "NA",
        "cpu_cores": "NA",
        "source": path if os.path.exists(path) else "NA",
    }
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                ln = line.strip()
                m = re.match(r"^Hostname:\s*(.+)$", ln)
                if m:
                    out["hostname"] = m.group(1).strip()
                    continue
                m = re.match(r"^CPU Model:\s*(.+)$", ln)
                if m:
                    out["cpu_model"] = m.group(1).strip()
                    continue
                m = re.match(r"^CPU Cores:\s*([0-9]+)", ln)
                if m:
                    out["cpu_cores"] = m.group(1)
                    continue
    except FileNotFoundError:
        pass
    return out


# ------------------------------ sections builder ------------------------------


def build_sections(species_dir: str) -> Dict[str, Any]:
    """
    Build the full BLOCK JSON composed of sections and blocks as requested.
    species_dir: the species base directory (e.g., "Brassica_rapra")
    """
    species, sp_root = _species_and_root(species_dir)

    def S(rel: str) -> str:
        return _join_species_root(sp_root, rel)

    sections: List[Dict[str, Any]] = []

    # 1. Seed PT and MT reads — seqkit (3 columns)
    sec1_blocks: List[Dict[str, Any]] = []
    t1 = _read_tsv_select_columns(
        S("v6/0/polap-assemble/annotate-read-pt/pt0.fq.seqkit.stats.ta.txt"),
        ["num_seqs", "sum_len", "AvgQual"],
    )
    if t1:
        sec1_blocks.append(_table_block(t1, "PT reads (seqkit, pt0)"))
    sections.append({"id": "1", "title": "Seed PT and MT reads", "blocks": sec1_blocks})

    # 2. ptDNA assembly — figures + space table
    sec2_blocks: List[Dict[str, Any]] = []
    f0 = _figure_block_from_path(
        S("v6/0/polap-assemble/annotate-read-pt/pt/30-contigger/graph_final.png"),
        "PT assembly graph (pt0)",
        sp_root,
        species,
    )
    if f0:
        sec2_blocks.append(f0)
    f1 = _figure_block_from_path(
        S("v6/0/polap-assemble/annotate-read-pt/pt1/assembly_graph.png"),
        "PT assembly graph (pt1)",
        sp_root,
        species,
    )
    if f1:
        sec2_blocks.append(f1)
    t2 = _read_space_table(
        S(
            "v6/0/polap-assemble/annotate-read-pt/pt1/pt-contig-annotation-depth-table.txt"
        )
    )
    if t2:
        sec2_blocks.append(_table_block(t2, "PT annotation table (pt1)"))
    sections.append({"id": "2", "title": "ptDNA assembly", "blocks": sec2_blocks})

    # 3. Filter out PT reads — pt_thresh diag + nonPT (3 col)
    sec3_blocks: List[Dict[str, Any]] = []
    t3 = _read_tsv_generic(S("v6/0/polap-assemble/mtseed/pt_thresh.diag.tsv"))
    if t3:
        sec3_blocks.append(
            _table_block(t3, "PT read selection cutoff (pt_thresh.diag.tsv)")
        )
    t4 = _read_tsv_select_columns(
        S("v6/0/polap-assemble/mtseed/reads.nonpt.fq.gz.seqkit.stats.ta.tsv"),
        ["num_seqs", "sum_len", "AvgQual"],
    )
    if t4:
        sec3_blocks.append(_table_block(t4, "Non-PT reads (after PT filtering)"))
    sections.append({"id": "3", "title": "Filter out PT reads", "blocks": sec3_blocks})

    # 4. Compute read overlapness — PDFs
    sec4_blocks: List[Dict[str, Any]] = []
    p1 = _figure_block_from_path(
        S("v6/0/polap-assemble/mtseed/03-allvsall/04-qc/scan_wdegree_hist.pdf"),
        "Weighted degree distribution",
        sp_root,
        species,
    )
    if p1:
        sec4_blocks.append(p1)
    p2 = _figure_block_from_path(
        S("v6/0/polap-assemble/mtseed/03-allvsall/04-qc/scan_cum_wdegree.pdf"),
        "Cumulative weighted degree (threshold)",
        sp_root,
        species,
    )
    if p2:
        sec4_blocks.append(p2)
    sections.append(
        {"id": "4", "title": "Compute read overlapness", "blocks": sec4_blocks}
    )

    # 5. Filter out NT reads — TSV table
    sec5_blocks: List[Dict[str, Any]] = []
    t5 = _read_tsv_generic(
        S("v6/0/polap-assemble/mtseed/05-round/threshold_from_nuclear.tsv")
    )
    if t5:
        sec5_blocks.append(
            _table_block(t5, "Weighted-degree threshold to filter NT reads")
        )
    sections.append({"id": "5", "title": "Filter out NT reads", "blocks": sec5_blocks})

    # 6. Assemble mt0 — figure
    sec6_blocks: List[Dict[str, Any]] = []
    f_mt0 = _figure_block_from_path(
        S("v6/0/polap-assemble/mtseed/07-flye/30-contigger/graph_final.png"),
        "Assembly mt0 (graph)",
        sp_root,
        species,
    )
    if f_mt0:
        sec6_blocks.append(f_mt0)
    sections.append(
        {
            "id": "6",
            "title": "Assemble seed contigs using miniasm or mtDNA assembly mt0",
            "blocks": sec6_blocks,
        }
    )

    # 7. mtDNA assembly mt1..mt3
    subsecs: List[Dict[str, Any]] = []
    for name, idx in (("mt1", "7.1"), ("mt2", "7.2"), ("mt3", "7.3")):
        base_rel = f"v6/0/polap-assemble/mtseed/{name}"
        blocks: List[Dict[str, Any]] = []

        tlen = _text_block_from_file(
            S(f"{base_rel}/01-contig/contig_total_length.txt"), "Contig length (bp)"
        )
        if tlen:
            blocks.append(tlen)
        fig = _figure_block_from_path(
            S(f"{base_rel}/assembly_graph.png"),
            f"Assembly graph ({name})",
            sp_root,
            species,
        )
        if fig:
            blocks.append(fig)
        t_mt = _read_space_table(S(f"{base_rel}/contig-annotation-depth-table.txt"))
        if t_mt:
            blocks.append(_table_block(t_mt, f"MT annotation table ({name})"))

        next_map = {"mt1": "mt2", "mt2": "mt3"}
        if name in next_map:
            sn = _text_block_from_file(
                S(f"{base_rel}/mt.contig.name-{next_map[name]}"),
                f"Seed contig names → {next_map[name]}",
            )
            if sn:
                blocks.append(sn)

        if blocks:
            subsecs.append(
                {"id": idx, "title": f"mtDNA assembly {name}", "blocks": blocks}
            )

    sections.append(
        {"id": "7", "title": "mtDNA assembly", "blocks": [], "subsections": subsecs}
    )

    # 8. Oatk pathfinder — FASTA paths as text
    sec8_blocks: List[Dict[str, Any]] = []
    pt_fa = S("v6/0/polap-assemble/extract/oatk.pltd.ctg.fasta")
    if _exists(pt_fa):
        sec8_blocks.append(
            {
                "type": "text",
                "text": {
                    "title": "Extracted ptDNA sequence (FASTA)",
                    "content": _species_href(pt_fa, sp_root, species),
                    "as_markdown": False,
                },
            }
        )
    mt_fa = S("v6/0/polap-assemble/extract/oatk.mito.ctg.fasta")
    if _exists(mt_fa):
        sec8_blocks.append(
            {
                "type": "text",
                "text": {
                    "title": "Extracted mtDNA sequence (FASTA)",
                    "content": _species_href(mt_fa, sp_root, species),
                    "as_markdown": False,
                },
            }
        )
    sections.append({"id": "8", "title": "Oatk's pathfinder", "blocks": sec8_blocks})

    # 9. Polishing — final FASTA as text
    sec9_blocks: List[Dict[str, Any]] = []
    polished = S("v6/0/polap-assemble/polish-longshort/polished.fa")
    if _exists(polished):
        sec9_blocks.append(
            {
                "type": "text",
                "text": {
                    "title": "Polished organelle sequence (FASTA)",
                    "content": _species_href(polished, sp_root, species),
                    "as_markdown": False,
                },
            }
        )
    sections.append({"id": "9", "title": "Polishing", "blocks": sec9_blocks})

    # 10. Coverage (mt) — three PDF plots
    cov_blocks: List[Dict[str, Any]] = []

    f_cov_A = _figure_block_from_path(
        S("v6/0/polap-coverage/coverage/mt/plots/A_flatness_lines.pdf"),
        "Coverage: Flatness (lines)", sp_root, species
    )
    if f_cov_A:
        cov_blocks.append(f_cov_A)

    f_cov_B = _figure_block_from_path(
        S("v6/0/polap-coverage/coverage/mt/plots/B_uniformity_ecdf.pdf"),
        "Coverage: Uniformity ECDF", sp_root, species
    )
    if f_cov_B:
        cov_blocks.append(f_cov_B)

    f_cov_C = _figure_block_from_path(
        S("v6/0/polap-coverage/coverage/mt/plots/C_lorenz_gini.pdf"),
        "Coverage: Lorenz–Gini", sp_root, species
    )
    if f_cov_C:
        cov_blocks.append(f_cov_C)

    if cov_blocks:
        sections.append({
            "id": "10",
            "title": "Coverage (mt)",
            "blocks": cov_blocks
        })


    # 11. Performance section (summary-polap-assemble + timing-polap-assemble)
    run_root = os.path.join(sp_root, "v6", "0")
    p_summary_polap = os.path.join(run_root, "summary-polap-assemble.txt")
    p_timing_polap = os.path.join(run_root, "timing-polap-assemble.txt")

    perf_summary = _parse_summary_polap(p_summary_polap)
    perf_timing = _parse_timing_polap(p_timing_polap)

    perf_rows = [
        ["Elapsed time", perf_summary.get("elapsed_hms", "NA")],
        ["Net memory increase (KB)", perf_summary.get("net_increase_kb", "NA")],
        ["Disk used (GB)", perf_summary.get("disk_used_gb", "NA")],
        ["CPU model", perf_timing.get("cpu_model", "NA")],
        ["CPU cores", perf_timing.get("cpu_cores", "NA")],
        ["Hostname", perf_timing.get("hostname", "NA")],
    ]

    sections.append(
        {
            "id": "11",
            "title": "System & Performance (polap-assemble)",
            "blocks": [
                {
                    "type": "table",
                    "table": {
                        "title": "Run metrics",
                        "headers": ["Metric", "Value"],
                        "rows": perf_rows,
                    },
                }
            ],
        }
    )

    # Page info
    page = {
        "species": species,
        "title": f"polap-assemble report ({species})",
        "generated_at": _iso_now(),
        "base_dir": sp_root,  # path to species dir
    }

    # Top-level perf for tooling
    perf = {"summary": perf_summary, "timing": perf_timing}

    return {"page": page, "sections": sections, "perf": perf}


# ------------------------------ main ------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Create BLOCK JSON for polap-assemble report."
    )
    ap.add_argument(
        "--base-dir",
        required=True,
        help="Species base directory (e.g., Brassica_rapra)",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON (…/v6/0/polap-assemble/report/assemble-report.json)",
    )
    args = ap.parse_args()

    species_dir = os.path.abspath(args.base_dir)
    if not os.path.isdir(species_dir):
        sys.exit(f"[ERR] base dir not found: {species_dir}")

    doc = build_sections(species_dir)
    model = {
        "meta": {
            "generated_at": doc["page"]["generated_at"],
            "version": VERSION,
            "species_dir": doc["page"]["base_dir"],
        },
        "perf": doc["perf"],
        "page": doc["page"],
        "sections": doc["sections"],
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(model, f, indent=2)
    print(f"[OK] wrote BLOCK JSON: {os.path.abspath(args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
