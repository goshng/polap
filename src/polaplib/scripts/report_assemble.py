#!/usr/bin/env python3
# report_assemble.py
# Version: v1.0.0
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
# - This script reads files from a conventional run directory under <Species>/v6/0/polap-assemble
#   and produces a single JSON with the exact blocks you requested.
# - Hrefs are species-rooted (e.g., "Brassica_rapa/v6/0/...") so the renderer can convert to
#   links relative to the final HTML location.
#
import os
import sys
import csv
import json
import argparse
import datetime
from typing import Dict, Any, List, Optional, Tuple

# ------------------------------ helpers ------------------------------


def species_and_root(base_dir: str) -> Tuple[str, str]:
    """Return (species_name, species_abs_root)."""
    sp = os.path.basename(os.path.abspath(base_dir))
    return sp, os.path.abspath(base_dir)


def P(sp_root: str, rel: str) -> str:
    """Join species root and a species-relative path."""
    return os.path.join(sp_root, rel)


def species_href(abs_path: str, species_root: str, species: str) -> str:
    """
    Produce a href like '<Species>/...' for any absolute file path under the species root.
    Fallback to basename if trimming fails.
    """
    ap = os.path.abspath(abs_path).replace("\\", "/")
    root = os.path.abspath(species_root).replace("\\", "/")
    if ap.startswith(root + "/"):
        tail = ap[len(root) + 1 :]
        return f"{species}/{tail}"
    # Fallback: search for '/<Species>/' fragment
    i = ap.find("/" + species + "/")
    if i != -1:
        return ap[i + 1 :]
    return os.path.basename(ap)


def exists(path: str) -> bool:
    return path and os.path.exists(path)


def read_tsv_select_columns(
    path: str, keep_cols: List[str]
) -> Optional[Dict[str, Any]]:
    """
    Generic TSV reader selecting columns by header names (case-sensitive).
    Returns table dict: {title, headers, rows}.
    If the file has one data row (e.g., seqkit -T/-Ta), we output that row only.
    """
    if not exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd, None)
        if not hdr:
            return None
        hpos = {h: i for i, h in enumerate(hdr)}
        rows = []
        for r in rd:
            row = []
            for col in keep_cols:
                i = hpos.get(col)
                row.append(r[i] if i is not None and i < len(r) else "")
            rows.append(row)
    if not rows:
        # still produce empty with headers
        rows = [[]]
    return {"title": os.path.basename(path), "headers": keep_cols, "rows": rows}


def read_tsv_generic(path: str, max_rows: int = 2000) -> Optional[Dict[str, Any]]:
    """Read generic TSV into headers, rows."""
    if not exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd, None)
        if hdr is None:
            # Attempt fallback: split by whitespace, build generic headers
            f.seek(0)
            lines = [ln.strip() for ln in f if ln.strip()]
            if not lines:
                return None
            rows = [ln.split() for ln in lines[:max_rows]]
            width = max(len(r) for r in rows)
            headers = [f"col{i+1}" for i in range(width)]
            return {"title": os.path.basename(path), "headers": headers, "rows": rows}
        rows = []
        for r in rd:
            rows.append(r)
            if len(rows) >= max_rows:
                break
    return {"title": os.path.basename(path), "headers": hdr, "rows": rows}


def read_space_table(path: str, max_rows: int = 2000) -> Optional[Dict[str, Any]]:
    """
    Read a whitespace-delimited table (space/tabs), inferring headers from first row if it looks textual.
    """
    if not exists(path):
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
    # Heuristic: if first row contains alpha, treat it as header
    if any(any(ch.isalpha() for ch in cell) for cell in rows[0]):
        headers = rows[0]
        data = rows[1:]
    else:
        width = max(len(r) for r in rows)
        headers = [f"col{i+1}" for i in range(width)]
        data = rows
    return {"title": os.path.basename(path), "headers": headers, "rows": data}


def text_block_from_file(
    path: str, title: str, limit_kb: int = 256
) -> Optional[Dict[str, Any]]:
    """Read small text file into a text block (truncated message if too large)."""
    if not exists(path):
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


def figure_block_from_path(
    path: str, label: str, species_root: str, species: str
) -> Optional[Dict[str, Any]]:
    if not exists(path):
        return None
    href = species_href(path, species_root, species)
    return {"type": "figure", "label": label, "href": href}


def table_block(table: Dict[str, Any], title: Optional[str] = None) -> Dict[str, Any]:
    if title:
        table = dict(table)
        table["title"] = title
    return {"type": "table", "table": table}


# ------------------------------ sections builder ------------------------------


def build_sections(species_dir: str) -> Dict[str, Any]:
    """
    Build the full BLOCK JSON composed of sections and blocks as requested.
    """
    species, sp_root = species_and_root(species_dir)

    # Convenience helper to join species-root relative strings
    def S(rel: str) -> str:
        return P(sp_root, rel)

    sections: List[Dict[str, Any]] = []

    # --- 1. Seed PT and MT reads: seqkit table (3 columns) ---
    sec1_blocks: List[Dict[str, Any]] = []
    pt0_seqkit = S("v6/0/polap-assemble/annotate-read-pt/pt0.fq.seqkit.stats.ta.txt")
    t1 = read_tsv_select_columns(pt0_seqkit, ["num_seqs", "sum_len", "AvgQual"])
    if t1:
        sec1_blocks.append(table_block(t1, "PT reads (seqkit, pt0)"))
    sections.append({"id": "1", "title": "Seed PT and MT reads", "blocks": sec1_blocks})

    # --- 2. ptDNA assembly: figures + space-delim table ---
    sec2_blocks: List[Dict[str, Any]] = []
    fig_pt0 = figure_block_from_path(
        S("v6/0/polap-assemble/annotate-read-pt/pt/30-contigger/graph_final.png"),
        "PT assembly graph (pt0)",
        sp_root,
        species,
    )
    if fig_pt0:
        sec2_blocks.append(fig_pt0)
    fig_pt1 = figure_block_from_path(
        S("v6/0/polap-assemble/annotate-read-pt/pt1/assembly_graph.png"),
        "PT assembly graph (pt1)",
        sp_root,
        species,
    )
    if fig_pt1:
        sec2_blocks.append(fig_pt1)

    pt1_anno = S(
        "v6/0/polap-assemble/annotate-read-pt/pt1/pt-contig-annotation-depth-table.txt"
    )
    t2 = read_space_table(pt1_anno)
    if t2:
        sec2_blocks.append(table_block(t2, "PT annotation table (pt1)"))
    sections.append({"id": "2", "title": "ptDNA assembly", "blocks": sec2_blocks})

    # --- 3. Filter out PT reads: pt_thresh diag table + seqkit 3-col ---
    sec3_blocks: List[Dict[str, Any]] = []
    diag_tsv = S("v6/0/polap-assemble/mtseed/pt_thresh.diag.tsv")
    t3 = read_tsv_generic(diag_tsv)
    if t3:
        sec3_blocks.append(
            table_block(t3, "PT read selection cutoff (pt_thresh.diag.tsv)")
        )

    nonpt_tsv = S("v6/0/polap-assemble/mtseed/reads.nonpt.fq.gz.seqkit.stats.ta.tsv")
    t4 = read_tsv_select_columns(nonpt_tsv, ["num_seqs", "sum_len", "AvgQual"])
    if t4:
        sec3_blocks.append(table_block(t4, "Non-PT reads (after PT filtering)"))
    sections.append({"id": "3", "title": "Filter out PT reads", "blocks": sec3_blocks})

    # --- 4. Compute read overlapness: two PDFs ---
    sec4_blocks: List[Dict[str, Any]] = []
    pdf_wdeg = figure_block_from_path(
        S("v6/0/polap-assemble/mtseed/03-allvsall/04-qc/scan_wdegree_hist.pdf"),
        "Weighted degree distribution",
        sp_root,
        species,
    )
    if pdf_wdeg:
        sec4_blocks.append(pdf_wdeg)

    pdf_cum = figure_block_from_path(
        S("v6/0/polap-assemble/mtseed/03-allvsall/04-qc/scan_cum_wdegree.pdf"),
        "Cumulative weighted degree (threshold)",
        sp_root,
        species,
    )
    if pdf_cum:
        sec4_blocks.append(pdf_cum)
    sections.append(
        {"id": "4", "title": "Compute read overlapness", "blocks": sec4_blocks}
    )

    # --- 5. Filter out NT reads: TSV table ---
    sec5_blocks: List[Dict[str, Any]] = []
    t5 = read_tsv_generic(
        S("v6/0/polap-assemble/mtseed/05-round/threshold_from_nuclear.tsv")
    )
    if t5:
        sec5_blocks.append(
            table_block(t5, "Weighted-degree threshold to filter NT reads")
        )
    sections.append({"id": "5", "title": "Filter out NT reads", "blocks": sec5_blocks})

    # --- 6. Assemble seed contigs using miniasm or mtDNA assembly mt0: figure ---
    sec6_blocks: List[Dict[str, Any]] = []
    f_mt0 = figure_block_from_path(
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

    # --- 7. mtDNA assembly (7.1..7.3) ---
    subsecs: List[Dict[str, Any]] = []
    for name, idx in (("mt1", "7.1"), ("mt2", "7.2"), ("mt3", "7.3")):
        base_rel = f"v6/0/polap-assemble/mtseed/{name}"
        blocks: List[Dict[str, Any]] = []

        # Contig total length (text)
        tlen = text_block_from_file(
            S(f"{base_rel}/01-contig/contig_total_length.txt"), "Contig length (bp)"
        )
        if tlen:
            blocks.append(tlen)

        # Assembly graph PNG
        f_asm = figure_block_from_path(
            S(f"{base_rel}/assembly_graph.png"),
            f"Assembly graph ({name})",
            sp_root,
            species,
        )
        if f_asm:
            blocks.append(f_asm)

        # MT annotation table (space-delimited)
        t_mt = read_space_table(S(f"{base_rel}/contig-annotation-depth-table.txt"))
        if t_mt:
            blocks.append(table_block(t_mt, f"MT annotation table ({name})"))

        # Seed contig names text (for mt1 → mt2; for mt2 → mt3 if present)
        next_map = {"mt1": "mt2", "mt2": "mt3"}
        if name in next_map:
            seed_names_path = S(f"{base_rel}/mt.contig.name-{next_map[name]}")
            t_seed = text_block_from_file(
                seed_names_path, f"Seed contig names → {next_map[name]}"
            )
            if t_seed:
                blocks.append(t_seed)

        if blocks:
            subsecs.append(
                {"id": idx, "title": f"mtDNA assembly {name}", "blocks": blocks}
            )

    sections.append(
        {"id": "7", "title": "mtDNA assembly", "blocks": [], "subsections": subsecs}
    )

    # --- 8. Oatk's pathfinder (text links to FASTA) ---
    sec8_blocks: List[Dict[str, Any]] = []
    pt_fa = S("v6/0/polap-assemble/extract/oatk.pltd.ctg.fasta")
    if exists(pt_fa):
        sec8_blocks.append(
            {
                "type": "text",
                "text": {
                    "title": "Extracted ptDNA sequence (FASTA)",
                    "content": species_href(pt_fa, sp_root, species),
                    "as_markdown": False,
                },
            }
        )
    mt_fa = S("v6/0/polap-assemble/extract/oatk.mito.ctg.fasta")
    if exists(mt_fa):
        sec8_blocks.append(
            {
                "type": "text",
                "text": {
                    "title": "Extracted mtDNA sequence (FASTA)",
                    "content": species_href(mt_fa, sp_root, species),
                    "as_markdown": False,
                },
            }
        )
    sections.append({"id": "8", "title": "Oatk's pathfinder", "blocks": sec8_blocks})

    # --- 9. Polishing (final polished.fa path as text) ---
    sec9_blocks: List[Dict[str, Any]] = []
    polished = S("v6/0/polap-assemble/polish-longshort/polished.fa")
    if exists(polished):
        sec9_blocks.append(
            {
                "type": "text",
                "text": {
                    "title": "Polished organelle sequence (FASTA)",
                    "content": species_href(polished, sp_root, species),
                    "as_markdown": False,
                },
            }
        )
    sections.append({"id": "9", "title": "Polishing", "blocks": sec9_blocks})

    # Page info
    page = {
        "species": species,
        "title": f"polap-assemble report ({species})",
        "generated_at": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "base_dir": species,
    }

    return {"page": page, "sections": sections}


# ------------------------------ main ------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Create BLOCK JSON for polap-assemble report."
    )
    ap.add_argument(
        "--base-dir", required=True, help="Species base directory (e.g., Brassica_rapa)"
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output JSON path (e.g., .../report/assemble-report.json)",
    )
    args = ap.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    if not os.path.isdir(base_dir):
        sys.exit(f"[ERR] base dir not found: {base_dir}")

    doc = build_sections(base_dir)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fo:
        json.dump(doc, fo, indent=2)
    print(f"[OK] wrote BLOCK JSON: {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
