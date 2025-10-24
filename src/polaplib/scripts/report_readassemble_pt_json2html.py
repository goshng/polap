#!/usr/bin/env python3
# scripts/report_readassemble_pt_json2html.py
# Version: v0.3.1
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Converts pt1-report.json -> compact HTML.
# Displays Total bases in both Mb and Gb units.
#
import os, re, sys, json, html, argparse


def esc(x):
    return html.escape("" if x is None else str(x))


def fmt_gb(x):
    try:
        return f"{float(x)/1e9:.2f} Gb"
    except Exception:
        return ""


def fmt_mb(x):
    try:
        return f"{float(x)/1e6:.1f} Mb"
    except Exception:
        return ""


def fmt_kb(x):
    try:
        return f"{float(x)/1e3:.1f} kb"
    except Exception:
        return ""


def find_species_name(json_path):
    """Guess species folder name, e.g. /.../Anthoceros_angustus/v6/0/..."""
    m = re.search(r"/([^/]+)/v\d+/", json_path)
    return m.group(1) if m else None


def rel_from_species(path, species):
    if not path or not species:
        return path
    i = path.find(species)
    return path[i:] if i != -1 else path


def file_exists_from_json_dir(json_path, relpath):
    if not relpath:
        return False
    species = find_species_name(json_path)
    if not species:
        return os.path.exists(relpath)
    cur = os.path.dirname(os.path.abspath(json_path))
    found = None
    for _ in range(10):
        if os.path.basename(cur) == species:
            found = cur
            break
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    if not found:
        return os.path.exists(relpath)
    abs_path = (
        os.path.join(found, relpath.split(species + os.sep, 1)[-1])
        if relpath.startswith(species + os.sep)
        else os.path.join(found, relpath)
    )
    return os.path.exists(abs_path)


def build_html(data, json_path):
    species = find_species_name(json_path)
    title = f"PT-1 Report ({species})" if species else "PT-1 Report"
    base_dir = data.get("base_dir", "")

    kpis = data.get("kpis", {})
    pt_reads = kpis.get("pt_reads", {})
    seqkit = kpis.get("pt_input_seqkit", {})
    qc = data.get("qc", {})

    # Input PT reads
    pt_count = pt_reads.get("count", "")
    pt_sum = pt_reads.get("total_bases", "")
    pt_bases_mb = fmt_mb(pt_sum)
    pt_bases_gb = fmt_gb(pt_sum)

    # Filtered PT reads (seqkit)
    in_num = seqkit.get("num_seqs", "")
    in_sum = seqkit.get("sum_len", "")
    in_sum_mb = fmt_mb(in_sum)
    in_sum_gb = fmt_gb(in_sum)
    in_avg = fmt_kb(seqkit.get("avg_len", ""))
    in_n50 = fmt_kb(seqkit.get("N50", ""))
    in_aq = seqkit.get("AvgQual", "")

    # QC images
    graph_png = rel_from_species(qc.get("graph_png", ""), species)
    pre_png = rel_from_species(qc.get("pre_cov_png", ""), species)
    post_png = rel_from_species(qc.get("post_cov_png", ""), species)

    css = """
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,Arial,sans-serif;margin:20px;line-height:1.45}
h1{font-size:1.5rem;margin-bottom:0.5rem}
h2{font-size:1.2rem;margin-top:1.4rem}
.card{border:1px solid #ddd;border-radius:8px;padding:12px;margin:10px 0;background:#fafafa}
.kv{display:grid;grid-template-columns:240px 1fr;gap:6px 12px;max-width:760px}
.kv div.k{font-weight:600}
img{max-width:100%;height:auto;border:1px solid #eee;border-radius:6px;margin:4px 0;background:#fff}
small{color:#666}
"""
    lines = []
    lines.append("<!DOCTYPE html>")
    lines.append("<html lang='en'><head><meta charset='utf-8'>")
    lines.append(f"<title>{esc(title)}</title>")
    lines.append(f"<style>{css}</style></head><body>")
    lines.append(f"<h1>{esc(title)}</h1>")
    if base_dir:
        lines.append(
            f"<div class='card'><small>Base dir: {esc(rel_from_species(base_dir, species))}</small></div>"
        )

    # 1) Input PT reads
    lines.append("<h2>1. Input PT reads</h2>")
    lines.append("<div class='card kv'>")
    lines.append("<div class='k'>Reads labeled as PT</div>")
    lines.append(f"<div>{esc(pt_count)}</div>")
    lines.append("<div class='k'>Total PT bases</div>")
    lines.append(f"<div>{esc(pt_bases_mb)} ({esc(pt_bases_gb)})</div>")
    lines.append("</div>")

    # 2) Filtered PT reads (seqkit)
    lines.append("<h2>2. Filtered PT reads (seqkit)</h2>")
    lines.append("<div class='card kv'>")
    lines.append(f"<div class='k'>Number of reads</div><div>{esc(in_num)}</div>")
    lines.append(
        f"<div class='k'>Total bases</div><div>{esc(in_sum_mb)} ({esc(in_sum_gb)})</div>"
    )
    lines.append(f"<div class='k'>Average read length</div><div>{esc(in_avg)}</div>")
    lines.append(f"<div class='k'>N50</div><div>{esc(in_n50)}</div>")
    lines.append(f"<div class='k'>Average quality</div><div>{esc(in_aq)}</div>")
    lines.append("</div>")

    # 3) QC images
    lines.append("<h2>3. QC</h2>")
    lines.append("<div class='card'>")

    def embed(label, relp):
        if not relp:
            return
        if file_exists_from_json_dir(json_path, relp):
            lines.append(f"<div><b>{esc(label)}</b></div>")
            lines.append(f"<img src='{esc(relp)}' alt='{esc(label)}'>")

    embed("Final assembly graph", graph_png)
    embed("Coverage (pre-recruit)", pre_png)
    embed("Coverage (post-recruit)", post_png)

    if not any((graph_png, pre_png, post_png)):
        lines.append("<small>No QC images available.</small>")

    lines.append("</div></body></html>")
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
