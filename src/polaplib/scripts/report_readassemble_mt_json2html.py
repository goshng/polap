#!/usr/bin/env python3
# scripts/report_readassemble_mt_json2html.py
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Convert mt-report.json into a compact, species-rooted HTML summary.
# Sections:
#   - PT-2 panel (IR/isomers/doubles sanity)
#   - PT-3 mapping (PT removal, non-PT pool, thresholds)
#   - MT-1 all-vs-all overlaps (embed QC PDFs)
#   - MT-2 busco (selected read count)
#   - MT-3 miniasm (seed stats)
#   - MT-4 iterations (mt0..mt3 table with stats & gene hints)
#
import os, re, sys, json, html, argparse

def esc(x): return html.escape("" if x is None else str(x))
def fmt_gb(x):
    try: return f"{float(x)/1e9:.2f} Gb"
    except: return ""
def fmt_mb(x):
    try: return f"{float(x)/1e6:.1f} Mb"
    except: return ""
def fmt_num(x):
    try:
        v=int(x)
        return f"{v:,d}"
    except Exception:
        return esc(x)

def find_species(json_path):
    m=re.search(r"/([^/]+)/v\d+/", json_path)
    return m.group(1) if m else None

def rel_from_species(path, species):
    if not path or not species: return path
    i=path.find(species)
    return path[i:] if i!=-1 else path

def file_exists_from_json_dir(json_path, relpath):
    if not relpath: return False
    species = find_species(json_path)
    if not species: return os.path.exists(relpath)
    cur = os.path.dirname(os.path.abspath(json_path))
    found=None
    for _ in range(12):
        if os.path.basename(cur)==species:
            found=cur; break
        parent=os.path.dirname(cur)
        if parent==cur: break
        cur=parent
    if not found: return os.path.exists(relpath)
    tail = relpath.split(species+os.sep,1)[-1] if relpath.startswith(species+os.sep) else relpath
    abs_path = os.path.join(found, tail)
    return os.path.exists(abs_path)

def build_html(data, json_path):
    species = find_species(json_path)
    title = f"MT Report ({species})" if species else "MT Report"
    steps = data.get("steps", {})

    css = """
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,Arial,sans-serif;margin:20px;line-height:1.45}
h1{font-size:1.6rem;margin-bottom:0.6rem}
h2{font-size:1.2rem;margin-top:1.4rem}
.card{border:1px solid #ddd;border-radius:8px;padding:12px;margin:10px 0;background:#fafafa}
.kv{display:grid;grid-template-columns:260px 1fr;gap:6px 14px;max-width:860px}
.kv div.k{font-weight:600}
table{border-collapse:collapse;width:100%;max-width:980px}
th,td{border:1px solid #e1e1e1;padding:6px 8px;text-align:left}
th{background:#f5f5f5}
object.pdf{width:100%; height:520px; border:1px solid #eee; border-radius:6px; background:#fff}
small{color:#666}
"""

    out=[]
    out.append("<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>")
    out.append(f"<title>{esc(title)}</title><style>{css}</style></head><body>")
    out.append(f"<h1>{esc(title)}</h1>")

    # PT-2 panel
    panel=steps.get("panel",{})
    irsum=panel.get("ir",{}).get("summary",[])
    iso=panel.get("isomers",{})
    out.append("<h2>PT-2: plastid reference preparation (isomers & doubles)</h2>")
    out.append("<div class='card kv'>")
    out.append("<div class='k'>Isomer A length</div><div>{}</div>".format(fmt_num(iso.get("A_len",""))))
    out.append("<div class='k'>Isomer A double</div><div>{} (2× check: {})</div>".format(fmt_num(iso.get("doubleA_len","")), esc(iso.get("doubleA_ok",""))))
    out.append("<div class='k'>Isomer B length</div><div>{}</div>".format(fmt_num(iso.get("B_len",""))))
    out.append("<div class='k'>Isomer B double</div><div>{} (2× check: {})</div>".format(fmt_num(iso.get("doubleB_len","")), esc(iso.get("doubleB_ok",""))))
    out.append("</div>")
    if irsum:
        out.append("<div class='card'><small>IR summary (first rows from ir.pick.tsv)</small><table><thead><tr>")
        hdr=list(irsum[0].keys())
        out.append("".join(f"<th>{esc(h)}</th>" for h in hdr))
        out.append("</tr></thead><tbody>")
        for row in irsum:
            out.append("<tr>"+"".join(f"<td>{esc(row.get(h,''))}</td>" for h in hdr)+"</tr>")
        out.append("</tbody></table></div>")

    # PT-3 mapping
    mapstep=steps.get("map",{})
    out.append("<h2>PT-3: map to doubled PT isomers; PT removal -> non-PT pool</h2>")
    out.append("<div class='card kv'>")
    out.append("<div class='k'>All read IDs</div><div>{} ({})</div>".format(fmt_num(mapstep.get("all_ids",{}).get("count","")), esc(mapstep.get("all_ids",{}).get("file",""))))
    out.append("<div class='k'>PT read IDs</div><div>{} ({})</div>".format(fmt_num(mapstep.get("pt_ids",{}).get("count","")), esc(mapstep.get("pt_ids",{}).get("file",""))))
    out.append("<div class='k'>PT removal fraction</div><div>{}</div>".format(esc(mapstep.get("removed_pt_frac",""))))
    out.append("<div class='k'>Kept non-PT IDs</div><div>{} ({})</div>".format(fmt_num(mapstep.get("nonpt_ids",{}).get("count","")), esc(mapstep.get("nonpt_ids",{}).get("file",""))))
    np = mapstep.get("nonpt_seqkit",{})
    out.append("<div class='k'>Non-PT sum_len</div><div>{} ({})</div>".format(fmt_mb(np.get("sum_len","")), fmt_gb(np.get("sum_len",""))))
    out.append("<div class='k'>Non-PT num_seqs</div><div>{}</div>".format(fmt_num(np.get("num_seqs",""))))
    out.append("<div class='k'>Threshold vars</div><div>{}</div>".format(esc(mapstep.get("thresholds",{}).get("vars",""))))
    out.append("</div>")

    # MT-1 all-vs-all
    allvsall=steps.get("allvsall",{})
    out.append("<h2>MT-1: shard-local all-vs-all overlaps</h2>")
    out.append("<div class='card'>")
    def embed_pdf(label, relp):
        if not relp: return
        # rewrite to species-rooted for portability
        rel = rel_from_species(relp, species)
        if file_exists_from_json_dir(json_path, rel):
            out.append(f"<div><b>{esc(label)}</b></div>")
            out.append(f"<object class='pdf' data='{esc(rel)}' type='application/pdf'>")
            out.append(f"<p><a href='{esc(rel)}'>{esc(rel)}</a></p></object>")
    embed_pdf("Degree histogram", allvsall.get("qc",{}).get("degree_hist_pdf",""))
    embed_pdf("Weighted degree histogram", allvsall.get("qc",{}).get("wdegree_hist_pdf",""))
    embed_pdf("Cumulative weighted degree", allvsall.get("qc",{}).get("cum_wdegree_pdf",""))
    if not any(allvsall.get("qc",{}).values()):
        out.append("<small>No QC PDFs available.</small>")
    out.append("</div>")

    # MT-2 busco
    busco=steps.get("busco",{})
    out.append("<h2>MT-2: BUSCO-calibrated selection (optional)</h2>")
    out.append("<div class='card kv'>")
    out.append("<div class='k'>Selected IDs</div><div>{} ({})</div>".format(
        fmt_num(busco.get("selected_ids",{}).get("count","")),
        esc(busco.get("selected_ids",{}).get("file",""))
    ))
    out.append("</div>")

    # MT-3 miniasm
    miniasm=steps.get("miniasm",{})
    s=miniasm.get("stats",{})
    out.append("<h2>MT-3: miniasm seeding</h2>")
    out.append("<div class='card kv'>")
    out.append("<div class='k'>Seeds fasta</div><div>{}</div>".format(esc(miniasm.get("seeds_fa",""))))
    out.append("<div class='k'>Miniasm GFA</div><div>{}</div>".format(esc(miniasm.get("gfa",""))))
    out.append("<div class='k'>Unitigs (S)</div><div>{}</div>".format(fmt_num(s.get("S",""))))
    out.append("<div class='k'>Total bp</div><div>{} ({})</div>".format(fmt_mb(s.get("total_bp","")), fmt_gb(s.get("total_bp",""))))
    out.append("<div class='k'>N50 / Max</div><div>{} / {}</div>".format(fmt_num(s.get("N50","")), fmt_num(s.get("max_bp",""))))
    out.append("</div>")

    # MT-4 iterations table
    its = steps.get("iterations",[])
    out.append("<h2>MT iterations (mt0 -> mt3)</h2>")
    if its:
        out.append("<div class='card'><table><thead><tr>")
        out.append("<th>Iter</th><th>#contigs</th><th>Total</th><th>N50</th><th>Max</th><th>MT genes</th><th>Annotated MT length</th><th>Graph</th></tr></thead><tbody>")
        for it in its:
            name = it.get("name","")
            stats= it.get("stats",{})
            genes= it.get("genes",{})
            png  = it.get("paths",{}).get("graph_png","")
            relpng = rel_from_species(png, species) if png else ""
            pngcell = f"<a href='{esc(relpng)}'>{esc(os.path.basename(relpng))}</a>" if relpng and file_exists_from_json_dir(json_path, relpng) else ""
            out.append("<tr>" +
                f"<td>{esc(name)}</td>" +
                f"<td>{fmt_num(stats.get('n',''))}</td>" +
                f"<td>{fmt_mb(stats.get('total_bp',''))} ({fmt_gb(stats.get('total_bp',''))})</td>" +
                f"<td>{fmt_num(stats.get('N50',''))}</td>" +
                f"<td>{fmt_num(stats.get('max_bp',''))}</td>" +
                f"<td>{fmt_num(genes.get('mt_gene_count',''))}</td>" +
                f"<td>{fmt_mb(genes.get('annotated_mt_length',''))}</td>" +
                f"<td>{pngcell}</td>" +
                "</tr>")
        out.append("</tbody></table></div>")
    else:
        out.append("<div class='card'><small>No iterations found.</small></div>")

    out.append("</body></html>")
    return "\n".join(out)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--html", required=True)
    args=ap.parse_args()

    with open(args.json, encoding="utf-8") as f:
        data=json.load(f)
    html_doc = build_html(data, os.path.abspath(args.json))
    os.makedirs(os.path.dirname(os.path.abspath(args.html)), exist_ok=True)
    with open(args.html,"w",encoding="utf-8") as f:
        f.write(html_doc)

if __name__=="__main__":
    sys.exit(main())
