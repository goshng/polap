#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-ensemble-report.py  v0.3.2  (combined mt + pt + Bandage + sizes)

Adds:
  - Size panel for mt.fastq.gz / pt.fastq.gz (if present)
"""

import argparse, os, base64, html


def load_ids(path):
    S = set()
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                t = line.strip().split()
                if t:
                    S.add(t[0])
    return S


def read_text_head(path, n=200):
    if not path or not os.path.exists(path):
        return ""
    try:
        with open(path) as f:
            return (
                "<pre class='small'>"
                + html.escape("".join(f.readlines()[:n]))
                + "</pre>"
            )
    except Exception as e:
        return f"<pre class='small'>(error reading {html.escape(path)}: {html.escape(str(e))})</pre>"


def embed_png(path, title=None):
    if not path or not os.path.exists(path):
        return ""
    try:
        with open(path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode("ascii")
        caption = (
            f"<div class='small'>{html.escape(os.path.basename(path))}</div>"
            if not title
            else f"<div class='small'>{html.escape(title)}</div>"
        )
        return f'{caption}<img alt="{html.escape(os.path.basename(path))}" src="data:image/png;base64,{b64}" style="max-width:100%;height:auto;border:1px solid #ccc;padding:4px;">'
    except Exception as e:
        return f"<div class='small'>(error embedding {html.escape(path)}: {html.escape(str(e))})</div>"


def head(title):
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>{html.escape(title)}</title>
<style>
body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,Arial,sans-serif;margin:24px;}}
h1,h2{{margin-top:1.1em}}
table{{border-collapse:collapse;margin:10px 0}}
td,th{{border:1px solid #ddd;padding:6px 10px;font-size:90%}}
code{{background:#f7f7f7;padding:2px 4px;border:1px solid #eee}}
.small{{font-size:85%;color:#555}}
.block{{display:block;margin:8px 0}}
hr{{border:0;border-top:1px solid #eee;margin:20px 0}}
.section-title{{margin-top:24px}}
</style></head><body>
"""


def tail():
    return "</body></html>"


def count_if(path):
    return len(load_ids(path)) if path and os.path.exists(path) else 0


def human(n):
    # n in bytes -> human-readable
    if n is None:
        return "(missing)"
    for unit in ["B", "KB", "MB", "GB", "TB", "PB"]:
        if n < 1024.0:
            return f"{n:3.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} PB"


def organelle_block(prefix_cov, prefix_ppr, prefix_ens, label="MT"):
    lab = label.lower()
    ids_cov = f"{prefix_cov}.{lab}.ids"
    ids_ppr = f"{prefix_ppr}.ids"
    ids_ppr_nuc = f"{prefix_ppr}.nuclear.ids"
    ppr_scores = f"{prefix_ppr}.scores.tsv"
    ppr_gstats = f"{prefix_ppr}.graph.stats.txt"

    ids_ens = f"{prefix_ens}.ids"
    ids_ens_nuc = f"{prefix_ens}.nuclear.ids"
    ens_report = f"{prefix_ens}.report.tsv"

    rows = []
    rows.append(f"<h2 class='section-title'>{html.escape(label)} summary</h2>")
    rows.append("<table><tr><th>Set</th><th>Count</th><th>Path</th></tr>")
    rows.append(
        f"<tr><td>Coverage {lab}</td><td>{count_if(ids_cov)}</td><td><code>{html.escape(ids_cov) if os.path.exists(ids_cov) else '(missing)'}</code></td></tr>"
    )
    rows.append(
        f"<tr><td>PPR {lab}</td><td>{count_if(ids_ppr)}</td><td><code>{html.escape(ids_ppr) if os.path.exists(ids_ppr) else '(missing)'}</code></td></tr>"
    )
    rows.append(
        f"<tr><td>Ensemble {lab}</td><td>{count_if(ids_ens)}</td><td><code>{html.escape(ids_ens) if os.path.exists(ids_ens) else '(missing)'}</code></td></tr>"
    )
    rows.append("</table>")

    # PPR details
    gtxt = read_text_head(ppr_gstats, n=200)
    stxt = read_text_head(ppr_scores, n=50)
    if gtxt or stxt:
        rows.append("<h3>PPR details</h3>")
        if gtxt:
            rows.append(gtxt)
        if stxt:
            rows.append("<p class='small'><b>Top scores (head)</b></p>" + stxt)

    # Ensemble report
    etxt = read_text_head(ens_report, n=200)
    if etxt:
        rows.append("<h3>Ensemble overlap</h3>")
        rows.append(etxt)

    # Bandage organelle image (optional)
    band_img = embed_png(
        f"{args.prefix}.bandage.{lab}.png", title=f"Bandage {label} graph"
    )
    if band_img:
        rows.append("<h3>Bandage</h3>")
        rows.append(band_img)

    return "\n".join(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--quickview", required=True)
    ap.add_argument("--xlda-prefix", required=True, dest="cov_prefix")
    ap.add_argument("--ppr-prefix", required=True, dest="ppr_mt_prefix")
    ap.add_argument("--ensemble-prefix", required=True, dest="ens_mt_prefix")
    ap.add_argument("--ppr-pt-prefix", dest="ppr_pt_prefix")
    ap.add_argument("--ensemble-pt-prefix", dest="ens_pt_prefix")
    ap.add_argument("--title", default="POLAP organelle selection report (mt + pt)")
    ap.add_argument("-o", "--out", required=True, dest="out_html")
    global args
    args = ap.parse_args()

    parts = [head(args.title), f"<h1>{html.escape(args.title)}</h1>"]

    # Coverage / inputs
    parts.append("<h2>Coverage & inputs</h2>")
    parts.append("<ul>")
    parts.append(
        f"<li>Quickview TSV: <code>{html.escape(args.quickview)}</code> {'(missing)' if not os.path.exists(args.quickview) else ''}</li>"
    )
    parts.append(
        f"<li>Coverage prefix: <code>{html.escape(args.cov_prefix)}</code></li>"
    )
    parts.append(
        f"<li>MT PPR prefix: <code>{html.escape(args.ppr_mt_prefix)}</code></li>"
    )
    parts.append(
        f"<li>MT ensemble prefix: <code>{html.escape(args.ens_mt_prefix)}</code></li>"
    )
    if args.ppr_pt_prefix or args.ens_pt_prefix:
        parts.append(
            f"<li>PT PPR prefix: <code>{html.escape(args.ppr_pt_prefix or '(none)')}</code></li>"
        )
        parts.append(
            f"<li>PT ensemble prefix: <code>{html.escape(args.ens_pt_prefix or '(none)')}</code></li>"
        )
    parts.append("</ul>")

    # Coverage histogram & overview bandage
    cov_hist = embed_png(f"{args.cov_prefix}.x_hist.png", title="Coverage X histogram")
    band_over = embed_png(
        f"{args.prefix}.bandage.overview.png", title="Bandage overview"
    )
    if cov_hist or band_over:
        parts.append("<h3>Graph overview</h3>")
        if cov_hist:
            parts.append(cov_hist)
        if band_over:
            parts.append(band_over)

    # Size panel (mt/pt fastqs)
    sizes_rows = []
    mt_fastq = f"{args.prefix}.mt.fastq.gz"
    pt_fastq = f"{args.prefix}.pt.fastq.gz"
    mt_size = os.path.getsize(mt_fastq) if os.path.exists(mt_fastq) else None
    pt_size = os.path.getsize(pt_fastq) if os.path.exists(pt_fastq) else None
    if mt_size is not None or pt_size is not None:
        sizes_rows.append("<h2>FASTQ sizes</h2>")
        sizes_rows.append("<table><tr><th>File</th><th>Size</th><th>Path</th></tr>")
        sizes_rows.append(
            f"<tr><td>mt.fastq.gz</td><td>{human(mt_size) if mt_size is not None else '(missing)'}</td><td><code>{html.escape(mt_fastq) if mt_size is not None else ''}</code></td></tr>"
        )
        sizes_rows.append(
            f"<tr><td>pt.fastq.gz</td><td>{human(pt_size) if pt_size is not None else '(missing)'}</td><td><code>{html.escape(pt_fastq) if pt_size is not None else ''}</code></td></tr>"
        )
        sizes_rows.append("</table>")
        parts.extend(sizes_rows)

    # MT block
    parts.append(
        organelle_block(
            args.cov_prefix, args.ppr_mt_prefix, args.ens_mt_prefix, label="MT"
        )
    )

    # PT block (optional)
    if args.ppr_pt_prefix and args.ens_pt_prefix:
        if (
            os.path.exists(f"{args.cov_prefix}.pt.ids")
            or os.path.exists(f"{args.ppr_pt_prefix}.ids")
            or os.path.exists(f"{args.ens_pt_prefix}.ids")
        ):
            parts.append("<hr/>")
            parts.append(
                organelle_block(
                    args.cov_prefix, args.ppr_pt_prefix, args.ens_pt_prefix, label="PT"
                )
            )

    parts.append(tail())
    os.makedirs(os.path.dirname(args.out_html) or ".", exist_ok=True)
    with open(args.out_html, "w") as w:
        w.write("".join(parts))
    print(f"[report] {args.out_html}")


if __name__ == "__main__":
    import os

    main()
