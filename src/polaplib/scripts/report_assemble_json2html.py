#!/usr/bin/env python3
# report_assemble_json2html.py
# Version: v0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Render JSON from report_assemble_from_base.py to a navigable HTML.
# - TOC for sections and 7.x subsections
# - Displays item count, total size, time range
# - Shows a table of entries; links/inline previews (PNG/JPG/SVG/PDF)
# - Path display always begins at '<Species>/...'
#
import os
import re
import json
import html
import argparse
from typing import Optional, Dict, Any, List

IMG_EXTS = {".png", ".jpg", ".jpeg", ".svg"}
PDF_EXTS = {".pdf"}


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def human_bytes(n) -> str:
    try:
        n = float(n)
    except Exception:
        return ""
    if n >= 1e12:
        return f"{n/1e12:.2f} TB"
    if n >= 1e9:
        return f"{n/1e9:.2f} GB"
    if n >= 1e6:
        return f"{n/1e6:.2f} MB"
    if n >= 1e3:
        return f"{n/1e3:.2f} KB"
    return f"{int(n)} B"


def species_root_from_output(html_out: str, species: Optional[str]) -> Optional[str]:
    """
    Walk upwards from the HTML output dir to locate a folder literally named <species>.
    """
    if not species:
        return None
    cur = os.path.dirname(os.path.abspath(html_out))
    seen = set()
    for _ in range(64):
        if os.path.basename(cur) == species:
            return cur
        nxt = os.path.dirname(cur)
        if nxt == cur or nxt in seen:
            break
        seen.add(nxt)
        cur = nxt
    return None


def make_href_and_exists(
    display_path: str, html_out: str, species: Optional[str]
) -> (str, bool, Optional[str]):
    """
    The JSON 'path' is '<Species>/...'. Resolve that to a filesystem path by
    finding <Species> above the HTML, then compute a relative href.
    """
    if not display_path:
        return "", False, None
    html_dir = os.path.dirname(os.path.abspath(html_out))
    species_root = species_root_from_output(html_out, species)
    abs_target = None
    if species_root and species and display_path.startswith(species + "/"):
        tail = display_path[len(species) + 1 :]
        abs_target = os.path.join(species_root, tail)
    else:
        # Fallback: treat as relative to HTML dir
        abs_target = os.path.join(html_dir, display_path)

    href = os.path.relpath(abs_target, start=html_dir)
    exists = os.path.exists(abs_target)
    return href.replace("\\", "/"), exists, abs_target


def section_html(sec: Dict[str, Any], html_out: str, species: Optional[str]) -> str:
    sid = sec.get("id", "")
    title = sec.get("title", "")
    stats = sec.get("stats", {}) or {}
    entries = sec.get("entries", []) or []
    subs = sec.get("subsections", []) or []

    parts: List[str] = []
    anchor = f"sec-{sid.replace('.', '-')}"
    parts.append(f"<h2 id='{esc(anchor)}'>{esc(sid)}. {esc(title)}</h2>")
    parts.append("<div class='card kv'>")
    parts.append(f"<div class='k'>Items</div><div>{len(entries):,d}</div>")
    parts.append(
        f"<div class='k'>Total size</div><div>{esc(human_bytes(stats.get('total_bytes', 0)))}</div>"
    )
    if stats.get("time_start") or stats.get("time_end"):
        parts.append(
            f"<div class='k'>Time range</div><div>{esc(stats.get('time_start',''))} → {esc(stats.get('time_end',''))}</div>"
        )
    parts.append("</div>")

    if entries:
        parts.append("<div class='card'>")
        parts.append(
            "<table><thead><tr>"
            "<th>Datetime (UTC)</th><th>Size</th><th>Path</th><th>Comment</th><th>Preview</th>"
            "</tr></thead><tbody>"
        )
        for e in entries:
            dt = e.get("datetime", "")
            size = human_bytes(e.get("size", 0))
            path_disp = e.get("path", "")
            comment = e.get("comment", "")
            ext = (e.get("ext") or "").lower()

            href, exists, abs_target = make_href_and_exists(
                path_disp, html_out, species
            )

            preview = ""
            if exists and ext in IMG_EXTS:
                preview = (
                    f"<details><summary>show</summary>"
                    f"<img src='{esc(href)}' alt='{esc(os.path.basename(href))}' "
                    f"style='max-width:100%;height:auto;border:1px solid #eee;border-radius:6px;margin-top:6px;'>"
                    f"</details>"
                )
            elif exists and ext in PDF_EXTS:
                preview = (
                    f"<details><summary>show</summary>"
                    f"<object data='{esc(href)}' type='application/pdf' "
                    f"style='width:100%;height:520px;border:1px solid #eee;border-radius:6px;margin-top:6px;'>"
                    f"<a href='{esc(href)}'>{esc(href)}</a></object>"
                    f"</details>"
                )

            parts.append("<tr>")
            parts.append(f"<td>{esc(dt)}</td>")
            parts.append(f"<td>{esc(size)}</td>")
            parts.append(f"<td><a href='{esc(href)}'>{esc(path_disp)}</a></td>")
            parts.append(f"<td>{esc(comment)}</td>")
            parts.append(f"<td>{preview}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table></div>")

    # subsections
    for ss in subs:
        parts.append(section_html(ss, html_out, species))

    return "\n".join(parts)


def build_html(doc: Dict[str, Any], html_out: str) -> str:
    species = doc.get("species")
    title = f"polap-assemble report ({species})" if species else "polap-assemble report"

    css = """
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,Arial,sans-serif;margin:20px;line-height:1.45}
h1{font-size:1.6rem;margin-bottom:0.6rem}
h2{font-size:1.25rem;margin-top:1.6rem}
.card{border:1px solid #ddd;border-radius:8px;padding:12px;margin:12px 0;background:#fafafa}
.kv{display:grid;grid-template-columns:200px 1fr;gap:6px 12px;max-width:980px}
.kv .k{font-weight:600}
table{border-collapse:collapse;width:100%;max-width:1200px}
th,td{border:1px solid #e5e5e5;padding:6px 8px;text-align:left;vertical-align:top}
th{background:#f7f7f7}
nav{background:#f3f6fb;border:1px solid #dce6f3;padding:10px;border-radius:8px;margin:10px 0}
nav a{margin-right:10px}
small{color:#666}
details summary{cursor:pointer}
"""

    # TOC
    toc: List[str] = []
    toc.append("<nav><b>Sections:</b> ")
    first = True
    for s in doc.get("sections", []):
        sid = s.get("id", "")
        if not sid:
            continue
        if not first:
            toc.append(" | ")
        anchor = f"#sec-{sid.replace('.', '-')}"
        toc.append(f"<a href='{esc(anchor)}'>{esc(sid)}. {esc(s.get('title',''))}</a>")
        subs = s.get("subsections") or []
        if subs:
            toc.append("<small> (")
            toc.append(
                ", ".join(
                    f"<a href='#sec-{esc(ss.get('id','').replace('.','-'))}'>{esc(ss.get('id',''))}</a>"
                    for ss in subs
                )
            )
            toc.append(")</small>")
        first = False
    toc.append("</nav>")

    header_card = (
        "<div class='card'><div class='kv'>"
        f"<div class='k'>Base dir</div><div>{esc(doc.get('base_dir',''))}</div>"
        f"<div class='k'>Species</div><div>{esc(doc.get('species',''))}</div>"
        f"<div class='k'>Generated at (UTC)</div><div>{esc(doc.get('generated_at',''))}</div>"
        "</div></div>"
    )

    body_parts: List[str] = []
    for s in doc.get("sections", []):
        body_parts.append(section_html(s, html_out, species))

    html_doc = [
        "<!doctype html>",
        "<html lang='en'><head><meta charset='utf-8'>",
        f"<title>{esc(title)}</title>",
        f"<style>{css}</style>",
        "</head><body>",
        f"<h1>{esc(title)}</h1>",
        header_card,
        "".join(toc),
        "\n".join(body_parts),
        "</body></html>",
    ]
    return "\n".join(html_doc)


def main():
    ap = argparse.ArgumentParser(description="Render polap-assemble JSON to HTML.")
    ap.add_argument(
        "--json", required=True, help="JSON from report_assemble_from_base.py"
    )
    ap.add_argument("--html", required=True, help="Output HTML file")
    args = ap.parse_args()

    with open(args.json, "r", encoding="utf-8") as f:
        doc = json.load(f)

    out_html = build_html(doc, os.path.abspath(args.html))
    os.makedirs(os.path.dirname(os.path.abspath(args.html)), exist_ok=True)
    with open(args.html, "w", encoding="utf-8") as fo:
        fo.write(out_html)

    print(f"[OK] Wrote HTML: {os.path.abspath(args.html)}")


if __name__ == "__main__":
    raise SystemExit(main())
