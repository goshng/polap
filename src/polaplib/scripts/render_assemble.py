#!/usr/bin/env python3
# render_assemble.py
# Version: v0.4.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Render a block-based polap-assemble report to HTML using Jinja2.
#
# Usage:
#   python3 render_assemble.py \
#       --base-dir Brassica_rapa \
#       --json Brassica_rapa/v6/0/polap-assemble/report/assemble-report.json \
#       --template Brassica_rapa/_templates \
#       --html-out Brassica_rapa/report-assemble.html
#
# Notes
# - This renderer DOES NOT parse TSV/PNG/PDF. It only renders what the JSON provides.
# - JSON should contain "page" and "sections" with block items (figure/table/text/kpi).
# - Hrefs in JSON should be species-rooted (e.g., 'Brassica_rapa/v6/0/...').
# - We convert species-rooted hrefs to paths relative to --html-out so local previews work.
#
import os
import re
import sys
import json
import argparse
from typing import Optional, Dict, Any, List, Tuple
from jinja2 import Environment, FileSystemLoader, select_autoescape

IMG_EXTS = {".png", ".jpg", ".jpeg", ".svg"}
PDF_EXTS = {".pdf"}

# ----------------------- helpers -----------------------


def esc(x) -> str:
    return "" if x is None else str(x)


def detect_species(
    base_dir: str, json_page: Dict[str, Any], sections: List[Dict[str, Any]]
) -> str:
    """
    Species priority:
      1) JSON page.species
      2) basename of --base-dir
      3) sniff from first figure/table href: /<Species>/vN/...
    """
    # 1) from JSON
    sp = esc((json_page or {}).get("species"))
    if sp:
        return sp
    # 2) from base-dir
    sp = os.path.basename(os.path.abspath(base_dir))
    if sp:
        return sp

    # 3) sniff from first href
    def sniff():
        for sec in sections or []:
            for b in sec.get("blocks", []):
                href = esc(b.get("href"))
                if href:
                    m = re.search(r"/([^/]+)/v\d+/", href)
                    if m:
                        return m.group(1)
            for s in sec.get("subsections", []):
                for b in s.get("blocks", []):
                    href = esc(b.get("href"))
                    if href:
                        m = re.search(r"/([^/]+)/v\d+/", href)
                        if m:
                            return m.group(1)
        return ""

    sp = sniff()
    return sp or "Species"


def find_species_root(html_out: str, species: str) -> Optional[str]:
    """
    Walk upward from the HTML output dir to locate a folder literally named <species>.
    This lets us resolve 'Species/...' to a real absolute for anchoring hrefs.
    """
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


def to_rel_href(
    species_root: Optional[str], html_out: str, species: str, href: str
) -> str:
    """
    Convert a species-rooted href (e.g., 'Species/v6/0/...') into a link relative to html_out.
    If species_root is found, we rebuild the absolute path then make it relative.
    Otherwise, we leave the href as-is (best effort).
    """
    if not href:
        return href
    h = href.replace("\\", "/")
    out_dir = os.path.dirname(os.path.abspath(html_out))

    # If starts with <Species>/... we can resolve via species_root
    if species_root and (h.startswith(species + "/") or ("/" + species + "/") in h):
        # ensure it starts with '<Species>/'
        if not h.startswith(species + "/"):
            i = h.find("/" + species + "/")
            if i != -1:
                h = h[i + 1 :]  # drop leading '/'
        tail = h[len(species) + 1 :] if h.startswith(species + "/") else h
        abs_target = os.path.join(species_root, tail)
        if os.path.exists(abs_target):
            rel = os.path.relpath(abs_target, start=out_dir).replace("\\", "/")
            return rel
        # If file doesn't exist, still provide a relative pointing under species_root tail
        rel = os.path.relpath(os.path.join(species_root, tail), start=out_dir).replace(
            "\\", "/"
        )
        return rel

    # If not species-rooted, try relative from html_out directory as-is
    # (Useful if JSON already stored a path relative to html_out)
    tentative = os.path.abspath(os.path.join(out_dir, h))
    if os.path.exists(tentative):
        return os.path.relpath(tentative, start=out_dir).replace("\\", "/")

    return h  # fallback: leave it


def normalize_blocks(
    blocks: List[Dict[str, Any]],
    species_root: Optional[str],
    html_out: str,
    species: str,
) -> List[Dict[str, Any]]:
    """
    Normalize block hrefs to be relative to html_out. Do NOT change table/text data.
    """
    out = []
    for b in blocks or []:
        btype = esc(b.get("type")).lower()
        if btype == "figure":
            href = esc(b.get("href"))
            new_href = to_rel_href(species_root, html_out, species, href)
            nb = dict(b)
            nb["href"] = new_href
            out.append(nb)
        elif btype in ("table", "text", "kpi"):
            out.append(b)
        else:
            # unknown type → pass through
            out.append(b)
    return out


def normalize_sections(
    sections: List[Dict[str, Any]],
    species_root: Optional[str],
    html_out: str,
    species: str,
) -> List[Dict[str, Any]]:
    norm = []
    for s in sections or []:
        node = {
            "id": esc(s.get("id")),
            "title": esc(s.get("title")),
            "blocks": normalize_blocks(
                s.get("blocks", []), species_root, html_out, species
            ),
        }
        # subsections?
        subs_norm = []
        for ss in s.get("subsections", []) or []:
            subs_norm.append(
                {
                    "id": esc(ss.get("id")),
                    "title": esc(ss.get("title")),
                    "blocks": normalize_blocks(
                        ss.get("blocks", []), species_root, html_out, species
                    ),
                }
            )
        if subs_norm:
            node["subsections"] = subs_norm
        norm.append(node)
    return norm


# ----------------------- main -----------------------


def main():
    ap = argparse.ArgumentParser(
        description="Render polap-assemble block JSON to HTML (Jinja2)."
    )
    ap.add_argument(
        "--base-dir", required=True, help="Species base directory (e.g., Brassica_rapa)"
    )
    ap.add_argument(
        "--json",
        required=True,
        help="Path to prebuilt assemble-report.json (block schema)",
    )
    ap.add_argument(
        "--template",
        required=True,
        help="Path to templates directory (must contain 'page.html' or 'assemble-blocks.html')",
    )
    ap.add_argument(
        "--html-out",
        required=True,
        help="Output HTML file (e.g., Brassica_rapa/report-assemble.html)",
    )
    args = ap.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    json_path = os.path.abspath(args.json)
    tdir = os.path.abspath(args.template)
    html_out = os.path.abspath(args.html_out)

    if not os.path.isdir(base_dir):
        sys.exit(f"[ERR] --base-dir not found: {base_dir}")
    if not os.path.isfile(json_path):
        sys.exit(f"[ERR] --json not found: {json_path}")
    if not os.path.isdir(tdir):
        sys.exit(f"[ERR] --template dir not found: {tdir}")

    # Load JSON (block schema expected)
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # Page info
    json_page = data.get("page", {})
    sections = data.get("sections", [])

    species = detect_species(base_dir, json_page, sections)
    species_root = find_species_root(html_out, species)

    # Prefer local Pico (assets/pico.min.css) under base-dir; else fallback to CDN
    local_pico = os.path.join(base_dir, "assets", "pico.min.css")
    pico_href = (
        os.path.relpath(local_pico, start=os.path.dirname(html_out)).replace("\\", "/")
        if os.path.isfile(local_pico)
        else "https://unpkg.com/@picocss/pico@2.0.6/css/pico.min.css"
    )

    # Normalize blocks' hrefs relative to html_out
    sections_vm = normalize_sections(sections, species_root, html_out, species)

    # Page VM
    page_vm = {
        "title": esc(json_page.get("title")) or f"polap-assemble report ({species})",
        "species": esc(json_page.get("species")) or species,
        "base_dir": esc(json_page.get("base_dir")) or base_dir,
        "generated_at": esc(json_page.get("generated_at")) or "",
        "pico_href": pico_href,
    }

    # Pick template file: allow either 'page.html' (preferred) or 'assemble-blocks.html'
    template_candidates = ("page.html", "assemble-blocks.html")
    env = Environment(
        loader=FileSystemLoader(tdir),
        autoescape=select_autoescape(enabled_extensions=("html",), default=True),
    )
    template_name = None
    for cand in template_candidates:
        if os.path.isfile(os.path.join(tdir, cand)):
            template_name = cand
            break
    if not template_name:
        sys.exit(
            f"[ERR] No page template found in {tdir}. Expected one of: {', '.join(template_candidates)}"
        )

    tmpl = env.get_template(template_name)
    html = tmpl.render(page=page_vm, sections=sections_vm)

    os.makedirs(os.path.dirname(html_out), exist_ok=True)
    with open(html_out, "w", encoding="utf-8") as fo:
        fo.write(html)

    print(f"[OK] wrote HTML: {html_out}")


if __name__ == "__main__":
    raise SystemExit(main())
