#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_assembly.py
Version: v0.2.1
SPDX-License-Identifier: GPL-3.0-or-later

Render a block-based assemble report JSON (from report_assemble.py) to HTML using Jinja2.

Usage:
  python3 render_assembly.py \
    --base-dir <Species_dir> \
    --json     <Species_dir>/v6/0/polap-assemble/report/assemble-report.json \
    --template <Species_dir>/_templates \
    --html-out <Species_dir>/report-assemble.html
"""

from __future__ import annotations
import argparse
import json
import os
import re
from typing import Any, Dict, List, Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

VERSION = "render_assembly.py v0.2.0"


def _esc(x: Any) -> str:
    return "" if x is None else str(x)


def _species_from_page_or_base(
    page: Dict[str, Any], base_dir: str, sections: List[Dict[str, Any]]
) -> str:
    sp = _esc((page or {}).get("species"))
    if sp:
        return sp
    sp = os.path.basename(os.path.abspath(base_dir))
    if sp:
        return sp
    # sniff from first href
    for s in sections or []:
        for b in s.get("blocks", []):
            href = _esc(b.get("href"))
            if href:
                m = re.search(r"/([^/]+)/v\d+/", href)
                if m:
                    return m.group(1)
        for ss in s.get("subsections", []):
            for b in ss.get("blocks", []):
                href = _esc(b.get("href"))
                if href:
                    m = re.search(r"/([^/]+)/v\d+/", href)
                    if m:
                        return m.group(1)
    return "Species"


def _find_species_root_up_from(html_out: str, species: str) -> Optional[str]:
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


def _to_rel_href(
    species_root: Optional[str], html_out: str, species: str, href: str
) -> str:
    if not href:
        return href
    h = href.replace("\\", "/")
    out_dir = os.path.dirname(os.path.abspath(html_out))
    if species_root and (h.startswith(species + "/") or ("/" + species + "/") in h):
        if not h.startswith(species + "/"):
            i = h.find("/" + species + "/")
            if i != -1:
                h = h[i + 1 :]
        tail = h[len(species) + 1 :] if h.startswith(species + "/") else h
        abs_target = os.path.join(species_root, tail)
        rel = os.path.relpath(abs_target, start=out_dir).replace("\\", "/")
        return rel
    abs_candidate = os.path.abspath(os.path.join(out_dir, h))
    if os.path.exists(abs_candidate):
        return os.path.relpath(abs_candidate, start=out_dir).replace("\\", "/")
    return h


def _normalize_figure_hrefs(
    sections: List[Dict[str, Any]],
    species_root: Optional[str],
    html_out: str,
    species: str,
) -> List[Dict[str, Any]]:
    norm: List[Dict[str, Any]] = []
    for s in sections or []:
        node = {"id": s.get("id"), "title": s.get("title"), "blocks": []}
        for b in s.get("blocks", []) or []:
            t = _esc(b.get("type")).lower()
            if t == "figure":
                nb = dict(b)
                nb["href"] = _to_rel_href(
                    species_root, html_out, species, _esc(b.get("href"))
                )
                node["blocks"].append(nb)
            else:
                node["blocks"].append(b)
        subs_out: List[Dict[str, Any]] = []
        for ss in s.get("subsections", []) or []:
            snode = {"id": ss.get("id"), "title": ss.get("title"), "blocks": []}
            for b in ss.get("blocks", []) or []:
                t = _esc(b.get("type")).lower()
                if t == "figure":
                    nb = dict(b)
                    nb["href"] = _to_rel_href(
                        species_root, html_out, species, _esc(b.get("href"))
                    )
                    snode["blocks"].append(nb)
                else:
                    snode["blocks"].append(b)
            subs_out.append(snode)
        if subs_out:
            node["subsections"] = subs_out
        norm.append(node)
    return norm


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Render assemble HTML from block JSON with Jinja2."
    )
    ap.add_argument(
        "--base-dir",
        required=True,
        help="Species base directory (e.g., Brassica_rapra)",
    )
    ap.add_argument(
        "--json", required=True, help="assemble-report.json from report_assemble.py"
    )
    ap.add_argument(
        "--template",
        required=True,
        help="Templates dir (must contain page.html or assemble.html)",
    )
    ap.add_argument("--html-out", required=True, help="Output HTML in species folder")
    args = ap.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    json_path = os.path.abspath(args.json)
    tmpl_dir = os.path.abspath(args.template)
    html_out = os.path.abspath(args.html_out)

    if not os.path.isdir(base_dir):
        sys.exit(f"[ERR] --base-dir not found: {base_dir}")
    if not os.path.isfile(json_path):
        sys.exit(f"[ERR] --json not found: {json_path}")
    if not os.path.isdir(tmpl_dir):
        sys.exit(f"[ERR] --template dir not found: {tmpl_dir}")

    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    page_in = data.get("page", {}) or data.get("meta", {})
    sections = data.get("sections", [])

    species = _species_from_page_or_base(page_in, base_dir, sections)
    species_root = _find_species_root_up_from(html_out, species)

    # prefer local Pico
    local_pico = os.path.join(base_dir, "assets", "pico.min.css")
    pico_href = (
        os.path.relpath(local_pico, start=os.path.dirname(html_out)).replace("\\", "/")
        if os.path.isfile(local_pico)
        else "https://unpkg.com/@picocss/pico@2.0.6/css/pico.min.css"
    )

    sections_norm = _normalize_figure_hrefs(sections, species_root, html_out, species)

    page = {
        "title": page_in.get("title") or f"polap-assemble report ({species})",
        "species": species,
        "base_dir": page_in.get("base_dir", base_dir),
        "generated_at": page_in.get("generated_at", ""),
        "pico_href": pico_href,
        "version": VERSION,
    }

    env = Environment(
        loader=FileSystemLoader(tmpl_dir),
        autoescape=select_autoescape(enabled_extensions=("html",), default=True),
    )
    template_name = (
        "page.html"
        if os.path.isfile(os.path.join(tmpl_dir, "page.html"))
        else "assemble.html"
    )
    tmpl = env.get_template(template_name)

    html = tmpl.render(page=page, sections=sections_norm)

    os.makedirs(os.path.dirname(html_out), exist_ok=True)
    with open(html_out, "w", encoding="utf-8") as out:
        out.write(html)
    print(f"[OK] wrote HTML: {html_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
