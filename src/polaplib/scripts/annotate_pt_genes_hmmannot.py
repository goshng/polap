#!/usr/bin/env python3
# FILE: scripts/annotate_pt_genes_hmmannot.py
# VERSION: 0.2.1
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Parse Oatk hmmannot output for PT assemblies and extract per-gene CDS/AA.

Usage:
  annotate_pt_genes_hmmannot.py --pt-fasta PT.fa --hmmannot pt_vs_hmms.annot.txt \
                                --outdir out/genes/Species --species Species
"""

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args():
    ap = ArgumentParser()
    ap.add_argument("--pt-fasta", required=True)
    ap.add_argument("--hmmannot", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--species", required=True)
    return ap.parse_args()


def default_indices_from_data():
    # Matches your example rows precisely
    return {
        "target": 0,
        "query": 2,
        "hmmfrom": 4,
        "hmmto": 5,
        "alifrom": 6,
        "alito": 7,
        "envfrom": 8,
        "envto": 9,
        "modlen": 10,
        "strand": 11,
        "evalue": 12,
    }


def indices_from_header(line):
    """Support both compact and split headers."""
    hdr = re.sub(r"^#\s*", "", line.strip().lower())
    cols = re.split(r"\s+", hdr)
    idx = {c: i for i, c in enumerate(cols)}

    def f(*names):
        for n in names:
            if n in idx:
                return idx[n]
        return None

    # Try compact first
    mm = {
        "target": f("target", "target_name", "targetname"),
        "query": f("query", "query_name", "queryname"),
        "hmmfrom": f("hmmfrom", "hmm_from"),
        "hmmto": f("hmmto", "hmm_to"),
        "alifrom": f("alifrom", "ali_from"),
        "alito": f("alito", "ali_to"),
        "envfrom": f("envfrom", "env_from"),
        "envto": f("envto", "env_to"),
        "modlen": f("modlen"),
        "strand": f("strand"),
        "evalue": f("e-value", "evalue"),
    }

    # If ali/env compact fields missing, attempt to stitch split headers.
    # Many hmmannot tables have: "ali from", "ali to", "env from", "env to"
    if mm["alifrom"] is None and "ali" in idx and "from" in idx:
        # The split words won’t give a single index; we’ll use data-position fallback later.
        pass
    if mm["envfrom"] is None and "env" in idx and "from" in idx:
        pass

    return mm


def safe_int(x):
    try:
        return int(x)
    except:
        return None


def choose_indices(idxmap):
    """Return an index map to use. If header-based map doesn’t include
    enough numeric fields, fall back to fixed data positions."""
    need_keys = ("target", "query")
    ok_basic = all(idxmap.get(k) is not None for k in need_keys)
    # require at least ali/env and evalue present
    ok_coords = (
        idxmap.get("envfrom") is not None and idxmap.get("envto") is not None
    ) or (idxmap.get("alifrom") is not None and idxmap.get("alito") is not None)
    ok_eval = idxmap.get("evalue") is not None
    if ok_basic and ok_coords and ok_eval:
        return idxmap
    return default_indices_from_data()


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    pt = {r.id: str(r.seq) for r in SeqIO.parse(a.pt_fasta, "fasta")}
    if not pt:
        sys.stderr.write(f"[annotate_pt_genes_hmmannot] No sequences in {a.pt_fasta}\n")
        sys.exit(1)

    header_idx = None
    use_idx = None
    hits = defaultdict(list)

    # First pass: find a header line if present
    with open(a.hmmannot, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if (
                raw.lstrip().startswith("#")
                and "target" in raw.lower()
                and "query" in raw.lower()
            ):
                header_idx = indices_from_header(raw)
                break

    if header_idx:
        use_idx = choose_indices(header_idx)
    else:
        use_idx = default_indices_from_data()

    # Second pass: parse data rows
    with open(a.hmmannot, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            # Guard against short lines
            if len(parts) < 13:  # up to e-value
                continue

            def get(k):
                i = use_idx.get(k)
                return parts[i] if (i is not None and i < len(parts)) else None

            gene = get("target")
            qname = get("query")
            if not gene or not qname or qname not in pt:
                continue
            gene = re.sub(r"\s+", "_", gene)

            # Prefer envfrom/envto else alifrom/alito
            s = safe_int(get("envfrom"))
            e = safe_int(get("envto"))
            if s is None or e is None:
                s = safe_int(get("alifrom"))
                e = safe_int(get("alito"))
            if s is None or e is None:
                continue

            strand = (get("strand") or "+").strip()
            ev_raw = get("evalue")
            try:
                evalue = float(ev_raw) if ev_raw is not None else 1.0
            except:
                evalue = 1.0

            start = min(s, e) - 1
            end = max(s, e)
            if start < 0 or end > len(pt[qname]) or end <= start:
                continue
            span = end - start

            hits[gene].append(
                {
                    "qname": qname,
                    "start": start,
                    "end": end,
                    "span": span,
                    "strand": strand,
                    "evalue": evalue,
                }
            )

    if not hits:
        sys.stderr.write(
            "[annotate_pt_genes_hmmannot] No usable hits parsed; "
            "check that query names match PT FASTA headers.\n"
        )
        sys.exit(2)

    # Best per gene
    for gene, arr in hits.items():
        arr.sort(key=lambda d: (d["evalue"], -d["span"]))
        h = arr[0]
        nt = pt[h["qname"]][h["start"] : h["end"]]
        if h["strand"] == "-":
            nt = str(Seq(nt).reverse_complement())

        # translate: pick longest frame
        best_aa = ""
        for frame in (0, 1, 2):
            aa = str(Seq(nt[frame:]).translate(to_stop=False))
            if len(aa) > len(best_aa):
                best_aa = aa

        with open(os.path.join(a.outdir, f"{gene}.prot.faa"), "a") as o:
            o.write(f">{a.species}\n{best_aa}\n")
        with open(os.path.join(a.outdir, f"{gene}.cds.fna"), "a") as o:
            o.write(f">{a.species}\n{nt}\n")

    print(f"[annotate_pt_genes_hmmannot] wrote per-gene prot/cds to {a.outdir}")


if __name__ == "__main__":
    main()
