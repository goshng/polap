#!/usr/bin/env python3
# Version: v0.9.0
"""
polap_py_mt_selfrepeats.py
Detect non-redundant intra-genome repeat pairs (≥ min length, ≥ min %identity)
from either:
  • MUMmer show-coords (-T -r -c -l)  [preferred]
  • BLASTN outfmt 6 (self vs self)

Robust to:
  • MUMmer preamble lines (file paths, 'NUCMER', blank lines)
  • Header present (tokens like '[S1] ... [% IDY] ...') or headerless (-H)

Outputs
-------
1) TSV (required):
   pair_id  r_name  r_s  r_e  q_name  q_s  q_e  orient  len  pid
   where orient ∈ {DIR, INV}, positions are 1-based inclusive

2) BED (required):
   two lines per pair (copy A and copy B), 0-based half-open

Examples
--------
# MUMmer path:
nucmer --maxmatch -p self asm.fa asm.fa
delta-filter -i 95 -l 1000 self.delta > self.filt.delta
show-coords -T -r -c -l self.filt.delta > self.coords.tsv

# BLAST path:
makeblastdb -in asm.fa -dbtype nucl
blastn -task megablast -db asm.fa -query asm.fa -perc_identity 95 -word_size 28 \
       -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen" \
       > self.blast6.tsv

Usage
-----
python3 polap_py_mt_selfrepeats.py \
  --mode mummer --coords self.coords.tsv \
  --assembly asm.fa --min-len 1000 --min-pid 95 \
  --out-tsv repeats.tsv --out-bed repeats.bed
"""

from __future__ import annotations
import argparse
import gzip
import itertools
import logging
from dataclasses import dataclass
from typing import Iterable, Dict, Tuple, List
import re as rx


@dataclass(frozen=True)
class Hit:
    rname: str
    rs: int
    re: int
    qname: str
    qs: int
    qe: int
    length: int
    pid: float
    orient: str  # 'DIR' or 'INV'


def open_auto(path: str):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path, "r")


def _tokenize_spaces(s: str) -> List[str]:
    return rx.split(r"\s+", s.strip())


def parse_coords(path: str) -> Iterable[Hit]:
    """
    Parse MUMmer show-coords (-T -r -c -l) with/without header and preamble.
    """
    with open_auto(path) as fh:
        lines = fh.readlines()

    # Find header row containing 'IDY' and starting with bracket tokens like "[S1]"
    header_i = None
    header_cols = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if not s:
            continue
        if "IDY" in s.upper() and s.startswith("["):
            header_i = i
            header_cols = _tokenize_spaces(s)
            break

    # Determine where data rows start and the %IDY column index
    if header_i is not None:
        data_start = header_i + 1
        try:
            i_pid = next(j for j, c in enumerate(header_cols) if "IDY" in c.upper())
        except StopIteration:
            raise SystemExit("ERROR: show-coords header found but lacks % IDY")
    else:
        # No header present -> find first non-preamble non-empty line
        data_start = None
        for i, ln in enumerate(lines):
            s = ln.strip()
            if not s:
                continue
            # skip obvious preamble lines (paths, 'NUCMER', bracket row if any)
            if s.startswith("/") or s.upper().startswith("NUCMER") or s.startswith("["):
                continue
            data_start = i
            break
        if data_start is None:
            return
        # Standard MUMmer4 (-T -r -c -l): %IDY is 7th column (1-based) -> index 6
        i_pid = 6

    # Fixed numeric fields for -T output:
    # [S1] [E1] [S2] [E2] [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]...
    i_s1, i_e1, i_s2, i_e2 = 0, 1, 2, 3
    # last two columns are reference and query tags (names)
    i_rname, i_qname = -2, -1

    for ln in lines[data_start:]:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("["):
            # safety: skip any repeated header blocks
            continue
        f = _tokenize_spaces(s)
        if len(f) < 12:
            # too short for -T -r -c -l, skip
            continue
        # Parse numeric coordinates; tolerate float pid
        try:
            s1 = int(f[i_s1])
            e1 = int(f[i_e1])
            s2 = int(f[i_s2])
            e2 = int(f[i_e2])
            pid = float(f[i_pid])
            rname = f[i_rname]
            qname = f[i_qname]
        except Exception:
            continue

        l1 = abs(e1 - s1) + 1
        l2 = abs(e2 - s2) + 1
        length = min(l1, l2)

        orient = "INV" if s2 > e2 else "DIR"
        rs, r_end = min(s1, e1), max(s1, e1)
        qs, q_end = min(s2, e2), max(s2, e2)

        yield Hit(rname, rs, r_end, qname, qs, q_end, length, pid, orient)


def parse_blast6(path: str) -> Iterable[Hit]:
    """
    Parse BLAST outfmt 6 rows with fields:
      qseqid sseqid pident length qstart qend sstart send qlen slen
    """
    with open_auto(path) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            parts = _tokenize_spaces(ln)
            if len(parts) < 10:
                # not the expected outfmt; skip
                continue
            qseqid, sseqid = parts[0], parts[1]
            try:
                pid = float(parts[2])
                L = int(parts[3])
                qs = int(parts[4])
                qe = int(parts[5])
                ss = int(parts[6])
                se = int(parts[7])
            except Exception:
                continue
            if L <= 0:
                continue

            orient = "INV" if ss > se else "DIR"
            rs, r_end = (se, ss) if ss > se else (ss, se)
            qs2, q_end = (qe, qs) if qs > qe else (qs, qe)

            yield Hit(sseqid, rs, r_end, qseqid, qs2, q_end, L, pid, orient)


def canonical_key(h: Hit) -> Tuple[Tuple, str]:
    """
    Symmetry-deduplication key (order independent), preserves orientation.
    """
    a = (h.rname, h.rs, h.re)
    b = (h.qname, h.qs, h.qe)
    key = (a + b) if a < b else (b + a)
    return (key, h.orient)


def main():
    ap = argparse.ArgumentParser(
        description="Detect non-redundant self repeats (≥len, ≥pid)."
    )
    ap.add_argument(
        "--mode",
        choices=["mummer", "blast"],
        required=True,
        help="Input type: MUMmer show-coords or BLAST outfmt 6",
    )
    ap.add_argument("--coords", help="show-coords (-T -r -c -l) output (MUMmer).")
    ap.add_argument(
        "--blast6",
        help="BLAST outfmt 6 with fields: qseqid sseqid pident length qstart qend sstart send qlen slen.",
    )
    ap.add_argument(
        "--assembly", required=True, help="Assembly FASTA (for provenance only)."
    )
    ap.add_argument(
        "--min-len", type=int, default=1000, help="Minimum repeat length (bp)."
    )
    ap.add_argument("--min-pid", type=float, default=95.0, help="Minimum %%identity.")
    ap.add_argument(
        "--out-tsv", required=True, help="Output TSV (non-redundant repeats)."
    )
    ap.add_argument("--out-bed", required=True, help="Output BED (two lines per pair).")
    ap.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    ap.add_argument("--version", action="store_true", help="Print version and exit.")
    args = ap.parse_args()

    if args.version:
        print("v0.9.0")
        return

    logging.basicConfig(
        level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s"
    )

    if args.mode == "mummer" and not args.coords:
        ap.error("--coords is required for --mode mummer")
    if args.mode == "blast" and not args.blast6:
        ap.error("--blast6 is required for --mode blast")

    # Parse hits
    if args.mode == "mummer":
        it = parse_coords(args.coords)
    else:
        it = parse_blast6(args.blast6)

    # Deduplicate symmetric pairs, keep the best by (length, pid)
    kept: Dict[Tuple[Tuple, str], Hit] = {}
    raw = 0
    for h in it:
        raw += 1
        # trivial diagonal (same locus) — skip if very close
        if h.rname == h.qname and abs(h.rs - h.qs) < 10 and abs(h.re - h.qe) < 10:
            continue
        if h.length < args.min_len or h.pid < args.min_pid:
            continue
        key = canonical_key(h)
        prev = kept.get(key)
        if prev is None or (h.length, h.pid) > (prev.length, prev.pid):
            kept[key] = h

    logging.info(f"Parsed hits: {raw}; kept non-redundant: {len(kept)}")

    # Write outputs
    with open(args.out_tsv, "w") as oh, open(args.out_bed, "w") as bed:
        oh.write("pair_id\tr_name\tr_s\tr_e\tq_name\tq_s\tq_e\torient\tlen\tpid\n")
        i = 0
        for (_, _), h in kept.items():
            i += 1
            oh.write(
                f"R{i}\t{h.rname}\t{h.rs}\t{h.re}\t{h.qname}\t{h.qs}\t{h.qe}\t{h.orient}\t{h.length}\t{h.pid:.3f}\n"
            )
            # BED is 0-based half-open
            bed.write(f"{h.rname}\t{h.rs-1}\t{h.re}\tR{i}|A|{h.orient}\n")
            bed.write(f"{h.qname}\t{h.qs-1}\t{h.qe}\tR{i}|B|{h.orient}\n")


if __name__ == "__main__":
    main()
