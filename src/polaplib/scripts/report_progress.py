#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_progress.py
Version: v1.3.0
SPDX-License-Identifier: GPL-3.0-or-later

Create an aggregate JSON of MT/PT assembly progress across species.

It reads:
  - species-codes.txt   (mapping code → Species folder name; code_full = code + 2-digit index)
  - For each <Species> under --root, scans:
      tmp/l.sra.txt  tmp/s.sra.txt
      v6/0/summary-data/l.fq.seqkit.stats.ta.txt
      v6/0/summary-data/s_1.fq.seqkit.stats.ta.txt
      v6/0/summary-data/s_2.fq.seqkit.stats.ta.txt
      v6/0/ncbi-ptdna/ptdna-reference.fa
      v6/0/ncbi-ptdna/00-bioproject/2-mtdna.accession
      v6/0/ncbi-mtdna/mtdna-reference.fa
      v6/0/ncbi-mtdna/00-bioproject/2-mtdna.accession
      v6/0/polap-assemble/pt.1.fasta
      v6/0/polap-assemble/mt.1.fasta
      v6/0/summary-polap-assemble.txt
      v6/0/timing-polap-assemble.txt
      v6/0/data-downsample-long/<long_sra>.fastq.seqkit.stats.ta.tsv
      v6/0/data-downsample-short/<short_sra>_1.fastq.seqkit.stats.ta.tsv
      v6/0/data-downsample-short/<short_sra>_2.fastq.seqkit.stats.ta.tsv

Conventions:
  - Missing files are tolerated; fields become "NA".
  - Converted fields:
      * long-read total bases (sum_len)           → sum_len_gb (decimal Gb, 2 dp; divisor = 1e9)
      * actually-used long-read bases (sum_len)   → sum_len_gb (decimal Gb, 2 dp; divisor = 1e9)
      * net memory increase (KB)                  → net_increase_gb (2 dp; divisor = 1024*1024)

Usage:
  python3 report_progress.py \
    --species-codes species-codes.txt \
    --root /PATH/TO/ROOT \
    --out progress-report.json
"""

from __future__ import annotations
import argparse
import csv
import io
import json
import os
import re
from typing import Dict, List, Tuple, Optional, Any


VERSION = "report_progress.py v1.3.0"


# ----------------------------- utilities -----------------------------


def iso_now() -> str:
    try:
        from datetime import datetime, timezone

        return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except Exception:
        return ""


def read_species_codes(path: str) -> List[Tuple[str, str]]:
    """
    Read 'species-codes.txt' which contains lines like:
      Aa Anthoceros_agrestis
      Ag Anthoceros_angustus
    Returns list of (code, species_dir_name) preserving order.
    """
    out: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 2:
                continue
            out.append((parts[0], parts[1]))
    return out


def slurp_first_line(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if ln:
                    return ln
        return "NA"
    except FileNotFoundError:
        return "NA"
    except Exception:
        return "NA"


def parse_seqkit_tsv_first_row(path: str) -> Dict[str, str]:
    """
    Parse a seqkit TSV (header + data rows) and return dict for the FIRST data row.
    If missing/invalid returns {}.
    """
    try:
        with open(path, "r", encoding="utf-8") as fh:
            content = fh.read()
        if not content.strip():
            return {}
        fh2 = io.StringIO(content)
        reader = csv.reader(fh2, delimiter="\t")
        rows = [r for r in reader if r]
        if not rows:
            return {}
        header = rows[0]
        data = None
        for r in rows[1:]:
            if len(r) >= len(header):
                data = r
                break
        if data is None:
            return {}
        return {
            header[i].strip(): data[i].strip() if i < len(data) else ""
            for i in range(len(header))
        }
    except FileNotFoundError:
        return {}
    except Exception:
        return {}


def parse_fasta_len_and_first_header(path: str) -> Tuple[int, int, str]:
    """
    Return (seq_count, total_bases, first_header_without_gt) or (0,0,"") on error.
    """
    try:
        seqs = 0
        total = 0
        first_header = ""
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith(">"):
                    seqs += 1
                    if not first_header:
                        first_header = ln[1:].strip()
                else:
                    total += len(ln.strip())
        return seqs, total, first_header
    except FileNotFoundError:
        return 0, 0, ""
    except Exception:
        return 0, 0, ""


def header_to_accession(header: str) -> str:
    """
    Extract an accession token from a FASTA header (best effort).
    - If pipe-delimited, use the second token (e.g., ref|NC_012345.1|...)
    - Else the first whitespace-delimited token.
    """
    if not header:
        return ""
    if "|" in header:
        toks = header.split("|")
        if len(toks) >= 2 and toks[1]:
            return toks[1]
    return header.split()[0]


def check_ncbi_reference(accession_file: str, fasta_path: str) -> Dict[str, Any]:
    """
    Compare accession in the file with the FASTA header; compute sequence length and seq_count.
    """
    acc = slurp_first_line(accession_file)
    seq_count, total_len, first_header = parse_fasta_len_and_first_header(fasta_path)
    fasta_acc = header_to_accession(first_header)
    if acc == "NA" and not fasta_acc:
        match = "NA"
    else:
        match = "true" if (acc != "NA" and acc in first_header) else "false"
    return {
        "accession": acc,
        "fasta_accession": fasta_acc if fasta_acc else "NA",
        "acc_match": match,
        "length": total_len if total_len > 0 else "NA",
        "seq_count": seq_count if seq_count > 0 else "NA",
    }


def parse_summary_polap(path: str) -> Dict[str, str]:
    """
    Parse summary-polap-assemble.txt for:
      - Elapsed time (HH:MM:SS)
      - Net increase KB (digits) -> plus computed net_increase_gb
      - Disk used GB (digits)
    Returns string dict with "NA" if not found.
    """
    out = {
        "elapsed_hms": "NA",
        "net_increase_kb": "NA",
        "net_increase_gb": "NA",
        "disk_used_gb": "NA",
    }
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                s = ln.strip()
                m = re.match(r"^Elapsed time:\s*([0-9]{2}:[0-9]{2}:[0-9]{2})", s)
                if m:
                    out["elapsed_hms"] = m.group(1)
                    continue
                m = re.match(r"^Net increase:\s*([0-9]+)\s*KB", s)
                if m:
                    out["net_increase_kb"] = m.group(1)
                    # compute GB here
                    try:
                        kb = float(m.group(1))
                        gb = kb / (1024.0 * 1024.0)  # 1 GB = 1024*1024 KB
                        out["net_increase_gb"] = f"{gb:.2f}"
                    except Exception:
                        out["net_increase_gb"] = "NA"
                    continue
                m = re.match(r"^Disk used:\s*([0-9]+)\s*GB", s)
                if m:
                    out["disk_used_gb"] = m.group(1)
                    continue
        return out
    except FileNotFoundError:
        return out
    except Exception:
        return out


def parse_timing_polap(path: str) -> Dict[str, str]:
    """
    Parse timing-polap-assemble.txt for Hostname, CPU Model, CPU Cores.
    Returns string dict with "NA" if not found.
    """
    out = {
        "hostname": "NA",
        "cpu_model": "NA",
        "cpu_cores": "NA",
    }
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                s = ln.strip()
                m = re.match(r"^Hostname:\s*(.+)$", s)
                if m:
                    out["hostname"] = m.group(1).strip()
                    continue
                m = re.match(r"^CPU Model:\s*(.+)$", s)
                if m:
                    out["cpu_model"] = m.group(1).strip()
                    continue
                m = re.match(r"^CPU Cores:\s*([0-9]+)", s)
                if m:
                    out["cpu_cores"] = m.group(1)
                    continue
        return out
    except FileNotFoundError:
        return out
    except Exception:
        return out


def select_fields(row: Dict[str, str], fields: List[str]) -> Dict[str, str]:
    """Pick fields (or NA) from a row dict."""
    if not row:
        return {f: "NA" for f in fields}
    out: Dict[str, str] = {}
    for f in fields:
        v = row.get(f, "")
        out[f] = v if v != "" else "NA"
    return out


def add_sum_len_gb(d: Dict[str, str]) -> Dict[str, str]:
    """
    From a dict containing 'sum_len' (bases), add 'sum_len_gb' (decimal GB, 2 dp; divisor = 1e9).
    If missing/invalid -> 'NA'.
    """
    val = d.get("sum_len", "NA")
    if val in ("NA", "", None):
        d["sum_len_gb"] = "NA"
        return d
    try:
        bases = float(val)
        gb = bases / 1e9  # decimal gigabases
        d["sum_len_gb"] = f"{gb:.2f}"
    except Exception:
        d["sum_len_gb"] = "NA"
    return d


def rel_or_na(path: str, root: str) -> str:
    if path == "NA":
        return "NA"
    try:
        if not os.path.exists(path):
            return "NA"
        return os.path.relpath(path, root)
    except Exception:
        return "NA"


# ----------------------------- per-species scan -----------------------------


def gather_species(root: str, code: str, idx1: int, species: str) -> Dict[str, Any]:
    """
    Build the record for one species. idx1 is 1-based index used to form code_full = code + 2-digit index.
    """
    species_dir = os.path.join(root, species)
    run_root = os.path.join(species_dir, "v6", "0")

    # SRA files
    p_l_sra = os.path.join(species_dir, "tmp", "l.sra.txt")
    p_s_sra = os.path.join(species_dir, "tmp", "s.sra.txt")
    long_sra = slurp_first_line(p_l_sra)
    short_sra = slurp_first_line(p_s_sra)

    # seqkit summaries
    p_l_sum = os.path.join(run_root, "summary-data", "l.fq.seqkit.stats.ta.txt")
    p_s1_sum = os.path.join(run_root, "summary-data", "s_1.fq.seqkit.stats.ta.txt")
    p_s2_sum = os.path.join(run_root, "summary-data", "s_2.fq.seqkit.stats.ta.txt")

    wanted = ["num_seqs", "sum_len", "AvgQual"]
    long_sum = add_sum_len_gb(
        select_fields(parse_seqkit_tsv_first_row(p_l_sum), wanted)
    )
    s1_sum = select_fields(parse_seqkit_tsv_first_row(p_s1_sum), wanted)
    s2_sum = select_fields(parse_seqkit_tsv_first_row(p_s2_sum), wanted)

    # used/actual data, derived from SRA IDs
    p_used_long = (
        os.path.join(
            run_root, "data-downsample-long", f"{long_sra}.fastq.seqkit.stats.ta.tsv"
        )
        if long_sra != "NA"
        else "NA"
    )
    p_used_s1 = (
        os.path.join(
            run_root,
            "data-downsample-short",
            f"{short_sra}_1.fastq.seqkit.stats.ta.tsv",
        )
        if short_sra != "NA"
        else "NA"
    )
    p_used_s2 = (
        os.path.join(
            run_root,
            "data-downsample-short",
            f"{short_sra}_2.fastq.seqkit.stats.ta.tsv",
        )
        if short_sra != "NA"
        else "NA"
    )

    used_long = add_sum_len_gb(
        select_fields(
            parse_seqkit_tsv_first_row(p_used_long) if p_used_long != "NA" else {},
            wanted,
        )
    )
    used_s1 = select_fields(
        parse_seqkit_tsv_first_row(p_used_s1) if p_used_s1 != "NA" else {}, wanted
    )
    used_s2 = select_fields(
        parse_seqkit_tsv_first_row(p_used_s2) if p_used_s2 != "NA" else {}, wanted
    )

    # NCBI FASTA + accession checks
    p_pt_fa = os.path.join(run_root, "ncbi-ptdna", "ptdna-reference.fa")
    p_pt_acc = os.path.join(
        run_root, "ncbi-ptdna", "00-bioproject", "2-mtdna.accession"
    )
    p_mt_fa = os.path.join(run_root, "ncbi-mtdna", "mtdna-reference.fa")
    p_mt_acc = os.path.join(
        run_root, "ncbi-mtdna", "00-bioproject", "2-mtdna.accession"
    )
    ncbi_pt = check_ncbi_reference(p_pt_acc, p_pt_fa)
    ncbi_mt = check_ncbi_reference(p_mt_acc, p_mt_fa)

    # polap fasta stats
    p_pt_polap = os.path.join(run_root, "polap-assemble", "pt.1.fasta")
    p_mt_polap = os.path.join(run_root, "polap-assemble", "mt.1.fasta")
    pt_seq_count, pt_total, _ = parse_fasta_len_and_first_header(p_pt_polap)
    mt_seq_count, mt_total, _ = parse_fasta_len_and_first_header(p_mt_polap)
    polap_pt = {
        "fasta": rel_or_na(p_pt_polap, root),
        "seq_count": pt_seq_count if pt_seq_count > 0 else "NA",
        "total_bases": pt_total if pt_total > 0 else "NA",
    }
    polap_mt = {
        "fasta": rel_or_na(p_mt_polap, root),
        "seq_count": mt_seq_count if mt_seq_count > 0 else "NA",
        "total_bases": mt_total if mt_total > 0 else "NA",
    }

    # perf summary + timing
    p_summary = os.path.join(run_root, "summary-polap-assemble.txt")
    p_timing = os.path.join(run_root, "timing-polap-assemble.txt")
    perf_summary = parse_summary_polap(p_summary)
    timing_info = parse_timing_polap(p_timing)

    record: Dict[str, Any] = {
        "code": code,
        "code_full": f"{code}{idx1:02d}",
        "order": idx1,
        "species": species,
        "paths": {
            "species_dir": rel_or_na(species_dir, root),
            "run_root": rel_or_na(run_root, root),
        },
        "sra": {
            "long": long_sra,
            "short": short_sra,
        },
        "reads_summary": {
            "long": {**long_sum, "source": rel_or_na(p_l_sum, root)},
            "short1": {**s1_sum, "source": rel_or_na(p_s1_sum, root)},
            "short2": {**s2_sum, "source": rel_or_na(p_s2_sum, root)},
        },
        "reads_used": {
            "long": {**used_long, "source": rel_or_na(p_used_long, root)},
            "short1": {**used_s1, "source": rel_or_na(p_used_s1, root)},
            "short2": {**used_s2, "source": rel_or_na(p_used_s2, root)},
        },
        "ncbi": {
            "pt": {
                **ncbi_pt,
                "accession_file": rel_or_na(p_pt_acc, root),
                "fasta_path": rel_or_na(p_pt_fa, root),
            },
            "mt": {
                **ncbi_mt,
                "accession_file": rel_or_na(p_mt_acc, root),
                "fasta_path": rel_or_na(p_mt_fa, root),
            },
        },
        "polap": {
            "pt": polap_pt,
            "mt": polap_mt,
        },
        "summaries": {"polap_assemble": perf_summary},  # polap-assemble summary
        "system": timing_info,  # timing-polap-assemble info
    }
    return record


# ----------------------------- main -----------------------------


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Create an aggregate JSON of MT/PT assembly progress across species."
    )
    ap.add_argument("--species-codes", required=True, help="Path to species-codes.txt")
    ap.add_argument(
        "--root",
        default=".",
        help="Root folder containing species directories (default: .)",
    )
    ap.add_argument("--out", required=True, help="Output JSON file path")
    args = ap.parse_args()

    pairs = read_species_codes(args.species_codes)
    if not pairs:
        print("[ERR] No species read from --species-codes", file=sys.stderr)
        return 2

    records: List[Dict[str, Any]] = []
    for i, (code, species) in enumerate(pairs, start=1):
        try:
            rec = gather_species(args.root, code, i, species)
        except Exception as e:
            print(f"[WARN] Failed on species {species}: {e}", file=sys.stderr)
            continue
        records.append(rec)

    model = {
        "meta": {
            "version": VERSION,
            "generated_at": iso_now(),
            "root": os.path.abspath(args.root),
            "species_count": len(records),
            "source_species_codes": os.path.abspath(args.species_codes),
        },
        "species": records,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out:
        json.dump(model, out, indent=2)
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
