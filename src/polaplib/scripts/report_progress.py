#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_progress.py
Create an aggregate JSON of MT/PT assembly progress across species.

Usage:
  python3 report_progress.py \
    --species-codes species-codes.txt \
    --root /path/to/root_with_species_dirs \
    --out progress-report.json

Notes:
- Expects each species directory at: <root>/<Species_Name>/
- Expects run paths inside each species at: <Species_Name>/v6/0/...
- Missing files are tolerated; fields become "NA".
"""

from __future__ import annotations
import argparse
import csv
import io
import json
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

VERSION = "report_progress.py v1.0.0"


def read_species_codes(path: str) -> List[Tuple[str, str]]:
    """
    Read species codes: each non-empty, non-comment line has:
      <code> <Species_Name>
    Returns a list of (code, species_name) in order.
    """
    pairs: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            code, species = parts[0], parts[1]
            pairs.append((code, species))
    return pairs


def slurp_first_line(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    return line
        return "NA"
    except FileNotFoundError:
        return "NA"
    except Exception:
        return "NA"


def parse_seqkit_tsv_first_row(path: str) -> Dict[str, str]:
    """
    Parse seqkit stats TSV with a header row + a single data row we care about.
    Returns a dict of all columns in that first data row (as strings).
    If file missing/empty, returns {}.
    """
    try:
        with open(path, "r", encoding="utf-8") as fh:
            content = fh.read()
        # Some files could have spaces; we expect *TSV* here.
        # Use csv with delimiter '\t'. If there is only whitespace, bail.
        if not content.strip():
            return {}
        fh2 = io.StringIO(content)
        reader = csv.reader(fh2, delimiter="\t")
        rows = [r for r in reader if r]
        if not rows:
            return {}
        header = rows[0]
        # find the first row that has same or more columns than header
        data = None
        for r in rows[1:]:
            if len(r) >= len(header):
                data = r
                break
        if data is None:
            return {}
        # Build dict
        rowdict = {
            h.strip(): data[i].strip() if i < len(data) else ""
            for i, h in enumerate(header)
        }
        return rowdict
    except FileNotFoundError:
        return {}
    except Exception:
        return {}


def parse_space_table_first_rows(path: str) -> Dict[str, List[List[str]]]:
    """
    For space-delimited files (not needed in this progress script, but here for parity).
    Returns {"header": [...], "rows": [[...], ...]} for the first 1-2 lines.
    (Not used in the three progress tables, kept for future extension.)
    """
    out = {"header": [], "rows": []}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        if not lines:
            return out
        # naive split on whitespace
        header = re.split(r"\s+", lines[0])
        out["header"] = header
        for ln in lines[1:3]:
            out["rows"].append(re.split(r"\s+", ln))
        return out
    except FileNotFoundError:
        return out
    except Exception:
        return out


def parse_fasta_stats_and_first_header(path: str) -> Tuple[int, int, str]:
    """
    Count sequences and total bases; also return the FIRST header line (without '>').
    Returns (seq_count, total_bases, first_header) or (0, 0, "") on failure.
    """
    try:
        seq_count = 0
        total = 0
        first_header = ""
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    seq_count += 1
                    if not first_header:
                        first_header = line[1:].strip()
                else:
                    # count sequence characters ignoring whitespace
                    total += len(line.strip())
        return seq_count, total, first_header
    except FileNotFoundError:
        return 0, 0, ""
    except Exception:
        return 0, 0, ""


def header_to_accession(header: str) -> str:
    """
    Try to extract an accession token from a FASTA header string.
    Strategy:
      - If bar-delimited: something like 'ref|NC_012345.1|...' -> return tokens[1]
      - Else: take the first whitespace-delimited token
    """
    if not header:
        return ""
    if "|" in header:
        toks = header.split("|")
        if len(toks) >= 2 and toks[1]:
            return toks[1]
    # fallback: first token
    return header.split()[0]


def check_ncbi_reference(accession_file: str, fasta_path: str) -> Dict[str, str]:
    """
    Read accession from file and compare with FASTA header.
    Return dict with: accession, fasta_accession, acc_match ("true"/"false"), length (int or "NA")
    """
    acc = slurp_first_line(accession_file)
    seq_count, tot_len, first_header = parse_fasta_stats_and_first_header(fasta_path)
    fasta_acc = header_to_accession(first_header)
    match = "NA"
    if acc == "NA" and not fasta_acc:
        match = "NA"
    else:
        # consider a 'match' if accession appears anywhere in header
        match = "true" if (acc != "NA" and acc in first_header) else "false"
    result = {
        "accession": acc,
        "fasta_accession": fasta_acc if fasta_acc else "NA",
        "acc_match": match,
        "length": tot_len if tot_len > 0 else "NA",
        "seq_count": seq_count if seq_count > 0 else "NA",
    }
    return result


def parse_summary_time_mem(path: str) -> Dict[str, str]:
    """
    Parse memory/time lines from a summary text file.
    We look for:
      Start used:   N KB
      Peak used:    N KB
      Net increase: N KB
      Elapsed time: HH:MM:SS (X.XX h)
    Returns dict with string values (or "NA").
    """
    out = {
        "start_used_kb": "NA",
        "peak_used_kb": "NA",
        "net_increase_kb": "NA",
        "elapsed_hms": "NA",
        "elapsed_hours": "NA",
    }
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                ln = line.strip()
                m1 = re.match(r"^Start used:\s+([0-9]+)\s*KB", ln)
                if m1:
                    out["start_used_kb"] = m1.group(1)
                    continue
                m2 = re.match(r"^Peak used:\s+([0-9]+)\s*KB", ln)
                if m2:
                    out["peak_used_kb"] = m2.group(1)
                    continue
                m3 = re.match(r"^Net increase:\s+([0-9]+)\s*KB", ln)
                if m3:
                    out["net_increase_kb"] = m3.group(1)
                    continue
                m4 = re.match(
                    r"^Elapsed time:\s*([0-9]{2}:[0-9]{2}:[0-9]{2})(?:\s*\(([\d\.]+)\s*h\))?",
                    ln,
                )
                if m4:
                    out["elapsed_hms"] = m4.group(1)
                    if m4.group(2):
                        out["elapsed_hours"] = m4.group(2)
                    continue
        return out
    except FileNotFoundError:
        return out
    except Exception:
        return out


def select_fields(row: Dict[str, str], fields: List[str]) -> Dict[str, str]:
    if not row:
        return {f: "NA" for f in fields}
    out = {}
    for f in fields:
        out[f] = row.get(f, "NA") if row.get(f, "") != "" else "NA"
    return out


def gather_species(root: str, code: str, idx1: int, species: str) -> Dict:
    """
    Build a record for one species.
    idx1: 1-based index in species list, used to create code_full = f"{code}{idx1:02d}"
    """
    species_dir = os.path.join(root, species)
    run_root = os.path.join(species_dir, "v6", "0")

    # Common paths
    p_tmp_l_sra = os.path.join(species_dir, "tmp", "l.sra.txt")
    p_tmp_s_sra = os.path.join(species_dir, "tmp", "s.sra.txt")

    p_sum_long = os.path.join(run_root, "summary-data", "l.fq.seqkit.stats.ta.txt")
    p_sum_s1 = os.path.join(run_root, "summary-data", "s_1.fq.seqkit.stats.ta.txt")
    p_sum_s2 = os.path.join(run_root, "summary-data", "s_2.fq.seqkit.stats.ta.txt")

    p_pt_ref_fa = os.path.join(run_root, "ncbi-ptdna", "ptdna-reference.fa")
    p_pt_acc_file = os.path.join(
        run_root, "ncbi-ptdna", "00-bioproject", "2-mtdna.accession"
    )
    p_mt_ref_fa = os.path.join(run_root, "ncbi-mtdna", "mtdna-reference.fa")
    p_mt_acc_file = os.path.join(
        run_root, "ncbi-mtdna", "00-bioproject", "2-mtdna.accession"
    )

    p_pt_polap_fa = os.path.join(run_root, "polap-assemble", "pt.1.fasta")
    p_mt_polap_fa = os.path.join(run_root, "polap-assemble", "mt.1.fasta")

    # Summaries (time/mem)
    p_sum_polap = os.path.join(run_root, "summary-polap-assemble.txt")
    p_sum_ptgaul = os.path.join(run_root, "summary-ptgaul.txt")
    p_sum_mtgaul = os.path.join(run_root, "summary-mtgaul.txt")
    p_sum_tippo = os.path.join(run_root, "summary-tippo.txt")
    p_sum_oatk = os.path.join(run_root, "summary-oatk.txt")

    # Downsample/used data
    long_sra = slurp_first_line(p_tmp_l_sra)
    short_sra = slurp_first_line(p_tmp_s_sra)

    # Derive used-data TSV paths from SRA IDs when possible
    p_used_long = (
        "NA"
        if long_sra == "NA"
        else os.path.join(
            run_root, "data-downsample-long", f"{long_sra}.fastq.seqkit.stats.ta.tsv"
        )
    )
    p_used_s1 = (
        "NA"
        if short_sra == "NA"
        else os.path.join(
            run_root,
            "data-downsample-short",
            f"{short_sra}_1.fastq.seqkit.stats.ta.tsv",
        )
    )
    p_used_s2 = (
        "NA"
        if short_sra == "NA"
        else os.path.join(
            run_root,
            "data-downsample-short",
            f"{short_sra}_2.fastq.seqkit.stats.ta.tsv",
        )
    )

    # Parse seqkit summaries (long + short)
    long_sum_row = parse_seqkit_tsv_first_row(p_sum_long)
    s1_sum_row = parse_seqkit_tsv_first_row(p_sum_s1)
    s2_sum_row = parse_seqkit_tsv_first_row(p_sum_s2)

    # Select required fields
    wanted = ["num_seqs", "sum_len", "AvgQual"]
    long_sum_sel = select_fields(long_sum_row, wanted)
    s1_sum_sel = select_fields(s1_sum_row, wanted)
    s2_sum_sel = select_fields(s2_sum_row, wanted)

    # Parse used/actual data (downsampled) – if path is "NA" or missing, fields become "NA"
    used_long_sel = select_fields(
        parse_seqkit_tsv_first_row(p_used_long) if p_used_long != "NA" else {}, wanted
    )
    used_s1_sel = select_fields(
        parse_seqkit_tsv_first_row(p_used_s1) if p_used_s1 != "NA" else {}, wanted
    )
    used_s2_sel = select_fields(
        parse_seqkit_tsv_first_row(p_used_s2) if p_used_s2 != "NA" else {}, wanted
    )

    # NCBI references & length checks
    ncbi_pt = check_ncbi_reference(p_pt_acc_file, p_pt_ref_fa)
    ncbi_mt = check_ncbi_reference(p_mt_acc_file, p_mt_ref_fa)

    # Polap assemblies (pt.1.fasta / mt.1.fasta)
    pt_seq_count, pt_total, _ = parse_fasta_stats_and_first_header(p_pt_polap_fa)
    mt_seq_count, mt_total, _ = parse_fasta_stats_and_first_header(p_mt_polap_fa)

    polap_pt = {
        "fasta": rel_or_na(p_pt_polap_fa, root),
        "seq_count": pt_seq_count if pt_seq_count > 0 else "NA",
        "total_bases": pt_total if pt_total > 0 else "NA",
    }
    polap_mt = {
        "fasta": rel_or_na(p_mt_polap_fa, root),
        "seq_count": mt_seq_count if mt_seq_count > 0 else "NA",
        "total_bases": mt_total if mt_total > 0 else "NA",
    }

    # Summary time/mem blocks
    summaries = {
        "polap_assemble": parse_summary_time_mem(p_sum_polap),
        "ptgaul": parse_summary_time_mem(p_sum_ptgaul),
        "mtgaul": parse_summary_time_mem(p_sum_mtgaul),
        "tippo": parse_summary_time_mem(p_sum_tippo),
        "oatk": parse_summary_time_mem(p_sum_oatk),
    }

    record = {
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
            "long": {
                **long_sum_sel,
                "source": rel_or_na(p_sum_long, root),
            },
            "short1": {
                **s1_sum_sel,
                "source": rel_or_na(p_sum_s1, root),
            },
            "short2": {
                **s2_sum_sel,
                "source": rel_or_na(p_sum_s2, root),
            },
        },
        "reads_used": {
            "long": {
                **used_long_sel,
                "source": rel_or_na(p_used_long, root),
            },
            "short1": {
                **used_s1_sel,
                "source": rel_or_na(p_used_s1, root),
            },
            "short2": {
                **used_s2_sel,
                "source": rel_or_na(p_used_s2, root),
            },
        },
        "ncbi": {
            "pt": {
                **ncbi_pt,
                "accession_file": rel_or_na(p_pt_acc_file, root),
                "fasta_path": rel_or_na(p_pt_ref_fa, root),
            },
            "mt": {
                **ncbi_mt,
                "accession_file": rel_or_na(p_mt_acc_file, root),
                "fasta_path": rel_or_na(p_mt_ref_fa, root),
            },
        },
        "polap": {
            "pt": polap_pt,
            "mt": polap_mt,
        },
        "summaries": summaries,
    }
    return record


def rel_or_na(path: str, root: str) -> str:
    if path == "NA":
        return "NA"
    try:
        if not os.path.exists(path):
            return "NA"
        return os.path.relpath(path, root)
    except Exception:
        return "NA"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Create an aggregate JSON of MT/PT assembly progress across species."
    )
    ap.add_argument("--species-codes", required=True, help="Path to species-codes.txt")
    ap.add_argument(
        "--root",
        default=".",
        help="Root folder containing species directories (default: current dir)",
    )
    ap.add_argument("--out", required=True, help="Output JSON file path")
    args = ap.parse_args()

    pairs = read_species_codes(args.species_codes)
    if not pairs:
        print("No species parsed from --species-codes", file=sys.stderr)
        return 2

    species_records: List[Dict] = []
    for i, (code, species) in enumerate(pairs, start=1):
        rec = gather_species(args.root, code, i, species)
        species_records.append(rec)

    model = {
        "meta": {
            "generated_at": _iso_now(),
            "version": VERSION,
            "root": os.path.abspath(args.root),
            "species_count": len(species_records),
            "source_species_codes": os.path.abspath(args.species_codes),
        },
        "species": species_records,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out:
        json.dump(model, out, indent=2)
    print(f"Wrote {args.out}")
    return 0


def _iso_now() -> str:
    try:
        from datetime import datetime, timezone

        return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except Exception:
        return ""


if __name__ == "__main__":
    raise SystemExit(main())
