#!/usr/bin/env python3
# ==============================================================================
# scripts/build_benchmark_tables-v0.1.0.py
#
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Prepare skeleton md/table-<set>-0.tsv if not present (or refresh when asked).
#   This is a light shim—your proper TSV builder can replace this later.
#
# Usage   :
#   python build_benchmark_tables-v0.1.0.py --csv POLAP.csv --out md --species-set some
# ==============================================================================
import os, sys, argparse, csv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True, help="CSV with at least 'species' column")
    p.add_argument("--out", required=True, help="output folder (e.g., md)")
    p.add_argument("--species-set", required=True, help="set label for naming")
    return p.parse_args()


def ensure_dir(d):
    os.makedirs(d, exist_ok=True)


def main():
    a = parse_args()
    ensure_dir(a.out)
    out_tsv = os.path.join(a.out, f"table-{a.species_set}-0.tsv")
    if not os.path.exists(out_tsv):
        with open(out_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                [
                    "_species_italic",
                    "_platform",
                    "_memory_gb_ptgaul",
                    "_total_hours_ptgaul",
                    "_memory_gb_polap_assemble_ont_pt",
                    "_total_hours_polap_assemble_ont_pt",
                    "_memory_gb_oatk_nextdenovo_30",
                    "_total_hours_oatk_nextdenovo_30",
                    "_memory_gb_tippo_nextdenovo_hifi",
                    "_total_hours_tippo_nextdenovo_hifi",
                ]
            )
    print(out_tsv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"[ERR] {e}", file=sys.stderr)
        sys.exit(1)
