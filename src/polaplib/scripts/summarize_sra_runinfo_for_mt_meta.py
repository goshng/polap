#!/usr/bin/env python3
# File: scripts/summarize_sra_runinfo_for_mt_meta.py
# Version: v0.1.0
"""
Summarize SRA RunInfo CSV (from Entrez Direct efetch -format runinfo)
into a simpler TSV for plant mitochondrial meta-analysis.

USAGE:
  python3 summarize_sra_runinfo_for_mt_meta.py plant_mt.sra_runinfo.csv out.tsv
"""
import csv
import sys


def main(in_csv, out_tsv):
    # Columns we care about (RunInfo has many fields; see NCBI docs)
    wanted = [
        "Run",
        "BioProject",
        "BioSample",
        "LibraryStrategy",
        "LibrarySource",
        "Platform",
        "Model",
        "TaxID",
        "ScientificName",
        "spots",
        "bases",
    ]

    with open(in_csv, newline="") as inf, open(out_tsv, "w", newline="") as outf:
        reader = csv.DictReader(inf)
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(wanted)
        for row in reader:
            writer.writerow([row.get(col, "") for col in wanted])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: summarize_sra_runinfo_for_mt_meta.py in.csv out.tsv\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
