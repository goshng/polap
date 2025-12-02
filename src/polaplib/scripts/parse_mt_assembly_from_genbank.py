#!/usr/bin/env python3
# File: scripts/parse_mt_assembly_from_genbank.py
# Version: v0.2.0
"""
Parse GenBank records for plant mitochondrial genomes and extract:
  - accession
  - organism
  - year
  - assembly_method (free text)
  - sequencing_technology (free text)
  - coarse seq_platform
  - coarse assembler
  - bioproject
  - biosample
  - sra_accessions (comma-separated SRR/ERR IDs)

USAGE:
  python3 parse_mt_assembly_from_genbank.py in.gbff out.tsv
"""

import sys
import re
from datetime import datetime
from Bio import SeqIO  # requires Biopython


def extract_year(date_str: str) -> str:
    if not date_str:
        return ""
    for fmt in ("%d-%b-%Y", "%d-%b-%y"):
        try:
            dt = datetime.strptime(date_str, fmt)
            return str(dt.year)
        except ValueError:
            continue
    m = re.search(r"(\d{4})", date_str)
    return m.group(1) if m else ""


def parse_structured_comment(comment: str):
    assembly_method = ""
    seq_tech = ""
    if not comment:
        return assembly_method, seq_tech
    text = re.sub(r"\s+", " ", comment)
    m1 = re.search(r"Assembly Method\s*::\s*([^;]+)", text, flags=re.IGNORECASE)
    if m1:
        assembly_method = m1.group(1).strip()
    m2 = re.search(r"Sequencing Technology\s*::\s*([^;]+)", text, flags=re.IGNORECASE)
    if m2:
        seq_tech = m2.group(1).strip()
    return assembly_method, seq_tech


def categorize_platform(seq_tech: str) -> str:
    s = (seq_tech or "").lower()
    if any(
        x in s for x in ["illumina", "miseq", "hiseq", "novaseq", "nextseq", "ga ii"]
    ):
        return "Illumina"
    if "pacbio" in s or "sequel" in s or "rs ii" in s:
        if "hifi" in s or "ccs" in s:
            return "PacBio HiFi"
        return "PacBio CLR"
    if any(
        x in s
        for x in ["nanopore", "minion", "promethion", "gridion", "oxford nanopore"]
    ):
        return "ONT"
    if "mgi" in s or "bgi" in s or "dnbseq" in s:
        return "MGI"
    if not seq_tech:
        return "Unknown"
    return "Other"


def categorize_assembler(assembly_method: str) -> str:
    s = (assembly_method or "").lower()
    if "flye" in s:
        return "Flye"
    if "canu" in s:
        return "Canu"
    if "hifiasm" in s:
        return "hifiasm"
    if "spades" in s:
        return "SPAdes"
    if "novoplasty" in s:
        return "NOVOPlasty"
    if "getorganelle" in s:
        return "GetOrganelle"
    if "mitoz" in s:
        return "MitoZ"
    if "velvet" in s:
        return "Velvet"
    if "megahit" in s:
        return "MEGAHIT"
    if not assembly_method:
        return "Unknown"
    return "Other"


def parse_links_from_dbxrefs(record):
    """Extract BioProject, BioSample, SRA accessions from record.dbxrefs and comments/DBLINK."""
    bioproject = ""
    biosample = ""
    sra_ids = set()

    # 1) dbxrefs is the most structured place
    for x in record.dbxrefs:
        if x.startswith("BioProject:"):
            bioproject = x.split(":", 1)[1]
        elif x.startswith("BioSample:"):
            biosample = x.split(":", 1)[1]
        elif x.startswith("SRA:"):
            sra_ids.add(x.split(":", 1)[1])

    # 2) fallback: scan comments for SRR/ERR
    comments = record.annotations.get("comment", "")
    if isinstance(comments, list):
        text = " ".join(comments)
    else:
        text = comments or ""

    # This will pick up "Sequence Read Archive: SRR30757341, SRR30757340"
    for m in re.finditer(r"(SRR\d+|ERR\d+)", text):
        sra_ids.add(m.group(1))

    return bioproject, biosample, ",".join(sorted(sra_ids))


def main(in_gb: str, out_tsv: str):
    with open(out_tsv, "w") as out:
        header = [
            "accession",
            "organism",
            "year",
            "assembly_method",
            "sequencing_technology",
            "seq_platform",
            "assembler",
            "bioproject",
            "biosample",
            "sra_accessions",
        ]
        out.write("\t".join(header) + "\n")

        for record in SeqIO.parse(in_gb, "genbank"):
            acc = record.id
            org = record.annotations.get("organism", "")
            date = record.annotations.get("date", "")
            year = extract_year(date)

            # COMMENT field(s)
            comments = record.annotations.get("comment", "")
            if isinstance(comments, list):
                comment_text = " ".join(comments)
            else:
                comment_text = comments or ""

            assembly_method, seq_tech = parse_structured_comment(comment_text)
            seq_platform = categorize_platform(seq_tech)
            assembler = categorize_assembler(assembly_method)

            bioproject, biosample, sra_ids_str = parse_links_from_dbxrefs(record)

            row = [
                acc,
                org,
                year,
                assembly_method,
                seq_tech,
                seq_platform,
                assembler,
                bioproject,
                biosample,
                sra_ids_str,
            ]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: parse_mt_assembly_from_genbank.py in.gbff out.tsv\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
