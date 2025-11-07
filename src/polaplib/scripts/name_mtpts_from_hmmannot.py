#!/usr/bin/env python3
# FILE: scripts/name_mtpts_from_hmmannot.py
# VERSION: 0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Parse Oatk's hmmannot (or hmm_annotation) tabular output and add gene-based names to MTPT tracts.

Inputs:
  1) tracts.tsv  : from merge_mt_hits_to_tracts.py (header expected)
     columns (first 11):
       tract_id  mt_contig  start  end  strand  n_hits  mean_pid  mean_alnlen  top_cp_qname  cp_min  cp_max
  2) annot.txt   : hmmannot output (table; comment lines start with '#')
     columns (tsv, per example):
       target_name  accession  query_name  accession  hmmfrom hmmto alifrom alito envfrom envto modlen strand E-value score bias description
     - target_name is the plastid gene (e.g., rbcL, rrn23, ndhB, ...)
     - query_name  is the FASTA sequence identifier of the query (should be the MTPT tract sequence ID)

Logic:
  - For each query_name (tract FASTA record), collect gene hits (target_name) and their alignment start (alifrom or ali to).
  - Sort genes by alignment coordinate (ascending along the tract).
  - For the corresponding tract_id row in tracts.tsv (tract_id should match FASTA record id), append:
      ordered_genes (semicolon-separated) and name (genes joined by '-'; fallback to cp[...] if none).

Usage:
  name_mtpts_from_hmmannot.py <tracts.tsv> <hmmannot.out> > tracts.named.tsv
"""
import sys, csv

if len(sys.argv) < 3:
    sys.stderr.write("Usage: name_mtpts_from_hmmannot.py <tracts.tsv> <hmmannot.out>\n")
    sys.exit(1)

tracts_tsv = sys.argv[1]
annot_txt = sys.argv[2]

# Parse hmmannot
# Map: query_name -> list of (pos_for_order, gene)
hits = {}
with open(annot_txt, "r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split()
        # Be defensive about column availability:
        if len(f) < 13:
            # not enough columns
            continue
        gene = f[0]  # target_name
        qname = f[2]  # query_name (should be MTPT tract FASTA ID)
        try:
            alifrom = int(f[6])
            alito = int(f[7])
        except Exception:
            # If alignment coords are not ints, skip ordering info but still store gene
            alifrom, alito = (0, 0)
        pos = min(alifrom, alito) if (alifrom and alito) else alifrom or alito or 0
        hits.setdefault(qname, []).append((pos, gene))

# Emit new TSV
with open(tracts_tsv, "r", encoding="utf-8", errors="ignore") as fh:
    header = fh.readline().rstrip("\n")
    cols = header.split("\t")
    # Ensure expected base columns count >= 11
    if len(cols) < 11 or cols[0] != "tract_id":
        sys.stderr.write(
            "[name_mtpts_from_hmmannot] Unexpected header in tracts TSV.\n"
        )
    print(header + "\tordered_genes\tname")
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        tract_id = f[0]
        cpname = f[8]
        cpmin = f[9]
        cpmax = f[10]
        genelist = []
        if tract_id in hits and len(hits[tract_id]) > 0:
            ordered = sorted(hits[tract_id], key=lambda x: x[0])
            genelist = [g for _, g in ordered]
        ordered_str = ";".join(genelist) if genelist else ""
        if genelist:
            name = "-".join(genelist)
        else:
            # fallback if no plastid gene detected in this tract
            name = f"cp[{cpname}:{cpmin}-{cpmax}]"
        print(line + "\t" + ordered_str + "\t" + name)
