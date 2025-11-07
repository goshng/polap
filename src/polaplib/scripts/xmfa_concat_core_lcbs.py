# FILE: scripts/xmfa_concat_core_lcbs.py
#!/usr/bin/env python3
# VERSION: 0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Concatenate core LCBs from a progressiveMauve XMFA, mask gap-heavy columns,
and emit a supermatrix FASTA + an IQ-TREE/NEXUS partitions file.

Usage
-----
  xmfa_concat_core_lcbs.py <in.xmfa> <min_len_bp> <gap_frac> <out.fna> <out.partitions.nex>

Parameters
----------
in.xmfa       : progressiveMauve alignment in XMFA format.
min_len_bp    : keep LCBs with raw alignment length >= this threshold (e.g., 2000).
gap_frac      : drop alignment columns whose gap fraction > this (e.g., 0.7).
out.fna       : output concatenated FASTA (all taxa).
out.partitions.nex : output NEXUS partitions file with one partition per kept LCB.

Notes
-----
- "Core" == LCB blocks that include *all taxa* present in the XMFA header mapping
  (or, if no mapping lines are present, all taxa observed across blocks).
- We *do not* reverse-complement block sequences; XMFA sequences are stored
  already in aligned orientation.
- We mask columns with gap_fraction > threshold (per block) before concatenation.
- Partitions are reported *after masking* with 1-based coordinates.

References
----------
- Darling AE, Mau B, Perna NT. 2010. progressiveMauve: multiple genome alignment
  with gene gain, loss and rearrangement. PLoS ONE 5(6): e11147.
"""

import sys
import os
import re
from collections import defaultdict, OrderedDict


def die(msg: str):
    sys.stderr.write(f"[xmfa] ERROR: {msg}\n")
    sys.exit(1)


if len(sys.argv) != 6:
    die(
        "Usage: xmfa_concat_core_lcbs.py <in.xmfa> <min_len_bp> <gap_frac> <out.fna> <out.partitions.nex>"
    )

xmfa_path = sys.argv[1]
min_len_bp = int(sys.argv[2])
gap_frac = float(sys.argv[3])
out_fna = sys.argv[4]
out_part = sys.argv[5]

if not os.path.isfile(xmfa_path):
    die(f"XMFA not found: {xmfa_path}")
if gap_frac < 0 or gap_frac > 1:
    die(f"gap_frac must be in [0,1], got {gap_frac}")

# ---------------------------------------------------------------------------
# 1) First pass: parse #Sequence header mapping (seqnum -> taxon name)
#    Fallback to header 'name' basename in block entries if mapping lines absent.
# ---------------------------------------------------------------------------
seqnum_to_taxon = {}
seq_header_re = re.compile(r"^\s*#\s*Sequence\s*([0-9]+)\s*[:=]\s*(\S.*)$", re.I)

with open(xmfa_path, "r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.startswith("#"):
            # stop after header area
            break
        m = seq_header_re.match(line)
        if m:
            idx = int(m.group(1))
            path = m.group(2).strip()
            # taxon name: basename without extension(s)
            b = os.path.basename(path)
            # keep last 2 suffixes for weird names, but strip common fasta ones
            taxon = re.sub(
                r"\.(fa|fna|fasta|fas|fa\.gz|fna\.gz|fasta\.gz)$", "", b, flags=re.I
            )
            seqnum_to_taxon[idx] = taxon

# We'll need a full file pass anyway; reopen for block parsing.
# ---------------------------------------------------------------------------
# 2) Parse XMFA blocks. Each block:
#      >seqnum:start-end strand name
#      <aligned sequence lines...>
#    Blocks separated by a single "=" line.
# ---------------------------------------------------------------------------
header_re = re.compile(r"^>(\d+):\s*([0-9]+)-([0-9]+)\s+([+-])\s+(.*)$")

blocks = []  # list of dicts: {"seq": {taxon: aligned_seq}, "raw_len": L}
taxa_seen_order = OrderedDict()  # preserve order of first appearance

with open(xmfa_path, "r", encoding="utf-8", errors="ignore") as fh:
    in_block = False
    cur_seqs = OrderedDict()  # taxon -> list of seq chunks (join at end)
    cur_seqtaxon = None

    def finish_block():
        # finalize current block (if any)
        nonlocal cur_seqs
        if not cur_seqs:
            return
        # join chunks
        seqs = {t: "".join(chunks) for t, chunks in cur_seqs.items()}
        # verify equal length
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            # malformed alignment block; skip
            cur_seqs = OrderedDict()
            return
        raw_len = next(iter(lengths))
        blocks.append({"seq": seqs, "raw_len": raw_len})
        cur_seqs = OrderedDict()

    for raw in fh:
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            # ignore comments here
            continue
        if line == "=":
            # end of a block
            finish_block()
            in_block = False
            cur_seqtaxon = None
            continue
        if line.startswith(">"):
            # new sequence within the current block
            in_block = True
            m = header_re.match(line)
            taxon = None
            if m:
                seqnum = int(m.group(1))
                name = m.group(5).strip()
                if seqnum in seqnum_to_taxon:
                    taxon = seqnum_to_taxon[seqnum]
                else:
                    # fallback: use provided name (basename)
                    b = os.path.basename(name)
                    taxon = re.sub(
                        r"\.(fa|fna|fasta|fas|fa\.gz|fna\.gz|fasta\.gz)$",
                        "",
                        b,
                        flags=re.I,
                    )
            else:
                # very old/odd header—try to salvage a name after '>'
                taxon = os.path.basename(line[1:].strip().split()[-1])
            taxa_seen_order.setdefault(taxon, True)
            if taxon not in cur_seqs:
                cur_seqs[taxon] = []
            cur_seqtaxon = taxon
            continue
        # sequence line
        if in_block and cur_seqtaxon is not None:
            cur_seqs[cur_seqtaxon].append(line.strip())

    # file end: flush last block if needed
    finish_block()

if not blocks:
    die("No alignment blocks parsed from XMFA.")

all_taxa = list(taxa_seen_order.keys())
if not all_taxa:
    die("Could not infer taxa from XMFA.")

# If header mapping existed, all taxa should match mapping taxa set; otherwise, use discovered taxa.
# We'll define "core" as "block contains *every* taxon in all_taxa".
n_taxa = len(all_taxa)


# ---------------------------------------------------------------------------
# 3) Filter to "core" blocks and by min_len_bp; mask gap-heavy columns; build concatenation.
# ---------------------------------------------------------------------------
def mask_block(block_seqs: dict, gap_thr: float):
    """
    Given dict taxon->aligned string (equal length), remove columns with gap_fraction > gap_thr.
    Returns:
      masked : dict taxon->string
      kept   : number of columns kept
    """
    taxa = sorted(block_seqs.keys())
    alns = [block_seqs[t] for t in taxa]
    L = len(alns[0])
    keep_cols = []
    for i in range(L):
        col = [s[i] for s in alns]
        gaps = sum(1 for c in col if c == "-" or c == "." or c == "N" or c == "n")
        # Here we treat only '-' and '.' as gaps typically; N can be kept. To be conservative, treat only '-' and '.' as gaps.
        gaps = sum(1 for c in col if c == "-" or c == ".")
        if gaps / float(len(col)) <= gap_thr:
            keep_cols.append(i)
    if not keep_cols:
        return {t: "" for t in taxa}, 0
    masked = {t: "".join(block_seqs[t][i] for i in keep_cols) for t in taxa}
    return masked, len(keep_cols)


kept_blocks = []  # tuples: (block_idx, masked_seqs_dict, kept_len)

for idx, blk in enumerate(blocks, start=1):
    seqmap = blk["seq"]
    # require presence of all taxa
    if len(seqmap) != n_taxa:
        continue
    # ensure all taxa keys match our all_taxa set
    if set(seqmap.keys()) != set(all_taxa):
        continue
    if blk["raw_len"] < min_len_bp:
        continue
    masked, kept_len = mask_block(seqmap, gap_frac)
    if kept_len <= 0:
        continue
    kept_blocks.append((idx, masked, kept_len))

if not kept_blocks:
    die(
        "No blocks passed core+length+masking filters. Consider lowering min_len_bp or increasing gap_frac."
    )

# ---------------------------------------------------------------------------
# 4) Emit concatenated FASTA and partitions (after masking)
# ---------------------------------------------------------------------------
# Build per-taxon concatenation in the order of all_taxa (stable across outputs)
concat = {t: [] for t in all_taxa}
partitions = []  # (label, start, end)
cur_start = 1

for idx, masked, kept_len in kept_blocks:
    for t in all_taxa:
        concat[t].append(masked[t])
    start = cur_start
    end = cur_start + kept_len - 1
    partitions.append((f"LCB{idx:04d}", start, end))
    cur_start = end + 1

# Write FASTA
with open(out_fna, "w") as fo:
    for t in all_taxa:
        seq = "".join(concat[t])
        fo.write(f">{t}\n")
        for i in range(0, len(seq), 80):
            fo.write(seq[i : i + 80] + "\n")

# Write NEXUS partitions
with open(out_part, "w") as po:
    po.write("#nexus\nbegin sets;\n")
    for label, s, e in partitions:
        po.write(f"  charset {label} = {s}-{e};\n")
    po.write("end;\n")

# ---------------------------------------------------------------------------
# 5) Summary to stderr
# ---------------------------------------------------------------------------
total_len = partitions[-1][2] if partitions else 0
sys.stderr.write(
    f"[xmfa] taxa: {n_taxa} | kept blocks: {len(kept_blocks)} | concat bp: {total_len}\n"
)
sys.stderr.write(f"[xmfa] out FASTA: {out_fna}\n")
sys.stderr.write(f"[xmfa] out partitions: {out_part}\n")
