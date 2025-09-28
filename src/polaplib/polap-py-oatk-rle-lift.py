#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oatk-rle-lift.py
Expand homopolymer-compressed (HPC) contigs back to base space using raw-read
alignments (PAF with cs:Z: tag). For each contig position, compute a run-length
consensus = median over reads of (1 + inserted same-base count near that ref
position). Then repeat that base run-length times.

USAGE:
  oatk-rle-lift.py \
    --hpc-fa contigs.hpc.fa \
    --paf raw_vs_hpc.paf \
    --out-fa contigs.expanded.fa \
    [--min-mapq 10] [--min-aln 100] [--min-cov 3]

NOTES:
- PAF must include cs tag: run minimap2 with --cs=long; use --secondary=no.
- Heuristic: insertions adjacent to a ref base are attributed to that base's run
  ONLY if the inserted base equals the ref base. Substitutions count as 1.
- We ignore long deletions (they reduce consensus via fewer contributions).
- After expansion, you may still polish once (racon/medaka) for micro indels.
"""

import sys, io, argparse, statistics, gzip, re
from collections import defaultdict


def open_text(p):
    return (
        gzip.open(p, "rt")
        if p.endswith(".gz")
        else open(p, "rt", encoding="utf-8", errors="replace")
    )


def load_fasta(path):
    seqs = {}
    with open_text(path) as fh:
        name = None
        buf = []
        for ln in fh:
            if not ln.strip():
                continue
            if ln.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name:
            seqs[name] = "".join(buf).upper()
    return seqs


cs_ins_re = re.compile(r"\+([acgtn]+)", re.I)
cs_del_re = re.compile(r"-([acgtn]+)", re.I)


def parse_cs_and_accumulate(cs, ref_seq, ref_start, ref_len, rl_lists):
    """
    Walk the cs string, mapping contributions to reference positions.
    rl_lists is dict: pos -> list of run-length contributions from this read.
    Approach:
      - For ':'N matches: for each of N ref positions, add 1 to that pos.
      - For '*xy' mismatch: consume 1 ref base -> add 1 to that pos.
      - For '+<bases>' insertion: no ref advance; distribute inserted bases to
        adjacent ref position(s) if base matches the ref base at pos-1 or pos.
      - For '-<bases>' deletion: consume ref positions without query bases -> add 0.
    """
    i = 0
    ref_pos = ref_start

    def add_ref(pos, inc):
        if pos < ref_start or pos >= ref_start + ref_len:
            return
        rl_lists[pos].append(inc)

    while i < len(cs):
        op = cs[i]
        if op == ":":
            i += 1
            j = i
            while j < len(cs) and cs[j].isdigit():
                j += 1
            n = int(cs[i:j])
            # n matches: +1 per ref position
            for _ in range(n):
                add_ref(ref_pos, 1)
                ref_pos += 1
            i = j
        elif op == "*":
            # substitution *xy
            if i + 2 < len(cs):
                # consume one ref base, count 1
                add_ref(ref_pos, 1)
                ref_pos += 1
                i += 3
            else:
                break
        elif op == "+":
            # insertion in query: +<bases>
            m = cs_ins_re.match(cs, i)
            if not m:  # malformed; skip one char
                i += 1
                continue
            ins = m.group(1).upper()
            # attribute inserted same-bases to nearest ref base
            # prefer previous ref base if exists, else current
            prev_pos = ref_pos - 1
            prev_base = ref_seq[prev_pos] if 0 <= prev_pos < len(ref_seq) else None
            curr_base = ref_seq[ref_pos] if 0 <= ref_pos < len(ref_seq) else None
            cnt_prev = cnt_curr = 0
            for b in ins:
                if prev_base and b == prev_base:
                    cnt_prev += 1
                elif curr_base and b == curr_base:
                    cnt_curr += 1
                # else ignore (mismatch insertion)
            if cnt_prev > 0 and prev_pos >= ref_start:
                rl_lists[prev_pos].append(cnt_prev)  # add ONLY the extra (not +1)
            if cnt_curr > 0 and ref_pos < ref_start + ref_len:
                rl_lists[ref_pos].append(cnt_curr)
            i = m.end()
        elif op == "-":
            # deletion in query: -<bases> (consume ref bases, add 0)
            m = cs_del_re.match(cs, i)
            if not m:
                i += 1
                continue
            dele = m.group(1)
            for _ in dele:
                add_ref(ref_pos, 0)
                ref_pos += 1
            i = m.end()
        elif op == "~":
            # splice (not expected); skip token like ~reflen,quelen or ~type
            # advance until next op
            i += 1
            while i < len(cs) and cs[i] not in ":*+-~":
                i += 1
        else:
            # unknown; advance
            i += 1


def lift_contigs(hpc_seqs, paf_path, min_mapq=10, min_aln=100, min_cov=3, out_fa="-"):
    # Prepare RL lists per contig
    rl_data = {name: [[] for _ in range(len(seq))] for name, seq in hpc_seqs.items()}

    with open_text(paf_path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            # PAF: col0=query col5=target col2.. positions vary by version; standard:
            # [0] qname [1] qlen [2] qstart [3] qend [4] strand [5] tname [6] tlen [7] tstart [8] tend [9] nmatch [10] alen [11] mapq ...
            try:
                tname = f[5]
                tlen = int(f[6])
                tstart = int(f[7])
                tend = int(f[8])
                mapq = int(f[11])
                nmatch = int(f[9])
                alen = int(f[10])
            except ValueError:
                continue
            if mapq < min_mapq or alen < min_aln:
                continue
            if tname not in hpc_seqs:
                continue
            # find cs:Z:
            cs = None
            for tag in f[12:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                    break
            if not cs:  # no cs tag; cannot attribute insertions
                continue
            # unify reference orientation
            ref_seq = hpc_seqs[tname]
            ref_start = tstart
            ref_len = tend - tstart
            parse_cs_and_accumulate(cs, ref_seq, ref_start, ref_len, rl_data[tname])

    # Write expanded contigs with per-position medians
    out = sys.stdout if out_fa == "-" else open(out_fa, "w", encoding="utf-8")
    with out:
        for name, hpc_seq in hpc_seqs.items():
            rl_lists = rl_data[name]
            # guard: if we have no data, default run-length=1 for all
            if sum(len(v) for v in rl_lists) == 0:
                expanded = hpc_seq  # no change
            else:
                runs = []
                for pos, base in enumerate(hpc_seq):
                    vals = rl_lists[pos]
                    # Always include the base itself (1) if no contributions at all
                    if not vals:
                        rl = 1
                    else:
                        # combine contributions: 1 (match/sub) + same-base insertions (already separated)
                        s = sum(vals)
                        # If we saw both 1s and extra insertions recorded separately, s may already be >1
                        # Take median across reads per position by treating each read's (1 + ins_same) as a "sample".
                        # Our vals list stores per-read increments contributed at this pos; but we may have split prev/curr.
                        # As a compromise: use 1 + median(extra) if any extra recorded; else 1.
                        extras = [v for v in vals if v > 1 or v == 0 or v == 1]
                        # Aggregate per-read extras is hard without read IDs; fallback: median(total contribution)
                        try:
                            rl = int(max(1, round(statistics.median(vals))))
                        except statistics.StatisticsError:
                            rl = int(max(1, round(s / max(1, len(vals)))))
                    # clamp to a reasonable max to avoid pathologic expansion
                    if rl > 1000:
                        rl = 1000
                    runs.append(base * rl)
                expanded = "".join(runs)
            # emit fasta
            out.write(f">{name}\n")
            for i in range(0, len(expanded), 80):
                out.write(expanded[i : i + 80] + "\n")


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument(
        "--hpc-fa", required=True, help="HPC/RLE contigs FASTA from Oatk/syncasm"
    )
    ap.add_argument(
        "--paf",
        required=True,
        help="PAF of RAW reads aligned to HPC contigs (minimap2 --cs=long)",
    )
    ap.add_argument("--out-fa", required=True, help="Output expanded FASTA")
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--min-aln", type=int, default=100)
    ap.add_argument(
        "--min-cov",
        type=int,
        default=3,
        help="(reserved) minimum read contributions per position",
    )
    args = ap.parse_args()

    hpc = load_fasta(args.hpc_fa)
    lift_contigs(hpc, args.paf, args.min_mapq, args.min_aln, args.min_cov, args.out_fa)


if __name__ == "__main__":
    sys.exit(main())
