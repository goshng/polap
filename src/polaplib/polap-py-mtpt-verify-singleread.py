#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-read "gold" verifier with interval-based fractional gold.

Outputs (tab-separated):
seqid  start  end  gold_reads  gold_decision  track_len  frac_threshold  frac_reads  frac_max
- track_len = interval length (passed in), matching mtpt.tsv
- frac_* computed over INTERVAL ONLY (the middle of JWHOLE)

Strict GOLD still uses JWHOLE margins by default.
"""
import argparse, os, re, csv
import pysam


def parse_edge_dir_name(basename: str):
    m = re.match(r"^(.*)_(\d+)_(\d+)$", basename)
    if not m:
        raise ValueError(f"Unrecognized edge directory name: {basename}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def jwhole_len(bam: pysam.AlignmentFile) -> int:
    return bam.get_reference_length("JWHOLE") if "JWHOLE" in bam.references else 0


def overlap_len(a1, a2, b1, b2):
    """Inclusive coords; returns non-negative overlap length."""
    lo = max(a1, b1)
    hi = min(a2, b2)
    return max(0, hi - lo + 1)


def main():
    ap = argparse.ArgumentParser(
        description="Gold verifier w/ interval-based fractional gold on JWHOLE."
    )
    ap.add_argument("--edge-dir", required=True, help="verify/edge_*/ directory")
    ap.add_argument(
        "--interval-len", type=int, required=True, help="interval length from mtpt.tsv"
    )
    ap.add_argument(
        "--margin", type=int, default=100, help="end margin for strict GOLD on JWHOLE"
    )
    ap.add_argument(
        "--gold-min", type=int, default=1, help="reads required to call GOLD"
    )
    ap.add_argument(
        "--frac-threshold",
        type=float,
        default=0.80,
        help="fractional-gold threshold (0-1)",
    )
    args = ap.parse_args()

    reads_dir = os.path.join(args.edge_dir, "reads")
    jwhole_bam = None
    for f in os.listdir(reads_dir):
        if f.startswith("jwhole.") and f.endswith(".bam"):
            jwhole_bam = os.path.join(reads_dir, f)
            break

    gold_reads = 0
    gold_dec = "NO_GOLD"
    frac_reads = 0
    frac_max = 0.0
    track_len = int(args.interval_len)

    if jwhole_bam and os.path.isfile(jwhole_bam):
        with pysam.AlignmentFile(jwhole_bam, "rb") as bam:
            if "JWHOLE" in bam.references:
                Lwhole = jwhole_len(bam)
                # derive flank and interval coords on JWHOLE
                if Lwhole >= track_len:
                    flank = (Lwhole - track_len) // 2
                else:
                    flank = 0
                I1 = flank + 1
                I2 = flank + track_len

                # Strict GOLD still uses JWHOLE margins
                left_lim = 1 + args.margin
                right_lim = Lwhole - args.margin

                names_gold = set()
                names_frac = set()
                for read in bam.fetch("JWHOLE"):
                    if read.is_unmapped:
                        continue
                    s1 = (read.reference_start or 0) + 1
                    e1 = read.reference_end or 0
                    # GOLD: near end-to-end on JWHOLE
                    if s1 <= left_lim and e1 >= right_lim:
                        names_gold.add(read.query_name)
                    # Fractional over INTERVAL
                    ov = overlap_len(s1, e1, I1, I2)
                    if track_len > 0:
                        frac = ov / track_len
                        if frac > frac_max:
                            frac_max = frac
                        if frac >= args.frac_threshold:
                            names_frac.add(read.query_name)

                gold_reads = len(names_gold)
                frac_reads = len(names_frac)

    if gold_reads >= args.gold_min:
        gold_dec = "GOLD"

    seqid, start, end = parse_edge_dir_name(os.path.basename(args.edge_dir))
    out_tsv = os.path.join(args.edge_dir, "singleread.tsv")
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "seqid",
                "start",
                "end",
                "gold_reads",
                "gold_decision",
                "track_len",
                "frac_threshold",
                "frac_reads",
                "frac_max",
            ]
        )
        w.writerow(
            [
                seqid,
                start,
                end,
                gold_reads,
                gold_dec,
                track_len,
                f"{args.frac_threshold:.2f}",
                frac_reads,
                f"{frac_max:.3f}",
            ]
        )

    print(
        f"{os.path.basename(args.edge_dir)}\t{gold_dec}\tgold_reads={gold_reads}\tfrac_reads={frac_reads}"
    )


if __name__ == "__main__":
    main()
