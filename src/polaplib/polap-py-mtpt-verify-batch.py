#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch mtpt verifier (no --verify-root).
Keeps CLI: --fasta, --reads, --mtpt-tsv, --out, --mode

Changes:
- track_len = length from mtpt.tsv (biological interval length)
- pass --interval-len to gold verifier so fractional-gold is interval-based
- track_len column placed after 'end'
"""

import argparse, csv, os, sys, shlex, subprocess
from typing import List, Tuple, Optional


def eprint(*a):
    print(*a, file=sys.stderr)


def run(cmd: str, cwd: Optional[str] = None) -> int:
    return subprocess.run(cmd, shell=True, cwd=cwd).returncode


def find_script(name: str) -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    cand = os.path.join(here, name)
    return cand if os.path.isfile(cand) else name


def read_tsv_one_row(path: str) -> Tuple[list, list]:
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        row = next(r)
        return header, row


def auto_detect_verify_root(anchors: List[Optional[str]]) -> str:
    candidates = []
    for a in anchors:
        if not a:
            continue
        cur = os.path.abspath(os.path.dirname(a))
        for _ in range(7):
            candidates.append(os.path.join(cur, "verify"))
            cur = os.path.dirname(cur)
    candidates += [
        os.path.join(os.getcwd(), "verify"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "verify"),
    ]
    for c in candidates:
        if os.path.isdir(c):
            eprint(f"[info] verify-root: {c}")
            return c
    raise FileNotFoundError(
        "Could not find 'verify' directory. Ensure your earlier step created it."
    )


def list_edge_dirs(verify_root: str) -> list:
    return [
        os.path.join(verify_root, b)
        for b in sorted(os.listdir(verify_root))
        if os.path.isdir(os.path.join(verify_root, b))
        and os.path.isdir(os.path.join(verify_root, b, "reads"))
    ]


def ensure_out_path(out_arg: str) -> str:
    if out_arg.endswith(os.sep) or os.path.isdir(out_arg):
        os.makedirs(out_arg, exist_ok=True)
        return os.path.join(out_arg, "mtpt_verify.summary.tsv")
    parent = os.path.dirname(out_arg) or "."
    os.makedirs(parent, exist_ok=True)
    return out_arg


def get_field(header: list, row: list, key: str, default: str = "") -> str:
    try:
        return row[header.index(key)]
    except Exception:
        return default


def load_mtpt_lengths(mtpt_path: Optional[str]) -> dict:
    """Return {(seqid,start,end): length} from mtpt.tsv; empty dict if missing."""
    m = {}
    if not mtpt_path or not os.path.isfile(mtpt_path):
        eprint(
            "[info] mtpt.tsv not provided or missing; fractional gold will still run (but track_len may be empty)."
        )
        return m
    with open(mtpt_path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        # expected columns: mt_contig, mt_start, mt_end, length, ...
        for row in r:
            try:
                key = (row["mt_contig"], int(row["mt_start"]), int(row["mt_end"]))
                m[key] = int(float(row["length"]))  # be forgiving
            except Exception:
                continue
    return m


def main():
    ap = argparse.ArgumentParser(description="Batch mtpt verifier (no --verify-root)")
    ap.add_argument("--fasta", required=False)
    ap.add_argument("--reads", required=False)
    ap.add_argument("--mtpt-tsv", required=False)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", required=False, default="both")
    # expert knobs
    ap.add_argument("--pass-any", type=int, default=1000)
    ap.add_argument("--pass-pri", type=int, default=50)
    ap.add_argument("--border-any", type=int, default=200)
    ap.add_argument("--border-pri", type=int, default=10)
    ap.add_argument("--W", type=int, default=100)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--gold-margin", type=int, default=100)
    ap.add_argument("--gold-min", type=int, default=1)
    ap.add_argument("--frac-threshold", type=float, default=0.80)
    args = ap.parse_args()

    verify_root = auto_detect_verify_root([args.mtpt_tsv, args.fasta, args.reads])
    out_summary = ensure_out_path(args.out)
    twob = find_script("polap-py-mtpt-verify-twoboundary.py")
    gold = find_script("polap-py-mtpt-verify-singleread.py")

    # Load mtpt lengths to drive track_len and pass to gold
    mtpt_lengths = load_mtpt_lengths(args.mtpt_tsv)

    edges = list_edge_dirs(verify_root)
    if not edges:
        eprint(f"[error] No edge directories under {verify_root}")
        sys.exit(3)

    rows, n_ok = [], 0
    for edge_dir in edges:
        base = os.path.basename(edge_dir)
        # Parse seqid,start,end from dirname
        try:
            parts = base.rsplit("_", 2)
            seqid, start, end = parts[0], int(parts[1]), int(parts[2])
        except Exception:
            eprint(f"[WARN] cannot parse edge_dir name: {base}")
            continue

        # Two-boundary (unchanged)
        cmd_tw = (
            f"{shlex.quote(sys.executable)} {shlex.quote(twob)} "
            f"--edge-dir {shlex.quote(edge_dir)} "
            f"--pass-any {args.pass_any} --pass-pri {args.pass_pri} "
            f"--border-any {args.border_any} --border-pri {args.border_pri} "
            f"--W {args.W} --step {args.step}"
        )
        if run(cmd_tw) != 0:
            eprint(f"[WARN] twoboundary failed for {base}")
            continue

        # Determine interval length for this edge (track_len)
        track_len = ""
        key = (seqid, start, end)
        if key in mtpt_lengths:
            track_len = str(mtpt_lengths[key])
            interval_len_arg = f"--interval-len {track_len}"
        else:
            # fallback: algebraic interval length
            track_len = str(end - start + 1)
            interval_len_arg = f"--interval-len {track_len}"

        # Single-read (interval-based fractional gold)
        cmd_gold = (
            f"{shlex.quote(sys.executable)} {shlex.quote(gold)} "
            f"--edge-dir {shlex.quote(edge_dir)} "
            f"{interval_len_arg} "
            f"--margin {args.gold_margin} --gold-min {args.gold_min} "
            f"--frac-threshold {args.frac_threshold}"
        )
        _ = run(cmd_gold)  # best-effort

        # Read outputs
        twob_tsv = os.path.join(edge_dir, "twoboundary.tsv")
        gold_tsv = os.path.join(edge_dir, "singleread.tsv")
        if not os.path.isfile(twob_tsv):
            eprint(f"[WARN] missing twoboundary.tsv in {base}; skipping")
            continue

        try:
            tw_h, tw_r = read_tsv_one_row(twob_tsv)
        except Exception as e:
            eprint(f"[WARN] cannot read twoboundary.tsv for {base}: {e}")
            continue

        # Defaults if singleread.tsv is missing
        gold_reads, gold_dec = "0", "NO_GOLD"
        frac_thr, frac_reads, frac_max = f"{args.frac_threshold:.2f}", "0", "0.000"
        if os.path.isfile(gold_tsv):
            try:
                g_h, g_r = read_tsv_one_row(gold_tsv)
                gold_reads = get_field(g_h, g_r, "gold_reads", "0")
                gold_dec = get_field(g_h, g_r, "gold_decision", "NO_GOLD")
                # singleread now outputs track_len (same as interval len)
                # but we already have it; keep the batch’s track_len as source of truth
                frac_thr = get_field(g_h, g_r, "frac_threshold", frac_thr)
                frac_reads = get_field(g_h, g_r, "frac_reads", "0")
                frac_max = get_field(g_h, g_r, "frac_max", "0.000")
            except Exception as e:
                eprint(f"[WARN] cannot read singleread.tsv for {base}: {e}")

        gi = lambda k: get_field(tw_h, tw_r, k, "")
        rows.append(
            [
                base,
                seqid,
                str(start),
                str(end),
                track_len,  # << after 'end'
                gi("pos_left"),
                gi("pos_right"),
                gi("w_left"),
                gi("w_right"),
                gi("twob_left"),
                gi("twob_left_primary"),
                gi("twob_right"),
                gi("twob_right_primary"),
                gi("twob_decision"),
                gi("twob_note"),
                gold_reads,
                gold_dec,
                f"{float(frac_thr):.2f}",
                frac_reads,
                frac_max,
            ]
        )
        n_ok += 1

    header = [
        "edge_dir",
        "seqid",
        "start",
        "end",
        "track_len",  # interval length from mtpt.tsv
        "pos_left",
        "pos_right",
        "w_left",
        "w_right",
        "twob_left",
        "twob_left_primary",
        "twob_right",
        "twob_right_primary",
        "twob_decision",
        "twob_note",
        "gold_reads",
        "gold_decision",
        "frac_threshold",
        "frac_reads",
        "frac_max",
    ]
    with open(out_summary, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)

    eprint(f"[OK] wrote {out_summary} with {n_ok} rows (from {len(edges)} edges)")


if __name__ == "__main__":
    main()
