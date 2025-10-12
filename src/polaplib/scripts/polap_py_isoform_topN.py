#!/usr/bin/env python3
# Version: v0.7.0
"""
Estimate top-N global isomer states from per-repeat recombinant fractions (no 2^k blow-up).
Assume independence across repeats as a first-order approximation.

Input
-----
repeat_fractions.tsv  (label pair_id orient P_raw R_raw P_norm R_norm f_hat CI_lo CI_hi ...)

Output
------
isoforms_top.tsv
  rank  logprob  isoform_bits  pairs
Where isoform_bits is a 0/1 vector in repeat order (0=parental-dominant, 1=recombined-dominant).

Usage
-----
python3 polap_py_isoform_topN.py --repeat-fractions repeat_fractions.tsv --topN 10 --beam 64 --out-tsv isoforms_top.tsv
"""
from __future__ import annotations
import argparse, math, logging
from typing import List, Tuple

__version__ = "0.7.0"


def main():
    ap = argparse.ArgumentParser(
        description="Top-N isoform estimator from per-repeat marginals."
    )
    ap.add_argument("--repeat-fractions", required=True)
    ap.add_argument("--topN", type=int, default=10)
    ap.add_argument("--beam", type=int, default=64)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument(
        "-v",
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    ap.add_argument("--version", action="store_true")
    a = ap.parse_args()
    if a.version:
        print(__version__)
        return
    logging.basicConfig(
        level=getattr(logging, a.log_level), format="%(levelname)s: %(message)s"
    )

    reps: List[Tuple[str, float]] = []  # (pair, f_hat)
    with open(a.repeat_fractions) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        if "pair_id" not in idx or "f_hat" not in idx:
            raise SystemExit(
                "ERROR: repeat_fractions.tsv must contain columns 'pair_id' and 'f_hat'."
            )
        for ln in fh:
            if not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            pair = f[idx["pair_id"]]
            try:
                fhv = max(1e-6, min(1 - 1e-6, float(f[idx["f_hat"]])))
            except:
                fhv = 0.5
            reps.append((pair, fhv))

    # Sort by uncertainty (closest to 0.5 first) improves beam pruning
    reps.sort(key=lambda x: abs(0.5 - x[1]))

    # Beam search over bit assignments; state = (logprob, [(pair,bit,f)])
    beam: List[Tuple[float, List[Tuple[str, int, float]]]] = [(0.0, [])]
    for pair, fhat in reps:
        new: List[Tuple[float, List[Tuple[str, int, float]]]] = []
        for logp, path in beam:
            # 0=parental with prob (1-f), 1=recombined with prob f
            new.append((logp + math.log(1.0 - fhat), path + [(pair, 0, fhat)]))
            new.append((logp + math.log(fhat), path + [(pair, 1, fhat)]))
        # keep best 'beam' states
        new.sort(key=lambda x: x[0], reverse=True)
        beam = new[: a.beam]

    # final topN
    beam.sort(key=lambda x: x[0], reverse=True)
    top = beam[: a.topN]

    with open(a.out_tsv, "w") as oh:
        oh.write("rank\tlogprob\tisoform_bits\tpairs\n")
        for i, (lp, path) in enumerate(top, start=1):
            bits = "".join("1" if b == 1 else "0" for _, b, _ in path)
            pairs = ";".join(f"{pair}:{b}" for pair, b, _ in path)
            oh.write(f"{i}\t{lp:.3f}\t{bits}\t{pairs}\n")


if __name__ == "__main__":
    main()
