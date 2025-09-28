#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-spike-numt-nupt.py
Inject NUMT/NUPT fragments into a nuclear reference.

Outputs:
  <out>.augmented.fa   : nuclear genome with inserts (or replacements)
  <out>.inserts.bed    : BED truth (chrom, start, end, name, 0, strand)
  <out>.inserts.tsv    : detailed table with provenance, lengths, divergence

Usage:
  polap-py-spike-numt-nupt.py \
    --nuc nuclear.fa --mt mt.fa --pt pt.fa \
    --numt 20 --numt-len 1000,3000 \
    --nupt 10 --nupt-len 500,1500 \
    --div-sub 0.02 --div-indel 0.002 \
    --inv-frac 0.2 \
    --len-mode uniform \
    --circular-mt --circular-pt \
    --mode insert \
    --seed 1 \
    --out spiked_nuc

Notes:
- Length as "min,max" (uniform) or "mean,sd" if --len-mode gaussian.
- Positions are chosen randomly per chromosome with non-overlap constraints (best-effort with backoff).
- Divergence model: iid per-base substitution + single-base indels; observed divergence reported post-mutation.
- In --mode insert, nuclear assembly length increases; in --mode replace, base count remains unchanged.
- Circular sampling allows NUMT/NUPT fragments to wrap around their references when near sequence ends.
"""
import sys, os, argparse, random, math
from typing import Dict, List, Tuple

DNA = ("A", "C", "G", "T")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_fa(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        name = None
        buf: List[str] = []
        for ln in fh:
            if not ln:
                continue
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper().replace("U", "T")
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            seqs[name] = "".join(buf).upper().replace("U", "T")
    if not seqs:
        raise SystemExit(f"[err] FASTA empty or unreadable: {path}")
    return seqs


def write_fa(seqs: Dict[str, str], path: str, width: int = 80) -> None:
    with open(path, "w") as out:
        for name, s in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(s), width):
                out.write(s[i : i + width] + "\n")


def parse_len(spec: str, mode: str, rng: random.Random) -> int:
    a, b = [float(x) for x in spec.split(",")]
    if mode == "uniform":
        return max(1, int(rng.uniform(a, b)))
    elif mode == "gaussian":
        x = max(1.0, rng.gauss(a, b))
        return int(x)
    else:
        raise ValueError("len-mode must be uniform or gaussian")


def sample_fragment(
    src_seq: str, L: int, circular: bool, rng: random.Random
) -> Tuple[int, int, str]:
    n = len(src_seq)
    if L >= n:
        # if requested > genome, use whole sequence (or tile)
        if not circular:
            return 0, n, src_seq
        # tile-wrapping if circular and L>n
        rep = (L // n) + 2
        s = (src_seq * rep)[:L]
        return 0, L, s
    if circular:
        start = rng.randrange(0, n)
        if start + L <= n:
            frag = src_seq[start : start + L]
        else:
            k = (start + L) - n
            frag = src_seq[start:] + src_seq[:k]
        # pseudo end = start+L mod n
        return start, (start + L) % n, frag
    else:
        start = rng.randrange(0, n - L + 1)
        return start, start + L, src_seq[start : start + L]


def mutate_iid(
    seq: str, p_sub: float, p_indel: float, rng: random.Random
) -> Tuple[str, float, float]:
    """Return mutated seq, observed_sub_rate, observed_indel_rate (per aligned length)."""
    out: List[str] = []
    subs = 0
    dels = 0
    ins = 0
    i = 0
    n = len(seq)
    while i < n:
        # indel?
        if rng.random() < p_indel:
            # 50/50 ins/del single-base
            if rng.random() < 0.5:
                out.append(rng.choice(DNA))  # insertion of random base
                ins += 1
                # do not advance reference (classic model: insertion relative to ref)
                continue
            else:
                # deletion of the current ref base
                i += 1
                dels += 1
                continue
        b = seq[i]
        if rng.random() < p_sub:
            out.append(rng.choice([x for x in DNA if x != b]))
            subs += 1
        else:
            out.append(b)
        i += 1
    mutated = "".join(out)
    # aligned length ~ ref length + insertions (approx), use ref length for rates to be conservative
    obs_sub = subs / float(n) if n else 0.0
    obs_indel = (ins + dels) / float(n) if n else 0.0
    return mutated, obs_sub, obs_indel


def reverse_complement(s: str) -> str:
    return "".join(COMP.get(b, "N") for b in s[::-1])


def try_place_nonoverlap(
    occupied: List[Tuple[int, int]],
    L: int,
    chrom_len: int,
    rng: random.Random,
    max_trials: int = 2000,
) -> int:
    """Return a start position with no overlap, or -1 if failed."""
    # occupied: list of (start,end), end-exclusive
    for _ in range(max_trials):
        pos = rng.randrange(0, chrom_len + 1)  # allow append at end (insert mode)
        ok = True
        for a, b in occupied:
            if max(pos, a) < min(pos + L, b):
                ok = False
                break
        if ok:
            return pos
    return -1


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--nuc", required=True, help="nuclear reference FASTA")
    ap.add_argument("--mt", required=True, help="mitochondrial FASTA")
    ap.add_argument("--pt", required=True, help="plastid FASTA")
    ap.add_argument("--numt", type=int, default=20, help="number of NUMT inserts")
    ap.add_argument("--nupt", type=int, default=10, help="number of NUPT inserts")
    ap.add_argument(
        "--numt-len",
        default="1000,3000",
        help='NUMT length: "min,max" or "mean,sd" (see --len-mode)',
    )
    ap.add_argument(
        "--nupt-len",
        default="500,1500",
        help='NUPT length: "min,max" or "mean,sd" (see --len-mode)',
    )
    ap.add_argument("--len-mode", choices=["uniform", "gaussian"], default="uniform")
    ap.add_argument(
        "--div-sub", type=float, default=0.02, help="substitution divergence per base"
    )
    ap.add_argument(
        "--div-indel",
        type=float,
        default=0.002,
        help="indel rate per base (single-base indels)",
    )
    ap.add_argument(
        "--inv-frac",
        type=float,
        default=0.2,
        help="fraction of inserts inverted (reverse-complement)",
    )
    ap.add_argument(
        "--circular-mt", action="store_true", help="treat mt as circular for sampling"
    )
    ap.add_argument(
        "--circular-pt", action="store_true", help="treat pt as circular for sampling"
    )
    ap.add_argument(
        "--mode",
        choices=["insert", "replace"],
        default="insert",
        help="insert increases length; replace keeps length constant",
    )
    ap.add_argument("--seed", type=int, default=1, help="RNG seed")
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    nuc = read_fa(args.nuc)
    mt = read_fa(args.mt)
    pt = read_fa(args.pt)

    # Prepare occupancy map per chromosome (end-exclusive intervals)
    occupied = {ch: [] for ch in nuc.keys()}

    # Build request list
    reqs: List[Tuple[str, str]] = [("NUMT", "mt")] * args.numt + [
        ("NUPT", "pt")
    ] * args.nupt

    inserts = []  # dicts with all fields and 'seq'
    for i, (tag, src_tag) in enumerate(reqs, start=1):
        Lspec = args.numt_len if tag == "NUMT" else args.nupt_len
        L = parse_len(Lspec, args.len_mode, rng)

        # sample fragment from source
        src_dict = mt if src_tag == "mt" else pt
        src_name = rng.choice(list(src_dict.keys()))
        sseq = src_dict[src_name]
        sst, sen, frag = sample_fragment(
            sseq,
            L,
            circular=(args.circular_mt if src_tag == "mt" else args.circular_pt),
            rng=rng,
        )

        # mutate
        frag_mut, obs_sub, obs_indel = mutate_iid(
            frag, args.div_sub, args.div_indel, rng
        )

        # optionally invert
        strand = "+"
        if rng.random() < args.inv_frac:
            strand = "-"
            frag_mut = reverse_complement(frag_mut)

        # choose chromosome
        ch = rng.choice(list(nuc.keys()))
        chrom_len = len(nuc[ch])

        # place
        if args.mode == "insert":
            pos = try_place_nonoverlap(occupied[ch], len(frag_mut), chrom_len, rng)
            if pos < 0:
                pos = rng.randrange(0, chrom_len + 1)
            occupied[ch].append((pos, pos + len(frag_mut)))
        else:  # replace
            if len(frag_mut) > chrom_len:
                frag_mut = frag_mut[:chrom_len]
            pos_max = max(0, chrom_len - len(frag_mut))
            pos = rng.randrange(0, pos_max + 1)
            occupied[ch].append((pos, pos + len(frag_mut)))

        inserts.append(
            {
                "chrom": ch,
                "start": pos,
                "end": pos + len(frag_mut),
                "name": f"{tag}_{i}",
                "strand": strand,
                "src": src_tag,
                "src_name": src_name,
                "src_start": sst,
                "src_end": sen,
                "req_len": L,
                "mut_len": len(frag_mut),
                "obs_sub": round(obs_sub, 6),
                "obs_indel": round(obs_indel, 6),
                "seq": frag_mut,
            }
        )

    # Apply edits per chromosome in descending order of start
    out_nuc = dict(nuc)  # strings (immutable), safe to reassign
    for ch in out_nuc.keys():
        edits = [x for x in inserts if x["chrom"] == ch]
        if not edits:
            continue
        s = out_nuc[ch]
        if args.mode == "insert":
            for ins in sorted(edits, key=lambda x: x["start"], reverse=True):
                s = s[: ins["start"]] + ins["seq"] + s[ins["start"] :]
        else:
            for ins in sorted(edits, key=lambda x: x["start"], reverse=True):
                a = ins["start"]
                b = ins["end"]
                s = s[:a] + ins["seq"] + s[b:]
        out_nuc[ch] = s

    # Write outputs
    fa_out = f"{args.out}.augmented.fa"
    bed_out = f"{args.out}.inserts.bed"
    tsv_out = f"{args.out}.inserts.tsv"

    # --- BEGIN PATCH: stable POSIX newlines for outputs ---
    fa_out = f"{args.out}.augmented.fa"
    bed_out = f"{args.out}.inserts.bed"
    tsv_out = f"{args.out}.inserts.tsv"

    # Write FASTA with explicit \n wrapping
    with open(fa_out, "w", encoding="utf-8", newline="\n") as out:
        for name, s in out_nuc.items():
            out.write(f">{name}\n")
            for i in range(0, len(s), 80):
                out.write(s[i : i + 80] + "\n")

    # Write BED and TSV with explicit \n per record
    with open(bed_out, "w", encoding="utf-8", newline="\n") as bed, open(
        tsv_out, "w", encoding="utf-8", newline="\n"
    ) as tsv:
        tsv.write(
            "chrom\tstart\tend\tname\tstrand\tsrc\tsrc_name\tsrc_start\tsrc_end\treq_len\tmut_len\tobs_sub\tobs_indel\n"
        )
        for ins in inserts:
            bed.write(
                f"{ins['chrom']}\t{ins['start']}\t{ins['end']}\t{ins['name']}\t0\t{ins['strand']}\n"
            )
            tsv.write(
                "{chrom}\t{start}\t{end}\t{name}\t{strand}\t{src}\t{src_name}\t{src_start}\t{src_end}\t{req_len}\t{mut_len}\t{obs_sub}\t{obs_indel}\n".format(
                    **ins
                )
            )
    # --- END PATCH ---
    #

    sys.stderr.write(f"[ok] augmented nuclear: {fa_out}\n")
    sys.stderr.write(f"[ok] truth BED        : {bed_out}\n")
    sys.stderr.write(f"[ok] details TSV      : {tsv_out}\n")


if __name__ == "__main__":
    main()
