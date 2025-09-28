#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-simulate-wgs-reads.py (clean regenerated version)

Simulate WGS reads for mixed plant samples:
 - PacBio HiFi (--platform hifi)
 - ONT simplex (--platform ont)
 - ONT high-quality (--platform onthq, R10.4.1 SUP/duplex-like)

Compartments:
 - Nuclear genome (random only)
 - Mitochondrial genome (--mt-fasta or random)
 - Plastid genome (--pt-fasta or random)

Outputs:
    <prefix>.fastq.gz
    <prefix>.nuclear.fastq.gz
    <prefix>.mito.fastq.gz
    <prefix>.plastid.fastq.gz
"""

import sys, gzip, random, argparse, math
from typing import List, Tuple, Optional, Dict, Callable

DNA = "ACGT"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# ------------------ I/O ------------------


def read_fasta(path: str) -> str:
    seqs, cur = [], []
    with open(path, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            if ln.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur = []
            else:
                cur.append(ln)
        if cur:
            seqs.append("".join(cur).upper())
    if not seqs:
        raise ValueError(f"No FASTA sequence found in {path}")
    return "".join(seqs) if len(seqs) > 1 else seqs[0]


def random_dna(n: int, seed: Optional[int] = None) -> str:
    rng = random.Random(seed) if seed is not None else random
    return "".join(rng.choice(DNA) for _ in range(n))


def write_fastq_gz(
    path: str, entries: List[Tuple[str, str, Optional[str]]], qchar: str
):
    if not entries:
        return
    with gzip.open(path, "wt") as gz:
        for name, seq, q in entries:
            if q is None:
                q = qchar * len(seq)
            gz.write(f"@{name}\n{seq}\n+\n{q}\n")


def parse_bp_with_suffix(x: str) -> int:
    s = x.strip().lower()
    if s.endswith("k"):
        return int(float(s[:-1]) * 1_000)
    if s.endswith("m"):
        return int(float(s[:-1]) * 1_000_000)
    if s.endswith("g"):
        return int(float(s[:-1]) * 1_000_000_000)
    return int(float(s))


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


# ------------------ Length models ------------------


def draw_length_gauss(rng, mean, sd, min_len, max_len):
    while True:
        L = int(rng.gauss(mean, sd))
        if min_len <= L <= max_len:
            return L


def _lognorm_mu_sigma_from_mean_sd(mean, sd):
    sigma2 = math.log(1.0 + (sd * sd) / (mean * mean))
    sigma = math.sqrt(sigma2)
    mu = math.log(mean) - 0.5 * sigma2
    return mu, sigma


def draw_length_lognorm(rng, mean, sd, min_len, max_len):
    mu, sigma = _lognorm_mu_sigma_from_mean_sd(mean, sd)
    while True:
        L = int(rng.lognormvariate(mu, sigma))
        if min_len <= L <= max_len:
            return L


# ------------------ Error models ------------------


def mutate_hifi(seq, sub, ins, dele, rng):
    counts = {"n_sub": 0, "n_ins": 0, "n_del": 0, "n_hp": 0}
    out = []
    i = 0
    L = len(seq)
    while i < L:
        r = rng.random()
        if r < dele:
            counts["n_del"] += 1
            i += 1
            continue
        if r < dele + ins:
            out.append(rng.choice(DNA))
            counts["n_ins"] += 1
        if r < dele + ins + sub:
            b = seq[i]
            out.append(rng.choice([x for x in DNA if x != b]))
            counts["n_sub"] += 1
        else:
            out.append(seq[i])
        i += 1
    return "".join(out), counts


def mutate_ont(seq, sub, ins, dele, hp_rate, rng):
    counts = {"n_sub": 0, "n_ins": 0, "n_del": 0, "n_hp": 0}
    out = []
    i = 0
    L = len(seq)
    while i < L:
        b = seq[i]
        j = i + 1
        while j < L and seq[j] == b:
            j += 1
        run_len = j - i
        stutter_prob = hp_rate * (run_len - 2) / 3.0 if run_len >= 3 else 0
        if rng.random() < stutter_prob:
            counts["n_hp"] += 1
            k = run_len + (1 if rng.random() < 0.5 else -1)
            out.extend(b * max(1, k))
        else:
            for _ in range(run_len):
                r = rng.random()
                if r < dele:
                    counts["n_del"] += 1
                    continue
                if r < dele + ins:
                    out.append(rng.choice(DNA))
                    counts["n_ins"] += 1
                if r < dele + ins + sub:
                    out.append(rng.choice([x for x in DNA if x != b]))
                    counts["n_sub"] += 1
                else:
                    out.append(b)
        i = j
    return "".join(out), counts


# ------------------ Helpers ------------------


def slice_chunk(ref, st, rl, circ):
    n = len(ref)
    if rl <= 0 or n == 0:
        return ""
    if circ:
        st %= n
        return ref[st : st + rl] if st + rl <= n else ref[st:] + ref[: (st + rl) - n]
    else:
        st = max(0, min(st, n - 1))
        en = min(n, st + rl)
        return ref[st:en]


# ------------------ Simulators ------------------


def simulate_long_reads_from_ref(
    ref, total_bases, draw_len, args, mut_fn, rng, prefix, circ, plat
):
    out = []
    acc = 0
    n = 0
    mean, sd, min_len, max_len = args
    while acc < total_bases:
        rl = draw_len(rng, mean, sd, min_len, max_len)
        st = rng.randrange(len(ref))
        chunk = slice_chunk(ref, st, rl, circ)
        mut, _ = mut_fn(chunk, rng)
        name = f"{prefix}_{n:06d}"
        out.append((name, mut, None))
        acc += len(mut)
        n += 1
    return out


# ------------------ CLI ------------------


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--platform", choices=["hifi", "ont", "onthq"], default="hifi")
    ap.add_argument(
        "--nuc-size", type=parse_bp_with_suffix, default=parse_bp_with_suffix("10m")
    )
    ap.add_argument(
        "--mt-size", type=parse_bp_with_suffix, default=parse_bp_with_suffix("500k")
    )
    ap.add_argument(
        "--pt-size", type=parse_bp_with_suffix, default=parse_bp_with_suffix("150k")
    )
    ap.add_argument(
        "--mt-fasta", type=str, default=None, help="Mitochondrial reference FASTA"
    )
    ap.add_argument(
        "--pt-fasta", type=str, default=None, help="Plastid reference FASTA"
    )
    ap.add_argument("--nuc-depth", type=float, default=5)
    ap.add_argument("--mt-depth", type=float, default=50)
    ap.add_argument("--pt-depth", type=float, default=250)
    ap.add_argument("-o", "--out-prefix", type=str, default="wgs.sim")
    ap.add_argument("--seed", type=int, default=43)
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)

    # nuclear always random
    nuc_ref = random_dna(args.nuc_size, args.seed + 11)
    nuc_target = int(args.nuc_size * args.nuc_depth)

    # mito: fasta or random
    if args.mt_fasta:
        mt_ref = read_fasta(args.mt_fasta)
        mt_size = len(mt_ref)
    else:
        mt_ref = random_dna(args.mt_size, args.seed + 22)
        mt_size = args.mt_size
    mt_target = int(mt_size * args.mt_depth)

    # plastid: fasta or random
    if args.pt_fasta:
        pt_ref = read_fasta(args.pt_fasta)
        pt_size = len(pt_ref)
    else:
        pt_ref = random_dna(args.pt_size, args.seed + 33)
        pt_size = args.pt_size
    pt_target = int(pt_size * args.pt_depth)

    # presets
    if args.platform == "hifi":
        draw = draw_length_gauss
        draw_args = (10000, 1500, 4000, 20000)
        mut = lambda s, r: mutate_hifi(s, 0.002, 0.0003, 0.0003, r)
        qchar = chr(33 + 40)
    elif args.platform == "ont":
        draw = draw_length_lognorm
        draw_args = (15000, 5000, 1000, 150000)
        mut = lambda s, r: mutate_ont(s, 0.004, 0.008, 0.008, 0.20, r)
        qchar = chr(33 + 20)
    else:  # onthq
        draw = draw_length_lognorm
        draw_args = (18000, 6000, 1000, 200000)
        mut = lambda s, r: mutate_ont(s, 0.001, 0.0015, 0.0015, 0.05, r)
        qchar = chr(33 + 30)

    nuc = simulate_long_reads_from_ref(
        nuc_ref, nuc_target, draw, draw_args, mut, rng, "nuc", False, args.platform
    )
    mt = simulate_long_reads_from_ref(
        mt_ref, mt_target, draw, draw_args, mut, rng, "mt", True, args.platform
    )
    pt = simulate_long_reads_from_ref(
        pt_ref, pt_target, draw, draw_args, mut, rng, "pt", True, args.platform
    )

    all_reads = nuc + mt + pt
    rng.shuffle(all_reads)
    out = args.out_prefix
    write_fastq_gz(f"{out}.fastq.gz", all_reads, qchar)
    write_fastq_gz(f"{out}.nuclear.fastq.gz", nuc, qchar)
    write_fastq_gz(f"{out}.mito.fastq.gz", mt, qchar)
    write_fastq_gz(f"{out}.plastid.fastq.gz", pt, qchar)
    return 0


if __name__ == "__main__":
    sys.exit(main())
