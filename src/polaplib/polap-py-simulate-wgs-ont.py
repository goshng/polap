#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate_wgs_ont.py
Simulate plant WGS reads (ONT-like) from nuclear/mt/pt with higher indel rates.

Example:
python3 simulate_wgs_ont.py \
  --nuclear-size 10m --nuclear-cov 5 \
  --mito-size 500k --mito-cov 50 \
  --plastid-size 150k --plastid-cov 250 \
  --rl-mean 10000 --rl-sd 3000 \
  --sub 0.01 --ins 0.01 --dele 0.01 \
  --emit-refs \
  -o wgs_mix.fq.gz \
  --summary wgs_mix.summary.tsv
"""
import argparse, os, sys, gzip, random
from collections import namedtuple

DNA = "ACGT"


def rand_dna(n):
    return "".join(random.choice(DNA) for _ in range(n))


def parse_num(s):
    s = s.strip().lower()
    mul = 1
    if s.endswith("k"):
        mul, s = 1000, s[:-1]
    elif s.endswith("m"):
        mul, s = 1000_000, s[:-1]
    elif s.endswith("g"):
        mul, s = 1000_000_000, s[:-1]
    return int(float(s) * mul)


def normal_length(mu, sd, min_len=200):
    import math, random

    x = int(random.gauss(mu, sd))
    return max(min_len, x)


def mutate_seq(seq, p_sub, p_ins, p_del):
    out = []
    i = 0
    while i < len(seq):
        r = random.random()
        if r < p_del:
            i += 1
            continue
        elif r < p_del + p_ins:
            out.append(random.choice(DNA))
            continue
        elif r < p_del + p_ins + p_sub:
            b = seq[i]
            choices = [c for c in DNA if c != b]
            out.append(random.choice(choices))
            i += 1
        else:
            out.append(seq[i])
            i += 1
    return "".join(out)


def open_gz_writer(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def main():
    ap = argparse.ArgumentParser(description="Simulate mixed WGS (ONT-like)")
    ap.add_argument("--nuclear-size", required=True)
    ap.add_argument("--nuclear-cov", type=float, required=True)
    ap.add_argument("--mito-size", required=True)
    ap.add_argument("--mito-cov", type=float, required=True)
    ap.add_argument("--plastid-size", required=True)
    ap.add_argument("--plastid-cov", type=float, required=True)
    ap.add_argument("--rl-mean", type=int, required=True)
    ap.add_argument("--rl-sd", type=int, required=True)
    ap.add_argument("--sub", type=float, default=0.01)
    ap.add_argument("--ins", type=float, default=0.01)
    ap.add_argument("--dele", type=float, default=0.01)
    ap.add_argument("--emit-refs", action="store_true")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    random.seed(19)
    sizes = {
        "nuclear": parse_num(args.nuclear_size),
        "mito": parse_num(args.mito_size),
        "plastid": parse_num(args.plastid_size),
    }
    covs = {
        "nuclear": args.nuclear_cov,
        "mito": args.mito_cov,
        "plastid": args.plastid_cov,
    }
    refs = {k: rand_dna(v) for k, v in sizes.items()}

    if args.emit_refs:
        for name, seq in refs.items():
            fa = os.path.splitext(args.out)[0] + f".{name}.fa"
            with open(fa, "w") as w:
                w.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    w.write(seq[i : i + 80] + "\n")

    def simulate_component(name, out_fq):
        L = sizes[name]
        cov = covs[name]
        target_bases = int(L * cov)
        made = 0
        reads = 0
        while made < target_bases:
            rl = normal_length(args.rl_mean, args.rl_sd, min_len=300)
            if rl > L:
                rl = L
            start = random.randrange(0, L)
            frag = (
                refs[name][start:] + refs[name][: max(0, (start + rl) - L)]
                if start + rl > L
                else refs[name][start : start + rl]
            )
            mfrag = mutate_seq(frag, args.sub, args.ins, args.dele)
            reads += 1
            out_fq.write(f"@{name}_read{reads}\n{mfrag}\n+\n{'I'*len(mfrag)}\n")
            made += rl
        return reads, made

    fq = open_gz_writer(args.out)
    S = []
    for name in ["nuclear", "mito", "plastid"]:
        r, b = simulate_component(name, fq)
        S.append((name, sizes[name], covs[name], r, b))
    fq.close()
    with open(args.summary, "w") as w:
        w.write("component\tgenome_size\tcov_target\treads\ttotal_bases\n")
        for row in S:
            w.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
