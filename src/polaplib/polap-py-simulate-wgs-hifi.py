#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate_wgs_hifi.py
Simulate plant WGS reads (HiFi-like) from three references: nuclear, mito, plastid.

Example:
python3 simulate_wgs_hifi.py \
  --nuclear-size 10m --nuclear-cov 5 \
  --mito-size 500k --mito-cov 50 \
  --plastid-size 150k --plastid-cov 250 \
  --rl-mean 10000 --rl-sd 1500 \
  --sub 0.002 --ins 0.0003 --dele 0.0003 \
  --emit-refs \
  -o wgs_mix.fq.gz \
  --summary wgs_mix.summary.tsv
"""
import argparse, os, sys, gzip, random, math
from math import ceil
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


def normal_length(mu, sd, min_len=100):
    x = int(random.gauss(mu, sd))
    return max(min_len, x)


def mutate_seq(seq, p_sub, p_ins, p_del):
    # simple i.i.d. per-base model, ONT/HiFi-ish; qualities fixed as 'I'
    out = []
    i = 0
    while i < len(seq):
        r = random.random()
        if r < p_del:
            # delete base
            i += 1
            continue
        elif r < p_del + p_ins:
            # insert a random base, do not consume template
            out.append(random.choice(DNA))
            continue
        elif r < p_del + p_ins + p_sub:
            # substitute
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
    ap = argparse.ArgumentParser(
        description="Simulate mixed WGS (HiFi-like) from nuclear/mt/pt"
    )
    ap.add_argument("--nuclear-size", required=True, help="e.g. 10m")
    ap.add_argument("--nuclear-cov", type=float, required=True)
    ap.add_argument("--mito-size", required=True, help="e.g. 500k")
    ap.add_argument("--mito-cov", type=float, required=True)
    ap.add_argument("--plastid-size", required=True, help="e.g. 150k")
    ap.add_argument("--plastid-cov", type=float, required=True)
    ap.add_argument("--rl-mean", type=int, required=True, help="read length mean")
    ap.add_argument("--rl-sd", type=int, required=True, help="read length sd")
    ap.add_argument("--sub", type=float, default=0.002)
    ap.add_argument("--ins", type=float, default=0.0003)
    ap.add_argument("--dele", type=float, default=0.0003)
    ap.add_argument("--emit-refs", action="store_true")
    ap.add_argument("-o", "--out", required=True, help="output fastq(.gz)")
    ap.add_argument("--summary", required=True, help="summary TSV")
    args = ap.parse_args()

    random.seed(13)

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
    refs = {name: rand_dna(sizes[name]) for name in sizes}

    if args.emit_refs:
        for name, seq in refs.items():
            fa = os.path.splitext(args.out)[0] + f".{name}.fa"
            with open(fa, "w") as w:
                w.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    w.write(seq[i : i + 80] + "\n")

    # expected total bases per component
    exp_bases = {name: covs[name] * sizes[name] for name in sizes}

    fq = open_gz_writer(args.out)
    Q = "I"

    summary_rows = []
    rid = 0
    for name in ["nuclear", "mito", "plastid"]:
        bases_target = exp_bases[name]
        bases_gen = 0
        n_reads = 0
        rseq = refs[name]
        L = len(rseq)
        while bases_gen < bases_target:
            rl = normal_length(args.rl_mean, args.rl_sd, min_len=200)
            if rl > L:
                rl = L
            start = random.randrange(0, L)
            # circular wrap
            if start + rl <= L:
                frag = rseq[start : start + rl]
            else:
                frag = rseq[start:] + rseq[: (start + rl) % L]
            # mutate
            mfrag = mutate_seq(frag, args.sub, args.ins, args.dele)
            rid += 1
            n_reads += 1
            fq.write(f"@{name}_read{rid}\n{mfrag}\n+\n{Q*len(mfrag)}\n")
            bases_gen += rl
        summary_rows.append((name, L, covs[name], n_reads, bases_gen))

    fq.close()
    with open(args.summary, "w") as w:
        w.write("component\tgenome_size\tcov_target\treads\ttotal_bases\n")
        for row in summary_rows:
            w.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
