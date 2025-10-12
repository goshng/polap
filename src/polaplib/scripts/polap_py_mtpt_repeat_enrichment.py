#!/usr/bin/env python3
# Quantify whether MTPT tracts are enriched near ≥1 kb repeats (within --window bp).
# Input : --mtpt-bed <BED>, --repeat-bed <BED>, --mt-fasta, --window, --permutations
# Output: TSV with observed overlap count/fraction and permutation-based p-value.
import argparse, random, sys
from collections import defaultdict
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--mtpt-bed", required=True)
ap.add_argument("--repeat-bed", required=True)
ap.add_argument("--mt-fasta", required=True)
ap.add_argument("--window", type=int, default=1000)
ap.add_argument("--permutations", type=int, default=1000)
ap.add_argument("--out-tsv", required=True)
args = ap.parse_args()


def read_bed(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            start = int(f[1])
            end = int(f[2])
            if start > end:
                start, end = end, start
            rows.append((chrom, start, end))
    return rows


# contig lengths
clen = {r.id: len(r.seq) for r in SeqIO.parse(args.mt_fasta, "fasta")}

mtpts = read_bed(args.mtpt_bed)
reps = read_bed(args.repeat_bed)

# expand repeats by window; merge naively
exp = defaultdict(list)
for c, s, e in reps:
    s2 = max(0, s - args.window)
    e2 = min(clen.get(c, 10**9), e + args.window)
    exp[c].append((s2, e2))


# merge intervals per contig
def merge(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


exp = {c: merge(v) for c, v in exp.items()}


def overlaps_any(c, s, e):
    if c not in exp:
        return False
    for a, b in exp[c]:
        if s < b and e > a:
            return True
    return False


obs = sum(1 for (c, s, e) in mtpts if overlaps_any(c, s, e))
tot = len(mtpts)
obs_frac = (obs / tot) if tot > 0 else 0.0

# permutation: keep lengths, randomize positions
lens_by_contig = []
for c, s, e in mtpts:
    lens_by_contig.append((c, e - s))

import math

rnd_counts = []
for it in range(args.permutations):
    cnt = 0
    for c, L in lens_by_contig:
        Lc = clen.get(c, 0)
        if Lc <= L:
            continue
        start = random.randint(0, Lc - L)
        end = start + L
        if overlaps_any(c, start, end):
            cnt += 1
    rnd_counts.append(cnt)

mean = sum(rnd_counts) / max(1, len(rnd_counts))
sd = (sum((x - mean) ** 2 for x in rnd_counts) / max(1, len(rnd_counts))) ** 0.5
z = (obs - mean) / sd if sd > 0 else 0.0
# empirical two-sided p
ge = sum(1 for x in rnd_counts if x >= obs)
le = sum(1 for x in rnd_counts if x <= obs)
p_emp = 2 * min(ge, le) / max(1, len(rnd_counts))
if p_emp > 1:
    p_emp = 1.0

with open(args.out_tsv, "w") as oh:
    oh.write(
        "sample\tobs_overlap\tmtpt_total\texp_mean\tz\tp_emp\toverlap_frac\twindow_bp\n"
    )
    oh.write(
        f"{args.sample}\t{obs}\t{tot}\t{mean:.2f}\t{z:.3f}\t{p_emp:.4f}\t{obs_frac:.3f}\t{args.window}\n"
    )
