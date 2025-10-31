#!/usr/bin/env python3
# Version: v0.2.0
import sys, argparse


def read_single(path):
    hdr = None
    seq = []
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if hdr is not None:
                    raise SystemExit("ERR: multiple FASTA records")
                hdr = ln.strip()[1:]
            else:
                seq.append(ln.strip())
    if hdr is None:
        raise SystemExit("ERR: no FASTA records")
    return hdr, "".join(seq).upper()


def best_overlap(seq, mn, mx):
    n = len(seq)
    mx = min(mx, n // 2)
    for k in range(mx, mn - 1, -1):
        if seq.endswith(seq[:k]):
            return k
    return 0


def rotate_to_seed(seq, seed):
    seed = seed.upper()
    i = seq.find(seed)
    return seq if i <= 0 else seq[i:] + seq[:i]


def write_fa(path, hdr, seq, w=80):
    with open(path, "w") as o:
        o.write(f">{hdr}\n")
        for i in range(0, len(seq), w):
            o.write(seq[i : i + w] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-ovl", type=int, default=2000)
    ap.add_argument("--max-ovl", type=int, default=20000)
    ap.add_argument("--rotate-seed", default="")
    a = ap.parse_args()
    hdr, seq = read_single(a.inp)
    k = best_overlap(seq, a.min_ovl, a.max_ovl)
    if k > 0:
        seq = seq[:-k]
        hdr = f"{hdr} circular ovl={k}"
    else:
        print(
            f"WARN: no exact overlap ≥{a.min_ovl}; emitting unchanged", file=sys.stderr
        )
    if a.rotate_seed:
        seq = rotate_to_seed(seq, a.rotate_seed)
        hdr = f"{hdr} rotated@{a.rotate_seed}"
    write_fa(a.out, hdr, seq)


if __name__ == "__main__":
    main()
