#!/usr/bin/env python3
# Version: v0.2.0
# Put ec:i first on every L line; optionally drop other tags.
import sys, argparse, gzip


def openg(p, mode="rt"):
    if p == "-":
        return sys.stdin
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


ap = argparse.ArgumentParser()
ap.add_argument("--gfa", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--drop-other-tags", action="store_true")
a = ap.parse_args()

with openg(a.gfa) as f, open(a.out, "w") as w:
    for ln in f:
        if not ln or ln[0] in "#\n":
            continue
        c = ln.rstrip("\n").split("\t")
        if c[0] != "L":
            w.write(ln)
            continue
        tags = c[6:] if len(c) > 6 else []
        idx = None
        for i, t in enumerate(tags):
            if t.startswith("ec:i:"):
                idx = i
                break
        if idx is None:
            sys.stderr.write("ERR: L without ec:i:\n" + ln)
            sys.exit(3)
        ec = tags[idx]
        if a.drop_other_tags:
            new_tags = [ec]
        else:
            new_tags = [ec] + tags[:idx] + tags[idx + 1 :]
        w.write("\t".join(c[:6] + new_tags) + "\n")
