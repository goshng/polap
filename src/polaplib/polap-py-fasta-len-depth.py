#!/usr/bin/env python3
# polap-py-fasta-len-depth.py
# Parse FASTA (plain or .gz), output TSV: id  Length  Depth
# Depth is taken from "SC:f:<float>" token in the header; default 0.0

import sys, gzip


def open_any(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def depth_from_header(h):
    # tokens separated by whitespace; look for SC:f:<float>
    for tok in h.split():
        if tok.startswith("SC:f:"):
            try:
                return float(tok.split(":", 2)[2])
            except Exception:
                return 0.0
    return 0.0


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: polap-py-fasta-len-depth.py IN.fa[.gz] OUT.tsv\n")
        sys.exit(2)
    infa, outtsv = sys.argv[1], sys.argv[2]

    with open_any(infa) as fin, open(outtsv, "w") as fout:
        fout.write("id\tLength\tDepth\n")
        sid, seq, dep = None, [], 0.0
        for line in fin:
            if line.startswith(">"):
                if sid is not None:
                    fout.write(f"{sid}\t{len(''.join(seq))}\t{dep}\n")
                hdr = line[1:].rstrip()
                sid = hdr.split()[0]
                dep = depth_from_header(hdr)
                seq = []
            else:
                seq.append(line.strip())
        if sid is not None:
            fout.write(f"{sid}\t{len(''.join(seq))}\t{dep}\n")


if __name__ == "__main__":
    main()
