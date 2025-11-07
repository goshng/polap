#!/usr/bin/env python3
# Version: v0.4.0
"""
Extract reads whose names appear in an ID list (one name per line).
Works on .fastq or .fastq.gz, single-end or paired (run twice).

USAGE
  python3 polap_py_fastq_filter_by_ids.py --fastq in.fq.gz --ids ids.txt --out out.fq.gz
"""
import sys, argparse, gzip

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()

def opener(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def main():
    a = parse_args()
    keep=set(x.strip().split()[0] for x in open(a.ids) if x.strip())
    fin = opener(a.fastq, "rt"); fout = opener(a.out, "wt")
    n_in=n_out=0
    while True:
        h = fin.readline()
        if not h: break
        s = fin.readline(); plus = fin.readline(); q = fin.readline()
        if not s or not plus or not q: break
        n_in += 1
        rid = h[1:].strip().split()[0]
        if rid in keep:
            fout.write(h); fout.write(s); fout.write(plus); fout.write(q)
            n_out += 1
    fin.close(); fout.close()
    sys.stderr.write(f"[fastq_filter_by_ids] kept {n_out}/{n_in} reads\n")

if __name__ == "__main__":
    main()

