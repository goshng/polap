#!/usr/bin/env python3
# Version: v0.4.0
"""
Filter minimap2 PAF records by identity and aligned length.

Identity is computed as nmatch/alnLen using PAF cols 10 and 11 (0-based 9 and 10).
Keeps records with (nmatch/alnLen >= --min-ident) and (alnLen >= --min-alen).

USAGE
  # from stdin to stdout (compatible with your existing shell usage)
  python3 polap_py_paf_filter.py 0.91 3000 < in.paf > out.paf

  # or with flags
  python3 polap_py_paf_filter.py --min-ident 0.91 --min-alen 3000 --in in.paf --out out.paf
"""
import sys, argparse

def parse_args(argv):
    # support legacy positional form: script MIN_IDENT MIN_ALEN
    if len(argv) == 3 and argv[1].replace('.','',1).isdigit() and argv[2].isdigit():
        return argparse.Namespace(min_ident=float(argv[1]), min_alen=int(argv[2]), inp='-', out='-')
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-ident', type=float, required=True)
    ap.add_argument('--min-alen',  type=int,   required=True)
    ap.add_argument('--in',  dest='inp',  default='-')
    ap.add_argument('--out', dest='out',  default='-')
    return ap.parse_args(argv[1:])

def open_io(path, mode):
    return (sys.stdin if path=='-' else open(path, mode))

def main(argv):
    args = parse_args(argv)
    fin  = open_io(args.inp, 'r')
    fout = open_io(args.out, 'w')
    n_in=n_out=0
    for ln in fin:
        if not ln or ln[0]=='#': continue
        parts = ln.rstrip('\n').split('\t')
        if len(parts) < 12:
            continue
        try:
            nmatch = int(parts[9])
            alen   = int(parts[10])
        except Exception:
            continue
        if alen <= 0: continue
        ident = nmatch / alen
        n_in += 1
        if ident >= args.min_ident and alen >= args.min_alen:
            fout.write(ln)
            n_out += 1
    if fout is not sys.stdout:
        fout.close()
    if fin is not sys.stdin:
        fin.close()
    # optional stderr summary
    sys.stderr.write(f"[paf_filter] kept {n_out}/{n_in} records (min_ident={args.min_ident}, min_alen={args.min_alen})\n")

if __name__ == "__main__":
    main(sys.argv)

