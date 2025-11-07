#!/usr/bin/env python3
# Version: v0.1.0
"""
Compute binned coverage from a coordinate-sorted BAM.
Outputs a TSV: contig  bin_start  bin_end  mean_depth  source

USAGE
  python3 polap_py_bam_coverage_bins.py \
    --fasta asm.fa --bam aln.sorted.bam \
    --bin 200 --min-mapq 0 \
    --source allONT \
    --out-tsv coverage.allONT.bin200.tsv
"""
import sys, argparse, os
import pysam

def fasta_lengths(fa):
    lens={}
    faidx = fa + ".fai"
    if os.path.exists(faidx):
        with open(faidx) as fh:
            for ln in fh:
                f=ln.split("\t")
                lens[f[0]]=int(f[1])
    else:
        name=None; L=0
        with open(fa) as fh:
            for ln in fh:
                if ln.startswith(">"):
                    if name is not None: lens[name]=L
                    name=ln[1:].split()[0]; L=0
                else:
                    L += len(ln.strip())
        if name is not None: lens[name]=L
    return lens

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--bin", type=int, default=200)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--source", default="NA")
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    lens = fasta_lengths(a.fasta)
    bam = pysam.AlignmentFile(a.bam, "rb")

    # Pre-allocate simple accumulators per bin per contig
    # For organelles, per-base arrays are fine, but do bin-on-the-fly for less RAM.
    with open(a.out_tsv, "w") as oh:
        oh.write("contig\tstart\tend\tmean_depth\tsource\n")
        for r in bam.header.references:
            if r not in lens: continue
            L = lens[r]
            nb = (L + a.bin - 1)//a.bin
            sums = [0]*nb
            counts = [0]*nb  # number of positions in bin (we fill later)
            # Initialize counts = bin sizes (except last)
            for bi in range(nb-1):
                counts[bi] = a.bin
            counts[-1] = L - a.bin*(nb-1)

            # pileup and count coverage with MAPQ filter and primary-only
            for col in bam.pileup(r, truncate=True, min_base_quality=0, stepper="all"):
                pos = col.reference_pos  # 0-based
                # count only reads meeting MAPQ and primary/non-supplementary
                cov = 0
                for pr in col.pileups:
                    rd = pr.alignment
                    if rd.is_unmapped or rd.is_secondary or rd.is_supplementary: continue
                    if rd.mapping_quality < a.min_mapq: continue
                    if pr.is_del or pr.is_refskip: continue
                    cov += 1
                bi = pos // a.bin
                sums[bi] += cov

            for bi in range(nb):
                s = bi*a.bin
                e = min(L, s + a.bin)
                denom = counts[bi] if counts[bi] > 0 else 1
                md = sums[bi] / denom
                oh.write(f"{r}\t{s}\t{e}\t{md:.6f}\t{a.source}\n")

    bam.close()
    sys.stderr.write(f"[cov_bins] wrote {a.out_tsv}\n")

if __name__ == "__main__":
    main()
