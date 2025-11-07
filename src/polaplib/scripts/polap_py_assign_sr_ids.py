#!/usr/bin/env python3
# Version: v0.4.0
"""
From a SAM/BAM of PE reads mapped to a COMBINED reference, emit read names
whose BOTH ends map (primary, MAPQ>=) to references starting with --target-prefix.

USAGE
  # map without --all for assignment (primary only)
  bowtie2 --very-sensitive -x combined -1 R1.fq.gz -2 R2.fq.gz \
    | samtools view -bh - \
    | python3 polap_py_assign_sr_ids.py --bam - --target-prefix 'mt|' --mapq 30 \
    > sr.target.ids
"""
import sys, argparse, pysam

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, help="BAM or '-' (stdin)")
    ap.add_argument("--target-prefix", required=True)
    ap.add_argument("--mapq", type=int, default=30)
    return ap.parse_args()

def main():
    a = parse_args()
    path = "/dev/stdin" if a.bam == "-" else a.bam
    mode = "rb" if a.bam != "-" else "rb"  # expect BAM; pipe samtools view -bh -
    af = pysam.AlignmentFile(path, mode)
    seenL = set(); seenR = set()
    for aln in af.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: continue
        if aln.mapping_quality < a.mapq: continue
        rname = af.get_reference_name(aln.reference_id)
        if not rname.startswith(a.target_prefix): continue
        q = aln.query_name
        if aln.is_read1: seenL.add(q)
        elif aln.is_read2: seenR.add(q)
    af.close()
    both = sorted(seenL & seenR)
    for q in both:
        print(q)

if __name__ == "__main__":
    main()

