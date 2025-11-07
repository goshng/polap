#!/usr/bin/env python3
# Version: v0.4.0
"""
Read a SAM/BAM from mapping to a COMBINED reference and output ONT read IDs
that are confidently assigned to the target (reference name startswith --target-prefix).

Filters:
  - primary alignments only (not secondary/supplementary)
  - MAPQ >= --mapq
  - rname startswith --target-prefix (e.g., "mt|" or "cp|")

USAGE
  # SAM from stdin:
  minimap2 -ax map-ont combined.fa ont.fq.gz \
    | samtools view -h -F 0x900 - \
    | python3 polap_py_assign_ont_ids.py --sam - --target-prefix 'mt|' --mapq 30 \
    > ont.target.ids

  # or BAM file:
  python3 polap_py_assign_ont_ids.py --bam combined.bam --target-prefix 'mt|' --mapq 30 > ids.txt
"""
import sys, argparse, pysam

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sam", default=None, help="SAM from stdin or path ('-' for stdin)")
    ap.add_argument("--bam", default=None, help="BAM path (mutually exclusive with --sam)")
    ap.add_argument("--target-prefix", required=True)
    ap.add_argument("--mapq", type=int, default=30)
    return ap.parse_args()

def main():
    a = parse_args()
    if (a.sam is None) == (a.bam is None):
        sys.stderr.write("Provide exactly one of --sam or --bam\n"); sys.exit(2)
    mode = None
    path = None
    if a.sam is not None:
        path = "/dev/stdin" if a.sam == "-" else a.sam
        mode = "r"   # SAM
    else:
        path = a.bam
        mode = "rb"  # BAM
    ids=set()
    af = pysam.AlignmentFile(path, mode)
    for aln in af.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: continue
        if aln.mapping_quality < a.mapq: continue
        rname = af.get_reference_name(aln.reference_id)
        if rname.startswith(a.target_prefix):
            ids.add(aln.query_name)
    af.close()
    for q in sorted(ids):
        print(q)

if __name__ == "__main__":
    main()

