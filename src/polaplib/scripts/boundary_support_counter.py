# FILE: scripts/boundary_support_counter.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Count ONT reads that span each tract boundary with >=anchor bp on both sides, MAPQ>=mapq.

Usage: boundary_support_counter.py <mt.fa> <tracts.bed> <ont.bam> <anchor> <mapq>
Output: tract_id  left_span_reads  right_span_reads
"""
import sys, pysam

if len(sys.argv) < 6:
    sys.stderr.write(
        "Usage: boundary_support_counter.py <mt.fa> <tracts.bed> <ont.bam> <anchor> <mapq>\n"
    )
    sys.exit(1)

fa, bed, bam_path, anchor, mapq = sys.argv[1:6]
anchor = int(anchor)
mapq = int(mapq)
bam = pysam.AlignmentFile(bam_path, "rb")


def spans_boundary(read, chrom, boundary, anchor):
    if read.reference_name != chrom or read.is_unmapped or read.mapping_quality < mapq:
        return False
    # compute aligned blocks on reference
    pos = read.reference_start
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):  # M, =, X
            block_start = pos
            block_end = pos + length
            # left anchor: boundary - anchor .. boundary
            # right anchor: boundary .. boundary + anchor
            left_ok = (block_start <= boundary - anchor) and (block_end >= boundary)
            right_ok = (block_start <= boundary) and (block_end >= boundary + anchor)
            if left_ok and right_ok:
                return True
            pos += length
        elif op in (2, 3):  # D,N ref skip
            pos += length
        elif op in (1, 4, 5):  # I,S,H
            continue
    return False


print("\t".join(["tract_id", "left_span_reads", "right_span_reads"]))
with open(bed) as fh:
    for line in fh:
        if not line.strip():
            continue
        chrom, s, e, tid, score, strand = line.rstrip("\n").split("\t")[:6]
        s = int(s)
        e = int(e)
        left = 0
        right = 0
        for read in bam.fetch(chrom, s - anchor - 5, s + anchor + 5):
            if spans_boundary(read, chrom, s, anchor):
                left += 1
        for read in bam.fetch(chrom, e - anchor - 5, e + anchor + 5):
            if spans_boundary(read, chrom, e, anchor):
                right += 1
        print("\t".join(map(str, [tid, left, right])))
