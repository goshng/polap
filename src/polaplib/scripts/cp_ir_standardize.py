# FILE: scripts/cp_ir_standardize.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Detects largest inverted repeat (IR) via minimap2 self-PAF, keeps one IR copy (IRa),
orders genome as LSC|IRa|SSC, rotates so rbcL starts at pos 1, and writes cp.std.fa.
Also reports IR length and coordinates to stdout (parsed by the orchestrator).

Usage: cp_ir_standardize.py <cp.fa> <plastid_hmms.hmm> <out.fa> <min_ir_len>
"""
import sys, os, subprocess, tempfile
from Bio import SeqIO

if len(sys.argv) < 5:
    sys.stderr.write(
        "Usage: cp_ir_standardize.py <cp.fa> <hmms.hmm> <out.fa> <min_ir_len>\n"
    )
    sys.exit(1)

cpfa, hmmlib, outfa, min_ir = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

tmp = tempfile.mkdtemp(prefix="cpir_")
paf = os.path.join(tmp, "self.paf")

# minimap2 self alignment (asm5; report all)
subprocess.run(
    ["minimap2", "-x", "asm5", cpfa, cpfa],
    check=True,
    stdout=open(paf, "w"),
    stderr=subprocess.DEVNULL,
)

# parse PAF: look for longest inverted hit (qname==tname, strand == "-")
best = None
with open(paf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        qn, qlen, qs, qe, strand, tn, tlen, ts, te = (
            f[0],
            int(f[1]),
            int(f[2]),
            int(f[3]),
            f[4],
            f[5],
            int(f[6]),
            int(f[7]),
            int(f[8]),
        )
        if qn != tn or strand != "-":
            continue
        aln_len = abs(qe - qs)
        if aln_len >= min_ir:
            if best is None or aln_len > best[0]:
                best = (aln_len, qs, qe, ts, te, qlen)
if best is None:
    # No IR detected; copy as-is
    seq = next(SeqIO.parse(cpfa, "fasta"))
    SeqIO.write(seq, outfa, "fasta")
    print("IR_LEN\t0")
    print("KEEP_IR_START\tNA")
    print("KEEP_IR_END\tNA")
    # rbcL position below
else:
    aln_len, qs, qe, ts, te, L = best
    # define two IR copies as intervals [qs,qe) and [ts,te), choose keep_IR as the one with larger start (IRa by convention)
    s1, e1 = min(qs, qe), max(qs, qe)
    s2, e2 = min(ts, te), max(ts, te)
    # choose IRa = with larger midpoint
    mid1, mid2 = (s1 + e1) // 2, (s2 + e2) // 2
    if mid1 >= mid2:
        keep_s, keep_e = s1, e1
        drop_s, drop_e = s2, e2
    else:
        keep_s, keep_e = s2, e2
        drop_s, drop_e = s1, e1

    rec = next(SeqIO.parse(cpfa, "fasta"))
    seq = rec.seq

    # Define single-copy segments between IRs in circular space
    # segment A: (drop_e .. keep_s), segment B: (keep_e .. drop_s) in circular sense
    def slice_circ(s, e):
        if s <= e:
            return seq[s:e]
        else:
            return seq[s:] + seq[:e]

    # Normalize coordinates
    drop_s %= L
    drop_e %= L
    keep_s %= L
    keep_e %= L

    segA = slice_circ(drop_e, keep_s)  # one single-copy region
    segB = slice_circ(keep_e, drop_s)  # the other single-copy region
    IRa = slice_circ(keep_s, keep_e)  # kept IR

    # Identify LSC by rbcL hit location on original cp
    # run hmmsearch for rbcL only
    rbcL_hmm = os.path.join(tmp, "rbcL.hmm")
    # write a one-HMM file filtered from hmmlib if available; otherwise try full lib with --incE 1e-20 and pick rbcL rows
    # for portability, just run hmmsearch vs full lib and pick rbcL
    tbl = os.path.join(tmp, "cp_vs_hmms.tbl")
    subprocess.run(
        ["hmmsearch", "--tblout", tbl, "-E", "1e-5", hmmlib, cpfa],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    rbc_pos = None
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            target = f[0]  # seq id
            query = f[2]  # HMM name
            if query.lower().startswith("rbcl"):  # match rbcL
                s = int(f[17])
                e = int(f[18])
                rbc_pos = min(s, e) - 1  # to 0-based
                break

    # Determine which single-copy segment contains rbcL
    def in_seg(start, seg_s, seg_e):
        if seg_s <= seg_e:
            return (start >= seg_s) and (start < seg_e)
        else:
            return (start >= seg_s) or (start < seg_e)

    LSC = segA
    SSC = segB
    # compute their ranges in original coords
    # segA: drop_e..keep_s, segB: keep_e..drop_s
    if rbc_pos is not None:
        if in_seg(rbc_pos, drop_e, keep_s):
            LSC = segA
            SSC = segB
        elif in_seg(rbc_pos, keep_e, drop_s):
            LSC = segB
            SSC = segA
        else:
            # fallback: pick longer as LSC
            LSC = segA if len(segA) >= len(segB) else segB
            SSC = segB if LSC is segA else segA

    cp_std = LSC + IRa + SSC

    # Rotate so rbcL at position 1 if known
    if rbc_pos is not None:
        # compute rbcL pos after reordering; rebuild a mapping is heavy; approximate by searching motif from original around rbc_pos
        # For simplicity, search for the first 30bp of original sequence starting at rbc_pos in cp_std (case-insensitive)
        motif = str(seq[rbc_pos : rbc_pos + 30]).upper()
        pos = str(cp_std).upper().find(motif)
        if pos == -1:
            rot = cp_std
            rb_emit = "NA"
        else:
            rot = cp_std[pos:] + cp_std[:pos]
            rb_emit = str(1)
        cp_std = rot
    else:
        rb_emit = "NA"

    rec.seq = cp_std
    rec.id = rec.id + "|std"
    rec.description = "LSC|IRa|SSC rotated_to_rbcL"
    SeqIO.write(rec, outfa, "fasta")
    print("IR_LEN\t{}".format(aln_len))
    print("KEEP_IR_START\t{}".format(keep_s))
    print("KEEP_IR_END\t{}".format(keep_e))
    print("RBCL_POS\t{}".format(rb_emit))
