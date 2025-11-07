# FILE: scripts/annotate_cp_genes_hmmer.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Annotate plastid CDS by hmmsearch vs HMM library; extract best non-overlapping hit per gene.
Outputs:
  <outdir>/<gene>.prot.faa (AA)
  <outdir>/<gene>.cds.fna (NT)
"""
import sys, os, subprocess, tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if len(sys.argv) < 4:
    sys.stderr.write("Usage: annotate_cp_genes_hmmer.py <cp.fa> <hmms.hmm> <outdir>\n")
    sys.exit(1)
cp, hmmlib, outd = sys.argv[1:4]
os.makedirs(outd, exist_ok=True)

# run hmmsearch
tbl = os.path.join(outd, "cp_vs_hmms.tbl")
subprocess.run(
    ["hmmsearch", "--tblout", tbl, "-E", "1e-5", hmmlib, cp],
    check=True,
    stdout=subprocess.DEVNULL,
)

# parse sequences
seqs = {r.id: str(r.seq) for r in SeqIO.parse(cp, "fasta")}

# parse tbl, choose best per gene
import collections

hits = collections.defaultdict(list)  # gene -> [(seqid, start, end, evalue)]
with open(tbl) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        target = f[0]  # sequence id (cp contig)
        gene = f[2]  # HMM name ~ gene
        s = int(f[17])
        e = int(f[18])
        evalue = float(f[4])
        hits[gene].append((target, min(s, e), max(s, e), evalue))

for gene, arr in hits.items():
    # best by evalue
    arr.sort(key=lambda x: x[3])
    t, s, e, _ = arr[0]
    nt = seqs[t][s - 1 : e]  # hmmsearch coords are 1-based inclusive
    # translate naive (frame unknown) -> try 3 frames, choose longest ORF
    best = None
    for frame in range(3):
        aa = Seq(nt[frame:]).translate(to_stop=False)
        if best is None or len(aa) > len(best):
            best = aa
    # write
    with open(os.path.join(outd, f"{gene}.prot.faa"), "a") as o:
        o.write(f">{os.path.basename(os.path.dirname(outd))}\n{str(best)}\n")
    with open(os.path.join(outd, f"{gene}.cds.fna"), "a") as o:
        o.write(f">{os.path.basename(os.path.dirname(outd))}\n{nt}\n")
