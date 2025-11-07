# FILE: scripts/align_and_backtranslate.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
For each gene subdir under genes/, pool *.prot.faa across species, align with MAFFT L-INS-i,
and back-translate using per-species cds.fna. Outputs into aln/:
  <gene>.aa.aln.faa, <gene>.codon.aln.fna
"""
import sys, os, subprocess
from Bio import SeqIO

if len(sys.argv) < 4:
    sys.stderr.write(
        "Usage: align_and_backtranslate.py <genes_dir> <aln_dir> <threads>\n"
    )
    sys.exit(1)
gdir, outd, threads = sys.argv[1:4]
threads = int(threads)
os.makedirs(outd, exist_ok=True)

genes = sorted(
    {
        os.path.basename(f).replace(".prot.faa", "")
        for f in sum(
            [
                [os.path.join(dp, f) for f in files if f.endswith(".prot.faa")]
                for dp, _, files in os.walk(gdir)
            ],
            [],
        )
    }
)

for gene in genes:
    # gather proteins
    pooled_prot = os.path.join(outd, f"{gene}.prot.faa")
    with open(pooled_prot, "w") as o:
        for sp in sorted(os.listdir(gdir)):
            prot = os.path.join(gdir, sp, f"{gene}.prot.faa")
            if os.path.exists(prot):
                for rec in SeqIO.parse(prot, "fasta"):
                    rec.id = sp
                    rec.description = ""
                    o.write(f">{sp}\n{str(rec.seq)}\n")
    aln_aa = os.path.join(outd, f"{gene}.aa.aln.faa")
    subprocess.run(
        [
            "mafft",
            "--localpair",
            "--maxiterate",
            "1000",
            "--thread",
            str(threads),
            pooled_prot,
        ],
        check=True,
        stdout=open(aln_aa, "w"),
    )

    # back-translate AA aln using species-specific cds (pal2nal-like)
    cds_map = {}
    for sp in sorted(os.listdir(gdir)):
        cds = os.path.join(gdir, sp, f"{gene}.cds.fna")
        if os.path.exists(cds):
            for rec in SeqIO.parse(cds, "fasta"):
                cds_map[sp] = str(rec.seq)
    # read aligned aa
    aa_aln = {rec.id: str(rec.seq) for rec in SeqIO.parse(aln_aa, "fasta")}
    # build codon aln
    out_nt = os.path.join(outd, f"{gene}.codon.aln.fna")
    with open(out_nt, "w") as o:
        for sp, aa in aa_aln.items():
            if sp not in cds_map:
                continue
            nt = cds_map[sp]
            # naive back-translation by consuming 3 nts per non-gap AA; gap-> '---'
            i = 0
            codons = []
            for a in aa:
                if a == "-":
                    codons.append("---")
                else:
                    codons.append(nt[i : i + 3] if i + 3 <= len(nt) else "---")
                    i += 3
            o.write(f">{sp}\n{''.join(codons)}\n")
