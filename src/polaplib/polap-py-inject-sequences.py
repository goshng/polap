#!/usr/bin/env python3
from Bio import SeqIO
import sys

masked_fasta = sys.argv[1]
inserts_fasta = sys.argv[2]
output_fasta = sys.argv[3]

inserts = list(SeqIO.parse(inserts_fasta, "fasta"))
insert_index = 0


def inject(seq):
    global insert_index
    new_seq = []
    i = 0
    while i < len(seq):
        if seq[i] == "N":
            frag = str(inserts[insert_index].seq)
            new_seq.append(frag)
            i += len(frag)
            insert_index += 1
        else:
            new_seq.append(seq[i])
            i += 1
    return "".join(new_seq)


with open(output_fasta, "w") as out:
    for record in SeqIO.parse(masked_fasta, "fasta"):
        injected_seq = inject(str(record.seq))
        record.seq = injected_seq
        SeqIO.write(record, out, "fasta")
