#!/bin/bash

fasta="$1"
# output="$2"
acc2taxid="nucl_gb.accession2taxid"

# if [[ ! -f "$fasta" || ! -f "$acc2taxid" ]]; then
#   echo "Usage: $0 <input.fasta> <output.tsv>"
#   echo "Make sure '$acc2taxid' (uncompressed) is in the current directory."
#   exit 1
# fi

# Step 1: Extract accession + length

seqkit fx2tab -ni -l "$fasta" >acc_len.tsv
cut -f1 acc_len.tsv >accessions.txt

# Step 2: Get TaxIDs
grep -F -f accessions.txt "$acc2taxid" >acc2taxid.tsv
cut -f1,3 acc2taxid.tsv >acc_taxid.tsv
cut -f3 acc2taxid.tsv | sort -u >taxids.txt

# Step 3: Get full lineage with taxon IDs and names
taxonkit lineage --show-lineage-taxids taxids.txt >lineage.tsv
