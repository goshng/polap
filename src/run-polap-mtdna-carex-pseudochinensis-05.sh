#!/usr/bin/env bash

rm -rf o/60-mt-1/o3000/seeds
mkdir -p o/60-mt-1/o3000/seeds
ln -s $PWD/o/n3k.fq.gz -t o/60-mt-1/o3000
seqkit grep -f o/50-annotation/mt.contig.name-1 o/30-contigger/graph_final.fasta -o o/60-mt-1/contig.fa
CONTIG_LENGTH=$(seqkit stats -Ta o/60-mt-1/contig.fa | csvtk cut -t -f sum_len | csvtk del-header)
echo "DATA: total number of bases in the contig data $CONTIG_LENGTH"

echo Step 6 is ready to run!
echo Your next command is:
echo $ src/run-polap-mtdna-carex-pseudochinensis-06.sh $CONTIG_LENGTH
