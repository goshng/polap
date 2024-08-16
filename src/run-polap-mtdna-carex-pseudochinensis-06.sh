#!/usr/bin/env bash

CONTIG_LENGTH=$1

minimap2 -cx map-ont o/60-mt-1/contig.fa o/n3k.fq.gz -t 56 -o o/60-mt-1/contig.paf >/dev/null 2>&1
cut -f1-11 o/60-mt-1/contig.paf | awk -v minlength=3000 '{if ($2>=minlength) {print}}' >o/60-mt-1/contig.tab

src/run-polap-pairs.R o/50-annotation/mt.contig.name-1 o/60-mt-1/contig.tab o/60-mt-1/o3000/seeds 3000 3000 >/dev/null 2>&1
cat o/60-mt-1/o3000/seeds/*.name o/60-mt-1/o3000/seeds/single.names | sort | uniq >o/60-mt-1/o3000/seeds/1.names
seqkit grep -f o/60-mt-1/o3000/seeds/1.names o/n3k.fq.gz -o o/60-mt-1/o3000/seeds/1.fq.gz >/dev/null 2>&1

TOTAL_LENGTH=$(seqkit stats -Ta o/60-mt-1/o3000/seeds/1.fq.gz | csvtk cut -t -f sum_len | csvtk del-header)
echo "DATA: total number of bases in the selected long-read data $TOTAL_LENGTH"
EXPECTED_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
echo "Compute: divide this number by the the number from step 05: $TOTAL_LENGTH / $CONTIG_LENGTH -> $EXPECTED_COVERAGE"
echo Your next command is:
echo $ src/run-polap-mtdna-carex-pseudochinensis-07.sh $EXPECTED_COVERAGE $CONTIG_LENGTH
