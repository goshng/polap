#!/usr/bin/env bash

# usage:
#
# src/run-step2.sh <step2 iteration number> <number of CPU cores>
# src/run-step2.sh 2 28

# prints the usage:
#
if [ $# -eq 0 ]; then
	echo "src/run-step2.sh <step2 iteration number> <number of CPU cores>"
	echo "e.g., src/run-step2.sh 2 4"
	exit 0
fi

# command line input options
#
STEP4=$1
NT=$2

if test -z "$NT"; then
	NT=2
fi

MR=3000

WDIR=src
ADIR=o/50-annotation
MTDIR=o/60-mt-${STEP4}
MTODIR=o/60-mt-${STEP4}/o3000
MTSEEDSDIR=o/60-mt-${STEP4}/o3000/seeds
ASSEMBLY_FASTA=o/30-contigger/graph_final.fasta
SINGLE_MIN=3000
LR3K=n3k.fq.gz
echo "INFO: mt.contig.name-${STEP4}: edge only"
echo 'INFO: use edges not contigs: o/30-contigger/graph_final.fasta'

rm -rf $MTODIR
mkdir -p $MTSEEDSDIR

seqkit grep -f $ADIR/mt.contig.name-$STEP4 $ASSEMBLY_FASTA -o $MTDIR/contig.fa >/dev/null 2>&1
echo "DATA: $MTDIR/contig.fa is created."

minimap2 -cx map-ont "$MTDIR"/contig.fa "$LR3K" -t "$NT" -o "$MTDIR"/contig.paf >/dev/null 2>&1
echo "DATA: $LR3K is used to select reads."
echo "DATA: $MTDIR/contig.paf is created."

cut -f1-11 "$MTDIR"/contig.paf | awk -v minlength=$MR '{if ($2>=minlength) {print}}' >"$MTDIR"/contig.tab
echo 'DATA: mininmum length of long reads in the selection mapping: 3000'
echo "DATA: $MTDIR/contig.tab is created."

CONTIG_LENGTH=$(seqkit stats -Ta $MTDIR/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
echo "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"

"$WDIR"/run-polap-pairs.R $ADIR/mt.contig.name-$STEP4 $MTDIR/contig.tab $MTSEEDSDIR $SINGLE_MIN $SINGLE_MIN >/dev/null 2>&1
echo "DATA: pair contig names in $MTSEEDSDIR are created."
echo "DATA: single contig name in $MTSEEDSDIR is created."

cat "$MTSEEDSDIR/"*".name" $MTSEEDSDIR/single.names | sort | uniq >$MTSEEDSDIR/1.names
echo "INFO: creates long read names and the single name in $MTSEEDSDIR/1.names"

seqkit grep -f $MTSEEDSDIR/1.names $LR3K -o $MTSEEDSDIR/1.fq.gz >/dev/null 2>&1
echo "DATA: organelle reads in $MTSEEDSDIR/1.fq.gz"

TOTAL_LENGTH=$(seqkit stats -Ta $MTSEEDSDIR/1.fq.gz | csvtk cut -t -f "sum_len" | csvtk del-header)
EXPECTED_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
echo "INFO: expected coverage: ${EXPECTED_COVERAGE}x"

echo 'DATA: no data reduction: COV=100'
ln -s 1.fq.gz $MTSEEDSDIR/2.fq.gz

flye --nano-raw $MTSEEDSDIR/2.fq.gz --out-dir $MTODIR --asm-coverage 30 --genome-size $CONTIG_LENGTH --threads $NT >/dev/null 2>&1
echo INFO: assembly graph: $MTODIR/assembly_graph.gfa
