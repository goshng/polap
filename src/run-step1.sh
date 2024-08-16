#!/usr/bin/env bash

# usage:
#
# src/run-step1.sh <long-read in fastq> <short-read in fastq> <short-read in fastq> <number of CPU cores>
# src/run-step1.sh long.fq short1.fq short2.fq 4
# src/run-step1.sh cpn.fq 4C-pse_1.fastq 4C-pse_2.fastq 28

# required packages installation
#
# install miniconda and execute the following in the base conda environment.
# conda env create -f src/environment.yaml
# conda env create -f src/environment-fmlrc.yaml

# prints the usage:
#
if [ $# -eq 0 ]; then
	echo "$0 <long-read in fastq> <short-read in fastq> <short-read in fastq> <number of CPU cores>"
	echo "example: $0 long.fq short1.fq short2.fq 4"
	exit 0
fi

# command line input options
#
LR=$1
SR1=$2
SR2=$3
NT=$4

# required input file check
#
if test -z "$LR"; then
	echo "ERROR: No long read data set!"
	exit 1
else
	if ! test -s "$LR"; then
		echo "ERROR: long read data file $LR does not exist!"
	fi
fi

if test -z "$SR1"; then
	echo "ERROR: No short read data set!"
	exit 1
else
	if ! test -s "$SR1"; then
		echo "ERROR: short read data file $SR1 does not exist!"
	fi
fi

if test -z "$SR2"; then
	echo "ERROR: No short read data set!"
	exit 1
else
	if ! test -s "$SR2"; then
		echo "ERROR: short read data file $SR2 does not exist!"
	fi
fi

if test -z "$NT"; then
	NT=2
fi

# required packages check
#
command -v seqkit >/dev/null 2>&1 || {
	echo >&2 "seqkit: not installed"
	exit 1
}
command -v minimap2 >/dev/null 2>&1 || {
	echo >&2 "minimap2: not installed"
	exit 1
}
command -v flye >/dev/null 2>&1 || {
	echo >&2 "flye: not installed"
	exit 1
}
command -v makeblastdb >/dev/null 2>&1 || {
	echo >&2 "makeblastdb: not installed"
	exit 1
}
command -v tblastn >/dev/null 2>&1 || {
	echo >&2 "tblastn: not installed"
	exit 1
}
command -v bedtools >/dev/null 2>&1 || {
	echo >&2 "bedtools: not installed"
	exit 1
}
command -v prefetch >/dev/null 2>&1 || {
	echo >&2 "prefetch: not installed"
	exit 1
}
command -v jellyfish >/dev/null 2>&1 || {
	echo >&2 "jellyfish: not installed"
	exit 1
}
command -v csvtk >/dev/null 2>&1 || {
	echo >&2 "csvtk: not installed"
	exit 1
}

# variables
ODIR=o
WDIR=src
MR=3000
LR3K=n3k.fq.gz

# step0
mkdir -p o
echo 'INFO: o:' an output folder is created.

seqkit stats -Ta "$LR" | csvtk cut -t -f "sum_len" | csvtk del-header >"$ODIR"/long_total_length.txt
LONG_TOTAL_LENGTH=$(cat "$ODIR"/long_total_length.txt)
echo "DATA: long reads total (bases): $LONG_TOTAL_LENGTH"

jellyfish count -t "$NT" -C -m 19 -s 5G -o "$ODIR"/19mer_out --min-qual-char=? "$SR1" "$SR2"
jellyfish histo -o "$ODIR"/19mer_out.histo "$ODIR"/19mer_out
"$WDIR"/run-polap-jellyfish.R "$ODIR"/19mer_out.histo "$LONG_TOTAL_LENGTH" "$ODIR"/long_coverage.txt "$ODIR"/short_expected_genome_size.txt
EXPECTED_GENOME_SIZE=$(cat "$ODIR"/short_expected_genome_size.txt)
EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
EXPECTED_COVERAGE=$(cat "$ODIR"/long_coverage.txt)
EXPECTED_COVERAGE=${EXPECTED_COVERAGE%.*}
echo "DATA: short reads expected genome size (bases): $EXPECTED_GENOME_SIZE"
echo "DATA: long reads expected coverage: ${EXPECTED_COVERAGE}x"

# step1
seqkit seq --quiet -m "$MR" --threads 4 "$LR" -o $LR3K >/dev/null 2>&1
echo 'DATA: long-read minimum 3kb reads data n3k.fq.gz is created'

flye --nano-raw "$LR3K" --out-dir "$ODIR" \
	--threads "$NT" \
	--stop-after contigger \
	--asm-coverage 30 \
	--genome-size "$EXPECTED_GENOME_SIZE" \
	>/dev/null 2>&1
echo 'INFO: assembly graph: o/30-contigger/graph_final.gfa'

# step2
src/run-polap-annotation o "$NT" --selective-annotate --contigger
echo "INFO: contig annotation with mitochondrial and plastid genes"

# step3
src/run-polap-mtcontig.R o o/50-annotation/mt.contig.name o/assembly_info_organelle_annotation_count.txt --contigger
echo 'INFO: o/50-annotation/mt.contig.name-1 for mtDNA contig candidates'
echo 'INFO: annotation: column -t o/assembly_info_organelle_annotation_count.txt'

# step4
STEP4=1
ADIR=o/50-annotation
MTDIR=o/60-mt-${STEP4}
MTODIR=o/60-mt-${STEP4}/o3000
MTSEEDSDIR=o/60-mt-${STEP4}/o3000/seeds
ASSEMBLY_FASTA=o/30-contigger/contigs.fasta
SINGLE_MIN=3000
LR3K=n3k.fq.gz
echo "INFO: mt.contig.name-${STEP4}: contigs only"
echo 'INFO: use contigs: o/30-contigger/contigs.fasta'

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
echo "column -t $ODIR/assembly_info_organelle_annotation_count.txt"
