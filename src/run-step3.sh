#!/usr/bin/env bash

# usage:
#
# src/run-step3.sh <unpolished genome sequence from long read assembly>
# src/run-step3.sh mt.0.fasta

# prints the usage:
#
if [ $# -eq 0 ]; then
	echo "src/run-step3.sh <unpolished genome sequence in fasta format>"
	echo "use Bandage to extract one from long read assembly"
	echo "e.g., src/run-step3.sh mt.0.fasta"
	exit 0
fi

# command line input options
#
FA=$1
NT=2

if test -z "$FA"; then
	echo "ERROR: No fasta file set!"
	exit 1
else
	if ! test -s "$FA"; then
		echo "ERROR: fasta file $FA does not exist!"
	fi
fi

# STEP 5: organelle genome candidate validation
# samtools coverage
#
# bio conda environment
command -v minimap2 >/dev/null 2>&1 || {
	echo >&2 "minimap2: not installed"
	exit 1
}
command -v samtools >/dev/null 2>&1 || {
	echo >&2 "samtools: not installed"
	exit 1
}

if test -z "$FA"; then
	echo ERROR: No assembled organelle sequence in FASTA format: ["$FA"]
	exit 1
fi

if test -z "$LR3K"; then
	LR3K=n3k.fq.gz
	echo "DATA: $LR3K is used for coverage analysis"
	if [[ ! -s $LR3K ]]; then
		echo ERROR: No 3K long-read data!
		exit 1
	fi
fi

# minimap2 -ax map-ont $FA $LR3K -o 1.sam
# samtools view -u .1.to.cpsen3k.sam | samtools sort -o cpseudochinensis.mt1.1.to.cpsen3k.bam
# samtools coverage -A -w 32 cpseudochinensis.mt1.1.to.cpsen3k.bam
echo INFO: long-reads assembled organelle genome in FASTA format: $FA
echo INFO: long-reads: $LR3K
echo INFO: number of cores: $NT

minimap2 -t $NT -ax map-ont $FA $LR3K | samtools view -u | samtools sort -o 1.bam
samtools coverage -A -w 32 1.bam

echo INFO: try to execute: samtools coverage -A -w 32 1.bam
echo Change your conda environment before performing STEP 6.
