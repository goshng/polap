#!/usr/bin/env bash

# usage:
#
# src/run-step4.sh <unpolished sequence fasta file> <output fasta file to polish> <short-read in fastq> <short-read in fastq>
# src/run-step4.sh mt.0.fasta mt.1.fa short1.fq short2.fq

# required packages installation
#
# conda env create -f src/environment.yaml
# conda env create -f src/environment-fmlrc.yaml

if [ $# -eq 0 ]; then
	echo "src/run-step4.sh <unpolished sequence fasta file> <output fasta file to polish> <short-read in fastq> <short-read in fastq> <number of CPU cores>"
	echo src/run-step4.sh mt.0.fasta mt.1.fa short1.fq short2.fq 4
	exit 0
fi

# command line input options
#
PA=$1
FA=$2
SR1=$3
SR2=$4
NT=$5
ODIR=o

if test -z "$PA"; then
	echo "ERROR: No fasta file set!"
	exit 1
else
	if ! test -s "$PA"; then
		echo "ERROR: fasta file $PA does not exist!"
	fi
fi

if test -z "$FA"; then
	echo "ERROR: No fasta file set!"
	exit 1
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

# STEP 6: polishing
#
function step6() {
	command -v msbwt >/dev/null 2>&1 || {
		echo >&2 "msbwt: not installed"
		echo "ERROR: execute: conda env create -f src/environment-fmlrc.yaml"
		echo "conda activate polap-fmlrc"
		exit 1
	}
	command -v ropebwt2 >/dev/null 2>&1 || {
		echo >&2 "ropebwt2: not installed"
		echo "ERROR: execute: conda env create -f src/environment-fmlrc.yaml"
		echo "conda activate polap-fmlrc"
		exit 1
	}

	if test -z "$SR1"; then
		echo "ERROR: No short read data set!"
		exit 1
	fi

	if [[ $SR1 = *.fastq || $SR1 = *.fq ]]; then
		cat $SR1 $SR2 |
			awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN |
			msbwt convert $ODIR/msbwt
	elif [[ $SR1 = *.fq.gz || $SR1 = *.fastq.gz ]]; then
		zcat $SR1 $SR2 |
			awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN |
			msbwt convert $ODIR/msbwt
	fi
}

function step7() {
	if test -z "$PA"; then
		echo ERROR: No assembled organelle sequence in FASTA format: ["$PA"]
		exit 1
	fi

	if test -z "$FA"; then
		echo ERROR: No assembled organelle sequence in FASTA format: ["$FA"]
		exit 1
	fi

	if [[ -s $PA ]]; then
		fmlrc -p $NT $ODIR/msbwt/comp_msbwt.npy $PA $FA
	else
		echo "ERROR: no such file: $PA"
		exit 1
	fi
}

step6
step7
