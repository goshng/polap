#!/usr/bin/env bash

################################################################################
# polap: Plant organelle long-read assembly pipeline
# Copyright (C) 2024 Sungshin Women's University
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
################################################################################

# setup
# 1. install conda environments
# 2. edit data file names

# install miniconda and create conda environments
#
# conda environment: polap-dev
# (base) $ conda env create -f src/environment.yaml
#
# conda environment: polap-fmlrc
# (base) $ conda env create -f src/environment-fmlrc.yaml

# input/output data file names
#
# LR=cpn.fq          # long-read data file
# SR1=4C-pse_1.fastq # paired short-read data file 1
# SR2=4C-pse_2.fastq # paired short-read data file 2
# PA=cp.mt.0.fasta   # assembled draft sequence extracted from bandage
# FA=cp.mt.1.fa      # polished sequence
LR=long.fq    # long-read data file
SR1=short1.fq # paired short-read data file 1
SR2=short2.fq # paired short-read data file 2
PA=mt.0.fasta # assembled draft sequence extracted from bandage
FA=mt.1.fa    # polished sequence

# variables
#
ODIR=o
IDIR=0
JDIR=1
WDIR="${BASH_SOURCE%/*}"
if [[ ! -d "$WDIR" ]]; then WDIR="$PWD"; fi
LR3K=$ODIR/n3k.fq.gz
MR=3000
MPAIR=3000
MBRIDGE=3000
MSING=3000
COV=50
NT=2
if test -z "$DEBUG"; then
	DEBUG=0
fi
STEP4=1
CIRCULARIZE="" # "--circularize"

# Constants
EXIT_SUCCESS=1
EXIT_FAIL=0
EXIT_ERROR=2

# functions
#
function run_check1() {
	command -v bc >/dev/null 2>&1 || {
		echo >&2 "bc: not installed"
		return $EXIT_FAIL
	}
	command -v seqkit >/dev/null 2>&1 || {
		echo >&2 "seqkit: not installed"
		return $EXIT_FAIL
	}
	command -v minimap2 >/dev/null 2>&1 || {
		echo >&2 "minimap2: not installed"
		return $EXIT_FAIL
	}
	command -v flye >/dev/null 2>&1 || {
		echo >&2 "flye: not installed"
		return $EXIT_FAIL
	}
	command -v makeblastdb >/dev/null 2>&1 || {
		echo >&2 "makeblastdb: not installed"
		return $EXIT_FAIL
	}
	command -v tblastn >/dev/null 2>&1 || {
		echo >&2 "tblastn: not installed"
		return $EXIT_FAIL
	}
	command -v bedtools >/dev/null 2>&1 || {
		echo >&2 "bedtools: not installed"
		return $EXIT_FAIL
	}
	command -v prefetch >/dev/null 2>&1 || {
		echo >&2 "prefetch: not installed"
		return $EXIT_FAIL
	}
	command -v jellyfish >/dev/null 2>&1 || {
		echo >&2 "jellyfish: not installed"
		return $EXIT_FAIL
	}
	command -v csvtk >/dev/null 2>&1 || {
		echo >&2 "csvtk: not installed"
		return $EXIT_FAIL
	}
	return $EXIT_SUCCESS
}

function run_check2() {
	command -v msbwt >/dev/null 2>&1 || {
		echo >&2 "msbwt: not installed"
		return $EXIT_FAIL
	}
	command -v ropebwt2 >/dev/null 2>&1 || {
		echo >&2 "ropebwt2: not installed"
		return $EXIT_FAIL
	}
	command -v fmlrc >/dev/null 2>&1 || {
		echo >&2 "ropebwt2: not installed"
		return $EXIT_FAIL
	}
	return $EXIT_SUCCESS
}

function run_output1() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi
	mkdir -p "$ODIR"
	echo INFO: output folder $ODIR is created.
	echo INFO: number of CPUs is $NT
	echo INFO: long read data file: $LR
	echo INFO: short read data file 1: $SR1
	echo INFO: short read data file 2: $SR2
}

function run_jellyfish() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	seqkit stats -Ta "$LR" | csvtk cut -t -f "sum_len" | csvtk del-header >"$ODIR"/long_total_length.txt
	LONG_TOTAL_LENGTH=$(cat "$ODIR"/long_total_length.txt)
	echo "DATA: long reads total (bases): $LONG_TOTAL_LENGTH"

	# See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
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

	echo "INFO: executing the whole-genome assembly using flye ... be patient!"
	flye --nano-raw "$LR3K" --out-dir "$ODIR" \
		--threads "$NT" \
		--stop-after contigger \
		--asm-coverage 30 \
		--genome-size "$EXPECTED_GENOME_SIZE" \
		>/dev/null 2>&1
	echo "INFO: assembly graph: "$ODIR"/30-contigger/graph_final.gfa"

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

run_annotation() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi
	# src/run-polap-annotation o "$NT" --selective-annotate --contigger
	echo "INFO: contig annotation with mitochondrial and plastid genes"

	MTAA=$WDIR/mt.1.c70.3.faa
	PTAA=$WDIR/pt.2.c70.3.faa
	FLYEDIR=$ODIR
	NUMTHREADS=$NT

	ALL_ANNOTATE=1
	FLYE_CONTIGGER=1

	ASSEMBLYINFO="$ODIR"/30-contigger/contigs_stats.txt
	ASSEMBLYFILE="$ODIR"/30-contigger/contigs.fasta
	ADIR="$ODIR"/50-annotation
	CONTIGFILE="$ODIR"/50-annotation/contig.fasta
	CONTIGDB="$ODIR"/50-annotation/contig
	MTAABLAST="$ODIR"/50-annotation/mtaa.blast
	MTAABED="$ODIR"/50-annotation/mtaa.bed
	MTGENECOUNT="$ODIR"/50-annotation/mt.gene.count
	PTAABLAST="$ODIR"/50-annotation/ptaa.blast
	PTAABED="$ODIR"/50-annotation/ptaa.bed
	PTGENECOUNT="$ODIR"/50-annotation/pt.gene.count

	mkdir -p "$ODIR"/50-annotation
	echo "INFO: $ODIR/50-annotation is created"

	CONTIGNAME="$ODIR"/50-annotation/contig.name

	#src/run-polap-select.R o/30-contigger/contigs_stats.txt o/50-annotation/contig.name
	grep -v "#" $ASSEMBLYINFO | cut -f 1 >$CONTIGNAME
	echo "INFO: contig sequence names in file: $CONTIGNAME"

	seqkit grep -f "$ODIR"/50-annotation/contig.name \
		"$ODIR"/30-contigger/contigs.fasta \
		-o "$ODIR"/50-annotation/contig.fasta \
		>/dev/null 2>&1
	echo "INFO: contig sequence file: $ODIR/50-annotation/contig.fasta"

	makeblastdb -dbtype nucl \
		-in "$ODIR"/50-annotation/contig.fasta \
		-out "$ODIR"/50-annotation/contig \
		>/dev/null 2>&1
	echo "INFO: BLASTDB of the contig sequences: $ODIR/50-annotation/contig"

	# mtDNA gene annotation and counts
	echo "INFO: BLAST of the mitochondrial proteins againt $ODIR/50-annotation/contig"
	echo "INFO: executing the tblastn ... be patient!"
	tblastn -query "$WDIR"/mt.1.c70.3.faa \
		-db "$ODIR"/50-annotation/contig \
		-out "$ODIR"/50-annotation/mtaa.blast \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads $NT \
		>/dev/null 2>&1
	"$WDIR"/run-polap-genes.R "$ODIR"/50-annotation/mtaa.blast \
		"$ODIR"/50-annotation/mtaa.blast.bed \
		>/dev/null 2>&1
	sort -k1,1 -k2,2n "$ODIR"/50-annotation/mtaa.blast.bed >"$MTAABLAST".sorted.bed
	mkdir "$ODIR"/50-annotation/mtaa.bed

	# + IFS=
	# + read -r contig
	# + grep -w contig_2905 o/50-annotation/ptaa.blast.sorted.bed
	# + bedtools merge -i o/50-annotation/ptaa.bed/contig_2905.bed
	# ++ wc -l
	# + printf '%s\t%d\n' contig_2905 0
	# + IFS=
	# + read -r contig

	echo "INFO: counting mitochondrial genes in the contigs"
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w $contig $MTAABLAST.sorted.bed >"$MTAABED"/$contig.bed
		bedtools merge -i $MTAABED/$contig.bed >"$MTAABED"/$contig.bed.txt
		printf "%s\t%d\n" $contig $(wc -l <"$MTAABED"/$contig.bed.txt)
	done <"$CONTIGNAME" | sort -k2 -rn >$MTGENECOUNT
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	echo "INFO: compressing the BLAST results of mitochondrial gene annotation"
	tar zcf "$ODIR"/50-annotation/mtaa.bed.tar.gz "$ODIR"/50-annotation/mtaa.bed
	rm -rf "$ODIR"/50-annotation/mtaa.bed

	# Plastid gene annotation and counts
	echo "INFO: BLAST of the plastid proteins againt $ODIR/50-annotation/contig"
	echo "INFO: executing the tblastn ... be patient!"
	tblastn -query "$WDIR"/pt.2.c70.3.faa \
		-db "$ODIR"/50-annotation/contig \
		-out "$ODIR"/50-annotation/ptaa.blast \
		-evalue 1e-30 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads $NT \
		>/dev/null 2>&1
	"$WDIR"/run-polap-genes.R "$ODIR"/50-annotation/ptaa.blast \
		"$ODIR"/50-annotation/ptaa.blast.bed \
		>/dev/null 2>&1
	sort -k1,1 -k2,2n "$ODIR"/50-annotation/ptaa.blast.bed >"$PTAABLAST".sorted.bed
	mkdir "$ODIR"/50-annotation/ptaa.bed

	echo "INFO: counting plastid genes in the contigs"
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w $contig $PTAABLAST.sorted.bed >"$PTAABED"/$contig.bed
		bedtools merge -i $PTAABED/$contig.bed >"$PTAABED"/$contig.bed.txt
		printf "%s\t%d\n" $contig $(wc -l <"$PTAABED"/$contig.bed.txt)
	done <"$CONTIGNAME" | sort -k2 -rn >$PTGENECOUNT
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	echo "INFO: compressing the BLAST results of plastid gene annotation"
	tar zcf "$ODIR"/50-annotation/ptaa.bed.tar.gz "$ODIR"/50-annotation/ptaa.bed
	rm -rf "$ODIR"/50-annotation/ptaa.bed

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

run_count_genes() {
	# step3
	"$WDIR"/run-polap-mtcontig.R "$ODIR" \
		"$ODIR"/50-annotation/mt.contig.name \
		"$ODIR"/assembly_info_organelle_annotation_count.txt \
		--contigger \
		>/dev/null 2>&1

	echo "USE: assembly graph: "$PWD"/"$ODIR"/30-contigger/graph_final.gfa"
	echo "USE: execute $ column -t "$ODIR"/assembly_info_organelle_annotation_count.txt"
	echo "INFO: edit "$ODIR"/50-annotation/mt.contig.name-1 for mtDNA contig candidates"
}

run_flye2() {
	# step4
	if [ "$DEBUG" -eq 1 ]; then set -x; fi
	STEP4=$1
	MR=$2
	ADIR="$ODIR"/50-annotation
	MTDIR="$ODIR"/60-mt-${STEP4}
	MTODIR="$ODIR"/60-mt-${STEP4}/o${MR}
	MTSEEDSDIR="$ODIR"/60-mt-${STEP4}/o${MR}/seeds
	# for contigs
	#	ASSEMBLY_FASTA=o/30-contigger/contigs.fasta
	#	for edges
	ASSEMBLY_FASTA="$ODIR"/30-contigger/graph_final.fasta
	PAIR_MIN=$MPAIR
	SINGLE_MIN=$MSING
	#	echo "INFO: mt.contig.name-${STEP4}: contigs only"
	#	echo 'INFO: use contigs: o/30-contigger/contigs.fasta'

	if [[ -d $MTODIR ]]; then
		echo "INFO: $MTODIR is deleted."
		rm -rf $MTODIR
	fi

	mkdir -p $MTSEEDSDIR
	ln -s "$PWD"/"$ODIR"/n3k.fq.gz -t $MTODIR
	echo "INFO: uses mt.contig.name at $ADIR/mt.contig.name-$STEP4"
	ADIR="$ODIR"/50-annotation
	if [[ ! -d $ADIR ]]; then
		mkdir -p $ADIR
	fi

	# seqkit grep -f $ADIR/mt.contig.name-$STEP4 $ASSEMBLY_FASTA -o $MTDIR/contig.fa
	seqkit grep -f $ADIR/mt.contig.name-$STEP4 $ASSEMBLY_FASTA -o $MTDIR/contig.fa >/dev/null 2>&1
	if [[ $CIRCULARIZE == "--circularize" ]]; then
		contig_count=$(wc -l <$ADIR/mt.contig.name-$STEP4)
		if [ "$contig_count" -eq 1 ]; then
			seqkit fx2tab --length --name $MTDIR/contig.fa -o $MTDIR/contig.fa.len >/dev/null 2>&1
			A=$(cut -f2 $MTDIR/contig.fa.len)
			B=$(echo "scale=0; $A/2" | bc)
			C=$((B + 1))
			seqkit subseq -r 1:$B $MTDIR/contig.fa -o $MTDIR/c1.fa >/dev/null 2>&1
			seqkit subseq -r $C:$A $MTDIR/contig.fa -o $MTDIR/c2.fa >/dev/null 2>&1
			cat $MTDIR/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o $MTDIR/contig.fa >/dev/null 2>&1
			cp $ADIR/mt.contig.name-$STEP4 $ADIR/mt.contig.name-${STEP4}-backup
			echo -e "edge_1\nedge_2" >$ADIR/mt.contig.name-$STEP4
			echo "INFO: creates new $MTDIR/contig.fa and $ADIR/mt.contig.name-$STEP4"
		else
			echo "DEV: not implemented yet"
			exit EXIT_ERROR
			"$WDIR"/run-polap-single.R $MTSEEDSDIR/contig.tab $MTSEEDSDIR $SINGLE_MIN >/dev/null 2>&1
			cat $MTSEEDSDIR/single.names | sort | uniq >$MTSEEDSDIR/1.names
			echo "INFO: creates long read single name in $MTSEEDSDIR/1.names"
		fi
	fi
	echo "DATA: $MTDIR/contig.fa is created."

	CONTIG_LENGTH=$(seqkit stats -Ta $MTDIR/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	echo $CONTIG_LENGTH >"$MTDIR"/contig_total_length.txt
	echo "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"

	if [[ -s "$MTDIR"/contig.paf ]]; then
		echo "DATA: previously created $MTDIR/contig.paf is used without executing minimap2."
	else
		minimap2 -cx map-ont "$MTDIR"/contig.fa "$LR3K" -t "$NT" -o "$MTDIR"/contig.paf >/dev/null 2>&1
		echo "DATA: $LR3K is used to select reads."
		echo "DATA: $MTDIR/contig.paf is created."
	fi

	cut -f1-11 "$MTDIR"/contig.paf | awk -v minlength=$MR '{if ($2>=minlength) {print}}' >"$MTODIR"/contig.tab
	echo "DATA: mininmum length of long reads in the selection mapping: $MR"
	echo "DATA: $MTODIR/contig.tab is created."

	"$WDIR"/run-polap-pairs.R $ADIR/mt.contig.name-$STEP4 $MTODIR/contig.tab $MTSEEDSDIR $PAIR_MIN $MBRIDGE >/dev/null 2>&1
	#	"$WDIR"/run-polap-pairs.R $ADIR/mt.contig.name-$STEP4 $MTODIR/contig.tab $MTSEEDSDIR $PAIR_MIN >/dev/null 2>&1
	echo "OPTION polap pairs alignment minimum: $PAIR_MIN"
	echo "OPTION polap pairs bridge minimum: $MBRIDGE"
	echo "DATA: pair contig names in $MTSEEDSDIR are created."
	echo "DATA: single contig name in $MTSEEDSDIR is created."

	cat "$MTSEEDSDIR/"*".name" $MTSEEDSDIR/single.names | sort | uniq >$MTSEEDSDIR/1.names
	echo "INFO: creates long read names and the single name in $MTSEEDSDIR/1.names"

	seqkit grep -f $MTSEEDSDIR/1.names $LR3K -o $MTSEEDSDIR/1.fq.gz >/dev/null 2>&1
	echo "DATA: organelle reads in $MTSEEDSDIR/1.fq.gz"

	TOTAL_LENGTH=$(seqkit stats -Ta $MTSEEDSDIR/1.fq.gz | csvtk cut -t -f "sum_len" | csvtk del-header)
	EXPECTED_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
	echo "INFO: expected coverage: ${EXPECTED_COVERAGE}x"

	if [ "$EXPECTED_COVERAGE" -lt $COV ]; then
		echo "DATA: no data reduction: COV=$COV"
		ln -s 1.fq.gz $MTSEEDSDIR/2.fq.gz
	else
		echo "SUGGESTION: you might want to increase the minimum read lengths"
		echo "SUGGESTION: since you have enough long-read data."
		RATE=$(echo "scale=10; $COV/$EXPECTED_COVERAGE" | bc)
		echo "DATA: data reduction by rate of $RATE"
		seqkit sample -p $RATE $MTSEEDSDIR/1.fq.gz -o $MTSEEDSDIR/2.fq.gz >/dev/null 2>&1
	fi

	# put the backup to the original
	if [[ $CIRCULARIZE == "--circularize" ]]; then
		if [[ -s $ADIR/mt.contig.name-${STEP4}-backup ]]; then
			mv $ADIR/mt.contig.name-${STEP4}-backup $ADIR/mt.contig.name-$STEP4
		else
			echo "DEV: not implemented yet"
			exit EXIT_ERROR
		fi
	fi

	echo "INFO: executing the organelle-genome assembly using flye ... be patient!"
	flye --nano-raw $MTSEEDSDIR/2.fq.gz \
		--out-dir $MTODIR \
		--asm-coverage 30 \
		--genome-size $CONTIG_LENGTH \
		--threads $NT \
		>/dev/null 2>&1
	echo "INFO: assembly graph in the flye contigger stage: "$PWD/$MTODIR"/30-contigger/graph_final.gfa"
	echo "RESULT: assembly graph: $PWD/$MTODIR/assembly_graph.gfa"
	echo "INFO: extract a draft organelle genome sequence from the assembly graph"
	# echo "column -t $ODIR/assembly_info_organelle_annotation_count.txt"

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function run_samtools() {
	echo "INFO: executing minimap2 and samtools for checking the long reads coverage on the $PA ... be patient!"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi
	minimap2 -t $NT -ax map-ont $PA $LR3K 2>/dev/null |
		samtools view -u 2>/dev/null |
		samtools sort -o "$ODIR"/1.bam \
			>/dev/null 2>&1
	samtools coverage -A -w 32 "$ODIR"/1.bam
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

run_polish1() {
	echo "INFO: excuting ropebwt2 and msbwt on the short reads ... be patient!"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi
	if [[ $SR1 = *.fastq || $SR1 = *.fq ]]; then
		cat $SR1 $SR2 |
			awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN |
			msbwt convert $ODIR/msbwt \
				>/dev/null 2>&1
	elif [[ $SR1 = *.fq.gz || $SR1 = *.fastq.gz ]]; then
		zcat $SR1 $SR2 |
			awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN |
			msbwt convert $ODIR/msbwt \
				>/dev/null 2>&1
	fi
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

run_polish2() {
	echo "INFO: executing fmlrc on the draft sequence $PA ... be patient!"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi
	if [[ -s $PA ]]; then
		fmlrc -p $NT $ODIR/msbwt/comp_msbwt.npy $PA $FA >/dev/null 2>&1
	else
		echo "ERROR: no such file: $PA"
		exit 1
	fi
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function usage() {
	cat >&2 <<EOM
$(basename "$0") [options]

Decription:
'P'lant 'o'rganelle DNA 'l'ong-read 'a'ssembly 'p'ipeline.

steps subcommand:
    -l <fastq>  long-reads data file in fastq format [$LR]
    -o <dir>    output folder name [$ODIR]
    -a <fastq>  short-read fastq file 1 [$SR1]
    -b <fastq>  short-read fastq file 2 [$SR2]
    -i <number> $ODIR/50-annotation/mt.contig.name-<number> [$STEP4]
				flye source output [$IDIR]
    -j <number> flye destination output [$JDIR]
    -m <number> minimum length of long reads [$MR]
    -a <number> minimum length of alignment in mapping [$MPAIR]
    -b <number> minimum length of bridging reads in mapping [$MBRIDGE]
    -c <number> coverage for the 2nd assembly [$COV]
	-u			circularize a contig [FALSE]
    -p <fasta>  polishing sequence in fasta format [$PA]
    -f <fasta>  final assembly in fasta format [$FA]
    -t <number> number of CPUs [$NT]
Usage:
  $WDIR/run-polap1.sh test
  $WDIR/run-polap1.sh assemble1 -l long.fq -a short1.fq -b short2.fq
  $WDIR/run-polap1.sh assemble2 -i 1 -t 4
  $WDIR/run-polap1.sh assemble2 -s 1 -i 1 -t 4
  $WDIR/run-polap1.sh assemble2 -m 7000 -i 2 -t 4
  $WDIR/run-polap1.sh assemble2 -a 1000 -m 5000 -i 2 -t 4
  $WDIR/run-polap1.sh assemble3 -a short1.fq -b short2.fq -p mt.0.fasta -f mt.1.fa
  $WDIR/run-polap1.sh coverage -p mt.0.fasta
EOM
	exit 1
}

function assemble1() {
	while getopts "l:a:b:o:t:" option; do
		case $option in
		l) LR=$OPTARG ;;
		a) SR1=$OPTARG ;;
		b) SR2=$OPTARG ;;
		o) ODIR=$OPTARG ;;
		t) NT=$OPTARG ;;
		\?) # incorrect option
			echo "Error: Invalid option; try -h option"
			exit 1
			;;
		esac
	done

	if run_check1; then
		echo "ERROR: change your conda environment to polap-dev."
		echo "INFO: (base) $ conda env create -f src/environment.yaml"
		echo "INFO: (base) $ conda activate polap-dev"
		exit $EXIT_ERROR
	fi
	run_output1

	if [[ ! -s $LR ]]; then
		echo "ERROR: no long-read data file found: $LR"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $SR1 ]]; then
		echo "ERROR: no short-read data file found: $SR1"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $SR2 ]]; then
		echo "ERROR: no short read-data file found: $SR2"
		exit $EXIT_ERROR
	fi
	LR3K=$ODIR/n3k.fq.gz

	run_jellyfish
	run_annotation
	run_count_genes
}

function assemble2() {
	MR=3000
	while getopts "i:j:o:t:m:a:b:c:u" option; do
		case $option in
		i) STEP4=$OPTARG ;;
		j) JDIR=$OPTARG ;;
		o) ODIR=$OPTARG ;;
		t) NT=$OPTARG ;;
		m) MR=$OPTARG ;;
		a) MPAIR=$OPTARG ;;
		b) MBRIDGE=$OPTARG ;;
		c) COV=$OPTARG ;;
		u) CIRCULARIZE="--circularize" ;;
		\?) # incorrect option
			echo "Error: Invalid option; try -h option"
			exit
			;;
		esac
	done
	LR3K=$ODIR/n3k.fq.gz
	if [[ "$MR" -lt 3000 ]]; then
		echo "ERROR: your minimum read length should be at least 3000 bp"
		exit $EXIT_ERROR
	fi

	# if [[ "$MR" -eq 3000 ]]; then
	# 	LR3K=$ODIR/n3k.fq.gz
	# elif [[ "$MR" -gt 3000 ]]; then
	# 	LR3K=$ODIR/nk.fq.gz
	# 	if [ -s "$LR" ]; then
	# 		seqkit seq --quiet -m "$MR" --threads 4 "$LR" -o $LR3K >/dev/null 2>&1
	# 		echo "DATA: long-read minimum $MR reads data nk.fq.gz is created"
	# 	else
	# 		echo "ERROR: no long-read data file while your minimum read length is less than 3000 bp"
	# 		exit $EXIT_ERROR
	# 	fi
	# else
	# 	echo "ERROR: your minimum read length should be at least 3000 bp"
	# 	exit $EXIT_ERROR
	# fi

	if run_check1; then
		echo "ERROR: change your conda environment to polap-dev."
		echo "INFO: (base) $ conda env create -f src/environment.yaml"
		echo "INFO: (base) $ conda activate polap-dev"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $LR3K ]]; then
		echo "ERROR: no $LR3K long read-data file found: $LR3K"
		exit $EXIT_ERROR
	fi

	run_flye2 $STEP4 $MR
}

function test() {
	STEP4=1
	ODIR=otest
	NT=56
	MR=3000
	while getopts "i:o:t:m:c:u" option; do
		case $option in
		i) STEP4=$OPTARG ;;
		o) ODIR=$OPTARG ;;
		t) NT=$OPTARG ;;
		m) MR=$OPTARG ;;
		c) COV=$OPTARG ;;
		u) CIRCULARIZE="--circularize" ;;
		\?) # incorrect option
			echo "Error: Invalid option; try -h option"
			exit
			;;
		esac
	done
	LR3K=$ODIR/n3k.fq.gz
	if [[ "$MR" -lt 3000 ]]; then
		echo "ERROR: your minimum read length should be at least 3000 bp"
		exit $EXIT_ERROR
	fi

	if run_check1; then
		echo "ERROR: change your conda environment to polap-dev."
		echo "INFO: (base) $ conda env create -f src/environment.yaml"
		echo "INFO: (base) $ conda activate polap-dev"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $LR3K ]]; then
		echo "ERROR: no $LR3K long read-data file found: $LR3K"
		exit $EXIT_ERROR
	fi

	run_annotation
	# run_count_genes
	run_flye2 $STEP4 $MR
}

function assemble3() {
	while getopts "a:b:p:f:o:t:" option; do
		case $option in
		a) SR1=$OPTARG ;;
		b) SR2=$OPTARG ;;
		p) PA=$OPTARG ;;
		f) FA=$OPTARG ;;
		o) ODIR=$OPTARG ;;
		t) NT=$OPTARG ;;
		\?) # incorrect option
			echo "Error: Invalid option; try -h option"
			exit
			;;
		esac
	done

	if run_check2; then
		echo "ERROR: change your conda environment to polap-fmlrc."
		echo "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
		echo "INFO: (base) $ conda activate polap-fmlrc"
		exit $EXIT_ERROR
	fi

	if [[ ! -s $SR1 ]]; then
		echo "ERROR: no short read-data file found: $SR1"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $SR2 ]]; then
		echo "ERROR: no short read-data file found: $SR2"
		exit $EXIT_ERROR
	fi

	LR3K=$ODIR/n3k.fq.gz
	run_polish1
	run_polish2
}

function coverage() {
	while getopts "o:p:t:" option; do
		case $option in
		o) ODIR=$OPTARG ;;
		p) PA=$OPTARG ;;
		t) NT=$OPTARG ;;
		\?) # incorrect option
			echo "Error: Invalid option; try -h option"
			exit
			;;
		esac
	done
	LR3K=$ODIR/n3k.fq.gz

	if run_check1; then
		echo "ERROR: change your conda environment to polap-dev."
		echo "INFO: (base) $ conda env create -f src/environment.yaml"
		echo "INFO: (base) $ conda activate polap-dev"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $PA ]]; then
		echo "ERROR: no draft sequence: $PA"
		exit $EXIT_ERROR
	fi
	if [[ ! -s $LR3K ]]; then
		echo "ERROR: no $LR3K long read-data file found: $LR3K"
		exit $EXIT_ERROR
	fi

	run_samtools
}
### conda activate polap-dev
# run_jellyfish
# run_annotation
# run_count_genes
# run_flye2 $1
#
# run_samtools
### conda activate polap-fmlrc
# run_polish1
# run_polish2

################################################################
# MAIN
if [ $# -eq 0 ]; then
	usage
fi

CMD="src/run-polap1.sh $*"
echo "CMD: $CMD"

# subcommand function call
if declare -f "$1" >/dev/null 2>&1; then
	# invoke that function, passing arguments through
	"$@" # same as "$1" "$2" "$3" ... for full argument list
else
	echo ERROR: no such subcommand "$1"
	exit $EXIT_ERROR
fi

ELAPSED="Time: $((SECONDS / 3600))hrs $(((SECONDS / 60) % 60))min $((SECONDS % 60))sec - $CMD"
echo "$ELAPSED"
