#!/usr/bin/env bash

JNUM=$1
EXPECTED_COVERAGE=$2
CONTIG_LENGTH=$3
COV=30

if [ "$EXPECTED_COVERAGE" -lt $COV ]; then
	echo "DATA: no data reduction: COV=$COV"
	ln -s 1.fq.gz o/60-mt-$JNUM/o3000/seeds/2.fq.gz
else
	echo "SUGGESTION: you might want to increase the minimum read lengths"
	echo "SUGGESTION: since you have enough long-read data."
	RATE=$(echo "scale=10; $COV/$EXPECTED_COVERAGE" | bc)
	echo "DATA: data reduction by rate of $RATE"
	seqkit sample -p $RATE o/60-mt-$JNUM/o3000/seeds/1.fq.gz -o o/60-mt-$JNUM/o3000/seeds/2.fq.gz >/dev/null 2>&1
fi

echo Step 6 is ready to run!
echo Your next command is:
echo "$ src/run-polap-mtdna-carex-pseudochinensis-08.sh <MT.CONTIG.NAME-NUMBER> $CONTIG_LENGTH"
echo MT.CONTIG.NAME-NUMBER is an integer: e.g., 1.
