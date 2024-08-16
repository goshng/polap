#!/usr/bin/env bash

EXPECTED_GENOME_SIZE=$1

seqkit seq --quiet -m 3000 --threads 4 l.fq -o o/n3k.fq.gz >/dev/null 2>&1
flye --nano-raw o/n3k.fq.gz --out-dir o --threads 56 --stop-after contigger --asm-coverage 30 --genome-size $EXPECTED_GENOME_SIZE
echo Step 2 has been finished!
echo Now, scan the contigs with organelle genes.
echo Your next command is:
echo "$ src/run-polap-mtdna-carex-pseudochinensis-04.sh"
