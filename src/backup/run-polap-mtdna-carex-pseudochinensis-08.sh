#!/usr/bin/env bash

JNUM=$1
CONTIG_LENGTH=$2

flye --nano-raw o/60-mt-$JNUM/o3000/seeds/2.fq.gz --out-dir o/60-mt-$JNUM/o3000 --asm-coverage 30 --genome-size $CONTIG_LENGTH --threads 56

echo RESULT: assembly graph: $PWD/o/60-mt-$JNUM/o3000/assembly_graph.gfa
echo INFO: extract a draft organelle genome sequence from the assembly graph
echo "  See the assembly graph by using Bandage availabe at https://rrwick.github.io/Bandage/"
echo "You could watch YouTube: https://youtu.be/UF-0UIc2ZDY"
echo "  for extracting a DNA sequence through a circular path."

echo Step 6 has been finished!
echo Your next command is:
echo "$ src/run-polap-mtdna-carex-pseudochinensis-09.sh"