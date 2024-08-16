#!/usr/bin/env bash

LONG_TOTAL_LENGTH=$1

echo -n "Please, wait ... "
jellyfish count -t 56 -C -m 19 -s 5G -o o/jellyfish_out '--min-qual-char=?' s1.fq s2.fq
jellyfish histo -o o/jellyfish_out.histo o/jellyfish_out
src/run-polap-jellyfish.R o/jellyfish_out.histo $LONG_TOTAL_LENGTH o/long_coverage.txt o/short_expected_genome_size.txt
echo "done!"
EXPECTED_GENOME_SIZE=$(cat o/short_expected_genome_size.txt)
EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
echo "INFO: short-read expected genome size in o/short_expected_genome_size.txt $EXPECTED_GENOME_SIZE"
echo Step 1 has been finished!
echo Your next command is:
echo "$ src/run-polap-mtdna-carex-pseudochinensis-03.sh $EXPECTED_GENOME_SIZE"
