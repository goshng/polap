#!/usr/bin/env bash
rm -rf o
echo INFO: deletes the output directory: o
mkdir o
echo INFO: creates a new output directory: o

LONG_TOTAL_LENGTH=$(seqkit stats -Ta l.fq | csvtk cut -t -f sum_len | csvtk del-header)
echo "DATA: total number of bases in the long-read data $LONG_TOTAL_LENGTH"
echo Step 1 is ready to run!
echo Your next command is:
echo $ src/run-polap-mtdna-carex-pseudochinensis-02.sh
