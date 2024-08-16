#!/usr/bin/env bash

################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with 
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

function usage() {
    cat <<EOM
$0 [options]

Usage:
    run-polap-confirm-of/genes <N0.tsv> 

A mitochondarial genome assembly pipeline.

EOM
    exit 1
}

# MAIN
if [ $# -eq 0 ]; then
    usage
fi

MDIR=$(dirname "$0")

N0=$1

N=$(cat $N0 | awk -F'\t' '{print NF}' | sort -nu | tail -n 1)

cut -f$N $N0 > of/genes.target.txt

rm -rf of/genes
mkdir of/genes

for (( i=4; i<$N; i++ )); do
  echo -n '-'
  cut -f1,$i $N0 > 1
  A=$(head -1 1 | cut -f2)
  seqkit seq -n of/$A.fa | cut -d' ' -f1-2 > 2

  k=$((i - 1))
  printf -v j "%02d" $i
  printf -v jk "%02d" $k

  "$MDIR"/run-polap-confirm-genes.R 1 2 $A of/genes/$j
  cut -f1 of/genes/$j > of/genes/$j.og
  cut -f2 of/genes/$j > of/genes/$j.gene
  diff of/genes/$j.og of/genes/$jk.og

done

paste of/genes/*.gene > of/genes.csv

echo You need to:
echo open of/genes.csv
echo append the content of of/genes.target.txt

"$MDIR"/run-polap-confirm-genes2.R of/genes.csv of/genes.target.txt phylogenies/genes.tsv

echo Check: phylogenies/genes.tsv
