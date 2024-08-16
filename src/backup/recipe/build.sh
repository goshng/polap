#!/usr/bin/env bash

mkdir -p $PREFIX/bin

for i in polap2.sh \
         run-polap-jellyfish.R \
         run-polap-mtcontig.R \
         run-polap-genes.R \
         run-polap-pairs.R \
         mt.1.c70.3.faa \
         pt.2.c70.3.faa
do
  cp $RECIPE_DIR/../src/$i $PREFIX/bin
done


