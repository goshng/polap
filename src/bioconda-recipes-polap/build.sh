#!/usr/bin/env bash

mkdir -p $PREFIX/bin

files=(polap
	polap-conda-environment-fmlrc.yaml
	polap-conda-environment.yaml
	polap-mt.1.c70.3.faa
	polap-parsing.sh
	polap-pt.2.c70.3.faa
	polap.sh
	run-polap-genes.R
	run-polap-jellyfish.R
	run-polap-mtcontig.R
	run-polap-pairs.R)

for i in "${files[@]}"; do
	cp src/$i $PREFIX/bin
done

chmod +x $PREFIX/bin/polap
chmod +x $PREFIX/bin/polap.sh
<<<<<<< HEAD
chmod +x $PREFIX/bin/run-polap-*.R
=======
>>>>>>> 587e5ac5efb38ad0aa0f8a5d6047d3eda435275c
