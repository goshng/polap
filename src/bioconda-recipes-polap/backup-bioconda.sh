#!/bin/bash

for i in $HOME/all/polap/bioconda-recipes-polap/recipes/polap/*; do
	cp -pu "$i" polap-$(basename "$i")
done

for i in $HOME/all/polap/bioconda-recipes-polap/recipes/cflye/*; do
	cp -pu "$i" cflye-$(basename "$i")
done

for i in $HOME/all/polap/bioconda-recipes-polap/recipes/dflye/*; do
	cp -pu "$i" dflye-$(basename "$i")
done

echo We have the bioconda-recipes-polap recipes at $PWD
