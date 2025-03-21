#!/bin/bash

package=${1:-polap}

for i in $HOME/all/polap/bioconda-recipes-polap/recipes/$package/*; do
	cp -pu "$i" $package-$(basename "$i")
done

echo We have $package recipes at $PWD
