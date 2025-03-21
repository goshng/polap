#!/bin/bash

package=${1:-polap}

A=$HOME/all/polap/bioconda-recipes-polap/recipes/$package/
cp -pu $package-build.sh $A/build.sh
cp -pu $package-meta.yaml $A/meta.yaml

echo We have copied the $package-meta.yaml and $package-build.sh to bioconda-recipes-polap recipe folder: $A
