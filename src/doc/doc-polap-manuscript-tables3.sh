#!/bin/bash

# mkdir t
# find polap-test -type d -links 2 -exec mkdir -p "t/{}" \;
# mv t/polap-test polap-github

S=('Spirodela_polyrhiza' 'Taraxacum_mongolicum' 'Trifolium_pratense' 'Salix_dunnii' 'Anthoceros_agrestis' 'Anthoceros_angustus' 'Brassica_rapa' 'Vigna_radiata' 'Macadamia_tetraphylla' 'Punica_granatum' 'Lolium_perenne')

for i in "${S[@]}"; do
	echo Copying ... $i

	src/contig-seeds.R run/$i/o/0/contig-annotation-table.txt run/$i/o/0/mt.contig.name-1 | pandoc -f tsv -t markdown_mmd >doc/polap-manuscript/tables3/$i.md
done
