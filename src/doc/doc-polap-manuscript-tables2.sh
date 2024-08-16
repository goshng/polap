#!/bin/bash

Species=('Spirodela_polyrhiza' 'Taraxacum_mongolicum' 'Trifolium_pratense' 'Salix_dunnii' 'Anthoceros_agrestis' 'Anthoceros_angustus' 'Brassica_rapa' 'Vigna_radiata' 'Macadamia_tetraphylla' 'Punica_granatum' 'Lolium_perenne')

for S in "${Species[@]}"; do
	A=$HOME/all/polap/proj/run/$S/o/0/assembly_info_organelle_annotation_count.txt
	T=$HOME/all/polap/proj/doc/polap-manuscript/tables2/$S-0.tsv
	sed 's/ /\t/g' $A | cut -f1,2,6,9,10,11 >$T
done
