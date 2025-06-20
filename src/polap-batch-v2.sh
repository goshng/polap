#!/bin/bash

_polap_script_bin_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
	echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
	exit 2
}

# Data28
S4=(
	Anthoceros_agrestis
	Codonopsis_lanceolata
	Juncus_roemerianus
)

_polap_data_cmd="${_polap_script_bin_dir}/polap-data-v2.sh"

for i in 0 1 2 3 4; do
	# ${_polap_data_cmd} local-batch all $i 0 on # in release mode only

	# ${_polap_data_cmd} get all $i add off

	# ${_polap_data_cmd} report all $i infer-1 # not use anymore

	${_polap_data_cmd} maintable1 all $i infer-1
	${_polap_data_cmd} supptable1 Eucalyptus_pauciflora $i infer-1 x
done

for i in "${S4[@]}"; do
	${_polap_data_cmd} supptable1 $i 2 infer-1 x
	# ${_polap_data_cmd} suppfigure1 $i 2 infer-1 1 1 yes
done

# ${_polap_data_cmd} supptable1 all 2 infer-1 x
# ${_polap_data_cmd} suppfigure1 all 2 infer-1 1 1 yes
${_polap_data_cmd} suppfigure3 2 infer-1 on yes
${_polap_data_cmd} suppfigure3 2 infer-1 off yes
${_polap_data_cmd} copy-figures

exit

${_polap_data_cmd} local-batch all 2 0 on
${_polap_data_cmd} local-batch all 0 0 on
${_polap_data_cmd} local-batch all 3 0 on
${_polap_data_cmd} get all 2 add off
${_polap_data_cmd} get all 0 add off
${_polap_data_cmd} get all 3 add off
${_polap_data_cmd} report all 2 infer-1
${_polap_data_cmd} report all 0 infer-1
${_polap_data_cmd} report all 3 infer-1
${_polap_data_cmd} maintable1 all 2 infer-1
${_polap_data_cmd} maintable1 all 0 infer-1
${_polap_data_cmd} maintable1 all 3 infer-1

${_polap_data_cmd} supptable1 all 2 infer-1 x

${_polap_data_cmd} supptable1 Eucalyptus_pauciflora 2 infer-1 x
${_polap_data_cmd} supptable1 Eucalyptus_pauciflora 0 infer-1 x
${_polap_data_cmd} supptable1 Eucalyptus_pauciflora 3 infer-1 x

${_polap_data_cmd} suppfigure1 all 2 infer-1 1 1 yes
${_polap_data_cmd} suppfigure1 Eucalyptus_pauciflora 2 infer-1 1 1 yes
${_polap_data_cmd} suppfigure3 2 infer-1 on yes
${_polap_data_cmd} suppfigure3 2 infer-1 off yes
${_polap_data_cmd} copy-figures

exit

# Data28
S28=(
	Anthoceros_agrestis
	Arabidopsis_thaliana
	Canavalia_ensiformis
	Cinchona_pubescens
	Codonopsis_lanceolata
	Cucumis_sativus_var_hardwickii
	Dioscorea_japonica
	Dunaliella_tertiolecta
	Eucalyptus_pauciflora
	Euonymus_alatus
	Gossypium_herbaceum
	Juncus_effusus
	Juncus_inflexus
	Juncus_roemerianus
	Juncus_validus
	Leiosporoceros_dussii
	Macadamia_jansenii
	Musa_acuminata_subsp_malaccensis
	Notothylas_orbicularis
	Ophrys_lutea
	Oryza_rufipogon
	Phaeomegaceros_chiloensis
	Populus_x_sibirica
	Prunus_mandshurica
	Solanum_lycopersicum
	Spirodela_polyrhiza
	Vaccinium_vitis-idaea
	Vitis_vinifera
)
