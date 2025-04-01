#!/bin/bash

bash src/polap-data-v2.sh get all 2 add off
# bash src/polap-data-v2.sh get all 0 add off
bash src/polap-data-v2.sh get all 3 add off
bash src/polap-data-v2.sh report all 2 infer-1
bash src/polap-data-v2.sh report all 0 infer-1
bash src/polap-data-v2.sh report all 3 infer-1
bash src/polap-data-v2.sh maintable1 all 2 infer-1
bash src/polap-data-v2.sh maintable1 all 0 infer-1
bash src/polap-data-v2.sh maintable1 all 3 infer-1
bash src/polap-data-v2.sh supptable1 all 2 infer-1 x
bash src/polap-data-v2.sh supptable1 Eucalyptus_pauciflora 2 infer-1 x
bash src/polap-data-v2.sh supptable1 Eucalyptus_pauciflora 0 infer-1 x
bash src/polap-data-v2.sh supptable1 Eucalyptus_pauciflora 3 infer-1 x
bash src/polap-data-v2.sh suppfigure1 all 2 infer-1 1 1 yes
bash src/polap-data-v2.sh suppfigure1 Eucalyptus_pauciflora 2 infer-1 1 1 yes
bash src/polap-data-v2.sh suppfigure3 2 infer-1 on yes
bash src/polap-data-v2.sh suppfigure3 2 infer-1 off yes
bash src/polap-data-v2.sh copy-figures

exit

Anthoceros_agrestis-a-2.tar.gz
Arabidopsis_thaliana-a-2.tar.gz
Canavalia_ensiformis-a-2.tar.gz
Cinchona_pubescens-a-2.tar.gz
Codonopsis_lanceolata-a-2.tar.gz
Cucumis_sativus_var_hardwickii-a-2.tar.gz
Dioscorea_japonica-a-2.tar.gz
Dunaliella_tertiolecta-a-2.tar.gz
Eucalyptus_pauciflora-a-2.tar.gz
Euonymus_alatus-a-2.tar.gz
Gossypium_herbaceum-a-2.tar.gz
Juncus_effusus-a-2.tar.gz
Juncus_inflexus-a-2.tar.gz
Juncus_roemerianus-a-2.tar.gz
Juncus_validus-a-2.tar.gz
Leiosporoceros_dussii-a-2.tar.gz
Macadamia_jansenii-a-2.tar.gz
Musa_acuminata_subsp_malaccensis-a-2.tar.gz
Notothylas_orbicularis-a-2.tar.gz
Ophrys_lutea-a-2.tar.gz
Oryza_rufipogon-a-2.tar.gz
Phaeomegaceros_chiloensis-a-2.tar.gz
Populus_x_sibirica-a-2.tar.gz
Prunus_mandshurica-a-2.tar.gz
Spirodela_polyrhiza-a-2.tar.gz

Vaccinium_vitis-idaea-a-2.tar.gz
