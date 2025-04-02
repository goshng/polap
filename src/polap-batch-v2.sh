#!/bin/bash

bash src/polap-data-v2.sh local-batch all 2 0 on
bash src/polap-data-v2.sh local-batch all 0 0 on
bash src/polap-data-v2.sh local-batch all 3 0 on
bash src/polap-data-v2.sh get all 2 add off
bash src/polap-data-v2.sh get all 0 add off
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
