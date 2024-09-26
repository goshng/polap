#/usr/bin/bash

source src/run-polap-function-utilities.sh

# 1000000000
#
v1=1000000000
v2=$(_polap_utility_convert_bp $v1)
echo ${v2}
