# 1
_run_polap_miniassemble

# 1
_polap_readassemble-pt

# 2
_polap_lib_readassemble-assemble-annotated-read-pt

bash "${_POLAPLIB_DIR}"/polap-bash-filter-pt-reads.sh

_polap_lib_assemble-sub
