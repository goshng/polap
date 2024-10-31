source run-polap-function-include.sh

################################################################################
# Ensure that the current script is sourced only once
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u; [[ -n "${!_POLAP_INCLUDE_}" ]] && return 0; set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

# The rest of your script content here
echo "Script is being sourced for the first time: $SCRIPT_NAME & $_POLAP_INCLUDE_"
