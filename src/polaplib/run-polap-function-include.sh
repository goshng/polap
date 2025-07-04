################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# This script emulates the include in C.
# It includes or sources a bash shell script only once.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
SCRIPT_NAME="${BASH_SOURCE[0]}"
SCRIPT_NAME=$(basename "${SCRIPT_NAME}")
SCRIPT_NAME="${SCRIPT_NAME//-/_}"
SCRIPT_NAME="${SCRIPT_NAME//./_}"
SCRIPT_NAME=$(echo "$SCRIPT_NAME" | tr '[:lower:]' '[:upper:]')
_POLAP_INCLUDE_="_POLAP_INCLUDE_${SCRIPT_NAME}"
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
  set -u
  return 0
fi 
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
  return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

function _polap_include {
	local input_string="$1"
	# Extract the basename
	local script_name=$(basename "${input_string}")
	# Replace hyphens and dots with underscores
	script_name="${script_name//-/_}"
	script_name="${script_name//./_}"
	# Convert to uppercase
	script_name=$(echo "$script_name" | tr '[:lower:]' '[:upper:]')
	# Create the _POLAP_INCLUDE_ variable
	local _POLAP_INCLUDE_="_POLAP_INCLUDE_${script_name}"

	echo "${_POLAP_INCLUDE_}"
}

################################################################################
# Ensure that the current script is sourced only once
# _POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
# set +u; [[ -n "${!_POLAP_INCLUDE_}" ]] && return 0; set -u
# declare "$_POLAP_INCLUDE_=1"
#
################################################################################
