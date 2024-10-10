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
# Ensure that the current script is sourced only once
SCRIPT_NAME="${BASH_SOURCE[0]}"
SCRIPT_NAME=$(basename "${SCRIPT_NAME}")
SCRIPT_NAME="${SCRIPT_NAME//-/_}"
SCRIPT_NAME="${SCRIPT_NAME//./_}"
SCRIPT_NAME=$(echo "$SCRIPT_NAME" | tr '[:lower:]' '[:upper:]')
_POLAP_INCLUDE_="_POLAP_INCLUDE_${SCRIPT_NAME}"
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

function _polap_include() {
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
# [[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
# declare "$_POLAP_INCLUDE_=1"
#
################################################################################
