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
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
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

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
function _polap_lib_unit-convert_bp {
	local _bp="$1"
	if [[ "$_bp" == "NA" ]]; then
		echo "NA"
		return
	fi

	local bp=${1%.*}
	if ((bp >= 1000000000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000000000") Tbp"
	elif ((bp >= 1000000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000000") Gbp"
	elif ((bp >= 1000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000") Mbp"
	elif ((bp >= 1000)); then
		echo "$(bc <<<"scale=1; $bp/1000") kbp"
	else
		echo "$bp bp"
	fi
}
