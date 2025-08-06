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
# Convert numbers between different units.
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
# Function to convert gfa to bandage png
# input1: gfa infile
# output: png outfile
# type: png
# colar: label0, label, label2, and other
################################################################################
_polap_lib_bandage() {
	local _infile="${1}"
	local _outfile="${2}"
	local _type="${3:-png}"
	local _color="${4:-label0}"

	if [[ "${_type}" == "png" ]]; then

		if [[ "${_color}" == "label0" ]]; then
			# gray color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour uniform \
				--unicolpos \#EEEEEE \
				--singlearr \
				--toutline 1
		elif [[ "${_color}" == "label" ]]; then
			# gray color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour uniform \
				--unicolpos \#EEEEEE \
				--singlearr \
				--toutline 1 \
				--names \
				--lengths \
				--depth
		elif [[ "${_color}" == "label2" ]]; then
			# gray color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour uniform \
				--unicolpos \#EEEEEE \
				--singlearr \
				--toutline 1 \
				--names \
				--lengths \
				--depth \
				--fontsize 3
		else
			# random color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour random \
				--singlearr \
				--toutline 1 \
				--names \
				--lengths \
				--depth \
				--fontsize 3
		fi
	fi

}
