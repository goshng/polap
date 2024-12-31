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
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

_polap_error_message() {
	local _polap_errno=$1

	case "${_polap_errno}" in
	${_POLAP_ERR_NO_EDGES_GFA})
		_polap_log0 "ERROR: No contigger edges gfa: ${_polap_var_ga_contigger_edges_gfa}"
		;;
	${_POLAP_ERR_NO_GENOME_SIZE})
		_polap_log0 "ERROR: No genome size estimate"
		;;
	${_POLAP_ERR_NO_NK_FQ})
		_polap_log0 "ERROR: No long-read data for the whole-genome assembly"
		;;
	${_POLAP_ERR_NO_SEEDS})
		_polap_log0 "ERROR: no seed contig file of the contig selection types!"
		;;
	${_POLAP_ERR_NO_DISK_SPACE})
		_polap_log0 "ERROR: no disk space available for a POLAP analysis!"
		;;
	2 | 3)
		echo 2 or 3
		;;
	*)
		echo default
		;;
	esac
}
