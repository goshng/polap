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

# ERR constants
_POLAP_ERR_NO_EDGES_GFA=11
_POLAP_ERR_NO_GENOME_SIZE=12
_POLAP_ERR_NO_NK_FQ=13
_POLAP_ERR_NO_SEEDS=14
_POLAP_ERR_NO_DISK_SPACE=15
_POLAP_ERR_CMD_OPTION_STEPS=16
_POLAP_ERR_CMD_OPTION_GENOMESIZE=17
_POLAP_ERR_CMD_OPTION_LONGREAD=18
_POLAP_ERR_CMD_OPTION_SHORTREAD=19
_POLAP_ERR_TOO_SMALL_SUBSAMPLE_SIZE=20
_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES=21
_POLAP_ERR_ALREADY_EXIST_OUT=31
_POLAP_ERR_MENU_MAP_READS=51

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
	${_POLAP_ERR_CMD_OPTION_STEPS})
		_polap_log0 "ERROR: --steps-include or --steps-exclude use a range (e.g., 2-4) or list (e.g., 1,3,4)."
		;;
	${_POLAP_ERR_CMD_OPTION_GENOMESIZE})
		_polap_log0 "ERROR: --genomesize is missing?"
		;;
	${_POLAP_ERR_CMD_OPTION_LONGREAD})
		_polap_log0 "ERROR: no long-read fastq file"
		;;
	${_POLAP_ERR_CMD_OPTION_SHORTREAD})
		_polap_log0 "ERROR: no short-read fastq file"
		;;
	${_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES})
		_polap_log0 "ERROR: subsampling has too small number of candidate assemblies"
		;;
	${_POLAP_ERR_ALREADY_EXIST_OUT})
		_polap_log0 "ERROR: output folder ${_arg_outdir} already exists."
		;;
	*)
		echo default
		;;
	esac
}
