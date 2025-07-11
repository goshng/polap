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
# This script sets variables for short- and long-read data files.
# It is used by functions that need read data files.
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
# Set the following variables:
# ${SR1}
# ${_arg_short_read2}
# depending on the options provided.
################################################################################
function _polap_set-variables-short-read {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	_polap_log1 "  determining short-read data files ..."

	# if --bioproject is used, we use 0-bioproject.
	# otherwise, -l option is used.
	if [ "${_arg_short_read1_is}" = "off" ]; then
		if [ -z "${_arg_bioproject}" ]; then
			_polap_log1 "    we use the default short-read1 data filename: ${_arg_short_read1}"
		else
			_polap_log1 "    we use the short-reads data 1 info: ${_polap_var_project_sra_short_read}"
			check_file_existence "${_polap_var_project_sra_short_read}"
			local SRA=$(cut -f1 "${_polap_var_project_sra_short_read}")
			_arg_short_read1="${_arg_outdir}/${SRA}_1.fastq"
		fi
	else
		_polap_log1 "    we use the short-read given by the option -a"
	fi

	if [ "${_arg_short_read2_is}" = "off" ]; then
		if [ -z "${_arg_bioproject}" ]; then
			_polap_log1 "    we use the default short-read2 data filename: ${_arg_short_read2}"
		else
			_polap_log1 "    we use the short-reads data 2 info: ${_polap_var_project_sra_short_read}"
			check_file_existence "${_polap_var_project_sra_short_read}"
			local SRA=$(cut -f1 "${_polap_var_project_sra_short_read}")
			_arg_short_read2="${_arg_outdir}/${SRA}_2.fastq"
		fi
	else
		_polap_log1 "    we use the short-read given by the option -b"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _polap_set-variables-long-read {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	_polap_log1 "  determining long-read data file ..."

	# if --bioproject is used, we use 0-bioproject.
	# otherwise, -l option is used.
	if [ "${_arg_long_reads_is}" = "off" ]; then
		if [ -z "${_arg_bioproject}" ]; then
			_polap_log1 "    we use the default long-reads data filename: ${_arg_long_reads}"
		else
			_polap_log1 "    we use the long-reads data: ${_polap_var_project_sra_long_read}"
			check_file_existence "${_polap_var_project_sra_long_read}"
			local SRA=$(cut -f1 "${_polap_var_project_sra_long_read}")
			_arg_long_reads="${_arg_outdir}/${SRA}.fastq"
		fi
	else
		_polap_log1 "    we use the long-reads given by the option -l"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
