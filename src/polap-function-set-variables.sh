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
# Template for an external shell script
#
# You could use this function template to create a new menu.
# Rename template and delete x in the name. You could execute such menu
# but such menus are not created as empty files by make-menus menu.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Set the following variables:
# ${SR1}
# ${SR2}
# depending on the options provided.
################################################################################
function _polap_set-variables-short-read() {

	# Set paths for bioproject data
	source "$script_dir/polap-variables-base.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'
	source "$script_dir/polap-variables-oga.sh"        # '.' means 'source'
	source "$script_dir/run-polap-function-utilities.sh"

	# if --bioproject is used, we use 0-bioproject.
	# otherwise, -l option is used.
	if [ "${_arg_short_read1_is}" = "off" ]; then
		if [ -z "${_arg_bioproject}" ]; then
			_polap_log1 "we use the default short-read1 data filename: ${_arg_short_read1}"
			SR1="${_arg_short_read1}"
		else
			_polap_log1 "we use the short-reads data 1 info: ${_polap_var_bioproject_sra_short_read}"
			check_file_existence "${_polap_var_bioproject_sra_short_read}"
			SRA=$(cut -f1 "${_polap_var_bioproject_sra_short_read}")
			SR1="${ODIR}/${SRA}_1.fastq"
		fi
	else
		_polap_log1 "we use the short-read given by the option -a"
		SR1="${_arg_short_read1}"
	fi

	if [ "${_arg_short_read2_is}" = "off" ]; then
		if [ -z "${_arg_bioproject}" ]; then
			_polap_log1 "we use the default short-read2 data filename: ${_arg_short_read2}"
			SR2="${_arg_short_read2}"
		else
			_polap_log1 "we use the short-reads data 2 info: ${_polap_var_bioproject_sra_short_read}"
			check_file_existence "${_polap_var_bioproject_sra_short_read}"
			SRA=$(cut -f1 "${_polap_var_bioproject_sra_short_read}")
			SR2="${ODIR}/${SRA}_2.fastq"
		fi
	else
		_polap_log1 "we use the short-read given by the option -b"
		SR2="${_arg_short_read2}"
	fi

}

function _polap_set-variables-long-read() {

	# Set paths for bioproject data
	source "$script_dir/polap-variables-base.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'
	source "$script_dir/polap-variables-oga.sh"        # '.' means 'source'
	source "$script_dir/run-polap-function-utilities.sh"

	# if --bioproject is used, we use 0-bioproject.
	# otherwise, -l option is used.
	if [ "${_arg_long_reads_is}" = "off" ]; then
		if [ -z "${_arg_bioproject}" ]; then
			_polap_log1 "we use the default long-reads data filename: ${_arg_long_reads}"
			LR="${_arg_long_reads}"
		else
			_polap_log1 "we use the long-reads data: ${_polap_var_bioproject_sra_long_read}"
			check_file_existence "${_polap_var_bioproject_sra_long_read}"
			SRA=$(cut -f1 "${_polap_var_bioproject_sra_long_read}")
			LR="${ODIR}/${SRA}.fastq"
		fi
	else
		_polap_log1 "we use the long-reads given by the option -l"
		LR="${_arg_long_reads}"
	fi
}
