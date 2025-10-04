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
# Arichve polap results.
# Polap creates many files with different sizes.
# Instead of saving particular result files, we save them by size limit.
# We also use a text file, which has file paths to save even if they are
# greater than the size limit.
#
# Functions:
# function _run_polap_archive-rsync-template { # archive a POLAP-cflye output folder for later use
# function _run_polap_archive { # archive a POLAP output folder for later use
# function _run_polap_cleanup { # cleanup an POLAP output folder
# function _run_polap_init { # initialize an output folder
# function _run_polap_log { # display the polap log
#
# See Also:
# polaplib/polap-template-aflye-archive-files.txt
# polaplib/polap-template-cflye-archive-files.txt
#
# TEST-SCC: not yet
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

function _run_polap_archive-rsync-template { # archive a POLAP-cflye output folder for later use
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	if [ "${_arg_short_read1_is}" = "on" ]; then
		_arg_archive="${_arg_short_read1}"
	fi

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Archive polap-cflye the ${_arg_outdir} folder to ${_arg_archive}
#
# First, use the following rsync command to copy a folder 
# to another new folder, including only files smaller than 2MB:
# rsync -aq --max-size=2M source_folder/ destination_folder/
# Second, use a template path file to copy more files.
#
# Arguments:
#   -o ${_arg_outdir}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
#   --template polap-archive-template.txt
#   --max-sizefile: 1M not implemented yet!
# Inputs:
#   ${_arg_outdir}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log1_log "archiving ${_arg_outdir} to ${_arg_archive} ... upto ${_arg_max_filesize}"

	# rsync -aq --max-size=2M source_folder/ destination_folder/
	_polap_lib_file-rsync "${_arg_outdir}" "${_arg_archive}" "${_arg_max_filesize}"
	_polap_lib_file-archive-folder \
		"${_arg_outdir}" "${_arg_archive}" "${_arg_template}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Archive the ${_arg_outdir} folder to ${_arg_archive}
################################################################################
function _run_polap_archive { # archive a POLAP output folder for later use
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	if [[ "${_arg_short_read1_is}" == "on" ]]; then
		_arg_archive="${_arg_short_read1}"
	elif [[ "${_arg_archive_is}" == "off" ]]; then
		_arg_archive="${_arg_outdir}-a"
	fi

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	source "${_POLAPLIB_DIR}/polap-package-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Archive the ${_arg_outdir} folder to ${_arg_archive}
#
# Step 1. package
# Step 2. do more
#
# Arguments:
#   -o ${_arg_outdir}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
# Inputs:
#   ${_arg_outdir}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
Example: $(basename "$0") ${_arg_menu[0]} cflye -o <folder1> -a <folder2>
Example: $(basename "$0") ${_arg_menu[0]} aflye -o <folder1> -a <folder2>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if [[ "${_arg_menu[1]}" == "cflye" ]]; then
		_arg_template="${_POLAPLIB_DIR}/polap-template-cflye-archive-files.txt"
		# _arg_max_filesize="1M"
		_run_polap_archive-rsync-template
	elif [[ "${_arg_menu[1]}" == "aflye" ]]; then
		_arg_template="${_POLAPLIB_DIR}/polap-template-aflye-archive-files.txt"
		# _arg_max_filesize="1M"
		_run_polap_archive-rsync-template
	else

		_polap_log0_log "archiving ${_arg_outdir} to ${_arg_archive} with file size upto ${_arg_max_filesize}"

		_run_polap_package

		cp -pr "${_polap_var_outdir_msbwt_dir}" "${_ppack_var_outdir}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Clean up the ${_arg_outdir}.
################################################################################
function _run_polap_cleanup { # cleanup an POLAP output folder
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Clean up the ${_arg_outdir}.
#
# Arguments:
#   -o ${_arg_outdir}: the output directory
#   -j ${_arg_jnum}: the target assembly
# Inputs:
#   ${_arg_outdir}
# Outputs:
#
Example: $(basename "$0") ${_arg_menu[0]} -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log3_cmd rm -f "${_polap_var_outdir_lk_fq_gz}"
	_polap_log3_cmd rm -rf "${_polap_var_oga}"/{00-assembly,10-consensus,20-repeat,40-polishing}
	_polap_log3_cmd rm -rf "${_polap_var_oga_contig}"/{contig.paf,contig.tab}
	for _pread_sel in ptgaul-intra-base-length single-intra-base-length combined-intra-base-length; do
		_polap_log3_cmd rm -rf "${_polap_var_oga_reads}/${_pread_sel}"
		_polap_log3_cmd rm -rf "${_polap_var_oga_seeds}/${_pread_sel}"
		_polap_log3_cmd rm -rf "${_polap_var_oga_subsample}/${_pread_sel}"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"
source "${_POLAPLIB_DIR}/polap-cmd-version.sh"
source "${_POLAPLIB_DIR}/polap-cmd-menus.sh"

################################################################################
# Initializes polap analysis in a starting folder,
#
# creating an output folder.
# Arguments:
#   -o ${_arg_outdir}
# Inputs: nothing
# Outputs:
#   ${_arg_outdir}
################################################################################
function _run_polap_init { # initialize an output folder
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
A polap analysis process is initiated in a specified starting folder
through the creation of an output folder.

1. Create an empty folder.
2. Create empty menu files to facilitate effortless command inputting.
3. Record all external software packages including their versions.

Arguments:
  -o ${_arg_outdir}: the output folder

Inputs: none

Outputs:
  ${_arg_outdir}

Example:
$(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if ! run_check1; then
		error_polap_conda
		exit $EXIT_ERROR
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -d "${_arg_outdir}" ]]; then
			ls -l "${_arg_outdir}" >&3
		else
			_polap_log0 "No such output folder: ${_arg_outdir}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return
	fi

	mkdir -p "${_arg_outdir}"
	_polap_log0 "create output folder [${_arg_outdir}] if no such folder exists ..."
	if [ "${_arg_outdir}" != "o" ]; then
		_polap_log1 "  Use -o ${_arg_outdir} option in all subsequent analysis"
		_polap_log1 "  because your output folder is not the default of 'o'."
	fi
	# _run_polap_make-menus
	_log_command_versions

	_polap_log1 "NEXT: $(basename "$0") summary-reads -o ${_arg_outdir} -l ${_arg_long_reads}"
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# View the polap log file.
################################################################################
function _run_polap_log { # display the polap log
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Display the polap log.
#
# Arguments:
#   -o ${_arg_outdir}: the output folder
#   --log <FILE>
# Inputs:
#   ${LOG_FILE}
# Outputs:
#   a page view of the log file: ${LOG_FILE}
Example: $(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -d "${_arg_outdir}" ]]; then
			ls -l "${_arg_outdir}" >&3
		else
			_polap_log0 "No such output folder: ${_arg_outdir}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	less "${LOG_FILE}" >&3

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
