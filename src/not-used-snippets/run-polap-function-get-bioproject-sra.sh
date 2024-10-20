################################################################################
# Gets BioProject ID of a given SRA ID.
#
# Arguments:
#   --sra ${_arg_sra}: SRA ID
#   -o ${ODIR}: output folder
# Inputs:
#   SRA ID
# Outputs:
#   BioProject ID
#   FILE: ${_polap_var_bioproject_txt}
################################################################################
function _run_polap_get-bioproject-sra() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-bioproject.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Gets BioProject ID of a given SRA ID.
#
# Arguments:
#   --sra ${_arg_sra}: SRA ID
#   -o ${ODIR}: output folder
# Inputs:
#   SRA ID
# Outputs:
#   BioProject ID
#   FILE: ${_polap_var_bioproject_txt}
Example: $(basename $0) ${_arg_menu[0]} --sra <SRA ID>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if [[ ${_arg_menu[1]} == "view" ]]; then
		if [[ -s "${_polap_var_bioproject_txt}" ]]; then
			_polap_log0_cat "${_polap_var_bioproject_txt}"
		else
			_polap_log0 "No such file: ${_polap_var_bioproject_txt}"
		fi
		exit $EXIT_SUCCESS
	fi

	if [[ ! -d "${ODIR}" ]]; then
		_polap_log1 "  no output folder, creating ${ODIR}"
		mkdir -p "${ODIR}"
	fi

	esearch -db sra -query "${_arg_sra}" |
		efetch -format runinfo |
		csvtk cut -f BioProject |
		csvtk del-header \
			>"${_polap_var_bioproject_txt}"

	_polap_log1_file "${ODIR}/bioproject.txt"
	_polap_log0 $(head -n 1 "${ODIR}/bioproject.txt")

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
