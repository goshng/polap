################################################################################
#
################################################################################
function _run_polap_get-bioproject-sra() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# the organelle-genome assembly.
# Arguments:
#   -b $SR2: bioproject ID
# Inputs:
#   SRA ID
# Outputs:
#   bioproject ID
Example: $(basename $0) ${_arg_menu[0]} --sra SRR11472010
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	esearch -db sra -query "${_arg_sra}" |
		efetch -format runinfo |
		csvtk cut -f BioProject |
		csvtk del-header \
			>"${ODIR}/bioproject.txt"

	_polap_log1 $(cat "${ODIR}/bioproject.txt")

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
