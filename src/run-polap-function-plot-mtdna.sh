################################################################################
# Gene table for importing to the Chloroplot R package.
################################################################################
function _run_polap_plot-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Set working directories and file paths
	local _polap_var_assembly="${ODIR}/${INUM}"
	local _polap_var_chloroplot="${_polap_var_assembly}/60-chloroplot"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Counts genes annotated on a genome assembly and plots the mtDNA genome.
# Arguments:
#   -f $FA
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   $FA: mtDNA sequence
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] [-f $FA]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	_polap_log1 "LOG: Plotting mitochondrial DNA genome for Flye assembly $INUM..."

	# Run the R script to generate the mtDNA plot
	"$WDIR/run-polap-plot-mtdna.R" \
		"${_polap_var_chloroplot}/annotation.bed" \
		"${FA}" \
		2>"$_polap_output_dest"

	# Output file information
	echoerr "FILE: mt.3.pdf has been created."

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
