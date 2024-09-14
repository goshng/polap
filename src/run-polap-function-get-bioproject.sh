################################################################################
#
################################################################################
function _run_polap_get-bioproject() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	local _polap_var_assembly="${ODIR}/${INUM}"
	local _polap_var_bioproject="${_polap_var_assembly}/70-bioproject"
	local _polap_var_bioproject_runinfo="${_polap_var_bioproject}/1-runinfo.tsv"
	local _polap_var_bioproject_sra_long_read="${_polap_var_bioproject}/1-sra-long-read.tsv"
	local _polap_var_bioproject_sra_short_read="${_polap_var_bioproject}/1-sra-short-read.tsv"
	local _polap_var_bioproject_species="${_polap_var_bioproject}/1-species.txt"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# the organelle-genome assembly.
# Arguments:
#   -o $ODIR
#   -i $INUM
#   -b $SR2: bioproject ID
# Inputs:
#   bioproject ID
# Outputs:
#   ${_polap_var_bioproject_runinfo}
Example: $(basename $0) ${_arg_menu[0]} -b PRJNA574453
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	# Clean and create the bioproject directory
	rm -rf "${_polap_var_bioproject}"
	mkdir -p "${_polap_var_bioproject}"

	# Set the BioProject ID
	local BIOPRJ="$SR2"

	# Fetch the run information for the given BioProject
	esearch -db bioproject -query "$BIOPRJ" |
		elink -target sra |
		efetch -format runinfo |
		csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName |
		csvtk csv2tab >"${_polap_var_bioproject_runinfo}"

	# Log the output file
	_polap_log1_file "downloaded: ${_polap_var_bioproject_runinfo}"

	# Select SRA long and short read data
	"$script_dir/run-polap-assemble-bioproject-1-select-sra.R" \
		"${_polap_var_bioproject_runinfo}" \
		"${_polap_var_bioproject_sra_long_read}" \
		"${_polap_var_bioproject_sra_short_read}" \
		"${_polap_var_bioproject_species}" \
		2>"$_polap_output_dest"

	# Log the long and short read SRA files
	_polap_log1_file "output1: ${_polap_var_bioproject_sra_long_read}"
	_polap_log1_file "output2: ${_polap_var_bioproject_sra_short_read}"
	_polap_log1_file "output3: ${_polap_var_bioproject_species}"

	# Check if the long-read dataset exists
	if [ ! -s "${_polap_var_bioproject_sra_long_read}" ]; then
		die "ERROR: no long-read dataset for the BioProject: $BIOPRJ"
	fi

	# Check if the short-read dataset exists
	if [ ! -s "${_polap_var_bioproject_sra_short_read}" ]; then
		die "ERROR: no short-read dataset for the BioProject: $BIOPRJ"
	fi

	if [ ! -s "${_polap_var_bioproject_species}" ]; then
		die "ERROR: no species name for the BioProject: $BIOPRJ"
	fi

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
