################################################################################
# Downloads the mtDNA sequence for a given species name.
################################################################################
function _run_polap_get-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Download the organelle-genome in FASTA format.
# Arguments:
#   --species "scientific name" highest priority
#   or
#   -o ${ODIR} next priority
#   or
#   --bioproject ${_arg_bioproject} the least priority
# Outputs:
#   $ODIR/bioproject/1-mtdna.fasta
Example: $(basename $0) ${_arg_menu[0]} --species "Anthoceros agrestis"
Example: $(basename $0) ${_arg_menu[0]} -o o
Example: $(basename $0) ${_arg_menu[0]} -b PRJNA574453
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if [[ ${_arg_menu[1]} == "view" ]]; then
		if [[ -s "${_polap_var_bioproject_mtdna_fasta2_accession}" ]]; then
			_polap_log0_cat "${_polap_var_bioproject_mtdna_fasta2_accession}"
			if [ -s "${_polap_var_bioproject_mtdna_fasta1}" ]; then
				seqkit stats "${_polap_var_bioproject_mtdna_fasta1}" >&2
				seqkit stats "${_polap_var_bioproject_mtdna_fasta2}" >&2
			fi
		else
			_polap_log0 "No result yet."
		fi
		exit $EXIT_SUCCESS
	fi

	# Ensure the bioproject directory exists
	# mkdir -p "${_polap_var_bioproject}"
	# _run_polap_get-bioproject

	# Determine species name
	local SPECIES=""
	if [ -n "${_arg_species}" ]; then
		_polap_log1 "LOG: species name is given as the option --species: $SPECIES"
		SPECIES="${_arg_species}"
		mkdir -p "${_polap_var_bioproject}"
	elif [ -s "${_polap_var_bioproject_species}" ]; then
		SPECIES=$(cut -f4 "${_polap_var_bioproject_species}")
		_polap_log1 "LOG: bioproject's species: $SPECIES"
	else
		_run_polap_get-bioproject
		if [ -s "${_polap_var_bioproject_species}" ]; then
			SPECIES=$(cut -f4 "${_polap_var_bioproject_species}")
			_polap_log1 "LOG: bioproject's species: $SPECIES"
		else
			die "No species name is provided."
		fi
	fi

	# Download the mitochondrial genome sequence for the given species
	if [ "${_arg_plastid}" = "off" ]; then
		esearch \
			-db nuccore \
			-query "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
			efetch -format fasta >"${_polap_var_bioproject_mtdna_fasta1}"
	else
		_polap_log1 "LOG: downloading chloroplast complete genomes of ${SPECIES} ..."
		esearch \
			-db nuccore \
			-query "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
			efetch -format fasta >"${_polap_var_bioproject_mtdna_fasta1}"
	fi

	# Check if the fasta file was successfully downloaded
	if [ -s "${_polap_var_bioproject_mtdna_fasta1}" ]; then
		_polap_log2_file "potentially multiple-sequence: ${_polap_var_bioproject_mtdna_fasta1}"
		seqkit stats -T "${_polap_var_bioproject_mtdna_fasta1}" \
			>"${_polap_var_bioproject_mtdna_fasta1_stats}"

		local n=$(seqkit stats -T "${_polap_var_bioproject_mtdna_fasta1}" |
			csvtk cut -t -f num_seqs |
			csvtk del-header)

		seqkit head -n 1 "${_polap_var_bioproject_mtdna_fasta1}" \
			-o "${_polap_var_bioproject_mtdna_fasta2}"

		seqkit seq -ni "${_polap_var_bioproject_mtdna_fasta2}" \
			-o "${_polap_var_bioproject_mtdna_fasta2_accession}"

		_polap_log1_file "mtDNA NCBI accession: ${_polap_var_bioproject_mtdna_fasta2_accession}"
	else
		echo "no mtDNA" >"${_polap_var_bioproject_mtdna_fasta2_accession}"
		_polap_log0 "No mtDNA sequence found for the species: ${SPECIES}"
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
