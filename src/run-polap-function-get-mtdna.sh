################################################################################
# Downloads the mtDNA sequence for a given species name.
################################################################################
function _run_polap_get-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Download the organelle-genome in FASTA format.
# Arguments:
#   -o $ODIR
#   -b $SR2: bioproject ID
#   --species "scientific name"
# Outputs:
#   $ODIR/bioproject/1-mtdna.fasta
Example: $(basename $0) ${_arg_menu[0]} -b PRJNA574453
Example: $(basename $0) ${_arg_menu[0]} --species "Anthoceros agrestis"
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	# Set variables for file paths
	local _polap_var_assembly="${ODIR}/${INUM}"
	local _polap_var_bioproject="${_polap_var_assembly}/70-bioproject"
	local _polap_var_bioproject_species="${_polap_var_bioproject}/1-species.txt"
	local _polap_var_bioproject_mtdna_fasta1="${_polap_var_bioproject}/1-mtdna.fasta"
	local _polap_var_bioproject_mtdna_fasta1_stats="${_polap_var_bioproject}/1-mtdna.fasta.stats"
	local _polap_var_bioproject_mtdna_fasta2="${_polap_var_bioproject}/2-mtdna.fasta"

	# Ensure the bioproject directory exists
	mkdir -p "${_polap_var_bioproject}"
	# _run_polap_get-bioproject

	# Determine species name
	local SPECIES=""
	if [ -n "${_arg_species}" ]; then
		SPECIES="${_arg_species}"
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
		_polap_log1_file "potentially multiple-sequence: ${_polap_var_bioproject_mtdna_fasta1}"
		seqkit stats -T "${_polap_var_bioproject_mtdna_fasta1}" \
			>"${_polap_var_bioproject_mtdna_fasta1_stats}"

		local n=$(seqkit stats -T "${_polap_var_bioproject_mtdna_fasta1}" |
			csvtk cut -t -f num_seqs |
			csvtk del-header)

		seqkit head -n 1 "${_polap_var_bioproject_mtdna_fasta1}" \
			-o "${_polap_var_bioproject_mtdna_fasta2}"

		_polap_log1_file "single-sequence fasta: ${_polap_var_bioproject_mtdna_fasta2}"
	else
		die "No mtDNA sequence found for the species: ${SPECIES}"
	fi

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
