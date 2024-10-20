################################################################################
# Reports the organelle-genome assembly results.
################################################################################
function _run_polap_report-assembly() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Reports the organelle-genome assembly results.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject
# Inputs:
#   ${ODIR}: output folder for BioProject
# Outputs:
Example: $(basename $0) ${_arg_menu[0]} [-o ${ODIR}]
Example: report-assembly -o PRJDB10540a 2>&1 | tr '\n' '\t' | sed 's/\t$/\n/'
HEREDOC
	)

	LRNK="${ODIR}/nk.fq.gz"

	# Set variables for file paths
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'
	source "$script_dir/polap-variables-base.sh"       # '.' means 'source'

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_log0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log1_log "reporting the organelle-genome assembly at ${ODIR} ..."

	if [ -d "${ODIR}" ]; then
		_polap_log2_file "the main output folder: ${ODIR}"
	else
		_polap_log2 "ERROR: no such output folder; use -o option"
		exit $EXIT_SUCCESS
	fi

	_polap_log0 $(cut -f1 "${_polap_var_bioproject_txt}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_sra_long_read}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_sra_short_read}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_species}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")

	# wc -l "${ODIR}/0/mt.contig.name-"? | awk '$2 != "total" {print $1}' | head -5 >&2

	# for i in "${_arg_select_contig_numbers[@]}"; do
	# 	# Call the function corresponding to the current number (index is i-1)
	# 	INUM="${i}"
	# 	source "$script_dir/polap-variables-oga.sh" # '.' means 'source'
	# 	_polap_log0 $(cat "${ODIR}/0/mt.contig.name-${i}" | wc -l)
	# 	_polap_log0_cat "${_polap_var_mtdna_compare}"
	# done

	# Array to store the names of the original files
	files=($(ls "${ODIR}/0/mt.contig.name-"?))

	# Temporary array to store the paths of unique files
	unique_files=()

	# Function to check if a file is unique
	is_unique() {
		for unique_file in "${unique_files[@]}"; do
			if cmp -s "$1" "$unique_file"; then
				echo "$unique_file"
				return 1 # Not unique
			fi
		done
		return 0 # Unique
	}

	_polap_log1 "Checking for unique files and their matches:"

	# Iterate over the files to find unique ones
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Call the function corresponding to the current number (index is i-1)
		FDIR="${ODIR}/0"
		JNUM="${i}"
		file="$FDIR"/mt.contig.name-$JNUM

		unique_file=$(is_unique "$file")
		if [ $? -eq 0 ]; then
			# If unique, add it to the unique_files array
			unique_files+=("$file")
			echo "$file is unique."

			MTCONTIGNAME="$file"
			INUM="${i}"
		else
			_polap_log1 "$file is the same as $unique_file."
			INUM="${unique_file##*-}"
		fi
		source "$script_dir/polap-variables-oga.sh" # '.' means 'source'
		_polap_log0 $(cat "${ODIR}/0/mt.contig.name-${i}" | wc -l)
		_polap_log0_cat "${_polap_var_mtdna_compare}"
	done

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
