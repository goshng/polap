################################################################################
# Compares the known mtDNA sequence and the assembled one.
################################################################################
function _run_polap_compare-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'
	source "$script_dir/polap-variables-oga.sh"        # '.' means 'source'

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Compares the known mtDNA sequence and the assembled one.
# Arguments:
#   -i $INUM
# Inputs:
#   ${_polap_var_mtdna3}
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [$ODIR]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	check_folder_existence "${ODIR}"

	local n1=$(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")
	local l1=$(seqkit stats -T "${_polap_var_bioproject_mtdna_fasta2}" |
		csvtk cut -t -f sum_len |
		csvtk del-header)

	# Run blastn between known mtDNA and assembled mtDNA (first round)
	blastn -query "${_polap_var_mtdna1}" \
		-subject "${_polap_var_bioproject_mtdna_fasta2}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand" \
		>"${_polap_var_oga_blastn1}"

	blastn -query "${_polap_var_mtdna1}" \
		-subject "${_polap_var_bioproject_mtdna_fasta2}" \
		>"${_polap_var_oga_blastn1}.full"

	_polap_log2_file "${_polap_var_oga_blastn1}"

	# Determine the strand orientation
	if [ -s "${_polap_var_oga_blastn1}" ]; then
		local _polap_var_bioproject_strand=$(cut -f5 "${_polap_var_oga_blastn1}" | head -n 1)
	else
		_polap_log1 "No hit in ${_polap_var_oga_blastn1}"
		local n1=$(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")
		printf "%s\t%d\t0\t0\n" ${n1} ${l1} >"${_polap_var_mtdna_compare}"
		return
	fi

	# Reverse sequence if the strand is negative
	if [ "${_polap_var_bioproject_strand}" = "plus" ]; then
		cp "${_polap_var_mtdna1}" "${_polap_var_mtdna2}"
	else
		seqkit seq -t dna -v -p -r "${_polap_var_mtdna1}" -o "${_polap_var_mtdna2}"
	fi

	# Run blastn (second round)
	blastn -query "${_polap_var_mtdna2}" -subject "${_polap_var_bioproject_mtdna_fasta2}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand" \
		>"${_polap_var_oga_blastn2}"

	_polap_log2_file "${_polap_var_oga_blastn2}"

	# Restart sequence alignment at the lowest start position
	local _polap_var_bioproject_restart_position=$(sort -n -k4 "${_polap_var_oga_blastn2}" | head -n 1 | cut -f2)
	_polap_log2 "LOG: restart position: ${_polap_var_bioproject_restart_position}"
	seqkit restart -i "${_polap_var_bioproject_restart_position}" "${_polap_var_mtdna2}" -o "${_polap_var_mtdna3}"

	# Run blastn (third round)
	blastn -query "${_polap_var_mtdna3}" -subject "${_polap_var_bioproject_mtdna_fasta2}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand length pident" \
		>"${_polap_var_oga_blastn3}"

	_polap_log2_file "${_polap_var_oga_blastn3}"

	# Analyze the length of the match using an R script
	"$WDIR"/run-polap-assemble-bioproject-3-length-match.R \
		"${_polap_var_oga_blastn3}" \
		"${_polap_var_oga_blastn3_length}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_oga_blastn3_length}"

	local l2=$(<"${_polap_var_oga_blastn3_length}")
	local c1=$(echo "scale=3; ${l2}/${l1}" | bc)
	_polap_log2 "Length of ${_polap_var_bioproject_mtdna_fasta2}: ${l1}"
	_polap_log2 "Length of match alignment: ${l2}"
	_polap_log1 "length coverage: ${c1}"
	printf "%s\t%d\t%d\t%f\n" ${n1} ${l1} ${l2} ${c1} >"${_polap_var_mtdna_compare}"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
