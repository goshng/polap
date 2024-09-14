################################################################################
# Compares the known mtDNA sequence and the assembled one.
################################################################################
function _run_polap_compare-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	local _polap_var_assembly="${ODIR}/${INUM}"
	local _polap_var_bioproject="${_polap_var_assembly}/70-bioproject"
	local _polap_var_bioproject_runinfo="${_polap_var_bioproject}/1-runinfo.tsv"
	local _polap_var_bioproject_sra_long_read="${_polap_var_bioproject}/1-sra-long-read.tsv"
	local _polap_var_bioproject_mtdna="${_polap_var_bioproject}/2-mtdna.fasta"
	local _polap_var_bioproject_blastn1="${_polap_var_bioproject}/3-blastn1.txt"
	local _polap_var_bioproject_blastn2="${_polap_var_bioproject}/3-blastn2.txt"
	local _polap_var_bioproject_blastn3="${_polap_var_bioproject}/3-blastn3.txt"
	local _polap_var_bioproject_blastn3_length="${_polap_var_bioproject}/3-blastn3.length.txt"

	# Paths for the mtDNA sequences to be compared
	local _polap_var_mtdna1="${_polap_var_assembly}/mt.1.fa"
	local _polap_var_mtdna2="${_polap_var_assembly}/mt.2.fa"
	local _polap_var_mtdna3="${_polap_var_assembly}/mt.3.fa"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Compares the known mtDNA sequence and the assembled one.
# Arguments:
#   -i $INUM
#   -f $FA
# Inputs:
#   $ODIR/$INUM/$FA
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [$ODIR]
HEREDOC
	)

	mkdir -p "${_polap_var_bioproject}"

	# Extract species name from bioproject's SRA long read
	# local SPECIES=$(cut -f4 "${_polap_var_bioproject_sra_long_read}")
	# echoerr "LOG: bioproject's species: $SPECIES"

	# Run blastn between known mtDNA and assembled mtDNA (first round)
	blastn -query "${_polap_var_mtdna1}" -subject "${_polap_var_bioproject_mtdna}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand" \
		>"${_polap_var_bioproject_blastn1}"

	blastn -query "${_polap_var_mtdna1}" -subject "${_polap_var_bioproject_mtdna}" \
		>"${_polap_var_bioproject_blastn1}.full"

	_polap_log2_file "${_polap_var_bioproject_blastn1}"

	# Determine the strand orientation
	local _polap_var_bioproject_strand=$(cut -f5 "${_polap_var_bioproject_blastn1}" | head -n 1)

	# Reverse sequence if the strand is negative
	if [ "${_polap_var_bioproject_strand}" = "plus" ]; then
		cp "${_polap_var_mtdna1}" "${_polap_var_mtdna2}"
	else
		seqkit seq -t dna -v -p -r "${_polap_var_mtdna1}" -o "${_polap_var_mtdna2}"
	fi

	# Run blastn (second round)
	blastn -query "${_polap_var_mtdna2}" -subject "${_polap_var_bioproject_mtdna}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand" \
		>"${_polap_var_bioproject_blastn2}"

	_polap_log2_file "${_polap_var_bioproject_blastn2}"

	# Restart sequence alignment at the lowest start position
	local _polap_var_bioproject_restart_position=$(sort -n -k4 "${_polap_var_bioproject_blastn2}" | head -n 1 | cut -f2)
	_polap_log2 "LOG: restart position: ${_polap_var_bioproject_restart_position}"
	seqkit restart -i "${_polap_var_bioproject_restart_position}" "${_polap_var_mtdna2}" -o "${_polap_var_mtdna3}"

	# Run blastn (third round)
	blastn -query "${_polap_var_mtdna3}" -subject "${_polap_var_bioproject_mtdna}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand length pident" \
		>"${_polap_var_bioproject_blastn3}"

	_polap_log2_file "${_polap_var_bioproject_blastn3}"

	# Analyze the length of the match using an R script
	"$WDIR"/run-polap-assemble-bioproject-3-length-match.R \
		"${_polap_var_bioproject_blastn3}" \
		"${_polap_var_bioproject_blastn3_length}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_bioproject_blastn3_length}"

	local l1=$(seqkit stats -T "${_polap_var_bioproject_mtdna}" |
		csvtk cut -t -f sum_len |
		csvtk del-header)

	local l2=$(cat "${_polap_var_bioproject_blastn3_length}")
	local c1=$(echo "scale=3; ${l2}/${l1}" | bc)
	_polap_log2 "Length of ${_polap_var_bioproject_mtdna}: ${l1}"
	_polap_log2 "Length of match alignment: ${l2}"
	_polap_log1 "length coverage: ${c1}"

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
