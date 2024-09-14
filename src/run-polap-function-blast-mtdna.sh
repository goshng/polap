################################################################################
# Blast the final mtDNA sequence against mitochondrial and plastid genes.
################################################################################
function _run_polap_blast-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	local MTAA="$WDIR/polap-mt.1.c70.3.faa"
	local PTAA="$WDIR/polap-pt.2.c70.3.faa"
	local _polap_var_assembly="${ODIR}/${INUM}"
	local _polap_var_chloroplot="${_polap_var_assembly}/60-chloroplot"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Blast the final mtDNA sequence against the POLAP organelle gene set.
# Arguments:
#   -f $FA
# Inputs:
#   $FA: mtDNA sequence
# Outputs:
#
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [$ODIR]
HEREDOC
	)

	# Display help message if needed
	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	# Check if the input FASTA file exists
	check_file_existence "${FA}"

	# Clean up and create the chloroplot directory
	if [ -d "${_polap_var_chloroplot}" ]; then
		rm -rf "${_polap_var_chloroplot}"
		_polap_log2 "LOG: ${_polap_var_chloroplot} is deleted."
	fi
	mkdir -p "${_polap_var_chloroplot}"
	_polap_log2 "LOG: ${_polap_var_chloroplot} is created."

	# Create the BLAST database for the mtDNA sequences
	makeblastdb -dbtype nucl \
		-in "${FA}" \
		-out "${_polap_var_chloroplot}/dna" \
		>/dev/null 2>&1
	_polap_log2 "LOG: BLASTDB of the contig sequences: ${_polap_var_chloroplot}/dna"

	# BLAST mtDNA gene annotation
	_polap_log1 "LOG: executing tblastn ... be patient!"
	_polap_log2 "LOG: BLAST of the mitochondrial proteins against ${_polap_var_chloroplot}/dna ..."

	tblastn -query "$MTAA" \
		-db "${_polap_var_chloroplot}/dna" \
		-out "${_polap_var_chloroplot}/mtaa.blast" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "$NT" \
		>/dev/null 2>&1

	_polap_log2 "INFO: tblastn complete"

	# Process the tblastn results
	"$WDIR"/run-polap-genes.R \
		"${_polap_var_chloroplot}/mtaa.blast" \
		"${_polap_var_chloroplot}/mtaa.blast.bed" \
		>/dev/null 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_chloroplot}/mtaa.blast.bed" \
		>"${_polap_var_chloroplot}/mtaa.blast.sorted.bed"

	# Create directory for gene bed files
	mkdir "${_polap_var_chloroplot}/mtaa.bed"

	"$WDIR"/run-polap-genes-bed4.R \
		"${_polap_var_chloroplot}/mtaa.blast" \
		"${_polap_var_chloroplot}/mtaa.blast.bed4"

	_polap_log2 "LOG: counting mitochondrial genes in the contigs ..."
	bedtools merge -i "${_polap_var_chloroplot}/mtaa.blast.sorted.bed" >"${_polap_var_chloroplot}/mtaa.blast.sorted.bed.txt"

	_polap_log2_file "${_polap_var_chloroplot}/mtaa.blast.sorted.bed.txt"

	# Process individual genes and their annotations
	local _polap_var_i=1
	while IFS= read -r gene; do
		printf "%s\n" "$gene" >"${_polap_var_chloroplot}/${_polap_var_i}.gene.bed"

		bedtools intersect -wa \
			-a "${_polap_var_chloroplot}/mtaa.blast.bed4" \
			-b "${_polap_var_chloroplot}/${_polap_var_i}.gene.bed" \
			>"${_polap_var_chloroplot}/${_polap_var_i}.bed4"

		"$WDIR"/run-polap-blast-mtdna-1-determine-gene.R \
			"${_polap_var_chloroplot}/${_polap_var_i}.bed4" \
			"$WDIR/polap-mt.1.c70.3.faa.name" \
			"$WDIR/polap-mtgenes.txt" \
			"${_polap_var_chloroplot}/${_polap_var_i}.gene" \
			"${_polap_var_chloroplot}/${_polap_var_i}.bed4.description" \
			"${_polap_var_chloroplot}/${_polap_var_i}.bed4.count" \
			>/dev/null 2>&1

		_polap_var_i=$((_polap_var_i + 1))

	done <"${_polap_var_chloroplot}/mtaa.blast.sorted.bed.txt"

	# Combine the gene annotations into one file
	paste \
		<(cat "${_polap_var_chloroplot}"/*.gene.bed) \
		<(cat "${_polap_var_chloroplot}"/*.gene) \
		>"${_polap_var_chloroplot}/annotation.bed"

	_polap_log1_file "Annotation: ${_polap_var_chloroplot}/annotation.bed"
	_polap_log1_file "Ambiguous gene annotations: ${_polap_var_chloroplot}/*.check"
	ls "${_polap_var_chloroplot}"/*.check \
		>"$_polap_output_dest"

	echoerr "NEXT: $(basename $0) gene-table-mtdna -o $ODIR [-i $INUM]"

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
