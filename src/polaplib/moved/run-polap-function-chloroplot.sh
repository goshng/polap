################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# These functions are moved from polap-cmd-mtdna.sh.
# function _run_polap_plot-mtdna {
# function _run_polap_blast-mtdna {
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

################################################################################
# BEGIN OF THE MOVE

################################################################################
# Blast the final mtDNA sequence against mitochondrial and plastid genes.
################################################################################
function _run_polap_blast-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Set variables for file paths
	local MTAA="${_POLAPLIB_DIR}/polap-mt.1.c70.3.faa"
	local PTAA="${_POLAPLIB_DIR}/polap-pt.2.c70.3.faa"
	local _polap_var_assembly="${_polap_var_ga}"
	local _polap_var_chloroplot="${_polap_var_assembly}/60-chloroplot"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Blast the final mtDNA sequence against the POLAP organelle gene set.
# Arguments:
#   -f ${_arg_final_assembly}
# Inputs:
#   ${_arg_final_assembly}: mtDNA sequence
# Outputs:
#
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [${_arg_outdir}]
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
		-num_threads "${_arg_threads}" \
		>/dev/null 2>&1

	_polap_log2 "INFO: tblastn complete"

	# Process the tblastn results
	"${_POLAPLIB_DIR}"/polap-r-genes.R \
		"${_polap_var_chloroplot}/mtaa.blast" \
		"${_polap_var_chloroplot}/mtaa.blast.bed" \
		>/dev/null 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_chloroplot}/mtaa.blast.bed" \
		>"${_polap_var_chloroplot}/mtaa.blast.sorted.bed"

	# Create directory for gene bed files
	mkdir "${_polap_var_chloroplot}/mtaa.bed"

	"${_POLAPLIB_DIR}"/run-polap-r-genes-bed4.R \
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

		"${_POLAPLIB_DIR}"/run-polap-r-blast-mtdna-1-determine-gene.R \
			"${_polap_var_chloroplot}/${_polap_var_i}.bed4" \
			"${_POLAPLIB_DIR}/polap-mt.1.c70.3.faa.name" \
			"${_POLAPLIB_DIR}/polap-mtgenes.txt" \
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

	echoerr "NEXT: $(basename $0) gene-table-mtdna -o ${_arg_outdir} [-i ${_arg_inum}]"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Gene table for importing to the Chloroplot R package.
################################################################################
function _run_polap_plot-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Set working directories and file paths
	local _polap_var_assembly="${_polap_var_ga}"
	local _polap_var_chloroplot="${_polap_var_assembly}/60-chloroplot"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Counts genes annotated on a genome assembly and plots the mtDNA genome.
# Arguments:
#   -f ${_arg_final_assembly}
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_arg_final_assembly}: mtDNA sequence
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] [-f ${_arg_final_assembly}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log1 "LOG: Plotting mitochondrial DNA genome for Flye assembly ${_arg_inum}..."

	# Run the R script to generate the mtDNA plot
	"${_POLAPLIB_DIR}/run-polap-r-plot-mtdna.R " \
		"${_polap_var_chloroplot}/annotation.bed" \
		"${FA}" \
		2>"$_polap_output_dest"

	# Output file information
	_polap_log0 "FILE: mt.3.pdf has been created."

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# END OF THE MOVE
################################################################################
