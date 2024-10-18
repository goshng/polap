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
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Blasts the genome assembly of a Flye run.
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   $assembly_contigs_fasta
# Outputs:
#   $CONTIGNAME
#   $MTGENECOUNT
#   $PTGENECOUNT
#   $ADIR/mtaa.blast
#   $ADIR/mtaa.blast.bed
#   $MTAABLAST.sorted.bed
#   $ADIR/mtaa.bed
#   $ADIR/ptaa.blast
#   $ADIR/ptaa.blast.bed
#   $PTAABLAST.sorted.bed
#   $ADIR/ptaa.bed
################################################################################
function _run_polap_blast-genome() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	# source "$script_dir/polap-variables-annotation.sh"
	source "$script_dir/polap-variables-ga.sh"
	local MTAA="$script_dir"/polap-mt.1.c70.3.faa
	local PTAA="$script_dir"/polap-pt.2.c70.3.faa

	help_message=$(
		cat <<HEREDOC
# Blasts the genome assembly of a Flye run againt the plant organelle genes.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_contigger_contigs_stats}
# Outputs:
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "blasting contigs with mitochondrial and plastid genes on the assembly number ${INUM} ..."

	# Checks the output files earlier than the input.
	if [[ -s "${_polap_var_ann_MTGENECOUNT}" ]] &&
		[[ -s "${_polap_var_ann_PTGENECOUNT}" ]] &&
		[ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found1: ${_polap_var_ann_MTGENECOUNT}"
		_polap_log0 "  found2: ${_polap_var_ann_PTGENECOUNT}"
		_polap_log0 "  so skipping the blast genome ..."
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	if [ ! -s "${_polap_var_contigger_contigs_stats}" ]; then
		_polap_log0 "ERROR: no contig_stats.txt file: ${_polap_var_contigger_contigs_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "  deleting and recreating the folder: ${_polap_var_ann}"
	if [ -d "${_polap_var_ann}" ]; then
		_polap_log3 rm -rf "${_polap_var_ann}"
		rm -rf "${_polap_var_ann}"
	fi
	_polap_log3 mkdir -p "${_polap_var_ann}"
	mkdir -p "${_polap_var_ann}"

	_polap_log1 "  input1: ${_polap_var_contigger_contigs_stats}"
	_polap_log1 "  input2: ${_polap_var_contigger_contigs_fasta}"

	#src/run-polap-select.R o/30-contigger/contigs_stats.txt o/50-annotation/contig.name
	_polap_log2 "  contig sequence names in file: ${_polap_var_ann_CONTIGNAME}"
	grep -v "#" "${_polap_var_contigger_contigs_stats}" |
		cut -f 1 >"${_polap_var_ann_CONTIGNAME}"

	# seqkit grep --threads $NT -f "$CONTIGNAME" \
	# 	"$assembly_contigs_fasta" \
	# 	-o "${_polap_var_ann}"/contig.fasta \
	# 	>/dev/null 2>&1
	_polap_log2 "  contig sequence file: ${_polap_var_ann_CONTIGFILE}"
	cp "${_polap_var_contigger_contigs_fasta}" "${_polap_var_ann_CONTIGFILE}"

	_polap_log2 "  making BLASTDB of the contig sequences: ${_polap_var_ann_CONTIGDB}"
	_polap_log3 makeblastdb -dbtype nucl \
		-in "${_polap_var_contigger_contigs_fasta}" \
		-out "${_polap_var_ann_CONTIGDB}" \
		>${_polap_output_dest} 2>&1
	makeblastdb -dbtype nucl \
		-in "${_polap_var_contigger_contigs_fasta}" \
		-out "${_polap_var_ann_CONTIGDB}" \
		>${_polap_output_dest} 2>&1

	# Mitochondrial gene annotation and counts
	_polap_log2 "  BLAST of the mitochondrial proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "    executing the tblastn ... be patient!"
	tblastn -query "${MTAA}" \
		-db "${_polap_var_ann_CONTIGDB}" \
		-out "${_polap_var_ann_MTAABLAST}" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "${NT}" \
		>${_polap_output_dest} 2>&1

	Rscript "$script_dir"/run-polap-genes.R \
		"${_polap_var_ann_MTAABLAST}" \
		"${_polap_var_ann_MTAABLASTBED}" \
		>${_polap_output_dest} 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_ann_MTAABLASTBED}" \
		>"${_polap_var_ann_MTAABLAST}".sorted.bed

	mkdir "${_polap_var_ann_MTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w "${contig}" "${_polap_var_ann_MTAABLAST}".sorted.bed \
			>"${_polap_var_ann_MTAABED}/${contig}".bed
		bedtools merge -i "${_polap_var_ann_MTAABED}/${contig}".bed \
			>"${_polap_var_ann_MTAABED}/${contig}".bed.txt
		printf "%s\t%d\n" "${contig}" \
			$(wc -l <"${_polap_var_ann_MTAABED}/${contig}".bed.txt)
	done <"${_polap_var_ann_CONTIGNAME}" |
		sort -k2 -rn >"${_polap_var_ann_MTGENECOUNT}"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	tar zcf "${_polap_var_ann_MTAABED}".tar.gz "${_polap_var_ann_MTAABED}"
	rm -rf "${_polap_var_ann_MTAABED}"

	_polap_log1 "  output1: ${_polap_var_ann_MTGENECOUNT}"

	# Plastid gene annotation and counts
	_polap_log2 "  BLAST of the plastid proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "    executing the tblastn ... be patient!"
	tblastn -query "${PTAA}" \
		-db "${_polap_var_ann_CONTIGDB}" \
		-out "${_polap_var_ann_PTAABLAST}" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "${NT}" \
		>${_polap_output_dest} 2>&1

	Rscript "$script_dir"/run-polap-genes.R \
		"${_polap_var_ann_PTAABLAST}" \
		"${_polap_var_ann_PTAABLASTBED}" \
		>${_polap_output_dest} 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_ann_PTAABLASTBED}" \
		>"${_polap_var_ann_PTAABLAST}".sorted.bed

	mkdir "${_polap_var_ann_PTAABED}"

	_polap_log2 "  counting plastid genes in the contigs ..."
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w "${contig}" "${_polap_var_ann_PTAABLAST}".sorted.bed \
			>"${_polap_var_ann_PTAABED}/${contig}".bed
		bedtools merge -i "${_polap_var_ann_PTAABED}/${contig}".bed \
			>"${_polap_var_ann_PTAABED}/${contig}".bed.txt
		printf "%s\t%d\n" "${contig}" \
			$(wc -l <"${_polap_var_ann_PTAABED}/${contig}".bed.txt)
	done <"${_polap_var_ann_CONTIGNAME}" |
		sort -k2 -rn >"${_polap_var_ann_PTGENECOUNT}"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	_polap_log2 "  compressing the BLAST results of plastid gene annotation"
	tar zcf "${_polap_var_ann_PTAABED}".tar.gz "${_polap_var_ann_PTAABED}"
	rm -rf "${_polap_var_ann_PTAABED}"

	_polap_log1 "  output2: ${_polap_var_ann_PTGENECOUNT}"

	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") count-gene --test"
	_polap_log1 "NEXT: $(basename $0) count-gene -o $ODIR [-i $INUM]"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Blasts the genome assembly of a Flye run.
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   graph_final.fasta not contigs.fasta
#
#   $assembly_contigs_fasta
# Outputs:
#   $CONTIGNAME
#   $MTGENECOUNT
#   $PTGENECOUNT
#   $ADIR/mtaa.blast
#   $ADIR/mtaa.blast.bed
#   $MTAABLAST.sorted.bed
#   $ADIR/mtaa.bed
#   $ADIR/ptaa.blast
#   $ADIR/ptaa.blast.bed
#   $PTAABLAST.sorted.bed
#   $ADIR/ptaa.bed
################################################################################
function _run_polap_blast-genome-edge() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	# source "$script_dir/polap-variables-annotation.sh"
	source "$script_dir/polap-variables-ga.sh"
	local MTAA="$script_dir"/polap-mt.1.c70.3.faa
	local PTAA="$script_dir"/polap-pt.2.c70.3.faa

	help_message=$(
		cat <<HEREDOC
# Blasts the genome assembly of a Flye run againt the plant organelle genes.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_contigger_edges_fasta}
#   ${_polap_var_contigger_edges_stats}
# Outputs:
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	_polap_log0 "blasting edge contigs with mitochondrial and plastid genes on the assembly number ${INUM} ..."

	# Checks the output files earlier than the input.
	if [[ -s "${_polap_var_ann_MTGENECOUNT}" ]] &&
		[[ -s "${_polap_var_ann_PTGENECOUNT}" ]] &&
		[ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found1: ${_polap_var_ann_MTGENECOUNT}"
		_polap_log0 "  found2: ${_polap_var_ann_PTGENECOUNT}"
		_polap_log0 "  so skipping the blast genome ..."
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	if [ ! -s "${_polap_var_contigger_edges_stats}" ]; then
		_polap_log0 "ERROR: no edges_stats.txt file: ${_polap_var_contigger_edges_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "  deleting and recreating the folder: ${_polap_var_ann}"
	_polap_log3_cmd rm -rf "${_polap_var_ann}"
	_polap_log3_cmd mkdir -p "${_polap_var_ann}"

	_polap_log1 "  input1: ${_polap_var_contigger_edges_stats}"
	_polap_log1 "  input2: ${_polap_var_contigger_edges_fasta}"

	#src/run-polap-select.R o/30-contigger/edges_stats.txt o/50-annotation/contig.name
	_polap_log2 "  contig sequence names in file: ${_polap_var_ann_CONTIGNAME}"
	grep -v "#" "${_polap_var_contigger_edges_stats}" |
		cut -f 1 >"${_polap_var_ann_CONTIGNAME}"

	# seqkit grep --threads $NT -f "$CONTIGNAME" \
	# 	"$assembly_edges_fasta" \
	# 	-o "${_polap_var_ann}"/contig.fasta \
	# 	>/dev/null 2>&1
	_polap_log2 "  contig sequence file: ${_polap_var_ann_CONTIGFILE}"
	cp "${_polap_var_contigger_edges_fasta}" "${_polap_var_ann_CONTIGFILE}"

	_polap_log2 "  making BLASTDB of the contig sequences: ${_polap_var_ann_CONTIGDB}"
	_polap_log3_pipe "makeblastdb -dbtype nucl \
		-in ${_polap_var_ann_CONTIGFILE} \
		-out ${_polap_var_ann_CONTIGDB} \
		2>${_polap_output_dest}"
	# makeblastdb -dbtype nucl \
	# 	-in "$_polap_var_ann_CONTIGFILE}" \
	# 	-out "${_polap_var_ann_CONTIGDB}" \
	# 	>${_polap_output_dest} 2>&1

	# Mitochondrial gene annotation and counts
	_polap_log2 "  BLAST of the mitochondrial proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "    executing the tblastn ... be patient!"
	tblastn -query "${MTAA}" \
		-db "${_polap_var_ann_CONTIGDB}" \
		-out "${_polap_var_ann_MTAABLAST}" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "${NT}" \
		>${_polap_output_dest} 2>&1

	Rscript "$script_dir"/run-polap-genes.R \
		"${_polap_var_ann_MTAABLAST}" \
		"${_polap_var_ann_MTAABLASTBED}" \
		>${_polap_output_dest} 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_ann_MTAABLASTBED}" \
		>"${_polap_var_ann_MTAABLAST}".sorted.bed

	mkdir "${_polap_var_ann_MTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w "${contig}" "${_polap_var_ann_MTAABLAST}".sorted.bed \
			>"${_polap_var_ann_MTAABED}/${contig}".bed
		bedtools merge -i "${_polap_var_ann_MTAABED}/${contig}".bed \
			>"${_polap_var_ann_MTAABED}/${contig}".bed.txt
		printf "%s\t%d\n" "${contig}" \
			$(wc -l <"${_polap_var_ann_MTAABED}/${contig}".bed.txt)
	done <"${_polap_var_ann_CONTIGNAME}" |
		sort -k2 -rn >"${_polap_var_ann_MTGENECOUNT}"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	tar zcf "${_polap_var_ann_MTAABED}".tar.gz "${_polap_var_ann_MTAABED}"
	rm -rf "${_polap_var_ann_MTAABED}"

	_polap_log1 "  output1: ${_polap_var_ann_MTGENECOUNT}"

	# Plastid gene annotation and counts
	_polap_log2 "  BLAST of the plastid proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "    executing the tblastn ... be patient!"
	tblastn -query "${PTAA}" \
		-db "${_polap_var_ann_CONTIGDB}" \
		-out "${_polap_var_ann_PTAABLAST}" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "${NT}" \
		>${_polap_output_dest} 2>&1

	Rscript "$script_dir"/run-polap-genes.R \
		"${_polap_var_ann_PTAABLAST}" \
		"${_polap_var_ann_PTAABLASTBED}" \
		>${_polap_output_dest} 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_ann_PTAABLASTBED}" \
		>"${_polap_var_ann_PTAABLAST}".sorted.bed

	mkdir "${_polap_var_ann_PTAABED}"

	_polap_log2 "  counting plastid genes in the contigs ..."
	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w "${contig}" "${_polap_var_ann_PTAABLAST}".sorted.bed \
			>"${_polap_var_ann_PTAABED}/${contig}".bed
		bedtools merge -i "${_polap_var_ann_PTAABED}/${contig}".bed \
			>"${_polap_var_ann_PTAABED}/${contig}".bed.txt
		printf "%s\t%d\n" "${contig}" \
			$(wc -l <"${_polap_var_ann_PTAABED}/${contig}".bed.txt)
	done <"${_polap_var_ann_CONTIGNAME}" |
		sort -k2 -rn >"${_polap_var_ann_PTGENECOUNT}"
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	_polap_log2 "  compressing the BLAST results of plastid gene annotation"
	tar zcf "${_polap_var_ann_PTAABED}".tar.gz "${_polap_var_ann_PTAABED}"
	rm -rf "${_polap_var_ann_PTAABED}"

	_polap_log1 "  output2: ${_polap_var_ann_PTGENECOUNT}"

	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") count-gene --test"
	_polap_log1 "NEXT: $(basename $0) count-gene -o $ODIR [-i $INUM]"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
