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
# Convert numbers between different units.
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
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
_polap_lib_annotate() {
	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	_polap_log1 "  annotate ${_arg_outdir}/${_arg_inum} using amino acid sequences"
	_polap_lib_annotate-edges-stats \
		-o "${_arg_outdir}" -i "${_arg_inum}" -j "${_arg_jnum}"
	_polap_lib_annotate-blast-genome \
		-o "${_arg_outdir}" -i "${_arg_inum}" -j "${_arg_jnum}"
	_polap_lib_annotate-count-gene \
		-o "${_arg_outdir}" -i "${_arg_inum}" -j "${_arg_jnum}"
}

# 2025-08-14
# Always redo or we do not use any previous output files.
# input: gfa
_polap_lib_annotate-edges-stats() {

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	_polap_log1 "  creating edges_stats.txt from graph_final.gfa ..."
	_polap_log1 "  creating GFA without sequence data: ${_polap_var_ga_gfa_all}"
	_polap_log2 "    input: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log2 "    output: ${_polap_var_ga_gfa_all}"

	if [[ ! -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
		_polap_log0 "ERROR: no such file: ${_polap_var_ga_contigger_edges_gfa}"
		return $RETURN_FAIL
	fi

	# if [ -s "${_polap_var_ga_gfa_all}" ] && [ "${_arg_redo}" = "off" ]; then
	# 	_polap_log1 "    found: ${_polap_var_ga_gfa_all}, so skipping ..."
	# else
	_polap_log3_pipe "gfatools view \
		  -S ${_polap_var_ga_contigger_edges_gfa} \
		  >${_polap_var_ga_gfa_all} \
		  2>$_polap_output_dest"
	# fi

	_polap_log2 "    extracting sequence part of GFA: ${_polap_var_ga_gfa_seq_part}"
	_polap_log3 "      input: ${_polap_var_ga_gfa_all}"
	_polap_log3 "      output: ${_polap_var_ga_gfa_seq_part}"
	# if [ -s "${_polap_var_ga_gfa_seq_part}" ] && [ "${_arg_redo}" = "off" ]; then
	# 	_polap_log1 "    found: ${_polap_var_ga_gfa_seq_part}, so skipping ..."
	# else
	_polap_log3_pipe "grep ^S ${_polap_var_ga_gfa_all} >${_polap_var_ga_gfa_seq_part}"
	# fi

	# Filter edges in GFA using depths.
	_polap_log2 "    filtering GFA sequence part using depth range"
	_polap_log3 "      input1: ${_polap_var_ga_gfa_seq_part}"
	_polap_log3 "      output1: ${_polap_var_ga_contigger_edges_stats}"
	# if [ -s "${_polap_var_ga_contigger_edges_stats}" ] && [ "${_arg_redo}" = "off" ]; then
	# 	_polap_log1 "    found: ${_polap_var_ga_contigger_edges_stats}, so skipping ..."
	# else
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/polap-r-edges-stats.R \
		--gfa ${_polap_var_ga_gfa_seq_part} \
		--out ${_polap_var_ga_contigger_edges_stats} \
		2>$_polap_output_dest"
	# fi

	if [[ -s "${_polap_var_ga_contigger_edges_stats}" ]]; then
		_polap_log2_column "${_polap_var_ga_contigger_edges_stats}"
	else
		_polap_log2 "ERROR: no such file: ${_polap_var_ga_contigger_edges_stats}"
	fi

	return 0
}

# 2025-08-14
# Always redo or we do not use any previous output files.
_polap_lib_annotate-count-gene() { # count MT and PT genes using edges_stats.txt

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	_polap_log1 "  counting mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."

	# Checks the output files earlier than the input.
	# if [[ -s "${_polap_var_ga_annotation_all}" ]] &&
	# 	[[ "${_arg_redo}" = "off" ]]; then
	# 	_polap_log0 "  found1: ${_polap_var_ga_annotation_all}, so skipping the blast genome ..."
	# 	# Disable debugging if previously enabled
	# 	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	# 	return 0
	# fi

	_polap_log2 "    creating annotation tables"
	_polap_log3 "      input1: ${_polap_var_ga_contigger_edges_stats}"
	_polap_log3 "      input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log3 "      input3: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log3 "      output1: ${_polap_var_ga_annotation}"
	_polap_log3 "      output2: ${_polap_var_ga_annotation_all}"
	_polap_log3 "      output3: ${_polap_var_ga_annotation_table}"
	_polap_log3 "      output4: ${_polap_var_ga_annotation_depth_table}"
	_polap_log3 "      output5: ${_polap_var_ga_pt_annotation_depth_table}"
	local _command1="Rscript ${_POLAPLIB_DIR}/polap-r-mtcontig.R \
		--flyeout-edges-stats ${_polap_var_ga_contigger_edges_stats} \
    --mt-gene-count ${_polap_var_ann_MTGENECOUNT} \
    --pt-gene-count ${_polap_var_ann_PTGENECOUNT} \
		--out-annotation ${_polap_var_ga_annotation} \
		--out-annotation-all ${_polap_var_ga_annotation_all} \
		--out-annotation-table ${_polap_var_ga_annotation_table} \
		--out-annotation-depth-table ${_polap_var_ga_annotation_depth_table} \
		--out-pt-annotation-depth-table ${_polap_var_ga_pt_annotation_depth_table} \
		--fixed \
		--contigger"
	if [[ "${_arg_plastid}" = "on" ]]; then
		_command1+=" \
      --plastid"
	fi
	_command1+=" \
      2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	if [[ -s "${_polap_var_ga_annotation_all}" ]]; then
		_polap_log3_column "${_polap_var_ga_annotation_all}"
	fi

	return 0
}

# 2025-08-14
# Always redo or we do not use any previous output files.
_polap_lib_annotate-blast-genome() {

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.faa
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.faa

	_polap_log1 "  blasting edge contigs with mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."
	_polap_log2 "    input1: ${_polap_var_ga_contigger_edges_stats}"
	_polap_log2 "    input2: ${_polap_var_ga_contigger_edges_fasta}"

	# Checks the output files earlier than the input.
	# if [[ -s "${_polap_var_ann_MTGENECOUNT}" ]] &&
	# 	[[ -s "${_polap_var_ann_PTGENECOUNT}" ]] &&
	# 	[[ "${_arg_redo}" = "off" ]]; then
	# 	_polap_log0 "  found1: ${_polap_var_ann_MTGENECOUNT}"
	# 	_polap_log0 "  found2: ${_polap_var_ann_PTGENECOUNT}"
	# 	_polap_log0 "  so skipping the blast genome ..."
	# 	# Disable debugging if previously enabled
	# 	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	# 	return 0
	# fi

	if [ ! -s "${_polap_var_ga_contigger_edges_stats}" ]; then
		_polap_log0 "ERROR: no edges_stats.txt file: ${_polap_var_ga_contigger_edges_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "    deleting and recreating the folder: ${_polap_var_ann}"
	_polap_log3_cmd rm -rf "${_polap_var_ann}"
	_polap_log3_cmd mkdir -p "${_polap_var_ann}"

	#src/run-polap-select.R o/30-contigger/edges_stats.txt o/50-annotation/contig.name
	_polap_log2 "    contig sequence names in file: ${_polap_var_ann_CONTIGNAME}"
	_polap_log3 "      input: ${_polap_var_ga_contigger_edges_stats}"
	_polap_log3 "      output: ${_polap_var_ann_CONTIGFILE}"
	_polap_log3_pipe "grep -v '#' ${_polap_var_ga_contigger_edges_stats} |
		cut -f 1 >${_polap_var_ann_CONTIGNAME}"

	_polap_log2 "  preparing edge contig sequence file: ${_polap_var_ann_CONTIGFILE}"
	# if [[ -s "${_polap_var_ga_contigger_edges_fasta}" ]]; then
	# 	_polap_log3 "    copying edge contig sequence file"
	# 	_polap_log3_cmd cp "${_polap_var_ga_contigger_edges_fasta}" "${_polap_var_ann_CONTIGFILE}"
	# else
	if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
		_polap_log3 "    extracting edge contig sequence from the gfa"
		_polap_log3_pipe "gfatools gfa2fa \
		    ${_polap_var_ga_contigger_edges_gfa} \
		    >${_polap_var_ga_contigger_edges_fasta} 2>${_polap_output_dest}"
		_polap_log3_cmd cp "${_polap_var_ga_contigger_edges_fasta}" "${_polap_var_ann_CONTIGFILE}"
	else
		_polap_log0 "ERROR: no such file: ${_polap_var_ga_contigger_edges_gfa}"
		return $RETURN_FAIL
	fi
	# fi

	_polap_log2 "    making BLASTDB of the contig sequences: ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "      input: ${_polap_var_ann_CONTIGFILE}"
	_polap_log3 "      output: ${_polap_var_ann_CONTIGDB}"
	_polap_log3_pipe "makeblastdb -dbtype nucl \
		-in ${_polap_var_ann_CONTIGFILE} \
		-out ${_polap_var_ann_CONTIGDB} \
		2>${_polap_output_dest}"

	# Mitochondrial gene annotation and counts
	_polap_log2 "    BLAST of the mitochondrial proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "      input query: ${MTAA}"
	_polap_log3 "      input BLAST DB: ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "      evalue cutoff: 1e-30"
	_polap_log3 "      number of CPUs: ${_arg_threads}"
	_polap_log3 "      executing the tblastn ... be patient!"
	_polap_log3_pipe "tblastn -query ${MTAA} \
		-db ${_polap_var_ann_CONTIGDB} \
		-out ${_polap_var_ann_MTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "    converting BLAST result in BED format for removing redundancy"
	_polap_log3 "      input1: ${_polap_var_ann_MTAABLAST}"
	_polap_log3 "      output1: ${_polap_var_ann_MTAABLASTBED}"
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/polap-r-genes.R \
		${_polap_var_ann_MTAABLAST} \
		${_polap_var_ann_MTAABLASTBED} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "    sorting BED format BLAST result"
	_polap_log3 "      input1: ${_polap_var_ann_MTAABLASTBED}"
	_polap_log3 "      output1: ${_polap_var_ann_MTAABLAST}.sorted.bed"
	_polap_log3_pipe "sort -k1,1 -k2,2n \
		${_polap_var_ann_MTAABLASTBED} \
		>${_polap_var_ann_MTAABLAST}.sorted.bed 2>${_polap_output_dest}"
	if [[ ! -s "${_polap_var_ann_MTAABLASTBED}" ]]; then
		_polap_log0 "  -> no mitochondrial genes"
	fi

	_polap_log2 "    creating folder: ${_polap_var_ann_MTAABED}"
	_polap_log3_cmd mkdir "${_polap_var_ann_MTAABED}"

	_polap_log2 "    counting mitochondrial genes in the contigs ..."
	_polap_log3 "      input1: ${_polap_var_ann_MTAABLAST}.sorted.bed"
	_polap_log3 "      output1: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log3_cmd rm -f "${_polap_var_ann_MTGENECOUNT_TMP}"
	while IFS= read -r contig; do
		_polap_log3_pipe "grep -w ${contig} ${_polap_var_ann_MTAABLAST}.sorted.bed \
			>${_polap_var_ann_MTAABED}/${contig}.bed"
		_polap_log3_pipe "bedtools merge -i ${_polap_var_ann_MTAABED}/${contig}.bed \
			>${_polap_var_ann_MTAABED}/${contig}.bed.txt"
		_polap_log3 commnad: printf \"%s\\t%d\\n\" "${contig}" $(wc -l <"${_polap_var_ann_MTAABED}/${contig}".bed.txt) ">>${_polap_var_ann_MTGENECOUNT_TMP}"
		printf "%s\t%d\n" ${contig} $(wc -l <${_polap_var_ann_MTAABED}/${contig}.bed.txt) >>"${_polap_var_ann_MTGENECOUNT_TMP}"
	done <"${_polap_var_ann_CONTIGNAME}"
	_polap_log3_pipe "sort -k2 -rn ${_polap_var_ann_MTGENECOUNT_TMP} >${_polap_var_ann_MTGENECOUNT}"

	_polap_log2 "    compressing the BLAST results of mitochondrial gene annotation"
	_polap_log3 "      input: ${_polap_var_ann_MTAABED}"
	_polap_log3 "      output: ${_polap_var_ann_MTAABED}.tar.gz"
	_polap_log3_cmd tar zcf "${_polap_var_ann_MTAABED}".tar.gz "${_polap_var_ann_MTAABED}"
	_polap_log2 "  deleting folder: ${_polap_var_ann_MTAABLASTBED}"
	_polap_log3_cmd rm -rf "${_polap_var_ann_MTAABED}"

	_polap_log1 "  output1: ${_polap_var_ann_MTGENECOUNT}"

	# Plastid gene annotation and counts
	_polap_log2 "    BLAST of the plastid proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "      input query: ${PTAA}"
	_polap_log3 "      input BLAST DB: ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "      evalue cutoff: 1e-30"
	_polap_log3 "      number of CPUs: ${_arg_threads}"
	_polap_log3 "      executing the tblastn ... be patient!"
	_polap_log3_pipe "tblastn -query ${PTAA} \
		-db ${_polap_var_ann_CONTIGDB} \
		-out ${_polap_var_ann_PTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "    converting BLAST result in BED format for removing redundancy"
	_polap_log3 "      input1: ${_polap_var_ann_PTAABLAST}"
	_polap_log3 "      output1: ${_polap_var_ann_PTAABLASTBED}"
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/polap-r-genes.R \
		${_polap_var_ann_PTAABLAST} \
		${_polap_var_ann_PTAABLASTBED} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "    sorting BED format BLAST result"
	_polap_log3 "      input1: ${_polap_var_ann_PTAABLASTBED}"
	_polap_log3 "      output1: ${_polap_var_ann_PTAABLAST}.sorted.bed"
	_polap_log3_pipe "sort -k1,1 -k2,2n \
		${_polap_var_ann_PTAABLASTBED} \
		>${_polap_var_ann_PTAABLAST}.sorted.bed 2>${_polap_output_dest}"
	if [[ ! -s "${_polap_var_ann_PTAABLASTBED}" ]]; then
		_polap_log0 "  -> no plastid genes"
	fi

	_polap_log2 "    creating folder: ${_polap_var_ann_PTAABED}"
	_polap_log3_cmd mkdir "${_polap_var_ann_PTAABED}"

	_polap_log2 "    counting mitochondrial genes in the contigs ..."
	_polap_log3 "      input1: ${_polap_var_ann_PTAABLAST}.sorted.bed"
	_polap_log3 "      output1: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log3_cmd rm -f "${_polap_var_ann_PTGENECOUNT_TMP}"
	while IFS= read -r contig; do
		_polap_log3_pipe "grep -w ${contig} ${_polap_var_ann_PTAABLAST}.sorted.bed \
			>${_polap_var_ann_PTAABED}/${contig}.bed"
		_polap_log3_pipe "bedtools merge -i ${_polap_var_ann_PTAABED}/${contig}.bed \
			>${_polap_var_ann_PTAABED}/${contig}.bed.txt"
		_polap_log3 commnad: printf \"%s\\t%d\\n\" "${contig}" $(wc -l <"${_polap_var_ann_PTAABED}/${contig}".bed.txt) ">>${_polap_var_ann_PTGENECOUNT_TMP}"
		printf "%s\t%d\n" ${contig} $(wc -l <${_polap_var_ann_PTAABED}/${contig}.bed.txt) >>"${_polap_var_ann_PTGENECOUNT_TMP}"
	done <"${_polap_var_ann_CONTIGNAME}"
	_polap_log3_pipe "sort -k2 -rn ${_polap_var_ann_PTGENECOUNT_TMP} >${_polap_var_ann_PTGENECOUNT}"

	_polap_log2 "    compressing the BLAST results of mitochondrial gene annotation"
	_polap_log3 "      input: ${_polap_var_ann_PTAABED}"
	_polap_log3 "      output: ${_polap_var_ann_PTAABED}.tar.gz"
	_polap_log3_cmd tar zcf "${_polap_var_ann_PTAABED}".tar.gz "${_polap_var_ann_PTAABED}"
	_polap_log2 "  deleting folder: ${_polap_var_ann_PTAABLASTBED}"
	_polap_log3_cmd rm -rf "${_polap_var_ann_PTAABED}"

	_polap_log1 "  output2: ${_polap_var_ann_PTGENECOUNT}"

	return 0
}
