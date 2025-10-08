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

_polap_lib_readassemble-common-variables() {
	local -n annotatedir_ref="$1"
	local -n pt_table_ref="$2"
	local -n mt_table_ref="$3"
	local -n at_table_ref="$4"
	local -n all_table_ref="$5"

	pt_table_ref="${annotatedir_ref}/pt-contig-annotation-depth-table.txt"
	mt_table_ref="${annotatedir_ref}/contig-annotation-depth-table.txt"
	at_table_ref="${annotatedir_ref}/at/contig-annotation-depth-table.txt"
	all_table_ref="${annotatedir_ref}/assembly_info_organelle_annotation_count-all.txt"
}

# delete or redo the annotate-read-pt or annotate-read-TYPE
_polap_lib_readassemble-annotate-read-pt() {
	local type="${1:-pt}"

	local annotatedir="${_arg_outdir}/annotate-read-${type}"
	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	rm -rf "${annotatedir}"
	mkdir -p "${annotatedir}/at"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.fna

	_polap_log2 "    map reads on mitochondrial genes using minimap2"
	_polap_log3_cmdout minimap2 -cx \
		${_arg_minimap2_data_type} \
		"${MTAA}" \
		"${_arg_long_reads}" \
		-t ${_arg_threads} \
		-o "${annotatedir}"/mt.paf
	if [[ ! -s "${annotatedir}"/mt.paf ]]; then
		_polap_log0 "ERROR: No minimap2 results: ${annotatedir}/mt.paf"
		return
	fi

	_polap_log2 "    map reads on plastid genes using minimap2"
	_polap_log3_cmdout minimap2 -cx \
		${_arg_minimap2_data_type} \
		"${PTAA}" \
		"${_arg_long_reads}" \
		-t ${_arg_threads} \
		-o "${annotatedir}"/pt.paf
	if [[ ! -s "${annotatedir}"/pt.paf ]]; then
		_polap_log0 "ERROR: No minimap2 results: ${annotatedir}/pt.paf"
		return
	fi

	_polap_log2 "    create the annotation table for mt and pt genes"

	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-reads.R \
		--mt "${annotatedir}"/mt.paf \
		--pt "${annotatedir}"/pt.paf \
		--output "${annotatedir}" \
		--min-mapq ${_arg_annotate_read_min_mapq} \
		--min-identity ${_arg_annotate_read_min_identity}

	# local annotatedir1="${_arg_outdir}/annotate-read-${type}-1"
	# rm -rf "${annotatedir1}"
	# mkdir -p "${annotatedir1}/at"
	#
	# _polap_log3_cmdout bash "${_POLAPLIB_DIR}"/polap-r-reads.sh \
	# 	--mt "${annotatedir}"/mt.paf \
	# 	--pt "${annotatedir}"/pt.paf \
	# 	--output "${annotatedir1}" \
	# 	--min-mapq ${_arg_annotate_read_min_mapq} \
	# 	--min-identity ${_arg_annotate_read_min_identity}

	# --min-pt 1
}

_polap_lib_readassemble-get-annotated-read() {
	local type="${1:-pt}"

	local annotatedir="${_arg_outdir}/annotate-read-${type}"

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	#######################################################################
	# BEGIN: function _run_polap_assemble-annotated-read
	#
	# _polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
	# 	--table "${pt_table}" \
	# 	--length 3e+7 \
	# 	--output "${annotatedir}"/pt.id.txt

	# tail -n +2 "${pt_table}" | cut -f1 >"${annotatedir}"/pt.id.txt
	bash "${_POLAPLIB_DIR}"/polap-bash-filter-pt-reads.sh \
		-t "${pt_table}" \
		-l 3e+7 \
		-o "${annotatedir}/pt"
	# 2025-09-10
	# I do not know why I used this all not 3e+7 size limited one.
	cp "${annotatedir}/pt/pt-contig-annotation-depth-table.txt.pt.filtered.all.txt" \
		"${annotatedir}"/pt.id.all.txt
	#
	# because we want to recruit more intergenic reads.
	# cp "${annotatedir}/pt/pt-contig-annotation-depth-table.txt.pt.id.txt" \
	# 	"${annotatedir}"/pt.id.all.txt

	# _polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
	# 	--table "${mt_table}" \
	# 	--length 1e+8 \
	# 	--output "${annotatedir}"/mt0.id.txt

	bash "${_POLAPLIB_DIR}"/polap-bash-filter-pt-reads.sh \
		-t "${mt_table}" \
		-l 1e+8 \
		-o "${annotatedir}/mt"
	# 2025-09-10
	# I do not know why I used this all not 3e+7 size limited one.
	# cp "${annotatedir}/pt/pt-contig-annotation-depth-table.txt.pt.filtered.all.txt" \
	cp "${annotatedir}/mt/contig-annotation-depth-table.txt.pt.filtered.all.txt" \
		"${annotatedir}"/mt.id.all.txt
	#
	# because we want to recruit more intergenic reads.
	# cp "${annotatedir}/mt/contig-annotation-depth-table.txt.pt.id.txt" \
	# 	"${annotatedir}"/mt.id.all.txt
}

# used for polap mtseed annotate command
_polap_lib_readassemble-select-organelle-reads() {
	local type="${1:-pt}"

	local annotatedir="${_arg_outdir}/annotate-read-${type}"

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	_polap_lib_readassemble-annotate-read-pt "${type}"

	bash "${_POLAPLIB_DIR}"/polap-bash-filter-pt-reads.sh \
		-t "${pt_table}" \
		-l 1e+12 \
		--disp-pt 2000 \
		-o "${annotatedir}/pt"

	cp "${annotatedir}/pt/pt-contig-annotation-depth-table.txt.pt.filtered.all.txt" \
		"${annotatedir}"/pt.id.all.txt

	bash "${_POLAPLIB_DIR}"/polap-bash-filter-pt-reads.sh \
		-t "${mt_table}" \
		-l 1e+12 \
		-o "${annotatedir}/mt"

	cp "${annotatedir}/mt/contig-annotation-depth-table.txt.pt.filtered.all.txt" \
		"${annotatedir}"/mt.id.all.txt
}

# input: _arg_long_reads
_polap_lib_readassemble-assemble-annotated-read-pt() {
	local type="${1:-pt}"

	local annotatedir="${_arg_outdir}/annotate-read-${type}"

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	#######################################################################
	# BEGIN: function _run_polap_assemble-annotated-read
	#
	# _polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
	# 	--table "${pt_table}" \
	# 	--length 3e+7 \
	# 	--output "${annotatedir}"/pt.id.txt

	# tail -n +2 "${pt_table}" | cut -f1 >"${annotatedir}"/pt.id.txt
	bash "${_POLAPLIB_DIR}"/polap-bash-filter-pt-reads.sh \
		-t "${pt_table}" \
		-l 3e+7 \
		-o "${annotatedir}/pt"
	cp "${annotatedir}/pt/pt-contig-annotation-depth-table.txt.pt.id.txt" \
		"${annotatedir}/pt.id.txt"

	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/pt.id.txt >"${annotatedir}"/pt.fq

	# subsample the data so that pt.fq is less than 30Mb
	rm -f "${annotatedir}"/pt0.fq
	_polap_lib_fastq-sample-to \
		"${annotatedir}"/pt.fq "${annotatedir}"/pt0.fq "10m"

	# flye v2.9.6
	if _polap_lib_version-check_flye_version; then
		rm -rf "${annotatedir}"/pt
		_polap_log1 "flye 2.9.6 assembly of ptDNA using selected reads: ${annotatedir}/pt0.fq"
		flye "${_arg_flye_data_type}" \
			"${annotatedir}"/pt0.fq \
			-t "${_arg_threads}" \
			--stop-after contigger \
			--out-dir "${annotatedir}"/pt \
			2>"${_polap_output_dest}"
		if [[ -s "${annotatedir}"/pt/30-contigger/graph_final.gfa ]]; then
			_polap_log2 "    output: PT assembly: ${annotatedir}/pt/30-contigger/graph_final.gfa"
		else
			_polap_log1 "output: no PT assembly"
			return
		fi
	else
		_polap_log0 "ERROR: Flye 2.9.6 is required. Aborting."
		exit 1
	fi
	#
	# END: function _run_polap_assemble-annotated-read
	#######################################################################

	#######################################################################
	# BEGIN: function readassemble-ont-pt-iterate_genus_species
	#
	local ptdna0_gfa="${annotatedir}/pt/30-contigger/graph_final.gfa"
	if [[ -s "${ptdna0_gfa}" ]]; then
		_polap_lib_pt-extract-dna \
			"${ptdna0_gfa}" \
			"${annotatedir}/pt/ptdna"

		_polap_lib_bandage \
			"${ptdna0_gfa}" \
			"${annotatedir}/pt/30-contigger/graph_final.png"

		ln -sf "pt/ptdna/pt.0.fa" \
			"${annotatedir}/pt.0.fa"

		ln -sf "pt/30-contigger/graph_final.gfa" \
			"${annotatedir}/pt.0.gfa"

		ln -sf "pt/30-contigger/graph_final.png" \
			"${annotatedir}/pt.0.png"
	else
		_polap_log0 "No ptDNA assembly stage 0"
		return
	fi

	ln -s pt "${annotatedir}"/pt0

	local i=0
	for ((i = 0; i < 8; i++)); do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the pt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i pt$i

		# NOTE: plastid seed
		_polap_lib_seed-plastid \
			-o "${annotatedir}" \
			-i pt$i -j pt$j

		if [[ ! -s "${annotatedir}/pt$i/mt.contig.name-pt$j" ]]; then
			_polap_log0 "No ptDNA seed for pt$j"
			if [[ "${_arg_data_type}" == "pacbio-raw" ]]; then
				_polap_log1 "use then input long reads with an adjusted omega for the final stage mtDNA assembly"
				_polap_lib_assemble-omega \
					-o "${annotatedir}" \
					-l "${resolved_fastq}" \
					-t pt \
					-i pt0 -j ptx
				_polap_lib_file-cleanup -d "${annotatedir}/ptx" -s 5M -a rm
			fi
			return
		fi

		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-w "${_arg_single_min}" \
			-i pt$i -j pt$j
		_polap_lib_file-cleanup -d "${annotatedir}/pt$j" -s 5M -a rm

		ptdna0_gfa="${annotatedir}/pt$j/30-contigger/graph_final.gfa"

		if [[ -s "${ptdna0_gfa}" ]]; then
			_polap_lib_pt-extract-dna \
				"${ptdna0_gfa}" \
				"${annotatedir}/pt$j/ptdna"

			_polap_lib_bandage \
				"${annotatedir}/pt$j/assembly_graph.gfa" \
				"${annotatedir}/pt$j/assembly_graph.png"

			ln -sf "pt$j/ptdna/pt.0.fa" \
				"${annotatedir}/pt.$j.fa"

			ln -sf "pt$j/assembly_graph.gfa" \
				"${annotatedir}/pt.$j.gfa"

			ln -sf "pt$j/assembly_graph.png" \
				"${annotatedir}/pt.$j.png"

			if [[ -s "${annotatedir}/pt.$j.fa" ]]; then
				_polap_log1 "ptDNA assembly: ${annotatedir}/pt.$j.gfa"
				break
			fi
		else
			_polap_log0 "No ptDNA assembly stage $j"
		fi
	done
	#
	# END: function readassemble-ont-pt-iterate_genus_species
	#######################################################################

	# use the generated reference to assemble a mitochondrial genome.
	local j=$((i + 1))
	if [[ "$i" == "7" ]]; then
		_polap_log1 "No ptDNA assembly stage $j"
	else
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i pt$j
	fi

	_polap_log0_column "${annotatedir}/pt$j/pt-contig-annotation-depth-table.txt"

	# final link
	ln -sf "annotate-read-${type}/pt.0.fa" \
		"${_arg_outdir}/${type}-pt.0.fa"

	ln -sf "annotate-read-${type}/pt.0.gfa" \
		"${_arg_outdir}/${type}-pt.0.gfa"

	ln -sf "annotate-read-${type}/pt.0.png" \
		"${_arg_outdir}/${type}-pt.0.png"

	ln -sf "annotate-read-${type}/pt.$j.fa" \
		"${_arg_outdir}/${type}-pt.1.fa"

	ln -sf "annotate-read-${type}/pt.$j.gfa" \
		"${_arg_outdir}/${type}-pt.1.gfa"

	ln -sf "annotate-read-${type}/pt.$j.png" \
		"${_arg_outdir}/${type}-pt.1.png"
}

################################################################################
# 2025-08-13
# We use hifi100k.sh or ont100k.sh
# for hifi data, we make longer unitig or reads.
#
# same as pt case
_polap_lib_readassemble-annotate-read-mt() {
	local type="${1:-mt}"

	local annotatedir="${_arg_outdir}/annotate-read-${type}"
	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	# rm -rf "${annotatedir}"
	# mkdir -p "${annotatedir}/at"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.fna
	# local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.3k.c80.fna

	if [[ ! -s "${annotatedir}"/mt.paf ]]; then
		_polap_log2 "    map reads on mitochondrial genes using minimap2"
		_polap_log3_cmdout minimap2 -cx \
			"${_arg_minimap2_data_type}" \
			"${MTAA}" \
			"${_arg_long_reads}" \
			-t ${_arg_threads} \
			-o "${annotatedir}"/mt.paf
		if [[ ! -s "${annotatedir}"/mt.paf ]]; then
			_polap_log0 "ERROR: No minimap2 results: ${annotatedir}/mt.paf"
			return
		fi
	fi

	if [[ ! -s "${annotatedir}"/pt.paf ]]; then
		_polap_log2 "    map reads on plastid genes using minimap2"
		_polap_log3_cmdout minimap2 -cx \
			"${_arg_minimap2_data_type}" \
			"${PTAA}" \
			"${_arg_long_reads}" \
			-t ${_arg_threads} \
			-o "${annotatedir}"/pt.paf
		if [[ ! -s "${annotatedir}"/pt.paf ]]; then
			_polap_log0 "ERROR: No minimap2 results: ${annotatedir}/pt.paf"
			return
		fi
	fi

	if [[ ! -s "${annotatedir}"/mt0.id.txt ]]; then
		_polap_log2 "    create the annotation table for mt and pt genes"

		_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-reads.R \
			--mt "${annotatedir}"/mt.paf \
			--pt "${annotatedir}"/pt.paf \
			--output "${annotatedir}" \
			--min-mapq ${_arg_annotate_read_min_mapq} \
			--min-identity ${_arg_annotate_read_min_identity} \
			--min-pt 1
	fi

	_polap_lib_lines-skip1 "${mt_table}" | cut -f1 >"${annotatedir}"/mt0.id.txt
	rm -f "${annotatedir}"/mt.fq
	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt.fq
}

# 2025-08-13
#
_polap_lib_readassemble-assemble-annotated-read-mt-100k() {

	local annotatedir="${_arg_outdir}"/annotate-read-mt
	local annotatedir_pt="${_arg_outdir}"/annotate-read-pt

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	# polap command: filter reference hifi
	local resolved_fastq="${_arg_long_reads}"
	local i=""
	local j=0
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-t mt \
			-i mt$i -j mt$j
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		_polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
		# -l "${annotatedir}"/mt.fq \
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-t mt \
			-i mt$i -j mt$j
	fi
	_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

	#######################################################################
	# BEGIN: function readassemble-ont-pt-iterate_genus_species
	#
	if [[ -s "${annotatedir}/mt0/assembly_graph.gfa" ]]; then
		_polap_lib_mt-extract-dna \
			"${annotatedir}/mt0/assembly_graph.gfa" \
			"${annotatedir}/mt0/mtdna"

		_polap_lib_bandage \
			"${annotatedir}/mt0/assembly_graph.gfa" \
			"${annotatedir}/mt0/assembly_graph.png"

		# ln -sf "mt0/mtdna/mt.0.fa" \
		# 	"${annotatedir}/mt.0.fa"

		ln -sf "mt0/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt0/assembly_graph.png" \
			"${annotatedir}/mt.0.png"

	else
		_polap_log0 "No mtDNA assembly stage 0"
	fi

	# use mt as mt0
	# ln -sf mt "${annotatedir}"/mt0

	# if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
	# 	if [[ -s "${_arg_outdir}/mt-pt.0.gfa" ]]; then
	# 		_polap_log3 "      polap command filter hifi data by ptDNA reference"
	# 		_polap_lib_filter-reads-by-reference \
	# 			-o "${annotatedir}" \
	# 			-l "${_arg_long_reads}" \
	# 			--reference "${_arg_outdir}"/mt-pt.0.gfa
	# 		resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	# 	else
	# 		_polap_log0 "No ptDNA for filtering: ${_arg_outdir}/mt-pt.0.gfa"
	# 	fi
	# fi

	local i
	for ((i = 0; i < 6; i++)); do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the mt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		_polap_lib_seed-mito-2 \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_polap_log0 "No mt seed for mt$j at ${annotatedir}/mt$i: ${annotatedir}/mt$i/mt.contig.name-mt$j"
			break
		fi

		# polap command: assemble-rate
		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
			_polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
			# -l "${annotatedir}"/mt.fq \
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		fi
		_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			_polap_lib_mt-extract-dna \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/mtdna"

			_polap_lib_bandage \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			# ln -sf "mt$j/mtdna/mt.0.fa" \
			# 	"${annotatedir}/mt.$j.fa"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
		else
			_polap_log0 "No mt assembly $j"
		fi

	done
	#
	# END: function readassemble-ont-pt-iterate_genus_species
	#######################################################################

	# use the generated reference to assemble a mitochondrial genome.
	local j=$((i + 1))
	# polap command: annotate
	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$i

	# polap command: seed-mito
	_polap_lib_seed-mito-2 \
		-o "${annotatedir}" \
		-i mt$i -j mt$j

	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		# --reference "${annotatedir_pt}"/pt.3.gfa \
		_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-t mt \
			-i mt$i -j mt$j
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# -l "${_arg_long_reads}" \
		# -l "${annotatedir}"/mt.fq \
		_polap_log1 "use then input long reads with an adjusted omega for the final stage mtDNA assembly"
		# -l "${_arg_long_reads}" \
		_polap_lib_assemble-omega \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-t mt \
			-i mt$i -j mt$j
	fi
	_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$j

	_polap_log0_column "${annotatedir}/mt$j/contig-annotation-depth-table.txt"

	if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt$j/assembly_graph.png"

		ln -sf "mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt.$j.gfa"

		ln -sf "mt$j/assembly_graph.png" \
			"${annotatedir}/mt.$j.png"

		_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
	else
		_polap_log0 "No MT assembly $j"
	fi

	# extract mtDNA sequence if you can
	#

	# final link
	i=$((j - 1))
	ln -sf "annotate-read-mt/mt.$i.gfa" \
		"${_arg_outdir}/mt.0.gfa"

	ln -sf "annotate-read-mt/mt.$i.png" \
		"${_arg_outdir}/mt.0.png"

	ln -sf "annotate-read-mt/mt.$j.gfa" \
		"${_arg_outdir}/mt.1.gfa"

	ln -sf "annotate-read-mt/mt.$j.png" \
		"${_arg_outdir}/mt.1.png"
}

# input: _arg_long_reads
_polap_lib_readassemble-assemble-annotated-read-mt() {

	local annotatedir="${_arg_outdir}"/annotate-read-mt
	local annotatedir_pt="${_arg_outdir}"/annotate-read-pt

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	if [[ ! -s "${mt_table}" ]]; then
		_polap_log0 "ERROR: no mt table: ${mt_table}"
		return 1
	fi

	#######################################################################
	# BEGIN: function _run_polap_assemble-annotated-read
	#
	# for animal mitochondrial genome assembly
	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${mt_table}" \
		--length 1e+8 \
		--output "${annotatedir}"/mt0.id.txt

	# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
	rm -f "${annotatedir}"/mt.fq
	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt.fq

	# subsample the data so that mt.fq is less than 100Mb
	rm -f "${annotatedir}"/mt0.fq
	_polap_lib_fastq-sample-to \
		"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "100m"

	# flye v2.9.6
	if _polap_lib_version-check_flye_version; then
		rm -rf "${annotatedir}"/mt
		_polap_log0 "flye 2.9.6 assembly of ptDNA using selected reads: ${annotatedir}/mt0.fq"
		flye "${_arg_flye_data_type}" \
			"${annotatedir}"/mt0.fq \
			-t "${_arg_threads}" \
			--out-dir "${annotatedir}"/mt \
			2>"${_polap_output_dest}"
		if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
			_polap_log2 "    output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
		else
			_polap_log0 "output: no MT assembly"
			return
		fi
	else
		echo "Flye 2.9.6 is required. Aborting."
		exit 1
	fi

	#
	# END: function _run_polap_assemble-annotated-read
	#######################################################################

	#######################################################################
	# BEGIN: function readassemble-ont-pt-iterate_genus_species
	#
	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		_polap_lib_mt-extract-dna \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/mtdna"

		_polap_lib_bandage \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/mtdna/mt.0.fa" \
			"${annotatedir}/mt.0.fa"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"

	else
		_polap_log0 "No mtDNA assembly stage 0"
	fi

	# use mt as mt0
	ln -sf mt "${annotatedir}"/mt0

	# polap command: filter reference hifi
	local resolved_fastq="${_arg_long_reads}"
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		if [[ -s "${_arg_outdir}/pt.0.gfa" ]]; then
			_polap_log3 "      polap command filter hifi data by ptDNA reference"
			_polap_lib_filter-reads-by-reference \
				-o "${annotatedir}" \
				-l "${_arg_long_reads}" \
				--reference "${_arg_outdir}"/pt.0.gfa
			resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
		else
			_polap_log0 "No ptDNA for filtering: ${_arg_outdir}/pt.0.gfa"
		fi
	fi

	local i
	for ((i = 0; i < 6; i++)); do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the mt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		_polap_lib_seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_polap_log0 "No mt seed for mt$j"
			break
		fi

		# polap command: assemble-rate
		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
			_polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${annotatedir}"/mt.fq \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		fi
		_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			_polap_lib_mt-extract-dna \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/mtdna"

			_polap_lib_bandage \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/mtdna/mt.0.fa" \
				"${annotatedir}/mt.$j.fa"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
		else
			_polap_log0 "No mt assembly $j"
		fi

	done
	#
	# END: function readassemble-ont-pt-iterate_genus_species
	#######################################################################

	# use the generated reference to assemble a mitochondrial genome.
	local j=$((i + 1))
	# polap command: annotate
	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$i

	# polap command: seed-mito
	_polap_lib_seed-mito \
		-o "${annotatedir}" \
		-i mt$i -j mt$j

	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		# --reference "${annotatedir_pt}"/pt.3.gfa \
		_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-i mt$i -j mt$j
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# -l "${_arg_long_reads}" \
		# -l "${annotatedir}"/mt.fq \
		_polap_log1 "use then input long reads with an adjusted omega for the final stage mtDNA assembly"
		_polap_lib_assemble-omega \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-i mt$i -j mt$j
	fi
	_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$j

	_polap_log0_column "${annotatedir}/mt$j/contig-annotation-depth-table.txt"

	if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt$j/assembly_graph.png"

		ln -sf "mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt.$j.gfa"

		ln -sf "mt$j/assembly_graph.png" \
			"${annotatedir}/mt.$j.png"

		_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
	else
		_polap_log0 "No MT assembly $j"
	fi

	# extract mtDNA sequence if you can
	#

	# final link
	i=$((j - 1))
	ln -sf "annotate-read-mt/mt.$i.gfa" \
		"${_arg_outdir}/mt.0.gfa"

	ln -sf "annotate-read-mt/mt.$i.png" \
		"${_arg_outdir}/mt.0.png"

	ln -sf "annotate-read-mt/mt.$j.gfa" \
		"${_arg_outdir}/mt.1.gfa"

	ln -sf "annotate-read-mt/mt.$j.png" \
		"${_arg_outdir}/mt.1.png"
}

# be called
_polap_lib_readassemble-annotated() {
	# local type="${1:-pt}"
	# local annotatedir="${_arg_outdir}/annotate-read-${type}"

	local mt="${_arg_inum}"
	local annotatedir="${_arg_outdir}"

	local annotatedir_pt="${_arg_outdir}"/annotate-read-pt

	#######################################################################
	# BEGIN: function readassemble-ont-pt-iterate_genus_species
	#
	if [[ -s "${annotatedir}/$mt/30-contigger/graph_final.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/$mt/30-contigger/graph_final.gfa" \
			"${annotatedir}/$mt/30-contigger/graph_final.png"

		ln -sf "$mt/30-contigger/graph_final.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "$mt/30-contigger/graph_final.png" \
			"${annotatedir}/mt.0.png"

	else
		_polap_log0 "No mtDNA assembly stage 0"
	fi

	# use mt as mt0
	_polap_log0 "Name the start i"
	ln -sf $mt "${annotatedir}"/mt0

	# polap command: filter reference hifi
	local resolved_fastq="${_arg_long_reads}"
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		if [[ -s "${_arg_outdir}/pt.0.gfa" ]]; then
			_polap_log3 "      polap command filter hifi data by ptDNA reference"
			_polap_lib_filter-reads-by-reference \
				-o "${annotatedir}" \
				-l "${_arg_long_reads}" \
				--reference "${_arg_outdir}"/pt.0.gfa
			resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
		else
			_polap_log0 "No ptDNA for filtering: ${_arg_outdir}/pt.0.gfa"
		fi
	fi

	local i
	for ((i = 0; i < ${_arg_readassemble_mtn}; i++)); do
		local j=$((i + 1))
		_polap_log0 "i: $i"

		# NOTE: annotate for seeding
		# select connected components of the mt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i mt$i

		_polap_log0 "i: $i"
		# NOTE: mito seed
		_polap_lib_seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		_polap_log0 "i: $i"
		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_polap_log0 "No mt seed for mt$j"
			break
		fi

		# polap command: assemble-rate
		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
			_polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
			_polap_log0 "i: $i"
			_polap_lib_assemble-omega \
				-o "${annotatedir}" \
				-l "${_arg_long_reads}" \
				-t mt \
				-i mt$i -j mt$j
			_polap_log0 "i: $i"

			# _polap_lib_assemble-rate \
			# 	-o "${annotatedir}" \
			# 	-l "${_arg_long_reads}" \
			# 	-w "${_arg_single_min}" \
			# 	-t mt \
			# 	-i mt$i -j mt$j
		fi
		_polap_log0 "i: $i"
		_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

		_polap_log0 "i: $i"
		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			# _polap_lib_mt-extract-dna \
			# 	"${annotatedir}/mt$j/assembly_graph.gfa" \
			# 	"${annotatedir}/mt$j/mtdna"

			_polap_lib_bandage \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"
			_polap_log0 "i: $i"

			ln -sf "mt$j/mtdna/mt.0.fa" \
				"${annotatedir}/mt.$j.fa"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			_polap_log0 "i: $i"
			_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
		else
			_polap_log0 "No mt assembly $j"
		fi

	done

	_polap_log0 "after the for-loop"
	_polap_log0 "i: $i"
	_polap_log0 "j: $j"

	#
	# END: function readassemble-ont-pt-iterate_genus_species
	#######################################################################

	# why mt2 not mt1?

	# use the generated reference to assemble a mitochondrial genome.
	local j=$((i + 1))
	_polap_log0 "after the j update"
	_polap_log0 "i: $i"
	_polap_log0 "j: $j"
	# polap command: annotate
	_polap_log0 _polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$i
	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$i

	# polap command: seed-mito
	_polap_lib_seed-mito \
		-o "${annotatedir}" \
		-i mt$i -j mt$j

	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		# --reference "${annotatedir_pt}"/pt.3.gfa \
		_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-i mt$i -j mt$j
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# -l "${_arg_long_reads}" \
		# -l "${annotatedir}"/mt.fq \
		_polap_log1 "use then input long reads with an adjusted omega for the final stage mtDNA assembly"
		_polap_lib_assemble-omega \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-t mt \
			-i mt$i -j mt$j
	fi
	_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$j

	_polap_log0_column "${annotatedir}/mt$j/contig-annotation-depth-table.txt"

	if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt$j/assembly_graph.png"

		ln -sf "mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt.$j.gfa"

		ln -sf "mt$j/assembly_graph.png" \
			"${annotatedir}/mt.$j.png"

		_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
	else
		_polap_log0 "No MT assembly $j"
	fi

	# extract mtDNA sequence if you can
	#

	# final link
	i=$((j - 1))
	ln -sf "mtseed/mt.$i.gfa" \
		"${_arg_outdir}/../mt.0.gfa"

	ln -sf "mtseed/mt.$i.png" \
		"${_arg_outdir}/../mt.0.png"

	ln -sf "mtseed/mt.$j.gfa" \
		"${_arg_outdir}/../mt.1.gfa"

	ln -sf "mtseed/mt.$j.png" \
		"${_arg_outdir}/../mt.1.png"
}

################################################################################
# This may work for some ONT mtDNA assembly but not all.
################################################################################
_polap_lib_readassemble-annotate-read-nt() {
	local type="nt"

	local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.c80.fna
	local annotatedir="${_arg_outdir}/annotate-read-${type}"

	_polap_lib_readassemble-annotate-read-pt nt

	_polap_log2 "    map reads on mitochondrial noncoding regions using minimap2"
	# v0.5.2.3
	# -k13 -w5 -m20 -p0.6 -N20 \
	_polap_log3_cmdout minimap2 -cx \
		"${_arg_minimap2_data_type}" \
		"${NTAA}" \
		"${_arg_long_reads}" \
		-t ${_arg_threads} \
		-o "${annotatedir}"/nt.paf
	if [[ ! -s "${annotatedir}"/nt.paf ]]; then
		_polap_log0 "ERROR: No minimap2 results: ${annotatedir}/pt.paf"
		return
	fi

	_polap_log2 "    create the annotation table for mito noncoding and pt genes"
	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-reads.R \
		--mt "${annotatedir}"/nt.paf \
		--pt "${annotatedir}"/pt.paf \
		--output "${annotatedir}"/at \
		--min-mapq ${_arg_annotate_read_min_mapq} \
		--min-identity ${_arg_annotate_read_min_identity} \
		--min-pt 1
}

_polap_lib_readassemble-assemble-annotated-read-nt() {

	local annotatedir="${_arg_outdir}"/annotate-read-nt
	local annotatedir_pt="${_arg_outdir}"/annotate-read-pt

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${pt_table}" \
		--length 3e+9 \
		--output "${annotatedir}"/pt0mt.id.txt
	_polap_log1 "  pt0mt.id: $(cat "${annotatedir}"/pt0mt.id.txt | wc -l)"

	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${mt_table}" \
		--length 1e+9 \
		--output "${annotatedir}"/mt0.id.txt
	_polap_log1 "  mt0.id: $(cat "${annotatedir}"/mt0.id.txt | wc -l)"

	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${at_table}" \
		--length 1e+9 \
		--output "${annotatedir}"/at0.id.txt
	_polap_log1 "  at0.id: $(cat "${annotatedir}"/at0.id.txt | wc -l)"

	# subtract pt.id.txt from at0.id.txt
	grep -vFf "${annotatedir}"/pt0mt.id.txt "${annotatedir}"/at0.id.txt \
		>"${annotatedir}"/at.id.txt
	_polap_log1 "  at.id (at0 - pt0mt): $(cat "${annotatedir}"/at.id.txt | wc -l)"

	if [[ "${_arg_noncoding}" == "on" ]]; then
		cat "${annotatedir}"/mt0.id.txt "${annotatedir}"/at.id.txt |
			sort | uniq >"${annotatedir}"/mt.id.txt
	else
		cat "${annotatedir}"/mt0.id.txt |
			sort | uniq >"${annotatedir}"/mt.id.txt
	fi
	_polap_log1 "mt.id (mt0 + at): $(cat "${annotatedir}"/mt.id.txt | wc -l)"

	# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
	rm -f "${annotatedir}"/mt.fq
	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt.id.txt >"${annotatedir}"/mt.fq

	# subsample the data so that mt.fq is less than 100Mb
	rm -f "${annotatedir}"/mt0.fq
	_polap_lib_fastq-sample-to \
		"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "100m"

	# flye v2.9.6
	if _polap_lib_version-check_flye_version; then
		rm -rf "${annotatedir}"/mt
		_polap_log1 "flye 2.9.6 assembly of ntDNA using selected reads: ${annotatedir}/mt0.fq"
		flye "${_arg_flye_data_type}" \
			"${annotatedir}"/mt0.fq \
			-t "${_arg_threads}" \
			--out-dir "${annotatedir}"/mt \
			2>"${_polap_output_dest}"
		if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
			_polap_log2 "    output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
		else
			_polap_log0 "output: no MT assembly"
			return
		fi
	else
		echo "Flye 2.9.6 is required. Aborting."
		exit 1
	fi

	# create lib for the following polap command.
	#
	# polap command: mt
	_polap_lib_mt-extract-dna \
		"${annotatedir}/mt/assembly_graph.gfa" \
		"${annotatedir}/mt/mtdna"

	# polap command: bandage png
	_polap_lib_bandage \
		"${annotatedir}/mt/assembly_graph.gfa" \
		"${annotatedir}/mt/assembly_graph.png"

	# create links to the gfa and png files
	ln -sf "mt/assembly_graph.gfa" \
		"${annotatedir}/mt.0.gfa"

	ln -sf "mt/assembly_graph.png" \
		"${annotatedir}/mt.0.png"

	# use mt as mt0
	ln -sf mt "${annotatedir}"/mt0

	# polap command: filter reference hifi
	local resolved_fastq="${_arg_long_reads}"
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		if [[ -s "${_arg_outdir}/pt-pt.1.gfa" ]]; then
			_polap_log3 "      polap command filter hifi data by ptDNA reference"
			_polap_lib_filter-reads-by-reference \
				-o "${annotatedir}" \
				-l "${_arg_long_reads}" \
				--reference "${_arg_outdir}"/pt-pt.1.gfa
			resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
		else
			_polap_log0 "No ptDNA for filtering: ${_arg_outdir}/pt-pt.1.gfa"
		fi
	fi

	# repeat several times
	# polap command: annotate
	# polap command: seed-mito
	# polap command: assemble-rate
	# polap command: mt
	# polap command: bandage png
	local i
	for ((i = 0; i < 8; i++)); do
		local j=$((i + 1))
		# polap command: annotate
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i mt$i

		# polap command: seed-mito
		_polap_lib_seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		# polap command: assemble-rate
		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
			_polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${annotatedir}"/mt.fq \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		fi
		_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

		# polap command: mt
		mtdna0_gfa="${annotatedir}/mt$j/30-contigger/graph_final.gfa"

		# polap command: bandage png
		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			_polap_lib_bandage \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
		else
			_polap_log0 "No MT assembly yet in stage $j"
		fi

	done

	# use the generated reference to assemble a mitochondrial genome.
	local j=$((i + 1))
	# polap command: annotate
	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$i

	# polap command: seed-mito
	_polap_lib_seed-mito \
		-o "${annotatedir}" \
		-i mt$i -j mt$j

	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		# --reference "${annotatedir_pt}"/pt.3.gfa \
		_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-t mt \
			-i mt$i -j mt$j
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# -l "${_arg_long_reads}" \
		# _polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
		#
		# -l "${annotatedir}"/mt.fq \
		_polap_log1 "use then input long reads with an adjusted omega for the final stage mtDNA assembly"
		_polap_lib_assemble-omega \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-i mt$i -j mt$j
	fi
	_polap_lib_file-cleanup -d "${annotatedir}/mt$j" -s 5M -a rm

	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$j

	_polap_log0_column "${annotatedir}/mt$j/contig-annotation-depth-table.txt"

	if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt$j/assembly_graph.png"

		ln -sf "mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt.$j.gfa"

		ln -sf "mt$j/assembly_graph.png" \
			"${annotatedir}/mt.$j.png"

		_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
	else
		_polap_log0 "No MT assembly $j"
	fi

	# final link
	ln -sf "annotate-read-nt/mt.0.gfa" \
		"${_arg_outdir}/nt.0.gfa"

	ln -sf "annotate-read-nt/mt.0.png" \
		"${_arg_outdir}/nt.0.png"

	ln -sf "annotate-read-nt/mt.$j.gfa" \
		"${_arg_outdir}/nt.1.gfa"

	ln -sf "annotate-read-nt/mt.$j.png" \
		"${_arg_outdir}/nt.1.png"
}
