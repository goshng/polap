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
	# local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.3k.c80.fna

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

	_polap_log2 "    create the annotation table for mt and pt genes"
	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-reads.R \
		--mt "${annotatedir}"/mt.paf \
		--pt "${annotatedir}"/pt.paf \
		--output "${annotatedir}" \
		--min-mapq ${_arg_annotate_read_min_mapq} \
		--min-identity ${_arg_annotate_read_min_identity}

}

# input: _arg_long_reads
_polap_lib_readassemble-assemble-annotated-read-pt() {

	local annotatedir="${_arg_outdir}"/annotate-read-pt

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	#######################################################################
	# BEGIN: function _run_polap_assemble-annotated-read
	#
	Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${pt_table}" \
		--length 3e+7 \
		--output "${annotatedir}"/pt.id.txt
	# tail -n +2 "${pt_table}" | cut -f1 >"${annotatedir}"/pt.id.txt

	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/pt.id.txt >"${annotatedir}"/pt.fq

	# subsample the data so that pt.fq is less than 30Mb
	rm -f "${annotatedir}"/pt0.fq
	_polap_lib_fastq-sample-to \
		"${annotatedir}"/pt.fq "${annotatedir}"/pt0.fq "30m"

	# flye v2.9.6
	if _polap_lib_version-check_flye_version; then
		rm -rf "${annotatedir}"/pt
		_polap_log0 "flye 2.9.6 assembly on ptDNA"
		flye "${_arg_flye_data_type}" "${annotatedir}"/pt0.fq -t "${_arg_threads}" \
			--out-dir "${annotatedir}"/pt 2>"${_polap_output_dest}"
		if [[ -s "${annotatedir}"/pt/assembly_graph.gfa ]]; then
			_polap_log1 "output: PT assembly: ${annotatedir}/pt/assembly_graph.gfa"
		else
			_polap_log0 "output: no PT assembly"
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
	if [[ -s "${annotatedir}/pt/assembly_graph.gfa" ]]; then
		_polap_lib_pt-extract-dna \
			"${annotatedir}/pt/assembly_graph.gfa" \
			"${annotatedir}/pt/ptdna"

		_polap_lib_bandage \
			"${annotatedir}/pt/assembly_graph.gfa" \
			"${annotatedir}/pt/assembly_graph.png"

		ln -sf "pt/ptdna/pt.0.fa" \
			"${annotatedir}/pt.0.fa"

		ln -sf "pt/assembly_graph.gfa" \
			"${annotatedir}/pt.0.gfa"

		ln -sf "pt/assembly_graph.png" \
			"${annotatedir}/pt.0.png"

	else
		_log_echo "No pt assembly 0"
	fi

	ln -s pt "${annotatedir}"/pt0

	local i
	for i in {0..1}; do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the pt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i pt$i

		# NOTE: mito seed
		_polap_lib_seed-plastid \
			-o "${annotatedir}" \
			-i pt$i -j pt$j

		if [[ ! -s "${annotatedir}/pt$i/mt.contig.name-pt$j" ]]; then
			_log_echo "No pt seed for pt$j"
			break
		fi

		# -l "${resolved_fastq}" \
		# -l "${annotatedir}"/pt.fq \
		#
		# ${_polap_cmd} assemble-rate \
		# 	"${option_data_type}" \
		# 	--plastid \
		# 	-o "${annotatedir}" \
		# 	-l "${resolved_fastq}" \
		# 	-i pt$i -j pt$j -w 1500

		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-i pt$i -j pt$j -w 1500

		if [[ -s "${annotatedir}/pt$j/assembly_graph.gfa" ]]; then
			_polap_lib_pt-extract-dna \
				"${annotatedir}/pt$j/assembly_graph.gfa" \
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

			_polap_log0 "ptDNA assembly: ${annotatedir}/pt.$j.gfa"
		else
			_polap_log0 "No pt assembly $j"
		fi

	done
	#
	# END: function readassemble-ont-pt-iterate_genus_species
	#######################################################################
}

# animal mtDNA
_polap_lib_readassemble-annotate-read-mt() {
	local type="${1:-pt}"

	local annotatedir="${_arg_outdir}/annotate-read-${type}"
	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	rm -rf "${annotatedir}"
	mkdir -p "${annotatedir}/at"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.fna
	# local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.3k.c80.fna

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

	_polap_log2 "    create the annotation table for mt and pt genes"
	_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-reads.R \
		--mt "${annotatedir}"/mt.paf \
		--pt "${annotatedir}"/pt.paf \
		--output "${annotatedir}" \
		--min-mapq ${_arg_annotate_read_min_mapq} \
		--min-identity ${_arg_annotate_read_min_identity}

}

# input: _arg_long_reads
_polap_lib_readassemble-assemble-annotated-read-mt() {

	local annotatedir="${_arg_outdir}"/annotate-read-pt

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	#######################################################################
	# BEGIN: function _run_polap_assemble-annotated-read
	#
	# for animal mitochondrial genome assembly
	Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${mt_table}" \
		--length 1e+6 \
		--output "${annotatedir}"/mt0.id.txt

	# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt.fq

	# subsample the data so that mt.fq is less than 100Mb
	rm -f "${annotatedir}"/mt0.fq
	_polap_lib_fastq-sample-to \
		"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "1m"

	# flye v2.9.6
	if _polap_lib_version-check_flye_version; then
		rm -rf "${annotatedir}"/mt
		_polap_log0 "flye 2.9.6 assembly on mtDNA"
		flye "${_arg_flye_data_type}" "${annotatedir}"/mt0.fq -t "${_arg_threads}" \
			--out-dir "${annotatedir}"/mt 2>"${_polap_output_dest}"
		if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
			_polap_log1 "output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
		else
			_polap_log0 "output: no MT assembly"
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
		_log_echo "No mt assembly 0"
	fi

	ln -s mt "${annotatedir}"/mt0

	local i
	for i in {0..1}; do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the mt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		_polap_lib_seed-plastid \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_log_echo "No mt seed for mt$j"
			break
		fi

		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-i mt$i -j mt$j -w 1500

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
}

################################################################################
#
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

	Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${pt_table}" \
		--length 3e+9 \
		--output "${annotatedir}"/pt0mt.id.txt
	_polap_log1 "pt0mt.id: $(cat "${annotatedir}"/pt0mt.id.txt | wc -l)"

	Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${mt_table}" \
		--length 1e+9 \
		--output "${annotatedir}"/mt0.id.txt
	_polap_log1 "mt0.id: $(cat "${annotatedir}"/mt0.id.txt | wc -l)"

	Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
		--table "${at_table}" \
		--length 1e+9 \
		--output "${annotatedir}"/at0.id.txt
	_polap_log1 "at0.id: $(cat "${annotatedir}"/at0.id.txt | wc -l)"

	# subtract pt.id.txt from at0.id.txt
	grep -vFf "${annotatedir}"/pt0mt.id.txt "${annotatedir}"/at0.id.txt \
		>"${annotatedir}"/at.id.txt
	_polap_log1 "at.id (at0 - pt0mt): $(cat "${annotatedir}"/at.id.txt | wc -l)"

	cat "${annotatedir}"/mt0.id.txt "${annotatedir}"/at.id.txt |
		sort | uniq >"${annotatedir}"/mt.id.txt
	_polap_log1 "mt.id (mt0 + at): $(cat "${annotatedir}"/mt.id.txt | wc -l)"

	# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
	seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt.id.txt >"${annotatedir}"/mt.fq

	# subsample the data so that mt.fq is less than 100Mb
	rm -f "${annotatedir}"/mt0.fq
	_polap_lib_fastq-sample-to \
		"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "100m"

	# flye v2.9.6
	if _polap_lib_version-check_flye_version; then
		rm -rf "${annotatedir}"/mt
		_polap_log0 "flye 2.9.6 assembly on mtDNA"
		flye "${_arg_flye_data_type}" "${annotatedir}"/mt0.fq -t "${_arg_threads}" \
			--out-dir "${annotatedir}"/mt 2>"${_polap_output_dest}"
		if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
			_polap_log1 "output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
		else
			_polap_log0 "output: no MT assembly"
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
		_polap_log3 "      polap command filter hifi data by ptDNA reference"
		_polap_lib_filter-reads-by-reference \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			--reference "${annotatedir_pt}"/pt.2.gfa
		resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	fi

	# repeat several times
	# polap command: annotate
	# polap command: seed-mito
	# polap command: assemble-rate
	# polap command: mt
	# polap command: bandage png
	local i
	for i in {0..5}; do
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
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				--reference "${annotatedir_pt}"/pt.2.gfa \
				-i mt$i -j mt$j -w 3000
		elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${annotatedir}"/mt.fq \
				-i mt$i -j mt$j -w 3000
		fi

		# polap command: mt

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
			_polap_log0 "No MT assembly $j"
		fi

	done
}

# 2025-08-06
# delete this.
_polap_lib_readassemble-annotate-nt2() {

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "  we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_arg_outdir}"/annotate-read-nt
	local annotatedir_pt="${_arg_outdir}"/annotate-read-pt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble mtDNA iteratively including noncoding seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble mtDNA iteratively including noncoding seed reads with ONT raw: ${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		_log_echo "found and use: ${annotatedir}/mt/assembly_graph.gfa"
		_log_echo "skip: assemble-annotated-read"
	else
		${_polap_cmd} assemble-annotated-read-nt \
			"${option_data_type}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		${_polap_cmd} mt \
			--infile "${annotatedir}/mt/assembly_graph.gfa" \
			-o "${annotatedir}/mt/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"

	else
		_log_echo "No MT assembly 0"
	fi

	ln -s mt "${annotatedir}"/mt0

	# TODO: remove pt reads from the input if not --plastid
	if [[ "${option_data_type}" == "--pacbio-hifi" ]]; then
		${_polap_cmd} filter reference hifi \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir_pt}"/pt.2.gfa
		local resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	fi

	local i
	for i in {0..5}; do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the MT contigs only
		${_polap_cmd} annotate \
			--quiet \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		${_polap_cmd} seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_log_echo "No MT seed for mt$j"
			break
		fi

		# -l "${resolved_fastq}" \
		# -l "${annotatedir}"/mt.fq \
		if [[ "${platform}" == "PACBIO_SMRT" ]]; then
			${_polap_cmd} assemble-rate \
				"${option_data_type}" \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				--reference "${annotatedir_pt}"/pt.2.gfa \
				-i mt$i -j mt$j -w 3000
		elif [[ "${platform}" == "ONT" ]]; then
			${_polap_cmd} assemble-rate \
				"${option_data_type}" \
				-o "${annotatedir}" \
				-l "${annotatedir}"/mt.fq \
				-i mt$i -j mt$j -w 3000
		fi

		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			${_polap_cmd} bandage png \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			echo mtDNA assembly: "${annotatedir}/mt.$j.gfa"
		else
			_log_echo "No MT assembly $j"
		fi

	done

}
