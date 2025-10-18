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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_polap_readassemble {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap readassemble - annotate reads before organelle genome assembly
Synopsis:
  polap readassemble [options]
Description:
  polap readassemble uses organelle genes to annotate reads before organelle genome assembly.
We have four-case workflows;
1) ptDNA from ONT: select plastid‑origin reads using protein‑to‑genome alignment (e.g., NCBI BLAST) and assemble with Flye v2.9.6. This leverages ONT read length to bridge repeats.
2) mtDNA from ONT: aggressively filter ptDNA/nuclear reads (e.g., protein markers; BUSCO to identify nuclear reads), bootstrap seed contigs with miniasm, then complete with Flye.
3) ptDNA from HiFi: either use the protein‑guided selection + Flye route (method 1) or a wrapper of Oatk with organelle gene–aware pathfinding.
4) mtDNA from HiFi: a wrapper of Oatk is generally preferred although (method 2).
We assemble ptDNA first before mtDNA assembly.

Options:
  -l FASTQ
    long reads data file
  --plastid
    assembly mtDNA instead of mtDNA
  --animal
    assemble animal mtDNA
  --nano-raw [default]
  --pacbio-hifi
  --use-oatk
Examples:
  Assemble plant mitochondrial sequences from ONT l.fq to produce l.mt.gfa:
    polap readassemble -l l.fq

  Assemble plant mitochondrial sequences from HiFi l.fq to produce l.mt.gfa:
    polap readassemble -l l.fq --pacbio-hifi

  Assemble plastid sequences from ONT l.fq to produce l.pt.fa and l.pt.gfa:
    polap readassemble -l l.fq --plastid

  Assemble plastid sequences from HiFi l.fq to produce l.pt.fa and l.pt.gfa:
    polap readassemble -l l.fq --plastid

  Oatk assembles plastid sequences from HiFi l.fq to produce l.pt.fa and l.pt.gfa:
    polap readassemble -l l.fq --plastid --use-oatk

  Oatk assembles mitochondrial sequences from HiFi l.fq to produce l.pt.fa and l.pt.gfa:
    polap readassemble -l l.fq --use-oatk

  (not tested) Oatk assembles animal mitochondrial sequences from HiFi l.fq to produce l.pt.fa and l.pt.gfa:
    polap readassemble -l l.fq --animal --use-oatk

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)
Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "annotated" ]]; then
		_polap_lib_conda-ensure_conda_env polap || exit 1
		_polap_lib_readassemble-annotated mtseed
		conda deactivate
		return 0
	fi

	if [[ "${_arg_use_oatk}" == "off" ]]; then
		_polap_lib_conda-ensure_conda_env polap || exit 1

		if [[ "${_arg_plastid}" == "on" ]]; then
			_polap_log0 "Assemble the plastid genome sequence from ${_arg_long_reads}"
			_polap_readassemble-pt
			local pt_fa="${_arg_long_reads%.*}.pt.fa"
			local pt_gfa="${_arg_long_reads%.*}.pt.gfa"
			local pt_png="${_arg_long_reads%.*}.pt.png"
			cp -p "${_arg_outdir}/pt.1.fa" "${pt_fa}"
			cp -p "${_arg_outdir}/pt.1.gfa" "${pt_gfa}"
			cp -p "${_arg_outdir}/pt.1.png" "${pt_png}"
			_polap_log0 "output fasta: ${pt_fa}"
			_polap_log0 "output assembly graph: ${pt_gfa}"
			_polap_log0 "output assembly graph figure: ${pt_png}"
		elif [[ "${_arg_animal}" == "on" ]]; then
			_polap_log1 "Read-assemble animal mtDNA (not tested yet)"
			_polap_readassemble-mt-animal
		else
			if [[ "${_arg_noncoding}" == "on" ]]; then
				local mt_gfa="${_arg_long_reads%.*}.mt.gfa"
				local mt_png="${_arg_long_reads%.*}.mt.png"
				_polap_log0 "Assemble the mitochondrial genome sequence from ${_arg_long_reads} in ${mt_gfa}"
				_polap_readassemble-nt
				cp -p "${_arg_outdir}/nt.1.gfa" "${mt_gfa}"
				cp -p "${_arg_outdir}/nt.1.png" "${mt_png}"
				_polap_log0 "output assembly graph: ${mt_gfa}"
				_polap_log0 "output assembly graph figure: ${mt_png}"
			else
				# mtDNA using ONT data
				_polap_log0 "Use miniasm to generate seed contigs"
				# _polap_log0 "Use miniasm to generate seed contigs ${_arg_verbose_str}"
				_polap_log1 "Read-assemble plant mtDNA without mitochondrial noncoding regions"
				_arg_plastid="on"
				if [[ "${_arg_readassemble_pt}" == "on" ]]; then
					_polap_readassemble-pt
				else
					_arg_long_reads="${_arg_outdir}/ld.fq"
				fi

				# use the original long-read input data
				# if [[ "${_arg_readassemble_use_all_long_read}" == "on" ]]; then
				# 	_arg_long_reads="${_arg_long_reads_original}"
				# fi
				_arg_long_reads="${_arg_long_reads_original}"

				_arg_plastid="off"
				_arg_menu[1]="annotate"
				if [[ "${_arg_readassemble_mtseed}" == "on" ]]; then
					_run_polap_mtseed
				fi

				_arg_menu[1]="fast"
				_arg_pt_ref="${_arg_outdir}/pt.1.fa"
				_polap_assert '[[ -s "${_arg_pt_ref}" ]]' \
					"pt ref must exist, '${_arg_pt_ref}'"
				# _arg_steps_include="1-9"
				# _arg_steps_include="3"
				if [[ -d "${_arg_outdir}/mtseed/mt1" ]]; then
					rm -rf "${_arg_outdir}/mtseed/mt1"
					rm -rf "${_arg_outdir}/mtseed/mt2"
					rm -rf "${_arg_outdir}/mtseed/mt3"
				fi
				_run_polap_mtseed
				local mt_gfa="${_arg_long_reads%.*}.mt.gfa"
				cp -p "${_arg_outdir}/mt.1.gfa" "${mt_gfa}"
				_polap_log0 "output assembly graph: ${mt_gfa}"
			fi
		fi
	else
		_polap_lib_conda-ensure_conda_env polap || exit 1

		_polap_assert '[[ "${_arg_data_type}" == "pacbio-hifi" ]]' \
			"_arg_data_type must be pacbio-hifi, got '${_arg_data_type}'"
	fi
	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

_polap_readassemble-nt() {
	# if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
	# 	_polap_log0 "Not Yet implemented!"
	# 	return
	# fi
	if [[ -s "${_arg_outdir}/pt-pt.1.gfa" ]]; then
		_polap_readassemble-pt
	fi
	_polap_lib_readassemble-annotate-read-nt
	_polap_lib_readassemble-assemble-annotated-read-nt
}

_polap_readassemble-pt-v1() {
	# downsampling
	# number of bases
	_polap_lib_readassemble-annotate-read-pt
	_polap_lib_readassemble-assemble-annotated-read-pt
}

_polap_readassemble-pt() {
	# downsampling
	# number of bases
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		if [[ "${_arg_reduction_reads}" == "on" ]]; then
			mkdir -p "${_arg_outdir}/genomesize"
			_polap_lib_genomesize-estimate -l "${_arg_long_reads}" \
				-o "${_arg_outdir}/genomesize"
			local genomesize=$(<"${_arg_outdir}/genomesize/genome_size.txt")

			if [[ "${_arg_downsample}" =~ [gGkKmM]$ ]]; then
				_polap_log0 "use --downsample 3 for HiFi reads"
				return
			fi
			_polap_lib_fastq-sample-to-coverage \
				"${_arg_long_reads}" \
				"${_arg_outdir}/ld.fq" \
				"${_arg_downsample}" \
				"${genomesize}"
		else
			ln -s "$(realpath "${_arg_long_reads}")" "${_arg_outdir}/ld.fq"
		fi
	elif [[ "${_arg_data_type}" == "pacbio-raw" ]]; then
		seqkit seq -m 3000 "${_arg_long_reads}" -o "${_arg_outdir}/ld.fq"
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		if [[ ! "${_arg_downsample}" =~ [gGkKmM]$ ]]; then
			_polap_log0 "use --downsample 1g for ONT reads"
			return
		fi
		_polap_log0 _polap_lib_fastq-sample-to \
			"${_arg_long_reads}" \
			"${_arg_outdir}/ld.fq" \
			"${_arg_downsample}"
		_polap_lib_fastq-sample-to \
			"${_arg_long_reads}" \
			"${_arg_outdir}/ld.fq" \
			"${_arg_downsample}"
	fi

	_arg_long_reads="${_arg_outdir}/ld.fq"

	_polap_lib_readassemble-annotate-read-pt
	_polap_lib_readassemble-assemble-annotated-read-pt
	_polap_lib_file-cleanup -d "${_arg_outdir}/annotate-read-pt" -s 5M -a rm
}

# 2025-08-13
# I want to add a version for hifi100k.sh and ont100k.sh
# I leave this as is.
_polap_readassemble-mt-v1() {
	_polap_lib_readassemble-annotate-read-mt
	_polap_lib_readassemble-assemble-annotated-read-mt
	# if not assemble
}

# We use hifi100k and ont100k.
# Tested only for hifi100k
_polap_readassemble-mt() {
	local type="mt"
	local annotatedir="${_arg_outdir}/annotate-read-${type}"

	rm -rf "$annotatedir"

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	local hifi1dir="$HOME/all/polap/hifi1"
	local speciesdir="${_arg_outdir%-*}"

	# Step 1
	# assemble pt first
	# we delete annotate-read-TYPE dir inside this function
	#
	if [[ ! -s "${_arg_outdir}/${type}-pt.0.fa" ]]; then
		if [[ -s "${_arg_outdir}/pt-pt.0.gfa" ]]; then
			ln -fs "pt-pt.0.gfa" "${_arg_outdir}/${type}-pt.0.gfa"
			ln -fs "pt-pt.0.fa" "${_arg_outdir}/${type}-pt.0.fa"
		fi
		if [[ -s "${_arg_outdir}/annotate-read-pt/pt.0.fa" ]]; then
			ln -fs "annotate-read-pt/pt.0.gfa" "${_arg_outdir}/${type}-pt.0.gfa"
			ln -fs "annotate-read-pt/pt.0.fa" "${_arg_outdir}/${type}-pt.0.fa"
		fi
		if [[ -s "$hifi1dir/${speciesdir}/t5/0/polap-readassemble-1-pt/pt.0.fa" ]]; then
			ln -fs $(realpath "${hifi1dir}/${speciesdir}/t5/0/polap-readassemble-1-pt/pt.0.gfa") \
				"${_arg_outdir}/${type}-pt.0.gfa"
			ln -fs $(realpath "${hifi1dir}/${speciesdir}/t5/0/polap-readassemble-1-pt/pt.0.fa") \
				"${_arg_outdir}/${type}-pt.0.fa"
		fi
	fi
	# ln -fs "pt-pt.0.fa" "${_arg_outdir}/${type}-pt.0.fa"

	if [[ ! -s "${_arg_outdir}/${type}-pt.0.gfa" ]]; then
		_polap_lib_readassemble-annotate-read-pt mt
		_polap_lib_readassemble-assemble-annotated-read-pt mt
	fi

	# _polap_lib_lines-skip1 "${mt_table}" | cut -f1 >"${annotatedir}"/mt0.id.txt
	# rm -f "${annotatedir}"/mt.fq
	# seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt.fq

	# Step 2
	# filter out ptDNA-origin reads
	# input: the long-read data
	# output: ${annotatedir}/kmer/ref-filtered.fastq
	#
	# -l "${_arg_long_reads}" \
	# -l "${annotatedir}/mt.fq" \
	#

	# "${annotatedir}/kmer/ref-filtered.fastq.gz" |
	if [[ "${_arg_data_type}" == "nano-raw" ]]; then
		mkdir -p "${annotatedir}/kmer"

		minimap2 -t 8 -x map-ont --secondary=no \
			"${_arg_outdir}/${type}-pt.0.fa" \
			"${_arg_long_reads}" >"${annotatedir}/kmer/ptdna-origin.paf"
		awk -v MINLEN=1000 -v MINID=0.3 \
			-f "${_POLAPLIB_DIR}/polap-awk-filter-paf-by-identity-len.awk" \
			"${annotatedir}/kmer/ptdna-origin.paf" \
			>"${annotatedir}/kmer/ptdna-origin.txt"

		seqkit grep -v -f "${annotatedir}/kmer/ptdna-origin.txt" \
			"${_arg_long_reads}" -o "${annotatedir}/kmer/ref-filtered.fastq.gz"
	elif [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		_polap_lib_filter-reads-by-reference \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			--reference "${_arg_outdir}/${type}-pt.0.gfa"

		# ln -sf "ref-filtered.fastq.gz" \
		# 	"${annotatedir}/kmer/ref-filtered2.fastq.gz"
	fi

	# Step 2-1
	# select reads with more MT genes than PT genes.

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	# Step 3
	# create 100k
	#
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		bash "${_POLAPLIB_DIR}/polap-bash-hifi100k.sh" \
			-r "${annotatedir}/kmer/ref-filtered.fastq.gz" \
			-g "${MTAA}" \
			-o "${annotatedir}/kmer/100k" \
			-t "${_arg_readassemble_t}" \
			-N "${_arg_readassemble_n}" \
			-T "${_arg_threads}"
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# echo "nano-raw"
		bash "${_POLAPLIB_DIR}/polap-bash-ont100k.sh" \
			-r "${annotatedir}/kmer/ref-filtered.fastq.gz" \
			-g "${MTAA}" \
			-o "${annotatedir}/kmer/100k" \
			-t "${_arg_readassemble_t}" \
			-N "${_arg_readassemble_n}" \
			-T "${_arg_threads}"
	else
		_polap_log0 "ERROR: no such data type available: ${_arg_data_type}"
		return
	fi

	# Step 4
	# Convert the seed contig fasta to gfa and mt.contig.name files.
	#
	mkdir -p "${annotatedir}/mt/30-contigger"
	bash "${_POLAPLIB_DIR}/polap-bash-fa2gfa.sh" \
		-o "${annotatedir}/mt/30-contigger/graph_final.gfa" \
		"${annotatedir}/kmer/100k/greedy_100k.fasta"
	#
	bash "${_POLAPLIB_DIR}/polap-bash-fa2mtcontigname.sh" \
		-o "${annotatedir}/mt/mt.contig.name-mt" \
		"${annotatedir}/kmer/100k/greedy_100k.fasta"

	# head -n "${_arg_readassemble_n}" \
	# 	"${annotatedir}/mt/mt.contig.name-mt" \
	# 	>"${annotatedir}/mt/mt.contig.name-mt0"

	gfatools gfa2fa \
		"${annotatedir}/mt/30-contigger/graph_final.gfa" \
		>"${annotatedir}/mt/30-contigger/graph_final.fasta" \
		2>${_polap_output_dest}

	# Step 5
	# select unitigs with more MT vs. PT genes.
	# FIXME: use DNA sequences not AA sequences.
	#
	# if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
	# 	_polap_lib_annotate \
	# 		-o "${annotatedir}" \
	# 		-i "mt"
	# elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
	# 	_arg_long_reads="${annotatedir}/mt/30-contigger/graph_final.fasta"
	# 	_polap_lib_readassemble-annotate-read-pt xt
	# fi
	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i "mt"

	_polap_lib_lines-skip1 "${annotatedir}/mt/contig-annotation-depth-table.txt" |
		cut -d' ' -f1 >"${annotatedir}/mt/mt.contig.name-mt0"

	# Step 5
	_arg_long_reads="${annotatedir}/kmer/ref-filtered.fastq.gz"
	_polap_lib_readassemble-assemble-annotated-read-mt-100k
}

_polap_readassemble-mt-animal() {
	_polap_log0 "Not Yet implemented!"
}
