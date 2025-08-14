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
Options:
  -l FASTQ
    long reads data file
  --plastid
    assembly mtDNA instead of mtDNA
  --animal
    assemble animal mtDNA
  --nano-raw [default]
  --pacbio-hifi
Examples:
  Assemble plant mitochondrial sequences:
    polap readassemble -l l.fq

  Assemble plastid sequences:
    polap readassemble -l l.fq --plastid
TODO:
  Dev.
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

	# Three cases
	# 1. plastid
	# 2. mitochondrial or mitochondrial + noncoding
	# 3. animial mitochondrial
	if [[ "${_arg_plastid}" == "on" ]]; then
		_polap_log1 "Read-assemble ptDNA"
		_polap_readassemble-pt
	elif [[ "${_arg_animal}" == "on" ]]; then
		_polap_log1 "Read-assemble animal mtDNA"
		_polap_readassemble-mt-animal
	else
		if [[ "${_arg_noncoding}" == "on" ]]; then
			_polap_log1 "Read-assemble plant mtDNA with mitochondrial noncoding regions"
			_polap_readassemble-nt
		else
			_polap_log1 "Read-assemble plant mtDNA without mitochondrial noncoding regions"
			# _polap_log0 "Before: _polap_readassemble-mt"
			# _polap_log0 "  -o ${_arg_outdir}"
			# _polap_log0 "  -i ${_arg_inum}"
			# _polap_log0 "  -j ${_arg_jnum}"
			_polap_readassemble-mt
			# _polap_log0 "After: _polap_readassemble-mt"
			# _polap_log0 "  -o ${_arg_outdir}"
			# _polap_log0 "  -i ${_arg_inum}"
			# _polap_log0 "  -j ${_arg_jnum}"

		fi
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

_polap_readassemble-nt() {
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		_polap_log0 "Not Yet implemented!"
		return
	fi
	_polap_lib_readassemble-annotate-read-nt
	_polap_lib_readassemble-assemble-annotated-read-nt
}

_polap_readassemble-pt() {
	_polap_lib_readassemble-annotate-read-pt
	_polap_lib_readassemble-assemble-annotated-read-pt
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

	local pt_table mt_table at_table all_table
	_polap_lib_readassemble-common-variables \
		annotatedir pt_table mt_table at_table all_table

	# Step 1
	# assemble pt first
	# we delete annotate-read-TYPE dir inside this function
	#
	_polap_lib_readassemble-annotate-read-pt mt
	_polap_lib_readassemble-assemble-annotated-read-pt mt

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
	_polap_lib_filter-reads-by-reference \
		-o "${annotatedir}" \
		-l "${_arg_long_reads}" \
		--reference "${_arg_outdir}/${type}-pt.0.gfa"

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
	#
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
