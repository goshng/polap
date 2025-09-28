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
# Subcommands for mtDNA or ptDNA sequences.
# 1. compare two sequences by a pairwise alignment.
# 2. fetch mtDNA or ptDNA sequences from the NCBI.
#
# Functions:
# function _run_polap_mafft-mtdna { -> we use this in v0.4.3
# function _run_polap_mauve-mtdna { -> we use this in v0.3.7.3
# function _run_polap_compare2ptdna { -> use the python code in disassemble
# function _run_polap_compare-mtdna { -> used by bioproject-postprocess
# function _run_polap_get-mtdna { -> used to download MT/PT DNA from the NCBI
# function _run_polap_get-dna-by-accession { -> not complicate, leave it
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
# mafft alignment
################################################################################
function _run_polap_mafft-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
Compare the known mitochondrial DNA (mtDNA) sequence with the newly assembled 
one for verification purposes using mafft.

Arguments:
  -a first.fa
  -b second.fa

Inputs:
  -a fasta: known mtDNA in fasta format
  -b fasta: another DNA sequence in fasta format

Outputs:
  summary of the pairwise alignment
  pident.txt

Example:
$(basename "$0") ${_arg_menu[0]} -a o/00-bioproject/2-mtdna.fasta -b o/1/assembly.fasta
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	local _mafft_dir="${_arg_outdir}"
	mkdir -p "${_mafft_dir}"
	cat "${_arg_short_read1}" >"${_mafft_dir}/in.fa"
	cat "${_arg_short_read2}" >>"${_mafft_dir}/in.fa"

	_polap_log3_pipe "mafft \
    --auto \
    ${_mafft_dir}/in.fa \
    >${_mafft_dir}/out.mafft \
    2>${_polap_output_dest}"

	_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-mafft.R \
    --input ${_mafft_dir}/out.mafft \
    --out ${_mafft_dir}/out.txt 2>${_polap_output_dest}"

	local pident_mafft=$(grep -oP '(?<=Percent Identity: )\S+' "${_mafft_dir}/out.txt")
	pident_mafft=${pident_mafft%\%}

	echo "$pident_mafft" >"${_mafft_dir}/pident.txt"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# mauve alignment
################################################################################
function _run_polap_mauve-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Compare the known mitochondrial DNA (mtDNA) sequence with the newly assembled 
# one for verification purposes.
#
# We utilize progressiveMauve to align two sequences and 
# calculate the total length of the LCB, which is then divided by the length 
# of a known mtDNA sequence. This approach provides a basic method for 
# comparing two sequences, albeit not highly sophisticated.
#
# Arguments:
#   -i ${_arg_inum}
# Inputs:
#   a.fasta: known mtDNA in fasta format
#   b.fasta: another DNA sequence in fasta format
# Outputs:
#   the ratio of LCB total length divided by the known mtDNA sequence
Example: $(basename "$0") ${_arg_menu[0]} -a o/00-bioproject/2-mtdna.fasta -b o/1/assembly.fasta
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log3_pipe "progressiveMauve \
    --output=${_arg_outdir}/mt.xmfa \
    ${_arg_short_read1} \
    ${_arg_short_read2} \
    >${_polap_output_dest} 2>&1"

	# awk 'NR > 1 {sum += $2 - $1} END {print sum}' "${_arg_outdir}/mt.xmfa.backbone" >"${_arg_outdir}/mt.identity.length.txt"
	awk 'NR > 1 {sum += ($1 > $2 ? $1 - $2 : $2 - $1)} END {print sum}' "${_arg_outdir}/mt.xmfa.backbone" >"${_arg_outdir}/mt.identity.length.txt"

	local _alen=$(<"${_arg_outdir}/mt.identity.length.txt")
	_polap_utility_get_contig_length \
		"${_arg_short_read1}" \
		"${_arg_outdir}/mt.reference.length.txt"
	local _blen=$(<"${_arg_outdir}/mt.reference.length.txt")
	local _percent_identity=$(echo "scale=3; ${_alen}/${_blen}" | bc)
	_polap_log0 "mauve_lcb_length_coverage: ${_percent_identity}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_compare2ptdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Compare the plastid DNA (ptDNA) sequence with the newly assembled 
# one for verification purposes.
#
# Inputs:
#   a.fasta: known mtDNA in fasta format
#   b.fasta: another DNA sequence in fasta format
# Outputs:
#   the ratio of LCB total length divided by the known mtDNA sequence
Example: $(basename "$0") ${_arg_menu[0]} -a o/00-bioproject/2-mtdna.fasta -b o/1/assembly.fasta
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	local _ptdir="${_arg_outdir}"
	_polap_log3_pipe "python ${_POLAPLIB_DIR}/polap-py-compare2ptdna.py \
    --seq1 ${_arg_short_read1} \
    --seq2 ${_arg_short_read2} \
		--out ${_ptdir} \
		2>${_polap_output_dest}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Compares the known mtDNA sequence and the assembled one.
################################################################################
function _run_polap_compare-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Compares the known mtDNA sequence and the assembled one.
#
# Arguments:
#   -i ${_arg_inum}
# Inputs:
#   ${_polap_var_compare_mtdna3}
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [${_arg_outdir}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	check_folder_existence "${_arg_outdir}"

	local n1=$(cut -f1 "${_polap_var_project_mtdna_fasta2_accession}")
	local l1=$(seqkit stats -T "${_polap_var_project_mtdna_fasta2}" |
		csvtk cut -t -f sum_len |
		csvtk del-header)

	# Run blastn between known mtDNA and assembled mtDNA (first round)
	blastn -query "${_polap_var_compare_mtdna1}" \
		-subject "${_polap_var_project_mtdna_fasta2}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand" \
		>"${_polap_var_compare_oga_blastn1}"

	blastn -query "${_polap_var_compare_mtdna1}" \
		-subject "${_polap_var_project_mtdna_fasta2}" \
		>"${_polap_var_compare_oga_blastn1}.full"

	_polap_log2_file "${_polap_var_compare_oga_blastn1}"

	# Determine the strand orientation
	if [ -s "${_polap_var_compare_oga_blastn1}" ]; then
		local _polap_var_project_strand=$(cut -f5 "${_polap_var_compare_oga_blastn1}" | head -n 1)
	else
		_polap_log1 "No hit in ${_polap_var_compare_oga_blastn1}"
		local n1=$(cut -f1 "${_polap_var_project_mtdna_fasta2_accession}")
		printf "%s\t%d\t0\t0\n" ${n1} ${l1} >"${_polap_var_compare_mtdna_compare}"
		return
	fi

	# Reverse sequence if the strand is negative
	if [ "${_polap_var_project_strand}" = "plus" ]; then
		cp "${_polap_var_compare_mtdna1}" "${_polap_var_compare_mtdna2}"
	else
		seqkit seq -t dna -v -p -r "${_polap_var_compare_mtdna1}" -o "${_polap_var_compare_mtdna2}"
	fi

	# Run blastn (second round)
	blastn -query "${_polap_var_compare_mtdna2}" -subject "${_polap_var_project_mtdna_fasta2}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand" \
		>"${_polap_var_compare_oga_blastn2}"

	_polap_log2_file "${_polap_var_compare_oga_blastn2}"

	# Restart sequence alignment at the lowest start position
	local _polap_var_project_restart_position=$(sort -n -k4 "${_polap_var_compare_oga_blastn2}" | head -n 1 | cut -f2)
	_polap_log2 "LOG: restart position: ${_polap_var_project_restart_position}"
	seqkit restart -i "${_polap_var_project_restart_position}" "${_polap_var_compare_mtdna2}" -o "${_polap_var_compare_mtdna3}"

	# Run blastn (third round)
	blastn -query "${_polap_var_compare_mtdna3}" -subject "${_polap_var_project_mtdna_fasta2}" \
		-outfmt "6 qseqid qstart sseqid sstart sstrand length pident" \
		>"${_polap_var_compare_oga_blastn3}"

	_polap_log2_file "${_polap_var_compare_oga_blastn3}"

	# Analyze the length of the match using an R script
	"${_POLAPLIB_DIR}"/polap-r-assemble-bioproject-3-length-match.R \
		"${_polap_var_compare_oga_blastn3}" \
		"${_polap_var_compare_oga_blastn3_length}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_compare_oga_blastn3_length}"

	local l2=$(<"${_polap_var_compare_oga_blastn3_length}")
	local c1=$(echo "scale=3; ${l2}/${l1}" | bc)
	_polap_log2 "Length of ${_polap_var_project_mtdna_fasta2}: ${l1}"
	_polap_log2 "Length of match alignment: ${l2}"
	_polap_log1 "length coverage: ${c1}"
	printf "%s\t%d\t%d\t%f\n" ${n1} ${l1} ${l2} ${c1} >"${_polap_var_compare_mtdna_compare}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

################################################################################
# Downloads the mtDNA sequence for a given species name.
################################################################################
function _run_polap_get-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local polap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - download the organelle-genome sequence in FASTA format from NCBI.

Synopsis:
  polap ${polap_cmd} [options]

Description:
  1. NCBI search for a mitochondrial genome:
    "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND <species>[Organism]"

  2. NCBI search for a plastid genome:
    "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND <species>[Organism]"

Options:
  --species
    "scientific name" (highest priority)

  -o ${_arg_outdir}
    (next priority)

Outputs:
  ${_polap_var_project_mtdna_fasta1}

  ${_polap_var_project_mtdna_fasta2}

  ${_polap_var_project_mtdna_fasta2_accession}

Examples:
  Get organelle genome sequences:
    polap ${polap_cmd} --species "Anthoceros agrestis"
    polap ${polap_cmd} without-genome --species "Anthoceros agrestis"
    polap ${polap_cmd} -o o
    polap ${polap_cmd} view

  Use preset species:
    add species: 'Trifolium pratense' to $HOME/.polap/profiles/tpra.yaml
    polap get-mtdna --plastid --preset tpra
    polap get-mtdna --plastid --preset tpra view

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

	if [[ ${_arg_menu[1]} == "file" ]]; then
		if [[ -s "${_polap_var_project_species}" ]]; then
			_polap_log0_file "${_polap_var_project_species}"
			_polap_log0_cat "${_polap_var_project_species}"
		else
			_polap_log0 "No species file."
		fi
		exit $EXIT_SUCCESS
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if ! run_check_ncbitools; then
		error_polap_conda
		exit $EXIT_ERROR
	fi

	if [[ ${_arg_menu[1]} == "view" ]]; then
		if [[ -s "${_polap_var_project_mtdna_fasta2_accession}" ]]; then
			_polap_log0_cat "${_polap_var_project_mtdna_fasta2_accession}"
			if [ -s "${_polap_var_project_mtdna_fasta1}" ]; then
				seqkit fx2tab -n -l "${_polap_var_project_mtdna_fasta1}" >&3
				seqkit stats "${_polap_var_project_mtdna_fasta1}" >&3
			else
				_polap_log1 "No such file: ${_polap_var_project_mtdna_fasta1}"
				_polap_log0 "No organelle genome sequence"
			fi
			if [ -s "${_polap_var_project_mtdna_fasta2}" ]; then
				seqkit stats "${_polap_var_project_mtdna_fasta2}" >&3
			else
				_polap_log1 "No such file: ${_polap_var_project_mtdna_fasta2}"
			fi

		else
			_polap_log0 "No result yet."
		fi
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_plastid}" == "off" ]]; then
		_polap_log1 "getting the mitochondrial organelle-genome sequence of a plant species ..."
	else
		_polap_log1 "getting the plastid organelle-genome sequence of a plant species ..."
	fi

	# Determine species name
	_polap_log1 "  step 1: determine the species name"
	local SPECIES=""
	if [ -n "${_arg_species}" ]; then
		_polap_log1 "  option --species: ${_arg_species}"
		SPECIES="${_arg_species}"
		_polap_log3_cmd mkdir -p "${_polap_var_project}"
	elif [ -s "${_polap_var_project_species}" ]; then
		if [[ -s "${_polap_var_project_txt}" ]]; then
			local bioproject_id=$(<"${_polap_var_project_txt}")
			local n=$(wc -l <"${_polap_var_project_species}")
			if [[ "${n}" -gt 1 ]]; then
				_polap_log0_file "${_polap_var_project_species}"
				_polap_log0_cat "${_polap_var_project_species}"
				die "ERROR: you have multiple species names in BioProject: ${bioproject_id}"
			fi
			SPECIES=$(<"${_polap_var_project_species}")
			_polap_log0 "  bioproject's species: $SPECIES"
		else
			die "ERROR: you have not run the menu: get-bioproject on the output folder [${_arg_outdir}], yet!"
		fi
	elif [ -n "${_arg_bioproject}" ]; then
		_run_polap_get-bioproject
		if [ -s "${_polap_var_project_species}" ]; then
			local bioproject_id=$(<"${_polap_var_project_txt}")
			local n=$(wc -l <"${_polap_var_project_species}")
			if [[ "${n}" -gt 1 ]]; then
				_polap_log0_file "${_polap_var_project_species}"
				_polap_log0_cat "${_polap_var_project_species}"
				die "ERROR: you have multiple species names in BioProject: ${bioproject_id}"
			fi
			SPECIES=$(<"${_polap_var_project_species}")
			_polap_log1 "  bioproject's species: $SPECIES"
		else
			die "  no species name is provided in the BioProject: ${_arg_bioproject}"
		fi
	else
		die "  use --species or --bioproject option."
	fi

	# Download the mitochondrial genome sequence for the given species
	_polap_log1 "  step 2: download the mitochondrial genome sequence for the given species"
	if [ "${_arg_plastid}" = "off" ]; then
		if [[ ${_arg_menu[1]} == "without-genome" ]]; then
			_polap_log1 "  downloading mitochondrial complete genomes of ${SPECIES} ..."
			esearch \
				-db nuccore \
				-query "(mitochondrion[Title] AND complete[Title]) AND ${SPECIES}[Organism]" |
				efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
		else
			_polap_log1 "  downloading mitochondrial complete sequence of ${SPECIES} ..."
			esearch \
				-db nuccore \
				-query "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
				efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
		fi
	else
		_polap_log1 "  downloading chloroplast complete genomes of ${SPECIES} ..."
		esearch \
			-db nuccore \
			-query "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
			efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
		if [[ ! -s "${_polap_var_project_mtdna_fasta1}" ]]; then
			esearch \
				-db nuccore \
				-query "(plastid[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
				efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
		fi
	fi
	_polap_log2_file "${_polap_var_project_mtdna_fasta1}"

	# Check if the fasta file was successfully downloaded
	_polap_log1 "  step 3: check if the fasta file was successfully downloaded"
	if [ -s "${_polap_var_project_mtdna_fasta1}" ]; then
		seqkit fx2tab --length --name --header-line \
			"${_polap_var_project_mtdna_fasta1}" \
			>"${_polap_var_project_mtdna_fasta1_stats}"

		local n=$(seqkit stats -T "${_polap_var_project_mtdna_fasta1}" |
			csvtk cut -t -f num_seqs |
			csvtk del-header)

		if [[ "${n}" -gt 1 ]]; then
			_polap_log2 "  multiple (${n}) accession numbers for the organelle genome sequences"
		fi

		seqkit head -n 1 "${_polap_var_project_mtdna_fasta1}" \
			-o "${_polap_var_project_mtdna_fasta2}"

		seqkit seq -ni "${_polap_var_project_mtdna_fasta2}" \
			-o "${_polap_var_project_mtdna_fasta2_accession}"

		local accession=$(<"${_polap_var_project_mtdna_fasta2_accession}")
		_polap_log2_cat "${_polap_var_project_mtdna_fasta1_stats}"
		_polap_log0 "Species: ${SPECIES}"
		_polap_log0 "NCBI accession: ${accession}"
	else
		echo "no mtDNA" >"${_polap_var_project_mtdna_fasta2_accession}"
		_polap_log0 "No organelle genome sequence found for the species: ${SPECIES}"
		rm -f \
			"${_polap_var_project_mtdna_fasta1_stats}" \
			"${_polap_var_project_mtdna_fasta1}" \
			"${_polap_var_project_mtdna_fasta2}"
	fi

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_get-dna-by-accession {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	if ! run_check_ncbitools; then
		error_polap_conda
		exit $EXIT_ERROR
	fi

	# Help message
	local help_message=$(
		cat <<HEREDOC
Obtain the DNA sequence by downloading it from NCBI in FASTA format 
using its corresponding accession ID.
Example: $(basename "$0") ${_arg_menu[0]} <accession> <output_file>
Example: $(basename "$0") ${_arg_menu[0]} NC_027276 ptdna-Podococcus_acaulis.fa
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0 "getting the DNA sequence using accesion ${_arg_menu[1]}"

	# Download the mitochondrial genome sequence for the given species
	esearch \
		-db nuccore \
		-query "${_arg_menu[1]}[ACCN]" |
		efetch -format fasta >"${_arg_menu[2]}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
