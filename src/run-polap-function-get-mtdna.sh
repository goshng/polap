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

source "$script_dir/run-polap-function-utilities.sh"

################################################################################
# Downloads the mtDNA sequence for a given species name.
################################################################################
function _run_polap_get-mtdna() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-bioproject.sh"

	if ! run_check_ncbitools; then
		error_polap_conda
		exit $EXIT_ERROR
	fi

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Download the organelle-genome sequence in FASTA format from NCBI.
#
# Note: A downloaded fasta file may have multiple sequences.
#
# 1. NCBI search for a mitochondrial genome:
#   "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]"
# 2. NCBI search for a plastid genome:
#   "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]"
#
# Arguments:
#   --species "scientific name" (highest priority)
#   or
#   -o ${ODIR} (next priority)
# Outputs:
#   ${_polap_var_bioproject_mtdna_fasta1}
#   ${_polap_var_bioproject_mtdna_fasta2}
#   ${_polap_var_bioproject_mtdna_fasta2_accession}
# Preconditions:
#   get-bioproject
Example: $(basename $0) ${_arg_menu[0]} --species "Anthoceros agrestis"
Example: $(basename $0) ${_arg_menu[0]} -o o
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if [[ ${_arg_menu[1]} == "view" ]]; then
		if [[ -s "${_polap_var_bioproject_mtdna_fasta2_accession}" ]]; then
			_polap_log0_cat "${_polap_var_bioproject_mtdna_fasta2_accession}"
			if [ -s "${_polap_var_bioproject_mtdna_fasta1}" ]; then
				seqkit stats "${_polap_var_bioproject_mtdna_fasta1}" >&2
			else
				_polap_log1 "No such file: ${_polap_var_bioproject_mtdna_fasta1}"
				_polap_log0 "No organelle genome sequence"
			fi
			if [ -s "${_polap_var_bioproject_mtdna_fasta2}" ]; then
				seqkit stats "${_polap_var_bioproject_mtdna_fasta2}" >&2
			else
				_polap_log1 "No such file: ${_polap_var_bioproject_mtdna_fasta2}"
			fi
		else
			_polap_log0 "No result yet."
		fi
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "getting the organelle-genome sequence of a plant species ..."

	# Determine species name
	_polap_log1 "  step 1: determine the species name"
	local SPECIES=""
	if [ -n "${_arg_species}" ]; then
		_polap_log1 "  option --species: ${_arg_species}"
		SPECIES="${_arg_species}"
		_polap_log3_cmd mkdir -p "${_polap_var_bioproject}"
	elif [ -s "${_polap_var_bioproject_species}" ]; then
		if [[ -s "${_polap_var_bioproject_txt}" ]]; then
			local bioproject_id=$(<"${_polap_var_bioproject_txt}")
			local n=$(wc -l <"${_polap_var_bioproject_species}")
			if [[ "${n}" -gt 1 ]]; then
				_polap_log0_file "${_polap_var_bioproject_species}"
				_polap_log0_cat "${_polap_var_bioproject_species}"
				die "ERROR: you have multiple species names in BioProject: ${bioproject_id}"
			fi
			SPECIES=$(<"${_polap_var_bioproject_species}")
			_polap_log0 "  bioproject's species: $SPECIES"
		else
			die "ERROR: you have not run the menu: get-bioproject on the output folder [${ODIR}], yet!"
		fi
	elif [ -n "${_arg_bioproject}" ]; then
		_run_polap_get-bioproject
		if [ -s "${_polap_var_bioproject_species}" ]; then
			local bioproject_id=$(<"${_polap_var_bioproject_txt}")
			local n=$(wc -l <"${_polap_var_bioproject_species}")
			if [[ "${n}" -gt 1 ]]; then
				_polap_log0_file "${_polap_var_bioproject_species}"
				_polap_log0_cat "${_polap_var_bioproject_species}"
				die "ERROR: you have multiple species names in BioProject: ${bioproject_id}"
			fi
			SPECIES=$(<"${_polap_var_bioproject_species}")
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
		_polap_log1 "  downloading mitochondrial complete genomes of ${SPECIES} ..."
		esearch \
			-db nuccore \
			-query "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
			efetch -format fasta >"${_polap_var_bioproject_mtdna_fasta1}"
	else
		_polap_log1 "  downloading chloroplast complete genomes of ${SPECIES} ..."
		esearch \
			-db nuccore \
			-query "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
			efetch -format fasta >"${_polap_var_bioproject_mtdna_fasta1}"
	fi
	_polap_log2_file "${_polap_var_bioproject_mtdna_fasta1}"

	# Check if the fasta file was successfully downloaded
	_polap_log1 "  step 3: check if the fasta file was successfully downloaded"
	if [ -s "${_polap_var_bioproject_mtdna_fasta1}" ]; then
		seqkit fx2tab --length --name --header-line \
			"${_polap_var_bioproject_mtdna_fasta1}" \
			>"${_polap_var_bioproject_mtdna_fasta1_stats}"

		local n=$(seqkit stats -T "${_polap_var_bioproject_mtdna_fasta1}" |
			csvtk cut -t -f num_seqs |
			csvtk del-header)

		if [[ "${n}" -gt 1 ]]; then
			_polap_log2 "  multiple (${n}) accession numbers for the organelle genome sequences"
		fi

		seqkit head -n 1 "${_polap_var_bioproject_mtdna_fasta1}" \
			-o "${_polap_var_bioproject_mtdna_fasta2}"

		seqkit seq -ni "${_polap_var_bioproject_mtdna_fasta2}" \
			-o "${_polap_var_bioproject_mtdna_fasta2_accession}"

		local accession=$(<"${_polap_var_bioproject_mtdna_fasta2_accession}")
		_polap_log2_cat "${_polap_var_bioproject_mtdna_fasta1_stats}"
		_polap_log0 "Species: ${SPECIES}"
		_polap_log0 "NCBI accession: ${accession}"
	else
		echo "no mtDNA" >"${_polap_var_bioproject_mtdna_fasta2_accession}"
		_polap_log0 "No mtDNA sequence found for the species: ${SPECIES}"
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_gm() {
	_run_polap_get-mtdna

}
