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

source "$script_dir/polap-constants.sh"

################################################################################
# Get BioProject information using BioProject accession.
################################################################################
function _run_polap_get-bioproject() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-bioproject.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Get BioProject information using BioProject accession.
#
# FIXME: PRJNA841235 no Illumina, PacBIO data: error run-polap-get-bioproject-1.R
#
# Arguments:
#   --bioproject <BioProject ID>
# Inputs:
#   bioproject ID
# Outputs:
#   ${_polap_var_bioproject_runinfo}
# Usage:
#   PRJNA557253 - multiple species
Example: $(basename $0) ${_arg_menu[0]} --bioproject PRJNA574453
Example: $(basename $0) ${_arg_menu[0]} view -o PRJNA5744532 2>&1 | grep DRR196916
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_bioproject_runinfo}" ]; then
			_polap_log0_file "${_polap_var_bioproject_runinfo}"
			_polap_log0_cat "${_polap_var_bioproject_runinfo}"
			_polap_log0 ""
		else
			_polap_log0 "No BioProject info."
		fi
		if [ -s "${_polap_var_bioproject_sra_long_read}" ]; then
			_polap_log0_file "${_polap_var_bioproject_sra_long_read}"
			_polap_log0_cat "${_polap_var_bioproject_sra_long_read}"
			_polap_log0 ""
			SRA=$(cut -f1 "${_polap_var_bioproject_sra_long_read}")
			echo $SRA >"${ODIR}/sra.txt"
		else
			_polap_log0 "No long-read info."
		fi
		if [ -s "${_polap_var_bioproject_sra_short_read}" ]; then
			_polap_log0_file "${_polap_var_bioproject_sra_short_read}"
			_polap_log0_cat "${_polap_var_bioproject_sra_short_read}"
			_polap_log0 ""
			SRA=$(cut -f1 "${_polap_var_bioproject_sra_short_read}")
			echo $SRA >>"${ODIR}/sra.txt"
		else
			_polap_log0 "No short-read info."
		fi

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	# Check if --bioproject option is provided.
	if [ -z "${_arg_bioproject}" ]; then
		_polap_log0 "ERROR: --bioproject option is required to download BioProject runinfo."
		_polap_log0 "       you might want to view the BioProject runinfo."
		_polap_log0 "       $0 ${_arg_menu[0]} view -o ${ODIR}"
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	else
		local BIOPRJ="${_arg_bioproject}"
	fi

	# Set the BioProject ID
	_polap_log0 "BioProject: ${BIOPRJ}"
	echo "${BIOPRJ}" >"${_polap_var_bioproject_txt}"

	# Delete the 0-bioproject folder.
	if [[ -d "${_polap_var_bioproject}" ]]; then
		_polap_log3_cmd rm -rf "${_polap_var_bioproject}"
	fi

	# Fetch the run information for the given BioProject
	if [[ ! -s "${_polap_var_bioproject_runinfo_all}" ]] || [[ "${_arg_redo}" = "on" ]]; then
		_polap_log0 "  downloading ${_polap_var_bioproject_runinfo_all} ..."
		esearch -db bioproject -query "$BIOPRJ" |
			elink -target sra |
			efetch -format runinfo \
				>${_polap_var_bioproject_runinfo_all}
	else
		_polap_log0 "  found: ${_polap_var_bioproject_runinfo_all}, so skipping downloading the runinfo ..."
	fi

	# create the 0-bioproject folder
	_polap_log1 "  creating output folder: ${_polap_var_bioproject}"
	_polap_log3_cmd mkdir -p "${_polap_var_bioproject}"

	check_file_existence "${_polap_var_bioproject_runinfo_all}"

	# Log the output file
	_polap_log1 "  extracting some columns in TSV: ${_polap_var_bioproject_runinfo}"
	_polap_log3_pipe "cat ${_polap_var_bioproject_runinfo_all} |\
    csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName |\
    csvtk csv2tab >${_polap_var_bioproject_runinfo}"

	# Select SRA long and short read data
	_polap_log1 "  run-polap-get-bioproject-1.R"
	_polap_log2 "    input1: ${_polap_var_bioproject_runinfo}"
	_polap_log2 "    output1: ${_polap_var_bioproject_sra_long_read}"
	_polap_log2 "    output2: ${_polap_var_bioproject_sra_short_read}"
	_polap_log2 "    output3: ${_polap_var_bioproject_species}"
	_polap_log2 "    output4: ${_polap_var_bioproject_sra_per_species}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-get-bioproject-1.R \
		${_polap_var_bioproject_runinfo} \
		${_polap_var_bioproject_sra_long_read} \
		${_polap_var_bioproject_sra_short_read} \
		${_polap_var_bioproject_species} \
		${_polap_var_bioproject_sra_per_species} \
		2>$_polap_output_dest"

	# Iterate through each line of the file, skipping the header
	_polap_log1 "  creating output folders per species ..."
	local input_file="${_polap_var_bioproject_sra_per_species}"
	tail -n +2 "$input_file" | while IFS=$'\t' read -r ScientificName Run_Nano Run_Illumina; do
		_polap_log2 "    taxon: ${ScientificName}"

		local bioproject_species_base="${ODIR}-$(echo $ScientificName | tr ' ' '_')"
		local bioproject_species_folder="${ODIR}-$(echo $ScientificName | tr ' ' '_')/0-bioproject"
		local taxon_id_folder="${ODIR}/0-bioproject/$(echo $ScientificName | tr ' ' '_')"

		# Create a directory with the taxon ID as the name
		_polap_log3_cmd mkdir -p "$taxon_id_folder"
		# grep "^${Run_Nano}" "${_polap_var_bioproject_runinfo}" >"$taxon_id_folder/1-sra-long-read.tsv"
		# grep "^${Run_Illumina}" "${_polap_var_bioproject_runinfo}" >"$taxon_id_folder/1-sra-short-read.tsv"
		_polap_log3_pipe "grep ^${Run_Nano} ${_polap_var_bioproject_runinfo} >$taxon_id_folder/1-sra-long-read.tsv"
		_polap_log3_pipe "grep ^${Run_Illumina} ${_polap_var_bioproject_runinfo} >$taxon_id_folder/1-sra-short-read.tsv"

		_polap_log0 "    creates a taxon directory: ${bioproject_species_base}"
		_polap_log3_cmd mkdir -p "${bioproject_species_folder}"
		_polap_log3_cmd touch "${bioproject_species_base}/log-need-to-fetch-data.txt"
		_polap_log3_cmd cp "${_polap_var_bioproject_txt}" "${bioproject_species_base}"
		_polap_log3_cmd cp "${_polap_var_bioproject}"/*.{tsv,txt} "${bioproject_species_folder}"
		_polap_log3_cmd cp "${taxon_id_folder}"/*.tsv "${bioproject_species_folder}"
		# _polap_log3_cmd cp "$taxon_id_folder/1-sra-long-read.tsv" "${bioproject_species_folder}"
		# _polap_log3_cmd cp "$taxon_id_folder/1-sra-short-read.tsv" "${bioproject_species_folder}"
		_polap_log3_pipe "echo ${ScientificName} >${bioproject_species_folder}/1-species.txt"
	done

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
