################################################################################
# Get BioProject information using BioProject accession.
#
#
################################################################################
function _run_polap_get-bioproject() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Get BioProject information using BioProject accession.
#
# Arguments:
#   --bioproject ${_arg_bioproject}: BioProject ID
#   -o $ODIR
# Inputs:
#   bioproject ID
# Outputs:
#   ${_polap_var_bioproject_runinfo}
Example: $(basename $0) ${_arg_menu[0]} --bioproject PRJNA574453
Example: $(basename $0) ${_arg_menu[0]} -o PRJDB8681-Dioscorea_cayenensis_subsp._rotundata --bioproject PRJDB8681 view 2>&1 | grep DRR196916
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Set variables for file paths
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'

	if [ -z "${_arg_bioproject}" ]; then
		_polap_echo0 "${help_message}"
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

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
		else
			_polap_log0 "No long-read info."
		fi
		if [ -s "${_polap_var_bioproject_sra_short_read}" ]; then
			_polap_log0_file "${_polap_var_bioproject_sra_short_read}"
			_polap_log0_cat "${_polap_var_bioproject_sra_short_read}"
			_polap_log0 ""
		else
			_polap_log0 "No short-read info."
		fi
		# _polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	# Clean and create the bioproject directory
	if [ -d "${_polap_var_bioproject}" ]; then
		_polap_log2 "  output folder: ${_polap_var_bioproject}"
	else
		_polap_log2 "  no output folder: ${_polap_var_bioproject} is created."
		mkdir -p "${_polap_var_bioproject}"
	fi

	# Set the BioProject ID
	_polap_log1 "BioProject: ${_arg_bioproject}"
	local BIOPRJ="${_arg_bioproject}"
	echo "${_arg_bioproject}" >"${_polap_var_bioproject_txt}"

	# Fetch the run information for the given BioProject
	if [ -s "${_polap_var_bioproject_runinfo_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "Note that we are using your previously downloaded NCBI runinfo."
		_polap_log0 "Use --redo if you download a new NCBI runinfo for the BioProject."
	else
		esearch -db bioproject -query "$BIOPRJ" |
			elink -target sra |
			efetch -format runinfo \
				>${_polap_var_bioproject_runinfo_all}
	fi

	check_file_existence "${_polap_var_bioproject_runinfo_all}"

	cat "${_polap_var_bioproject_runinfo_all}" |
		csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName |
		csvtk csv2tab >"${_polap_var_bioproject_runinfo}"

	# Log the output file
	_polap_log2_file "downloaded: ${_polap_var_bioproject_runinfo}"

	# Select SRA long and short read data
	"$script_dir/run-polap-assemble-bioproject-1-select-sra.R" \
		"${_polap_var_bioproject_runinfo}" \
		"${_polap_var_bioproject_sra_long_read}" \
		"${_polap_var_bioproject_sra_short_read}" \
		"${_polap_var_bioproject_species}" \
		"${_polap_var_bioproject_sra_per_species}" \
		2>"$_polap_output_dest"

	#!/bin/bash

	# Define the input file
	local input_file="${_polap_var_bioproject_sra_per_species}"

	# Iterate through each line of the file, skipping the header
	tail -n +2 "$input_file" | while IFS=$'\t' read -r ScientificName Run_Nano Run_Illumina; do

		_polap_log2 "Taxon: ${ScientificName}"

		# local taxon_id_folder="${ODIR}/0-bioproject/${taxon_id}"
		local bioproject_species_base="${ODIR}-$(echo $ScientificName | tr ' ' '_')"
		local bioproject_species_folder="${ODIR}-$(echo $ScientificName | tr ' ' '_')/0-bioproject"
		local taxon_id_folder="${ODIR}/0-bioproject/$(echo $ScientificName | tr ' ' '_')"

		# Create a directory with the taxon ID as the name
		mkdir -p "$taxon_id_folder"

		# Write the Run_Nano to 1-sra-long-read.tsv
		# echo "$Run_Nano" >"$taxon_id_folder/1-sra-long-read.tsv"
		grep "^${Run_Nano}" "${_polap_var_bioproject_runinfo}" >"$taxon_id_folder/1-sra-long-read.tsv"

		# Write the Run_Illumina to 1-sra-short-read.tsv
		# echo "$Run_Illumina" >"$taxon_id_folder/1-sra-short-read.tsv"
		grep "^${Run_Illumina}" "${_polap_var_bioproject_runinfo}" >"$taxon_id_folder/1-sra-short-read.tsv"

		# Inform the user
		_polap_log2 "Created files for $ScientificName with Taxon ID $taxon_id."

		# Create a new bioproject folder with the species name.
		mkdir -p "${bioproject_species_folder}"
		cp "${_polap_var_bioproject_txt}" "${bioproject_species_base}"
		cp "${_polap_var_bioproject}"/* "${bioproject_species_folder}"
		cp "$taxon_id_folder/1-sra-long-read.tsv" "${bioproject_species_folder}"
		cp "$taxon_id_folder/1-sra-short-read.tsv" "${bioproject_species_folder}"
	done

	if [ ! -s "${_polap_var_bioproject_sra_long_read}" ]; then
		_polap_log0 "No long-read data for BioProject: ${_arg_bioproject}"
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "${_polap_var_bioproject_sra_short_read}" ]; then
		_polap_log0 "No short-read data for BioProject: ${_arg_bioproject}"
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	# Check if the long-read dataset exists
	check_file_existence "${_polap_var_bioproject_sra_long_read}"
	check_file_existence "${_polap_var_bioproject_sra_short_read}"
	check_file_existence "${_polap_var_bioproject_species}"

	# Log the long and short read SRA files
	_polap_log2_file "output1: ${_polap_var_bioproject_sra_long_read}"
	_polap_log2_cat "${_polap_var_bioproject_sra_long_read}"
	_polap_log2_file "output2: ${_polap_var_bioproject_sra_short_read}"
	_polap_log2_cat "${_polap_var_bioproject_sra_short_read}"
	_polap_log2_file "output3: ${_polap_var_bioproject_species}"
	_polap_log2_cat "${_polap_var_bioproject_species}"

	# local taxon=$(cat "${_polap_var_bioproject_species}")
	#
	# esearch -db taxonomy -query "${taxon}" |
	# 	efetch -format docsum |
	# 	xtract -pattern DocumentSummary -element TaxId \
	# 		>"${_polap_var_bioproject_taxon_id}"
	#
	# _polap_log2_file "output4: ${_polap_var_bioproject_taxon_id}"
	# _polap_log2_cat "${_polap_var_bioproject_taxon_id}"
	#
	# local taxid=$(cat "${_polap_var_bioproject_taxon_id}")
	# datasets summary taxonomy taxon "${taxid}" \
	# 	>"${_polap_var_bioproject_taxonomy}" \
	# 	2>"$_polap_output_dest"
	#
	# _polap_log2_file "output5: ${_polap_var_bioproject_taxonomy}"
	# _polap_log2_cat "${_polap_var_bioproject_taxonomy}"

	_polap_log1 "LOG: $BIOPRJ - passed"
	touch "${_polap_var_bioproject_passed}"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
