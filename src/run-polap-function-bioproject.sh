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
set +u
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

source "$script_dir/polap-constants.sh"

################################################################################
# Gets BioProject ID of a given SRA ID.
#
# Arguments:
#   --sra ${_arg_sra}: SRA ID
#   -o ${ODIR}: output folder
# Inputs:
#   SRA ID
# Outputs:
#   BioProject ID
#   FILE: ${_polap_var_bioproject_txt}
################################################################################
function _run_polap_get-bioproject-sra() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Gets BioProject ID of a given SRA ID.
#
# Arguments:
#   --sra ${_arg_sra}: SRA ID
#   -o ${ODIR}: output folder
# Inputs:
#   SRA ID
# Outputs:
#   BioProject ID
#   FILE: ${_polap_var_bioproject_txt}
Example: $(basename $0) ${_arg_menu[0]} --sra <SRA ID>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if [[ ${_arg_menu[1]} == "view" ]]; then
		if [[ -s "${_polap_var_bioproject_txt}" ]]; then
			_polap_log0_cat "${_polap_var_bioproject_txt}"
		else
			_polap_log0 "No such file: ${_polap_var_bioproject_txt}"
		fi
		exit $EXIT_SUCCESS
	fi

	if [[ ! -d "${ODIR}" ]]; then
		_polap_log1 "  no output folder, creating ${ODIR}"
		mkdir -p "${ODIR}"
	fi

	esearch -db sra -query "${_arg_sra}" |
		efetch -format runinfo |
		csvtk cut -f BioProject |
		csvtk del-header \
			>"${_polap_var_bioproject_txt}"

	_polap_log1_file "${ODIR}/bioproject.txt"
	_polap_log0 $(head -n 1 "${ODIR}/bioproject.txt")

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Copy SRA fastq files for a given bioproject accession.
################################################################################
function _run_polap_copy-sra-bioproject() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	local _polap_var_assembly="${ODIR}/${INUM}"
	# local _polap_var_assembly="${ODIR}"
	local _polap_var_bioproject="${_polap_var_assembly}/70-bioproject"
	local _polap_var_bioproject_runinfo_all="${_polap_var_bioproject}/1-runinfo.all"
	local _polap_var_bioproject_runinfo="${_polap_var_bioproject}/1-runinfo.tsv"
	local _polap_var_bioproject_sra_long_read="${_polap_var_bioproject}/1-sra-long-read.tsv"
	local _polap_var_bioproject_sra_short_read="${_polap_var_bioproject}/1-sra-short-read.tsv"
	local _polap_var_bioproject_species="${_polap_var_bioproject}/1-species.txt"
	local _polap_var_bioproject_passed="${_polap_var_bioproject}/1-passed.txt"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Copy SRA fastq files for a given bioproject accession.
#
# first do the menu: get-bioproject
#
# 1. Do not delete $ODIR.
#
# Arguments:
#   -o $ODIR
#   -i $INUM
#   -b $SR2: bioproject ID not working
# Inputs:
#   bioproject ID
#   ${_polap_var_bioproject_runinfo}
# Outputs:
#   copy commands ...
Example: $(basename $0) ${_arg_menu[0]} -o PRJNA574453
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	# Clean and create the bioproject directory
	# rm -rf "${_polap_var_bioproject}"
	# mkdir -p "${_polap_var_bioproject}"

	# Set the BioProject ID
	local BIOPRJ="$SR2"
	#
	# Fetch the run information for the given BioProject
	# esearch -db bioproject -query "$BIOPRJ" |
	# 	elink -target sra |
	# 	efetch -format runinfo \
	# 		>${_polap_var_bioproject_runinfo_all}

	check_file_existence "${_polap_var_bioproject_runinfo_all}"
	check_file_existence "${_polap_var_bioproject_runinfo}"
	check_file_existence "${_polap_var_bioproject_sra_long_read}"
	check_file_existence "${_polap_var_bioproject_sra_short_read}"
	check_file_existence "${_polap_var_bioproject_species}"

	# cat "${_polap_var_bioproject_runinfo_all}" |
	# 	csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName |
	# 	csvtk csv2tab >"${_polap_var_bioproject_runinfo}"
	#
	# # Log the output file
	# _polap_log2_file "downloaded: ${_polap_var_bioproject_runinfo}"
	#
	# # Select SRA long and short read data
	# "$script_dir/run-polap-assemble-bioproject-1-select-sra.R" \
	# 	"${_polap_var_bioproject_runinfo}" \
	# 	"${_polap_var_bioproject_sra_long_read}" \
	# 	"${_polap_var_bioproject_sra_short_read}" \
	# 	"${_polap_var_bioproject_species}" \
	# 	2>"$_polap_output_dest"
	#
	# # Log the long and short read SRA files
	# _polap_log2_file "output1: ${_polap_var_bioproject_sra_long_read}"
	# _polap_log2_file "output2: ${_polap_var_bioproject_sra_short_read}"
	# _polap_log2_file "output3: ${_polap_var_bioproject_species}"

	# Check if the long-read dataset exists
	if [ ! -s "${_polap_var_bioproject_sra_long_read}" ]; then
		die "ERROR: no long-read dataset for the BioProject: $BIOPRJ"
	fi

	# Check if the short-read dataset exists
	if [ ! -s "${_polap_var_bioproject_sra_short_read}" ]; then
		die "ERROR: no short-read dataset for the BioProject: $BIOPRJ"
	fi

	if [ ! -s "${_polap_var_bioproject_species}" ]; then
		die "ERROR: no species name for the BioProject: $BIOPRJ"
	fi

	# _polap_log1 "LOG: $BIOPRJ - passed"
	# touch "${_polap_var_bioproject_passed}"

	local lsra=$(cut -f 1 "${_polap_var_bioproject_sra_long_read}")
	_polap_log1 "scp gingko:${lsra}.fastq.gz $BIOPRJ/"
	scp gingko:${lsra}.fastq.gz $BIOPRJ/

	local ssra=$(cut -f 1 "${_polap_var_bioproject_sra_short_read}")
	_polap_log1 "scp gingko:${ssra}_1.fastq.gz $BIOPRJ/"
	scp gingko:${ssra}_1.fastq.gz $BIOPRJ/
	_polap_log1 "scp gingko:${ssra}_2.fastq.gz $BIOPRJ/"
	scp gingko:${ssra}_2.fastq.gz $BIOPRJ/

	_polap_log1 "rm ${lsra}.fastq.gz"
	_polap_log1 "rm ${ssra}_1.fastq.gz"
	_polap_log1 "rm ${ssra}_2.fastq.gz"

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Get BioProject information using BioProject accession.
################################################################################
function _run_polap_get-bioproject() { # get BioProject info from NCBI
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Get BioProject information using BioProject accession.
#
# FIXME: PRJNA841235 no Illumina, PacBIO data: error run-polap-r-get-bioproject-1.R
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
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

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

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
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
	_polap_log1 "  run-polap-r-get-bioproject-1.R"
	_polap_log2 "    input1: ${_polap_var_bioproject_runinfo}"
	_polap_log2 "    output1: ${_polap_var_bioproject_sra_long_read}"
	_polap_log2 "    output2: ${_polap_var_bioproject_sra_short_read}"
	_polap_log2 "    output3: ${_polap_var_bioproject_species}"
	_polap_log2 "    output4: ${_polap_var_bioproject_sra_per_species}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-get-bioproject-1.R \
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

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# bioproject data preparation
################################################################################
function _run_polap_bioproject-prepare() { # bioproject data preparation
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	source "$script_dir/polap-variables-common.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# Set the output directory for the current job number
	if [[ -s "${_polap_var_bioproject_txt}" ]]; then
		_arg_bioproject=$(<"${_polap_var_bioproject_txt}")
		local BIOPRJ="${_arg_bioproject}"
	else
		_polap_log0 "NOTE: you have no ${_polap_var_bioproject_txt}"
		_polap_log0 "  save BioProject ID in ${_polap_var_bioproject_txt}"
		_polap_log0 "  if you already have the folder ${_polap_var_bioproject}"
		if [ "${_arg_short_read2_is}" = "on" ]; then
			_arg_bioproject="${SR2}"
		fi

		if [ -z "${_arg_bioproject}" ]; then
			_polap_log0 "ERROR: use --bioproject option"
			exit $EXIT_SUCCESS
		fi
	fi

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Runs the organelle-genome assembly for a given BioProject ID.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject with bioproject.txt
#   or
#   -b or --bioproject ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${ODIR}/bioproject.txt
#   or
#   -b or --bioproject <BioProject ID>
#
#   long-read input data in the order of priority
#   ${_polap_var_base_nk_fq_gz}
#   ${_polap_var_base_l_fq_gz}
#   ${_polap_var_bioproject_sra_long_read}
#   to download the data if no such file
#
#   short-read input data in the order of priority
#   ${_polap_var_base_msbwt_tar_gz}
# 	${_polap_var_base_msbwt}
#   ${_polap_var_bioproject_sra_short_read} 
#   to download the data if no such file
# Outputs:
#   Long- and short-read data files
Example: $(basename $0) ${_arg_menu[0]} -o ${ODIR}
Example: $(basename $0) ${_arg_menu[0]} -o ${ODIR} [-b ${BIOPRJ}]
HEREDOC
	)

	local _long_read_data=""
	local _short_read_data1=""
	local _short_read_data2=""

	LRNK="${ODIR}/nk.fq.gz"

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

	_polap_log0 "preparing for assembling organelle genomes of BioProject ${BIOPRJ} ..."

	if [ -d "${ODIR}" ]; then
		_polap_log0 "  the main output folder: ${ODIR}"
	else
		_polap_log0 "  no output folder; creating ${ODIR}"
		_polap_log3_cmd mkdir -p "${ODIR}"
	fi

	# Get the long-read information from NCBI
	if [ -s "${_polap_var_bioproject_sra_long_read}" ]; then
		_polap_log1 "  using SRA info at: ${_polap_var_bioproject_sra_long_read}"
		_polap_log1_cat "${_polap_var_bioproject_sra_long_read}"
	else
		_polap_log0 "ERROR: no such file: ${_polap_var_bioproject_sra_long_read}"
		return
		# _run_polap_get-bioproject
	fi

	# Fetch the long-read dataset
	if [ -s "${LRNK}" ]; then
		_polap_log1 "  we use the long-read data: ${LRNK}"
		_long_read_data="${_polap_var_base_nk_fq_gz}"
	elif [[ -s "${_polap_var_base_l_fq_gz}" ]]; then
		_polap_log1 "  we use the long-read data: ${_polap_var_base_l_fq_gz}"
		_long_read_data="${_polap_var_base_l_fq_gz}"
	elif [ -s "${_polap_var_bioproject_sra_long_read}" ]; then
		local SRA=$(cut -f1 "${_polap_var_bioproject_sra_long_read}")
		LR="${ODIR}/${SRA}.fastq"
		if [ -s "${LR}" ]; then
			_polap_log2_file "${LR}"
		elif [ -s "${LR}.gz" ]; then
			_polap_log2_file "${LR}.gz is being extracted ..."
			gunzip "${LR}.gz"
		else
			"$script_dir"/run-polap-ncbitools fetch sra "$SRA"
			mv "$SRA.fastq" "${ODIR}/"
		fi
		_long_read_data="${LR}"
	else
		die "ERROR: no long-read dataset for the BioProject: $BIOPRJ"
	fi

	# Fetch the short-read dataset
	if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
		_polap_log1 "  we use the short-read data in: ${_polap_var_base_msbwt_tar_gz}"
		_short_read_data1=${_polap_var_base_msbwt_tar_gz}
		# check_file_existence "${_polap_var_base_genome_size}"
	elif [[ -s "${_polap_var_base_msbwt}" ]]; then
		_polap_log1 "  we use the short-read data in: ${_polap_var_base_msbwt}"
		_short_read_data1=${_polap_var_base_msbwt}
		# check_file_existence "${_polap_var_base_genome_size}"
	elif [ -s "${_polap_var_bioproject_sra_short_read}" ]; then
		local SRA=$(cut -f1 "${_polap_var_bioproject_sra_short_read}")
		SR1="${ODIR}/${SRA}_1.fastq"
		SR2="${ODIR}/${SRA}_2.fastq"
		if [ -s "${SR1}" ]; then
			_polap_log2_file "${SR1}"
		elif [ -s "${SR1}.gz" ]; then
			_polap_log2_file "${SR1}.gz is being extracted ..."
			gunzip "${SR1}.gz"
		fi

		if [ -s "${SR2}" ]; then
			_polap_log2_file "${SR2}"
		elif [ -s "${SR2}.gz" ]; then
			_polap_log2_file "${SR2}.gz is being extracted ..."
			gunzip "${SR2}.gz"
		fi
		# if [ -s "${SR1}" ]; then
		# 	seqkit stats -T "${SR1}" >"${_polap_var_base_fq_stats}"
		# fi
		# if [ -s "${SR2}" ]; then
		# 	seqkit stats -T "${SR2}" >>"${_polap_var_base_fq_stats}"
		# fi

		if [ ! -s "${SR1}" ] || [ ! -s "${SR2}" ]; then
			_polap_log1 "  downloading the paired-end short-read data: $SRA"
			"$script_dir"/run-polap-ncbitools fetch sra "$SRA"
			mv "${SRA}_1.fastq" "${ODIR}/"
			mv "${SRA}_2.fastq" "${ODIR}/"
		fi

		_short_read_data1=${SR1}
		_short_read_data2=${SR2}
	else
		die "ERROR: no read short-read dataset for the BioProject: $BIOPRJ"
	fi

	# check_file_existence "${LR}"
	# check_file_existence "${SR1}"
	# check_file_existence "${SR2}"
	#
	# if [ -s "${_polap_var_base_fq_stats}" ]; then
	# 	_polap_log2_file "${_polap_var_base_fq_stats}"
	# 	_polap_log3_cat "${_polap_var_base_fq_stats}"
	# else
	# 	_run_polap_summary-reads
	# fi

	# _run_polap_total-length-long
	# _run_polap_find-genome-size
	# _run_polap_reduce-data

	if [ -s "${_polap_var_base_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_base_msbwt_tar_gz} ... later when we polish it with the short-read data."
			# tar zxf "${_polap_var_base_msbwt_tar_gz}"
		else
			_polap_log1 "  Do the preparation of short-read polishing ... early"
			check_file_existence "${SR1}"
			check_file_existence "${SR2}"
			_run_polap_prepare-polishing
		fi
	fi

	_polap_log0 "Long-read data: ${_long_read_data}"
	_polap_log0 "Short-read data: ${_short_read_data1}"
	if [[ -n "${_short_read_data2}" ]]; then
		_polap_log0 "Short-read data: ${_short_read_data2}"
	fi
	_polap_log0 "LR: ${LR}"
	_polap_log0 "SR1: ${SR1}"
	_polap_log0 "SR2: ${SR2}"
	_polap_log0 "NEXT: assemble-draft"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# called after bioproject-prepare
################################################################################
function _run_polap_x-assemble-draft() { # called after bioproject-prepare
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh"  # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Runs the POLAP organelle-genome assembly with sequencing data.
# 
# Arguments:
#   -o $ODIR
#   -l $LR: a long-read fastq data file
#   -a $SR1: a short-read fastq data file
#   -b $SR2: another short-read fastq data file
# Inputs:
#   $LR: a long-read fastq 
#   $SR1: a short-read fastq data file
#   $SR2: another short-read fastq data file
# Outputs:
#   $MTDIR/assembly_graph.gfa
Example: $(basename $0) ${_arg_menu[0]} --test
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Not delete the output directory.
	mkdir -p "$ODIR"

	# Run assembly, annotation, and contig selection steps
	if [ -s "${_polap_var_wga_contigger_gfa}" ]; then
		_polap_log1 "  skipping the whole-genome assembly"
	else
		_run_polap_assemble1
	fi

	if [ -s "${_polap_var_wga_annotation}" ]; then
		_polap_log1 "  skipping the organelle annotation on the whole-genome"
	else
		_run_polap_edges-stats
		_run_polap_annotate
	fi

	# Select seed contigs
	_run_polap_select-contigs

	# Array to store the names of the original files
	files=($(ls "${_polap_var_wga}/mt.contig.name-"?))

	# Temporary array to store the paths of unique files
	unique_files=()

	# Function to check if a file is unique
	is_unique() {
		for unique_file in "${unique_files[@]}"; do
			if cmp -s "$1" "$unique_file"; then
				echo "$unique_file"
				return 1 # Not unique
			fi
		done
		return 0 # Unique
	}

	_polap_log1 "Checking for unique files and their matches:"

	# Iterate over the files to find unique ones
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Call the function corresponding to the current number (index is i-1)
		INUM=0
		FDIR="${ODIR}/${INUM}"
		JNUM="${i}"
		file="$FDIR"/mt.contig.name-$JNUM
		mkdir -p "${ODIR}/${JNUM}"

		unique_file=$(is_unique "$file")
		if [ $? -eq 0 ]; then
			# If unique, add it to the unique_files array
			unique_files+=("$file")
			echo "$file is unique."

			MTCONTIGNAME="$file"
			# check the mt.contig.name-1
			if [ -s "$MTCONTIGNAME" ]; then
				# Run secondary assembly, polishing, and mtDNA selection steps
				_polap_log1_file "${MTCONTIGNAME}"

				if [ -s "${ODIR}/${JNUM}/30-contigger/graph_final.gfa" ] && [ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  skipping organelle-genome assembly -i 0 -j ${JNUM} ..."
				else
					_polap_log2 "  not skipping organelle-genome assembly -i 0 -j ${JNUM} ..."
					_arg_yes="on"
					_run_polap_assemble2
				fi

				if [ -s "${ODIR}/${JNUM}/contig-annotation-table.txt" ] && [ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  skipping organelle-genome annotation -i ${JNUM} ..."
				else
					INUM="${i}" _run_polap_annotate
				fi

				if [ -s "${ODIR}/${JNUM}/assembly_graph.gfa" ] && [ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  skipping organelle-genome long-read polishing -i ${JNUM} ..."
				else
					JNUM="${i}" _run_polap_flye-polishing
				fi

				if [ -s "${ODIR}/${JNUM}/mt.0.fasta" ] && [ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  skipping organelle-genome extraction -i ${JNUM} ..."
				else
					INUM="${i}" _run_polap_select-mtdna
				fi
			else
				_polap_log1 "LOG: $MTCONTIGNAME is empty for select-contig type $i ..."
			fi
		else
			_polap_log1 "$file is the same as $unique_file."
		fi
	done

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Postprocess the bioproject assembly
################################################################################
function _run_polap_bioproject-postprocess() { # postprocess the bioproject assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	source "$script_dir/polap-variables-common.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# Set the output directory for the current job number
	if [[ -s "${_polap_var_bioproject_txt}" ]]; then
		_arg_bioproject=$(<"${_polap_var_bioproject_txt}")
		local BIOPRJ="${_arg_bioproject}"
	else
		_polap_log0 "NOTE: you have no ${_polap_var_bioproject_txt}"
		_polap_log0 "  save BioProject ID in ${_polap_var_bioproject_txt}"
		_polap_log0 "  if you already have the folder ${_polap_var_bioproject}"
		if [ "${_arg_short_read2_is}" = "on" ]; then
			_arg_bioproject="${SR2}"
		fi

		if [ -z "${_arg_bioproject}" ]; then
			_polap_log0 "ERROR: use --bioproject option"
			exit $EXIT_SUCCESS
		fi
	fi

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Postprocess a bioproject assembly.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject with bioproject.txt
#   or
#   -b or --bioproject ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${ODIR}/bioproject.txt
#   or
#   -b or --bioproject <BioProject ID>
#
#   long-read input data in the order of priority
#   ${_polap_var_base_nk_fq_gz}
#   ${_polap_var_base_l_fq_gz}
#   ${_polap_var_bioproject_sra_long_read}
#   to download the data if no such file
#
#   short-read input data in the order of priority
#   ${_polap_var_base_msbwt_tar_gz}
# 	${_polap_var_base_msbwt}
#   ${_polap_var_bioproject_sra_short_read} 
#   to download the data if no such file
# Outputs:
#   Long- and short-read data files
Example: $(basename $0) ${_arg_menu[0]} -o ${ODIR}
Example: $(basename $0) ${_arg_menu[0]} -o ${ODIR} [-b ${BIOPRJ}]
HEREDOC
	)

	LRNK="${ODIR}/nk.fq.gz"

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

	_polap_log0 "protprocessing the organelle genome assembly of BioProject ${BIOPRJ} ..."

	if [ -s "${_polap_var_base_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_base_msbwt_tar_gz} ..."
			tar -zxf "${_polap_var_base_msbwt_tar_gz}" -C "${ODIR}"
		else
			_polap_log1 "  Do the preparation of short-read polishing ... early"
			check_file_existence "${SR1}"
			check_file_existence "${SR2}"
			_run_polap_prepare-polishing
		fi
	fi

	# Run the polishing step
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Define the paths for mtDNA sequences to be polished
		PA="${ODIR}/${i}/mt.0.fasta"
		FA="${ODIR}/${i}/mt.1.fa"
		check_file_existence "${_polap_var_base_msbwt}"

		if [ -s "${PA}" ]; then
			if [ -s "${FA}" ]; then
				_polap_log1 "  skipping the short-read polishing ..."
			else
				_polap_log1 "  mtDNA is being polished for select-contig type ${i}..."
				_run_polap_polish
			fi
		else
			_polap_log1 "  no mtDNA candidate for select-contig type ${i}..."
		fi
	done

	# Download the known reference mtDNA sequence in fasta format if available
	# -o PRJNA914763
	# INUM=0
	# if [ -s "${_polap_var_bioproject_mtdna_fasta2}" ]; then
	# 	_polap_log2 "  skipping downloading mtDNA from NCBI"
	# 	_polap_log2 "  we use the mtDNA: ${_polap_var_bioproject_mtdna_fasta2}"
	# else
	# 	_run_polap_get-mtdna
	#
	# fi

	if [[ -s "${_polap_var_bioproject_mtdna_fasta2}" ]]; then

		# Compare the known mtDNA and the assembled one.
		for i in "${_arg_select_contig_numbers[@]}"; do
			# Define the paths for mtDNA sequences to be polished
			FA="${ODIR}/${i}/mt.1.fa"

			if [ -s "${FA}" ]; then
				INUM="${i}"
				_run_polap_compare-mtdna
			else
				_polap_log1 "  skipping the short-read polishing ..."
				INUM="${i}"
				source "$script_dir/polap-variables-common.sh" # '.' means 'source'
				local n1=$(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")
				local l1="0"
				local l2="0"
				local c1="0"
				printf "%s\t%d\t%d\t%f\n" ${n1} ${l1} ${l2} ${c1} >"${_polap_var_mtdna_compare}"
			fi
		done
	else
		_polap_log0 "No known mtDNA for $(<${_polap_var_bioproject_species})"
	fi

	# Finalize the assemble-bioproject function.
	_polap_log1 "rsync back to thorne ..."
	touch "${ODIR}/log-assembled.txt"
	# touch "$BIOPRJ/log-copying-to-xxx.txt"
	# touch "$BIOPRJ/log-not-started-yet.txt"
	if [ "$(hostname)" = "thorne" ]; then
		if [[ "${PWD}" = "/home/goshng/run/polap" ]]; then
			cp -pr "${ODIR}/" "/home/goshng/all/polap/figshare/bioprojects/"
			touch "/home/goshng/all/polap/figshare/bioprojects/log-assembled-${ODIR}.txt"
		fi
	else
		touch "log-assembled-${ODIR}.txt"
		scp -pq "log-assembled-${ODIR}.txt" thorne:"$PWD/"
		rsync -a -e ssh "${ODIR}/" thorne:"$PWD/${ODIR}/"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Assembles organelle genome from a BioProject using long and short reads.
# Compares the known and assembled mtDNA using BLAST.
#
# bioproject.txt is necessary.
################################################################################
function _run_polap_assemble-bioproject() { # main function for this bioproject module
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	source "$script_dir/polap-variables-common.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# Set the output directory for the current job number
	if [[ -s "${_polap_var_bioproject_txt}" ]]; then
		_arg_bioproject=$(<"${_polap_var_bioproject_txt}")
		local BIOPRJ="${_arg_bioproject}"
	else
		_polap_log0 "NOTE: you have no ${_polap_var_bioproject_txt}"
		_polap_log0 "  save BioProject ID in ${_polap_var_bioproject_txt}"
		_polap_log0 "  if you already have the folder ${_polap_var_bioproject}"
		if [ "${_arg_short_read2_is}" = "on" ]; then
			_arg_bioproject="${SR2}"
		fi

		if [ -z "${_arg_bioproject}" ]; then
			_polap_log0 "ERROR: use --bioproject option"
			exit $EXIT_SUCCESS
		fi
	fi

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Runs the organelle-genome assembly for a given BioProject ID.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject with bioproject.txt
#   or
#   -b or --bioproject ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${ODIR}/bioproject.txt
#   or
#   -b or --bioproject <BioProject ID>
#
#   long-read input data in the order of priority
#   ${_polap_var_base_nk_fq_gz}
#   ${_polap_var_base_l_fq_gz}
#   ${_polap_var_bioproject_sra_long_read}
#   to download the data if no such file
#
#   short-read input data in the order of priority
#   ${_polap_var_base_msbwt_tar_gz}
# 	${_polap_var_base_msbwt}
#   ${_polap_var_bioproject_sra_short_read} 
#   to download the data if no such file
# Outputs:
#   mt.1.fa
Example: $(basename $0) ${_arg_menu[0]} -o ${ODIR} [-b ${BIOPRJ}]
HEREDOC
	)

	LRNK="${ODIR}/nk.fq.gz"

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

	_polap_log0 "assembling organelle genomes of BioProject ${BIOPRJ} ..."

	_run_polap_bioproject-prepare

	_run_polap_x-assemble-draft

	_run_polap_bioproject-postprocess

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}