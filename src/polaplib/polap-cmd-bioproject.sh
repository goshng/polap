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

# TODO: documentation

################################################################################
# BioProject subcommand collects data from the NCBI.
# This script is often changed to adapt to the different polap data analysis.
# It is used to prepare the data for organelle genome assembly.
# It fetches BioProject information, SRA data, and prepares the data for assembly.
# The get-bioproject command is used to fetch BioProject information from NCBI.
# It is the main entry point for the BioProject module.
# But, the start was from the summary file from the NCBI.
# One could search the NCBI with the following query:
# "Viridiplantae"[Organism] AND oxford[All Fields]
# to save bioproject_result.txt file. We want to iterate over the summary
# so that we could create a table for the BioProject ID, SRA ID, and species name.
#
# Collect animal data for animal mitochondrial genome assembly.
#
# NCBI search:
# NCBI:BioProject:"Viridiplantae"[Organism] AND oxford[All Fields]
# NCBI:BioProject:"Viridiplantae"[Organism] AND pacbio[All Fields]
# NCBI:BioProject:"animalia"[Organism] AND oxford[All Fields]
#
# function _run_polap_download-runinfo-bioproject {
# function _run_polap_get-bioproject-sra {
# function _run_polap_copy-sra-bioproject {
# function _run_polap_get-bioproject { # get BioProject info from NCBI
# function _run_polap_bioproject-prepare { # bioproject data preparation
# function _run_polap_x-assemble-draft { # called after bioproject-prepare
# function _run_polap_bioproject-postprocess { # postprocess the bioproject assembly
# function _run_polap_assemble-bioproject { # main function for this bioproject module
#
# TODO:
#   - menus need to be documented.
#   - Collect animal data for animal mitochondrial genome assembly.
#
# TEST-SCC: not yet
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

source "${_POLAPLIB_DIR}/polap-constants.sh"

################################################################################
# Gets BioProject ID of a given SRA ID.
#
# Arguments:
#   --sra ${_arg_sra}: SRA ID
#   -o ${_arg_outdir}: output folder
# Inputs:
#   SRA ID
# Outputs:
#   BioProject ID
#   FILE: ${_polap_var_project_txt}
################################################################################
function _run_polap_get-bioproject-sra {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Gets BioProject ID of a given SRA ID.
#
# Arguments:
#   --sra ${_arg_sra}: SRA ID
#   -o ${_arg_outdir}: output folder
# Inputs:
#   SRA ID
# Outputs:
#   BioProject ID
#   FILE: ${_polap_var_project_txt}
Example: $(basename $0) ${_arg_menu[0]} --sra <SRA ID>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  if [[ ${_arg_menu[1]} == "view" ]]; then
    if [[ -s "${_polap_var_project_txt}" ]]; then
      _polap_log0_cat "${_polap_var_project_txt}"
    else
      _polap_log0 "No such file: ${_polap_var_project_txt}"
    fi
    exit $EXIT_SUCCESS
  fi

  if [[ ! -d "${_arg_outdir}" ]]; then
    _polap_log1 "  no output folder, creating ${_arg_outdir}"
    mkdir -p "${_arg_outdir}"
  fi

  esearch -db sra -query "${_arg_sra}" |
    efetch -format runinfo |
    csvtk cut -f BioProject |
    csvtk del-header \
      >"${_polap_var_project_txt}"

  _polap_log1_file "${_arg_outdir}/bioproject.txt"
  _polap_log0 $(head -n 1 "${_arg_outdir}/bioproject.txt")

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Copy SRA fastq files for a given bioproject accession.
################################################################################
function _run_polap_copy-sra-bioproject {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

  # Set variables for file paths
  local _polap_var_assembly="${_arg_outdir}/${_arg_inum}"
  # local _polap_var_assembly="${_arg_outdir}"
  local _polap_var_project="${_polap_var_assembly}/70-bioproject"
  local _polap_var_project_runinfo_all="${_polap_var_project}/1-runinfo.all"
  local _polap_var_project_runinfo="${_polap_var_project}/1-runinfo.tsv"
  local _polap_var_project_sra_long_read="${_polap_var_project}/1-sra-long-read.tsv"
  local _polap_var_project_sra_short_read="${_polap_var_project}/1-sra-short-read.tsv"
  local _polap_var_project_species="${_polap_var_project}/1-species.txt"
  local _polap_var_project_passed="${_polap_var_project}/1-passed.txt"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Copy SRA fastq files for a given bioproject accession.
#
# first do the menu: get-bioproject
#
# 1. Do not delete ${_arg_outdir}.
#
# Arguments:
#   -o ${_arg_outdir}
#   -i ${_arg_inum}
#   -b ${_arg_short_read1}: bioproject ID not working
# Inputs:
#   bioproject ID
#   ${_polap_var_project_runinfo}
# Outputs:
#   copy commands ...
Example: $(basename $0) ${_arg_menu[0]} -o PRJNA574453
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

  # Clean and create the bioproject directory
  # rm -rf "${_polap_var_project}"
  # mkdir -p "${_polap_var_project}"

  # Set the BioProject ID
  local BIOPRJ="$SR2"
  #
  # Fetch the run information for the given BioProject
  # esearch -db bioproject -query "$BIOPRJ" |
  # 	elink -target sra |
  # 	efetch -format runinfo \
  # 		>${_polap_var_project_runinfo_all}

  check_file_existence "${_polap_var_project_runinfo_all}"
  check_file_existence "${_polap_var_project_runinfo}"
  check_file_existence "${_polap_var_project_sra_long_read}"
  check_file_existence "${_polap_var_project_sra_short_read}"
  check_file_existence "${_polap_var_project_species}"

  # cat "${_polap_var_project_runinfo_all}" |
  # 	csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName |
  # 	csvtk csv2tab >"${_polap_var_project_runinfo}"
  #
  # # Log the output file
  # _polap_log2_file "downloaded: ${_polap_var_project_runinfo}"
  #
  # # Select SRA long and short read data
  # "${_POLAPLIB_DIR}/run-polap-assemble-bioproject-1-select-sra.R" \
  # 	"${_polap_var_project_runinfo}" \
  # 	"${_polap_var_project_sra_long_read}" \
  # 	"${_polap_var_project_sra_short_read}" \
  # 	"${_polap_var_project_species}" \
  # 	2>"$_polap_output_dest"
  #
  # # Log the long and short read SRA files
  # _polap_log2_file "output1: ${_polap_var_project_sra_long_read}"
  # _polap_log2_file "output2: ${_polap_var_project_sra_short_read}"
  # _polap_log2_file "output3: ${_polap_var_project_species}"

  # Check if the long-read dataset exists
  if [ ! -s "${_polap_var_project_sra_long_read}" ]; then
    die "ERROR: no long-read dataset for the BioProject: $BIOPRJ"
  fi

  # Check if the short-read dataset exists
  if [ ! -s "${_polap_var_project_sra_short_read}" ]; then
    die "ERROR: no short-read dataset for the BioProject: $BIOPRJ"
  fi

  if [ ! -s "${_polap_var_project_species}" ]; then
    die "ERROR: no species name for the BioProject: $BIOPRJ"
  fi

  # _polap_log1 "LOG: $BIOPRJ - passed"
  # touch "${_polap_var_project_passed}"

  local lsra=$(cut -f 1 "${_polap_var_project_sra_long_read}")
  _polap_log1 "scp gingko:${lsra}.fastq.gz $BIOPRJ/"
  scp gingko:${lsra}.fastq.gz $BIOPRJ/

  local ssra=$(cut -f 1 "${_polap_var_project_sra_short_read}")
  _polap_log1 "scp gingko:${ssra}_1.fastq.gz $BIOPRJ/"
  scp gingko:${ssra}_1.fastq.gz $BIOPRJ/
  _polap_log1 "scp gingko:${ssra}_2.fastq.gz $BIOPRJ/"
  scp gingko:${ssra}_2.fastq.gz $BIOPRJ/

  _polap_log1 "rm ${lsra}.fastq.gz"
  _polap_log1 "rm ${ssra}_1.fastq.gz"
  _polap_log1 "rm ${ssra}_2.fastq.gz"

  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Get BioProject information using BioProject accession.
################################################################################
function _run_polap_get-bioproject { # get BioProject info from NCBI
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
Get BioProject information using BioProject accession.

1. search the NCBI for a list of BioProject accession IDs
  -> bioproject_result.accession.txt
2. download runinfo for a list of BioProject
  polap download-runinfo-bioproject bioproject_result.accession.txt -o runinfo
3. use get-bioproject polap command for each runinfo CSV file
  see polap-bash-get-bioproject.sh for reference

NCBI search:
  NCBI:BioProject:"Viridiplantae"[Organism] AND oxford[All Fields]
  NCBI:BioProject:"Viridiplantae"[Organism] AND pacbio[All Fields]
  NCBI:BioProject:"animalia"[Organism] AND oxford[All Fields]

FIXME: PRJNA841235 no Illumina, PacBIO data: error polap-r-get-bioproject.R

Arguments:
  --bioproject <BioProject ID>

Inputs:
  bioproject ID

Outputs:
  ${_polap_var_project_runinfo}

Usage:
  PRJNA557253 - multiple species

Example:
$(basename $0) ${_arg_menu[0]} --bioproject PRJNA574453
$(basename $0) ${_arg_menu[0]} view -o PRJNA5744532 2>&1 | grep DRR196916
$(basename $0) ${_arg_menu[0]} pacbio --bioproject PRJNA574453
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then
    if [ -s "${_polap_var_project_runinfo}" ]; then
      _polap_log0_file "${_polap_var_project_runinfo}"
      _polap_log0_cat "${_polap_var_project_runinfo}"
      _polap_log0 ""
    else
      _polap_log0 "No BioProject info."
    fi
    if [ -s "${_polap_var_project_sra_long_read}" ]; then
      _polap_log0_file "${_polap_var_project_sra_long_read}"
      _polap_log0_cat "${_polap_var_project_sra_long_read}"
      _polap_log0 ""
      SRA=$(cut -f1 "${_polap_var_project_sra_long_read}")
      echo $SRA >"${_arg_outdir}/sra.txt"
    else
      _polap_log0 "No long-read info."
    fi
    if [ -s "${_polap_var_project_sra_short_read}" ]; then
      _polap_log0_file "${_polap_var_project_sra_short_read}"
      _polap_log0_cat "${_polap_var_project_sra_short_read}"
      _polap_log0 ""
      SRA=$(cut -f1 "${_polap_var_project_sra_short_read}")
      echo $SRA >>"${_arg_outdir}/sra.txt"
    else
      _polap_log0 "No short-read info."
    fi

    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    exit $EXIT_SUCCESS
  fi

  if [[ -d "${_polap_var_project}" ]]; then
    return ${_POLAP_ERR_ALREADY_EXIST_OUT}
  fi

  # Check if --bioproject option is provided.
  if [ -z "${_arg_bioproject}" ]; then
    _polap_log0 "ERROR: --bioproject option is required to download BioProject runinfo."
    _polap_log0 "       you might want to view the BioProject runinfo."
    _polap_log0 "       $0 ${_arg_menu[0]} view -o ${_arg_outdir}"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    exit $EXIT_SUCCESS
  else
    local BIOPRJ="${_arg_bioproject}"
  fi

  # Set the BioProject ID
  _polap_log0 "BioProject: ${BIOPRJ}"
  echo "${BIOPRJ}" >"${_polap_var_project_txt}"

  # Delete the 00-bioproject folder.
  if [[ -d "${_polap_var_project}" ]]; then
    _polap_log3_cmd rm -rf "${_polap_var_project}"
  fi

  # Fetch the run information for the given BioProject
  if [[ ! -s "${_polap_var_project_runinfo_all}" ]] || [[ "${_arg_redo}" = "on" ]]; then
    _polap_log0 "  downloading ${_polap_var_project_runinfo_all} ..."
    esearch -db bioproject -query "$BIOPRJ" |
      elink -target sra |
      efetch -format runinfo \
        >${_polap_var_project_runinfo_all}
  else
    _polap_log0 "  found: ${_polap_var_project_runinfo_all}, so skipping downloading the runinfo ..."
  fi

  # create the 0-bioproject folder
  _polap_log1 "  creating output folder: ${_polap_var_project}"
  _polap_log3_cmd mkdir -p "${_polap_var_project}"

  check_file_existence "${_polap_var_project_runinfo_all}"

  # Log the output file
  _polap_log1 "  extracting some columns in TSV: ${_polap_var_project_runinfo}"
  _polap_log3_pipe "cat ${_polap_var_project_runinfo_all} |\
    csvtk cut -f Run,bases,avgLength,size_MB,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,Model,ScientificName |\
    csvtk csv2tab >${_polap_var_project_runinfo}"

  # Select SRA long and short read data
  _polap_log1 "  polap-r-get-bioproject.R"
  _polap_log2 "    input1: ${_polap_var_project_runinfo}"
  _polap_log2 "    output1: ${_polap_var_project_sra_long_read}"
  _polap_log2 "    output2: ${_polap_var_project_sra_short_read}"
  _polap_log2 "    output3: ${_polap_var_project_species}"
  _polap_log2 "    output4: ${_polap_var_project_sra_per_species}"

  if [[ "${_arg_menu[1]}" == "pacbio" ]]; then
    _polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/polap-r-get-bioproject.R \
    --input ${_polap_var_project_runinfo_all} \
    --out ${_arg_outdir}/00-bioproject \
    --pacbio
		2>$_polap_output_dest"
  else
    _polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/polap-r-get-bioproject.R \
    --input ${_polap_var_project_runinfo_all} \
    --out ${_arg_outdir}/00-bioproject \
		2>$_polap_output_dest"
  fi

  # Check if the long-read dataset exists
  if [[ -s "${_polap_var_project_sra_long_read}" ]]; then
    _polap_log1_file "${_polap_var_project_sra_long_read}"
    _polap_log1_cat "${_polap_var_project_sra_long_read}"
  else
    _polap_log0 "INFO: no long-read SRA data for the BioProject: $BIOPRJ"
  fi

  # Check if the short-read dataset exists
  if [[ -s "${_polap_var_project_sra_short_read}" ]]; then
    _polap_log1_file "${_polap_var_project_sra_short_read}"
    _polap_log1_cat "${_polap_var_project_sra_short_read}"
  else
    _polap_log0 "INFO: no short-read SRA data for the BioProject: $BIOPRJ"
  fi

  # Iterate through each line of the file, skipping the header
  _polap_log1 "  creating output folders per species ..."
  local input_file="${_polap_var_project_sra_per_species}"
  tail -n +2 "$input_file" | while IFS=$'\t' read -r ScientificName Run_Nano Run_Illumina; do
    _polap_log2 "    taxon: ${ScientificName}"

    local bioproject_species_base="${_arg_outdir}-$(echo $ScientificName | tr ' ' '_')"
    local bioproject_species_folder="${_arg_outdir}-$(echo $ScientificName | tr ' ' '_')/0-bioproject"
    local taxon_id_folder="${_arg_outdir}/0-bioproject/$(echo $ScientificName | tr ' ' '_')"

    # Create a directory with the taxon ID as the name
    _polap_log3_cmd mkdir -p "$taxon_id_folder"
    # grep "^${Run_Nano}" "${_polap_var_project_runinfo}" >"$taxon_id_folder/1-sra-long-read.tsv"
    # grep "^${Run_Illumina}" "${_polap_var_project_runinfo}" >"$taxon_id_folder/1-sra-short-read.tsv"
    _polap_log3_pipe "grep ^${Run_Nano} ${_polap_var_project_runinfo} >$taxon_id_folder/1-sra-long-read.tsv"
    _polap_log3_pipe "grep ^${Run_Illumina} ${_polap_var_project_runinfo} >$taxon_id_folder/1-sra-short-read.tsv"

    _polap_log0 "    creates a taxon directory: ${bioproject_species_base}"
    _polap_log3_cmd mkdir -p "${bioproject_species_folder}"
    _polap_log3_cmd touch "${bioproject_species_base}/log-need-to-fetch-data.txt"
    _polap_log3_cmd cp "${_polap_var_project_txt}" "${bioproject_species_base}"
    _polap_log3_cmd cp "${_polap_var_project}"/*.{tsv,txt} "${bioproject_species_folder}"
    _polap_log3_cmd cp "${taxon_id_folder}"/*.tsv "${bioproject_species_folder}"
    # _polap_log3_cmd cp "$taxon_id_folder/1-sra-long-read.tsv" "${bioproject_species_folder}"
    # _polap_log3_cmd cp "$taxon_id_folder/1-sra-short-read.tsv" "${bioproject_species_folder}"
    _polap_log3_pipe "echo ${ScientificName} >${bioproject_species_folder}/1-species.txt"
  done

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# bioproject data preparation
################################################################################
function _run_polap_bioproject-prepare { # bioproject data preparation
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set variables for file paths
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  # Set the output directory for the current job number
  if [[ -s "${_polap_var_project_txt}" ]]; then
    _arg_bioproject=$(<"${_polap_var_project_txt}")
    local BIOPRJ="${_arg_bioproject}"
  else
    _polap_log0 "NOTE: you have no ${_polap_var_project_txt}"
    _polap_log0 "  save BioProject ID in ${_polap_var_project_txt}"
    _polap_log0 "  if you already have the folder ${_polap_var_project}"
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
#   -o ${_arg_outdir}: output folder for BioProject with bioproject.txt
#   or
#   -b or --bioproject ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${_arg_outdir}/bioproject.txt
#   or
#   -b or --bioproject <BioProject ID>
#
#   long-read input data in the order of priority
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_polap_var_outdir_l_fq_gz}
#   ${_polap_var_project_sra_long_read}
#   to download the data if no such file
#
#   short-read input data in the order of priority
#   ${_polap_var_outdir_msbwt_tar_gz}
# 	${_polap_var_outdir_msbwt}
#   ${_polap_var_project_sra_short_read} 
#   to download the data if no such file
# Outputs:
#   Long- and short-read data files
Example: $(basename $0) ${_arg_menu[0]} -o ${_arg_outdir}
Example: $(basename $0) ${_arg_menu[0]} -o ${_arg_outdir} [-b ${BIOPRJ}]
HEREDOC
  )

  local _long_read_data=""
  local _short_read_data1=""
  local _short_read_data2=""

  LRNK="${_arg_outdir}/nk.fq.gz"

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

  _polap_log0 "preparing for assembling organelle genomes of BioProject ${BIOPRJ} ..."

  if [ -d "${_arg_outdir}" ]; then
    _polap_log0 "  the main output folder: ${_arg_outdir}"
  else
    _polap_log0 "  no output folder; creating ${_arg_outdir}"
    _polap_log3_cmd mkdir -p "${_arg_outdir}"
  fi

  # Get the long-read information from NCBI
  if [ -s "${_polap_var_project_sra_long_read}" ]; then
    _polap_log1 "  using SRA info at: ${_polap_var_project_sra_long_read}"
    _polap_log1_cat "${_polap_var_project_sra_long_read}"
  else
    _polap_log0 "ERROR: no such file: ${_polap_var_project_sra_long_read}"
    return
    # _run_polap_get-bioproject
  fi

  # Fetch the long-read dataset
  if [ -s "${LRNK}" ]; then
    _polap_log1 "  we use the long-read data: ${LRNK}"
    _long_read_data="${_polap_var_outdir_nk_fq_gz}"
  elif [[ -s "${_polap_var_outdir_l_fq_gz}" ]]; then
    _polap_log1 "  we use the long-read data: ${_polap_var_outdir_l_fq_gz}"
    _long_read_data="${_polap_var_outdir_l_fq_gz}"
  elif [ -s "${_polap_var_project_sra_long_read}" ]; then
    local SRA=$(cut -f1 "${_polap_var_project_sra_long_read}")
    LR="${_arg_outdir}/${SRA}.fastq"
    if [ -s "${LR}" ]; then
      _polap_log2_file "${LR}"
    elif [ -s "${LR}.gz" ]; then
      _polap_log2_file "${LR}.gz is being extracted ..."
      gunzip "${LR}.gz"
    else
      "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$SRA"
      mv "$SRA.fastq" "${_arg_outdir}/"
    fi
    _long_read_data="${LR}"
  else
    die "ERROR: no long-read dataset for the BioProject: $BIOPRJ"
  fi

  # Fetch the short-read dataset
  if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
    _polap_log1 "  we use the short-read data in: ${_polap_var_outdir_msbwt_tar_gz}"
    _short_read_data1=${_polap_var_outdir_msbwt_tar_gz}
    # check_file_existence "${_polap_var_outdir_genome_size}"
  elif [[ -s "${_polap_var_outdir_msbwt}" ]]; then
    _polap_log1 "  we use the short-read data in: ${_polap_var_outdir_msbwt}"
    _short_read_data1=${_polap_var_outdir_msbwt}
    # check_file_existence "${_polap_var_outdir_genome_size}"
  elif [ -s "${_polap_var_project_sra_short_read}" ]; then
    local SRA=$(cut -f1 "${_polap_var_project_sra_short_read}")
    SR1="${_arg_outdir}/${SRA}_1.fastq"
    SR2="${_arg_outdir}/${SRA}_2.fastq"
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
    # 	seqkit stats -T "${SR1}" >"${_polap_var_outdir_fq_stats}"
    # fi
    # if [ -s "${SR2}" ]; then
    # 	seqkit stats -T "${SR2}" >>"${_polap_var_outdir_fq_stats}"
    # fi

    if [ ! -s "${SR1}" ] || [ ! -s "${SR2}" ]; then
      _polap_log1 "  downloading the paired-end short-read data: $SRA"
      "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$SRA"
      mv "${SRA}_1.fastq" "${_arg_outdir}/"
      mv "${SRA}_2.fastq" "${_arg_outdir}/"
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
  # if [ -s "${_polap_var_outdir_fq_stats}" ]; then
  # 	_polap_log2_file "${_polap_var_outdir_fq_stats}"
  # 	_polap_log3_cat "${_polap_var_outdir_fq_stats}"
  # else
  # 	_run_polap_summary-reads
  # fi

  # _run_polap_total-length-long
  # _run_polap_find-genome-size
  # _run_polap_reduce-data

  if [ -s "${_polap_var_outdir_msbwt}" ]; then
    _polap_log1 "  skipping the preparation of short-read polishing ..."
  else
    if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
      _polap_log1 "  decompressing ${_polap_var_outdir_msbwt_tar_gz} ... later when we polish it with the short-read data."
      # tar zxf "${_polap_var_outdir_msbwt_tar_gz}"
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
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# called after bioproject-prepare
################################################################################
function _run_polap_x-assemble-draft { # called after bioproject-prepare
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  help_message=$(
    cat <<HEREDOC
# Runs the POLAP organelle-genome assembly with sequencing data.
# 
# Arguments:
#   -o ${_arg_outdir}
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
  mkdir -p "${_arg_outdir}"

  # Run assembly, annotation, and contig selection steps
  if [ -s "${_polap_var_wga_contigger_edges_gfa}" ]; then
    _polap_log1 "  skipping the whole-genome assembly"
  else
    _run_polap_assemble1
  fi

  if [ -s "${_polap_var_wga_annotation_all}" ]; then
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
    _arg_inum=0
    FDIR="${_arg_outdir}/${_arg_inum}"
    _arg_jnum="${i}"
    file="$FDIR"/mt.contig.name-${_arg_jnum}
    mkdir -p "${_arg_outdir}/${_arg_jnum}"

    unique_file=$(is_unique "$file")
    if [ $? -eq 0 ]; then
      # If unique, add it to the unique_files array
      unique_files+=("$file")
      echo "$file is unique."

      _polap_var_mtcontigname="$file"
      # check the mt.contig.name-1
      if [ -s "$_polap_var_mtcontigname" ]; then
        # Run secondary assembly, polishing, and mtDNA selection steps
        _polap_log1_file "${_polap_var_mtcontigname}"

        if [ -s "${_arg_outdir}/${_arg_jnum}/30-contigger/graph_final.gfa" ] && [ "${_arg_redo}" = "off" ]; then
          _polap_log2 "  skipping organelle-genome assembly -i 0 -j ${_arg_jnum} ..."
        else
          _polap_log2 "  not skipping organelle-genome assembly -i 0 -j ${_arg_jnum} ..."
          _arg_yes="on"
          _run_polap_assemble2
        fi

        if [ -s "${_arg_outdir}/${_arg_jnum}/contig-annotation-table.txt" ] && [ "${_arg_redo}" = "off" ]; then
          _polap_log2 "  skipping organelle-genome annotation -i ${_arg_jnum} ..."
        else
          _arg_inum="${i}" _run_polap_annotate
        fi

        if [ -s "${_arg_outdir}/${_arg_jnum}/assembly_graph.gfa" ] && [ "${_arg_redo}" = "off" ]; then
          _polap_log2 "  skipping organelle-genome long-read polishing -i ${_arg_jnum} ..."
        else
          _arg_jnum="${i}" _run_polap_flye-polishing
        fi

        if [ -s "${_arg_outdir}/${_arg_jnum}/mt.0.fasta" ] && [ "${_arg_redo}" = "off" ]; then
          _polap_log2 "  skipping organelle-genome extraction -i ${_arg_jnum} ..."
        else
          _polap_log2 "  not implemented yet; skipping organelle-genome extraction -i ${_arg_jnum} ..."
          # _arg_inum="${i}" _run_polap_select-mtdna
        fi
      else
        _polap_log1 "LOG: $_polap_var_mtcontigname is empty for select-contig type $i ..."
      fi
    else
      _polap_log1 "$file is the same as $unique_file."
    fi
  done

  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Postprocess the bioproject assembly
################################################################################
function _run_polap_bioproject-postprocess { # postprocess the bioproject assembly
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set variables for file paths
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  # Set the output directory for the current job number
  if [[ -s "${_polap_var_project_txt}" ]]; then
    _arg_bioproject=$(<"${_polap_var_project_txt}")
    local BIOPRJ="${_arg_bioproject}"
  else
    _polap_log0 "NOTE: you have no ${_polap_var_project_txt}"
    _polap_log0 "  save BioProject ID in ${_polap_var_project_txt}"
    _polap_log0 "  if you already have the folder ${_polap_var_project}"
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
#   -o ${_arg_outdir}: output folder for BioProject with bioproject.txt
#   or
#   -b or --bioproject ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${_arg_outdir}/bioproject.txt
#   or
#   -b or --bioproject <BioProject ID>
#
#   long-read input data in the order of priority
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_polap_var_outdir_l_fq_gz}
#   ${_polap_var_project_sra_long_read}
#   to download the data if no such file
#
#   short-read input data in the order of priority
#   ${_polap_var_outdir_msbwt_tar_gz}
# 	${_polap_var_outdir_msbwt}
#   ${_polap_var_project_sra_short_read} 
#   to download the data if no such file
# Outputs:
#   Long- and short-read data files
Example: $(basename $0) ${_arg_menu[0]} -o ${_arg_outdir}
Example: $(basename $0) ${_arg_menu[0]} -o ${_arg_outdir} [-b ${BIOPRJ}]
HEREDOC
  )

  LRNK="${_arg_outdir}/nk.fq.gz"

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

  _polap_log0 "protprocessing the organelle genome assembly of BioProject ${BIOPRJ} ..."

  if [ -s "${_polap_var_outdir_msbwt}" ]; then
    _polap_log1 "  skipping the preparation of short-read polishing ..."
  else
    if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
      _polap_log1 "  decompressing ${_polap_var_outdir_msbwt_tar_gz} ..."
      tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
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
    PA="${_arg_outdir}/${i}/mt.0.fasta"
    FA="${_arg_outdir}/${i}/mt.1.fa"
    check_file_existence "${_polap_var_outdir_msbwt}"

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
  # _arg_inum=0
  # if [ -s "${_polap_var_project_mtdna_fasta2}" ]; then
  # 	_polap_log2 "  skipping downloading mtDNA from NCBI"
  # 	_polap_log2 "  we use the mtDNA: ${_polap_var_project_mtdna_fasta2}"
  # else
  # 	_run_polap_get-mtdna
  #
  # fi

  if [[ -s "${_polap_var_project_mtdna_fasta2}" ]]; then

    # Compare the known mtDNA and the assembled one.
    for i in "${_arg_select_contig_numbers[@]}"; do
      # Define the paths for mtDNA sequences to be polished
      FA="${_arg_outdir}/${i}/mt.1.fa"

      if [ -s "${FA}" ]; then
        _arg_inum="${i}"
        _run_polap_compare-mtdna
      else
        _polap_log1 "  skipping the short-read polishing ..."
        _arg_inum="${i}"
        source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
        local n1=$(cut -f1 "${_polap_var_project_mtdna_fasta2_accession}")
        local l1="0"
        local l2="0"
        local c1="0"
        printf "%s\t%d\t%d\t%f\n" ${n1} ${l1} ${l2} ${c1} >"${_polap_var_mtdna_compare}"
      fi
    done
  else
    _polap_log0 "No known mtDNA for $(<${_polap_var_project_species})"
  fi

  # Finalize the assemble-bioproject function.
  _polap_log1 "rsync back to thorne ..."
  touch "${_arg_outdir}/log-assembled.txt"
  # touch "$BIOPRJ/log-copying-to-xxx.txt"
  # touch "$BIOPRJ/log-not-started-yet.txt"
  if [ "$(hostname)" = "thorne" ]; then
    if [[ "${PWD}" = "/home/goshng/run/polap" ]]; then
      cp -pr "${_arg_outdir}/" "/home/goshng/all/polap/figshare/bioprojects/"
      touch "/home/goshng/all/polap/figshare/bioprojects/log-assembled-${_arg_outdir}.txt"
    fi
  else
    touch "log-assembled-${_arg_outdir}.txt"
    scp -pq "log-assembled-${_arg_outdir}.txt" thorne:"$PWD/"
    rsync -a -e ssh "${_arg_outdir}/" thorne:"$PWD/${_arg_outdir}/"
  fi

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Assembles organelle genome from a BioProject using long and short reads.
# Compares the known and assembled mtDNA using BLAST.
#
# bioproject.txt is necessary.
################################################################################
function _run_polap_assemble-bioproject { # main function for this bioproject module
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set variables for file paths
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  # Set the output directory for the current job number
  if [[ -s "${_polap_var_project_txt}" ]]; then
    _arg_bioproject=$(<"${_polap_var_project_txt}")
    local BIOPRJ="${_arg_bioproject}"
  else
    _polap_log0 "NOTE: you have no ${_polap_var_project_txt}"
    _polap_log0 "  save BioProject ID in ${_polap_var_project_txt}"
    _polap_log0 "  if you already have the folder ${_polap_var_project}"
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
#   -o ${_arg_outdir}: output folder for BioProject with bioproject.txt
#   or
#   -b or --bioproject ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${_arg_outdir}/bioproject.txt
#   or
#   -b or --bioproject <BioProject ID>
#
#   long-read input data in the order of priority
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_polap_var_outdir_l_fq_gz}
#   ${_polap_var_project_sra_long_read}
#   to download the data if no such file
#
#   short-read input data in the order of priority
#   ${_polap_var_outdir_msbwt_tar_gz}
# 	${_polap_var_outdir_msbwt}
#   ${_polap_var_project_sra_short_read} 
#   to download the data if no such file
# Outputs:
#   mt.1.fa
Example: $(basename $0) ${_arg_menu[0]} -o ${_arg_outdir} [-b ${BIOPRJ}]
HEREDOC
  )

  LRNK="${_arg_outdir}/nk.fq.gz"

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

  _polap_log0 "assembling organelle genomes of BioProject ${BIOPRJ} ..."

  _run_polap_bioproject-prepare

  _run_polap_x-assemble-draft

  _run_polap_bioproject-postprocess

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Download runinfo for a list of BioProject accessions.
# This function is used to prepare the runinfo for the BioProject.
################################################################################
function _run_polap_download-runinfo-bioproject { # download runinfo for a list of BioProject
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Print help message if requested
  help_message=$(
    cat <<HEREDOC
Downloads runinfo for a list of BioProject accessions.

Inputs:
a text file with BioProject accessions, one per line.

Outputs:
files with runinfo for the BioProject accessions.

Example: $(basename $0) ${_arg_menu[0]} <bioproject_accessions.txt>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then

    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
    exit $EXIT_SUCCESS
  fi

  local input_file="${_arg_menu[1]}"
  if [[ ! -s "${input_file}" ]]; then
    _polap_log0 "ERROR: no such file: ${input_file}"
    _polap_log0 "  please provide a file with BioProject accessions, one per line."
    _polap_log0 "  Example: $(basename $0) ${_arg_menu[0]} <bioproject_accessions.txt>"
  fi

  # while IFS= read -r BIOPRJ; do
  # 	_polap_log0 "Downloading runinfo for BioProject: ${BIOPRJ}"
  # done < "${input_file}"
  # return 0

  # Some of the commands in polap-ncbitools seems to steal the stdin.
  # So we use a file descriptor to read the input file.
  # Because we use the file descriptor number 3 for log,
  # we use file descriptor number 4 for reading the input file.
  while IFS= read -r BIOPRJ <&4; do
    [[ -z "$BIOPRJ" ]] && continue # Skip blank lines
    local output_file="${_arg_outdir}/${BIOPRJ}.csv"
    if [[ -f "${output_file}" ]]; then
      _polap_log1 "  already exists: ${output_file}"
      continue
    else
      # Download the runinfo for the BioProject
      _polap_log0 "Downloading runinfo for BioProject: ${BIOPRJ}"
      "${_polap_script_bin_dir}"/polap-ncbitools fetch runinfo "$BIOPRJ" >"${output_file}"
      if [[ $? -ne 0 ]]; then
        _polap_log0 "INFO: failed to download runinfo for BioProject: ${BIOPRJ}"
        continue
      fi
      if [[ -s "${output_file}" ]]; then
        _polap_log1 "  runinfo for ${BIOPRJ} saved to ${output_file}"
      else
        _polap_log1 "INFO: runinfo file is empty for BioProject: ${BIOPRJ}"
      fi
    fi
    sleep 30 # Sleep for a second to avoid overwhelming the server
  done 4<"${input_file}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}
