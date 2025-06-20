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

function _run_polap_archive-rsync-template { # archive a POLAP-cflye output folder for later use
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  if [ "${_arg_short_read1_is}" = "on" ]; then
    _arg_archive="${_arg_short_read1}"
  fi

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Print help message if requested
  help_message=$(
    cat <<HEREDOC
# Archive polap-cflye the ${_arg_outdir} folder to ${_arg_archive}
#
# First, use the following rsync command to copy a folder 
# to another new folder, including only files smaller than 2MB:
# rsync -aq --max-size=2M source_folder/ destination_folder/
# Second, use a template path file to copy more files.
#
# Arguments:
#   -o ${_arg_outdir}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
#   --template polap-archive-template.txt
#   --max-sizefile: 1M not implemented yet!
# Inputs:
#   ${_arg_outdir}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log1_log "archiving ${_arg_outdir} to ${_arg_archive} ... upto ${_arg_max_filesize}"

  # rsync -aq --max-size=2M source_folder/ destination_folder/
  _polap_lib_file-rsync "${_arg_outdir}" "${_arg_archive}" "${_arg_max_filesize}"
  _polap_lib_file-archive-folder \
    "${_arg_outdir}" "${_arg_archive}" "${_arg_template}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Archive the ${_arg_outdir} folder to ${_arg_archive}
################################################################################
function _run_polap_archive { # archive a POLAP output folder for later use
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  if [[ "${_arg_short_read1_is}" == "on" ]]; then
    _arg_archive="${_arg_short_read1}"
  elif [[ "${_arg_archive_is}" == "off" ]]; then
    _arg_archive="${_arg_outdir}-a"
  fi

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  source "${_POLAPLIB_DIR}/polap-package-common.sh"

  # Print help message if requested
  help_message=$(
    cat <<HEREDOC
# Archive the ${_arg_outdir} folder to ${_arg_archive}
#
# Step 1. package
# Step 2. do more
#
# Arguments:
#   -o ${_arg_outdir}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
# Inputs:
#   ${_arg_outdir}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
Example: $(basename "$0") ${_arg_menu[0]} cflye -o <folder1> -a <folder2>
Example: $(basename "$0") ${_arg_menu[0]} aflye -o <folder1> -a <folder2>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  if [[ "${_arg_menu[1]}" == "cflye" ]]; then
    _arg_template="${_POLAPLIB_DIR}/polap-template-cflye-archive-files.txt"
    # _arg_max_filesize="1M"
    _run_polap_archive-rsync-template
  elif [[ "${_arg_menu[1]}" == "aflye" ]]; then
    _arg_template="${_POLAPLIB_DIR}/polap-template-aflye-archive-files.txt"
    # _arg_max_filesize="1M"
    _run_polap_archive-rsync-template
  else

    _polap_log0_log "archiving ${_arg_outdir} to ${_arg_archive} with file size upto ${_arg_max_filesize}"

    _run_polap_package

    cp -pr "${_polap_var_outdir_msbwt_dir}" "${_ppack_var_outdir}"
  fi

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# To filter a GFA file by a specific depth range
################################################################################
function _polap_archive_gfa-depth-filtered {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  local _gfa_infile=$1
  local _depth_range_string=$2
  local _gfa_outfile=$3

  # Convert the string to an array by changing the IFS to a comma
  IFS=',' read -r -a _depth_values <<<"${_depth_range_string}"

  if [[ "${_depth_values[0]}" -gt 0 ]]; then
    local depth_lower="${_depth_values[0]}"
    local depth_upper="${_depth_values[1]}"
    _polap_log2 "  depth range for graph filtering: $depth_lower ~ $depth_upper"
  else
    die "ERROR: no depth ranges"
  fi

  local _gfa_all=$(_polap_create_tempfile)

  _polap_log1 "    step 3-1: creating GFA without sequence data: ${_gfa_all}"
  _polap_log1 "      input1: ${_gfa_infile}"
  _polap_log2 "      output: ${_gfa_all}"
  _polap_log3_pipe "gfatools view \
		-S ${_gfa_infile} \
		>${_gfa_all} \
		2>$_polap_output_dest"

  local _gfa_seq_part=$(_polap_create_tempfile)

  _polap_log1 "    step 3-2: extracting sequence part of GFA: ${_gfa_seq_part}"
  _polap_log2 "      input1: ${_gfa_all}"
  _polap_log2 "      output: ${_gfa_seq_part}"
  _polap_log3_pipe "grep ^S ${_gfa_all} >${_gfa_seq_part}"

  local _gfa_seq_filtered=$(_polap_create_tempfile)

  # Filter edges in GFA using depths.
  _polap_log1 "    step 3-3: filtering GFA sequence part using depth range"
  _polap_log2 "      input1: ${_gfa_seq_part}"
  _polap_log2 "      input2: depth range: ${depth_lower} ~ ${depth_upper}"
  _polap_log2 "      output: ${_gfa_seq_filtered}"
  _polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-depthfilter-gfa.R \
			--gfa ${_gfa_seq_part} \
      --lower-bound-depth ${depth_lower} \
      --upper-bound-depth ${depth_upper} \
			--out ${_gfa_seq_filtered} \
			2>$_polap_output_dest"

  local _gfa_seq_filtered_edge=$(_polap_create_tempfile)

  _polap_log2 "    step 4-1: subsetting GFA using the depth-filtered GFA sequence part"
  _polap_log2 "      input1: ${_gfa_seq_filtered}"
  _polap_log2 "      output: ${_gfa_outfile}"

  _polap_log3_pipe "cut -f1 \
    ${_gfa_seq_filtered} \
    >${_gfa_seq_filtered_edge}"

  _polap_log3_pipe "gfatools view \
		-l @${_gfa_seq_filtered_edge} \
		${_gfa_infile} \
		2>$_polap_output_dest \
		>${_gfa_outfile}"

  _polap_delete_tempfiles

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

function _run_polap_filter-gfa-with-edge { # archive a POLAP output folder for later use
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  source "${_POLAPLIB_DIR}/polap-package-common.sh"

  if [ "${_arg_short_read1_is}" = "on" ]; then
    _arg_archive="${_arg_short_read1}"
  fi

  # Print help message if requested
  help_message=$(
    cat <<HEREDOC
# Create graph_final-<number>.gfa and graph_final.gfa.without.sequences
# for the v0.2.6 version.
#
# Arguments:
#   -o ${_arg_outdir}: the source output folder to archive
# Inputs:
#   ${_arg_outdir}
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log0 "reducing the whole-genome assembly with the annotated and selected MT contigs ..."
  _polap_log1 "  input1: ${_polap_var_ga_annotation_table}"
  _polap_log1 "  input2: ${_polap_var_mtcontigname}"

  _arg_inum="0"
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
    # Copy all contents from the "whole-genome" assembly folder into a new location.
    _polap_log1 "  step 3-1: creating GFA without sequence data"
    _polap_log1 "    input1: ${_polap_var_ga_contigger_edges_gfa}"
    _polap_log2 "    output: ${_polap_var_ga_contigger_edges_gfa}.without.sequences"
    _polap_log3_pipe "gfatools view \
		    -S ${_polap_var_ga_contigger_edges_gfa} \
		    >${_polap_var_ga_contigger_edges_gfa}.without.sequences \
		    2>$_polap_output_dest"
  fi

  local _edge_names=$(_polap_create_tempfile)

  # Filter edges in GFA using depths.
  _polap_log1 "  get all edges"
  _polap_log2 "    input1: ${_polap_var_ga_annotation_table}"
  _polap_log2 "    output: ${_edge_names}"
  _polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-contig2edge.R \
		--table ${_polap_var_ga_annotation_table} \
		--out ${_edge_names} \
		2>$_polap_output_dest"

  # Copy the contents of the seed contigs for further analysis or processing.
  for _file in "${_polap_var_ga}"/mt.contig.name-*; do
    if [[ -f "${_file}" ]]; then
      # Extract the number part
      local number="${_file##*-}"
      echo "Processing file: $_file with number: $number"
      # Your processing code here
      local _file_selected_edges=$(_polap_create_tempfile)
      cat "${_file}" "${_edge_names}" | sort | uniq >"${_file_selected_edges}"

      _polap_log3_pipe "gfatools view \
		      -l @${_file_selected_edges} \
    	  	${_polap_var_ga_contigger_edges_gfa} \
    	  	2>$_polap_output_dest \
    	  	>${_polap_var_ga_contigger}/graph_final-${number}.gfa"
    fi
  done

  _polap_delete_tempfiles

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

function _run_polap_package { # archive a POLAP output folder for later use
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  if [ "${_arg_short_read1_is}" = "on" ]; then
    _arg_archive="${_arg_short_read1}"
  fi

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  source "${_POLAPLIB_DIR}/polap-package-common.sh"

  # Print help message if requested
  help_message=$(
    cat <<HEREDOC
# Package the ${_arg_outdir} folder to ${_arg_archive}
#
# Arguments:
#   -o ${_arg_outdir}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
# Inputs:
#   ${_arg_outdir}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log0 "packaging ${_arg_outdir} to ${_arg_archive} ..."

  # Create a directory structure that mirrors the existing hierarchy by duplicating all folders and subfolders.
  rsync -a -f"+ */" -f"- *" "${_arg_outdir}"/ "${_arg_archive}"/

  # Copy all contents from the base directory into a new location.
  cp -fp "${_polap_var_outdir_genome_size}" \
    "${_ppack_var_outdir_genome_size}" 2>"${_polap_output_dest}"
  cp -fp "${_polap_var_outdir_long_total_length}" \
    "${_ppack_var_outdir_long_total_length}" 2>"${_polap_output_dest}"
  cp -fp "${_polap_var_outdir}"/*.stats \
    "${_ppack_var_outdir}" 2>"${_polap_output_dest}"
  cp -fp "${_polap_var_outdir}"/*.random.* \
    "${_ppack_var_outdir}" 2>"${_polap_output_dest}"
  rsync -a "${_polap_var_project}"/ "${_ppack_var_project}"/ 2>"${_polap_output_dest}"

  # copy mt.contig.name files
  for dir in $(find ${_arg_outdir} -maxdepth 1 -type d -regex '.*/[0-9]+$' | sort); do
    local _assembly_number=$(basename "$dir")

    _arg_inum="${_assembly_number}"
    _arg_jnum="${_assembly_number}"
    source "${_POLAPLIB_DIR}/polap-variables-common.sh"
    source "${_POLAPLIB_DIR}/polap-package-common.sh"

    if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
      # Copy all contents from the "whole-genome" assembly folder into a new location.
      _polap_log1 "    step 3-1: creating GFA without sequence data"
      _polap_log1 "      input1: ${_polap_var_ga_contigger_edges_gfa}"
      _polap_log2 "      output: ${_ppack_var_ga_contigger_edges_gfa}"
      _polap_log3_pipe "gfatools view \
		    -S ${_polap_var_ga_contigger_edges_gfa} \
		    >${_ppack_var_ga_contigger_edges_gfa}.without.sequences \
		    2>$_polap_output_dest"
    fi

    if [[ -s "${_polap_var_ga_annotation}" ]]; then
      cut -f 1 ${_polap_var_ga_annotation} | grep edge |
        cat - "${_polap_var_ga}"/mt.contig.name-* |
        sort | uniq >"${_polap_var_ga}"/all.mt.contig.name

      _polap_log3_pipe "gfatools view \
		    -l @${_polap_var_ga}/all.mt.contig.name \
    	  ${_polap_var_ga_contigger_edges_gfa} \
    	  2>${_polap_output_dest} \
    	  >${_ppack_var_ga_contigger_edges_gfa}"
    fi

    # gfatools gfa2fa ${_ppack_var_ga_contigger_edges_gfa} \
    # 	>${_ppack_var_ga_contigger_edges_fasta}

    # _polap_archive_gfa-depth-filtered \
    # 	${_polap_var_ga_contigger_edges_gfa} \
    # 	1000,1500 \
    # 	${_ppack_var_assembly_graph_final_gfa}

    # Copy the contents of the seed contigs for further analysis or processing.
    for _file in "${_polap_var_ga}"/mt.contig.name-*; do
      if [[ -f "${_file}" ]]; then
        # Extract the number part
        local number="${_file##*-}"
        echo "Processing file: $_file with number: $number"
        # Your processing code here
        cp -fp "${_file}" "${_ppack_var_ga}" 2>/dev/null
        local _file_selected_edges=$(_polap_create_tempfile)
        local _annotation_depth_table_seed_target="${_polap_var_ga}/contig-annotation-depth-table-seed-${number}.txt"
        if [[ -s "${_annotation_depth_table_seed_target}" ]]; then
          tail -n +2 "${_annotation_depth_table_seed_target}" |
            cut -f1 >"${_file_selected_edges}"
        else
          _polap_log0 "  no depth-table-seed-target: ${_annotation_depth_table_seed_target}"
          cp -fp "${_file}" "${_file_selected_edges}"
        fi

        # _polap_log3_pipe "gfatools view \
        #     -l @${_file_selected_edges} \
        # 	  	${_polap_var_ga_contigger_edges_gfa} \
        # 	  	2>$_polap_output_dest \
        # 	  	>${_ppack_var_ga_contigger}/graph_final-${number}.gfa"

        cp -fp -t "${_ppack_var_ga}" \
          "${_polap_var_ga_annotation_all}" \
          "${_polap_var_ga_annotation}" \
          "${_polap_var_ga_annotation_depth_table}" 2>/dev/null
        cp -fp -t "${_ppack_var_ga}" \
          "${_polap_var_ga}"/*table-seed-*.txt 2>/dev/null
      fi
    done

    # Check if basename is a positive number
    if [[ "${_assembly_number}" =~ ^[0-9]+$ ]] && [ "${_assembly_number}" -gt 0 ]; then
      _polap_log0 "Processing directory: $dir with basename: ${_assembly_number}"
      _arg_inum="${_assembly_number}"
      _arg_jnum="${_assembly_number}"
      source "${_POLAPLIB_DIR}/polap-variables-common.sh"
      source "${_POLAPLIB_DIR}/polap-package-common.sh"
      # Your processing code here
      rsync -a "${_polap_var_oga}/01-contig"/*.txt \
        "${_ppack_var_oga}/01-contig/" 2>/dev/null
      rsync -a "${_polap_var_oga}/02-reads/" \
        "${_ppack_var_oga}/02-reads/" 2>/dev/null
      rsync -a "${_polap_var_oga}/04-subsample/" \
        "${_ppack_var_oga}/04-subsample/" 2>/dev/null
      find ${_ppack_var_oga}/04-subsample -type f -name "*.gz" -delete
      rsync -a "${_polap_var_oga}/06-summary/" \
        "${_ppack_var_oga}/06-summary/" 2>/dev/null
      rsync -a "${_polap_var_oga}/07-plot/" \
        "${_ppack_var_oga}/07-plot/" 2>/dev/null
      cp -fp "${_polap_var_oga_assembly_graph_gfa}" \
        "${_ppack_var_oga}" 2>/dev/null
      cp -fp "${_polap_var_oga}"/*.png \
        "${_ppack_var_oga}" 2>/dev/null
      cp -fp "${_polap_var_ga_contigger_edges_gfa}" \
        "${_ppack_var_ga_contigger_edges_gfa}"
    fi
  done

  cp -pf "${LOG_FILE}" "${_arg_archive}"

  find ${_arg_archive} -type d -empty -delete

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Clean up the ${_arg_outdir}.
################################################################################
function _run_polap_cleanup { # cleanup an POLAP output folder
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Print help message if requested
  help_message=$(
    cat <<HEREDOC
# Clean up the ${_arg_outdir}.
#
# Arguments:
#   -o ${_arg_outdir}: the output directory
#   -j ${_arg_jnum}: the target assembly
# Inputs:
#   ${_arg_outdir}
# Outputs:
#
Example: $(basename "$0") ${_arg_menu[0]} -j <arg>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log3_cmd rm -f "${_polap_var_outdir_lk_fq_gz}"
  _polap_log3_cmd rm -rf "${_polap_var_oga}"/{00-assembly,10-consensus,20-repeat,40-polishing}
  _polap_log3_cmd rm -rf "${_polap_var_oga_contig}"/{contig.paf,contig.tab}
  for _pread_sel in ptgaul-intra-base-length single-intra-base-length combined-intra-base-length; do
    _polap_log3_cmd rm -rf "${_polap_var_oga_reads}/${_pread_sel}"
    _polap_log3_cmd rm -rf "${_polap_var_oga_seeds}/${_pread_sel}"
    _polap_log3_cmd rm -rf "${_polap_var_oga_subsample}/${_pread_sel}"
  done

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"
source "${_POLAPLIB_DIR}/run-polap-function-menus.sh"

################################################################################
# Initializes polap analysis in a starting folder,
#
# creating an output folder.
# Arguments:
#   -o ${_arg_outdir}
# Inputs: nothing
# Outputs:
#   ${_arg_outdir}
################################################################################
function _run_polap_init { # initialize an output folder
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  help_message=$(
    cat <<HEREDOC
# A polap analysis process is initiated in a specified starting folder
# through the creation of an output folder.
#
# 1. Create an empty folder.
# 2. Create empty menu files to facilitate effortless command inputting.
# 3. Record all external software packages including their versions.
#
# Arguments:
#   -o ${_arg_outdir}: the output folder
# Inputs: none
# Outputs:
#   ${_arg_outdir}
Example: $(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir}
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  if ! run_check1; then
    error_polap_conda
    exit $EXIT_ERROR
  fi

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then
    if [[ -d "${_arg_outdir}" ]]; then
      ls -l "${_arg_outdir}" >&3
    else
      _polap_log0 "No such output folder: ${_arg_outdir}"
    fi
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return
  fi

  mkdir -p "${_arg_outdir}"
  _polap_log0 "creating output folder [${_arg_outdir}] if no such folder exists ..."
  if [ "${_arg_outdir}" != "o" ]; then
    _polap_log1 "  Use -o ${_arg_outdir} option in all subsequent analysis"
    _polap_log1 "  because your output folder is not the default of 'o'."
  fi
  _run_polap_make-menus
  _log_command_versions

  _polap_log1 "NEXT: $(basename "$0") summary-reads -o ${_arg_outdir} -l ${_arg_long_reads}"
  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# View the polap log file.
################################################################################
function _run_polap_log { # display the polap log
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  help_message=$(
    cat <<HEREDOC
# Display the polap log.
#
# Arguments:
#   -o ${_arg_outdir}: the output folder
#   --log <FILE>
# Inputs:
#   ${LOG_FILE}
# Outputs:
#   a page view of the log file: ${LOG_FILE}
Example: $(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir}
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then
    if [[ -d "${_arg_outdir}" ]]; then
      ls -l "${_arg_outdir}" >&3
    else
      _polap_log0 "No such output folder: ${_arg_outdir}"
    fi
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  less "${LOG_FILE}" >&3

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}
