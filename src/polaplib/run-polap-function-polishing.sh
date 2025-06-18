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
# Polap uses msbwt and fmlrc for polishing a draft organelle genome sequence.
# This script has such subcommands.
#
# TODO: We need to explore different error correction or polihsing tools.
#
# function _run_polap_prepare-polishing { # prepare the polishing using FMLRC
# function _run_polap_prepare-polishing_v2 { # prepare the polishing using FMLRC
# function _run_polap_polish { # polish organelle genome sequences using FMLRC
# function _run_polap_polish_v2 { # polish organelle genome sequences using FMLRC
# function _run_polap_prepare-polishing_sort_efficient { # prepare the polishing using FMLRC
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

# NOTE: the following disassemble has its own polishing script, which is using
# Version 1.
# Look for these lines in run-polap-function-disassemble.sh
# command time -v bash "${_POLAPLIB_DIR}"/polap-build-msbwt.sh \
# command time -v fmlrc -p "${_arg_threads}" \
# The two lines above are for subsampling-polishing.
# The simple-polishing of disassemble uses _run_polap_polish the Version 1.
# Have a separate msbwt folder like msbwt2: polap-variables-common.sh
# local _polap_var_outdir_msbwt2_dir="${_polap_var_outdir}/msbwt2"
# local _polap_var_outdir_msbwt2="${_polap_var_outdir}/msbwt2/comp_msbwt2.npy"
# local _polap_var_outdir_msbwt2_tar_gz="${_polap_var_outdir}/msbwt2.tar.gz"

################################################################################
# 2025-05-27
#
# function _run_polap_prepare-polishing { # prepare the polishing using FMLRC
#   -> Version 1 has been used since the beginning of the polap development.
#
# function _run_polap_prepare-polishing_v2 { # prepare the polishing using FMLRC
#   -> Version 2 needs testing.
#
# function _run_polap_polish { # polish organelle genome sequences using FMLRC
#   -> Version 1 has been used since the beginning of the polap development.
#
# function _run_polap_polish_v2 { # polish organelle genome sequences using FMLRC
#   -> Version 2 needs testing.
#
# function _run_polap_prepare-polishing_sort_efficient { # prepare the polishing using FMLRC
#   -> Version 1 with sort changed.
#   -> The code is more original than Version 1 itself except the sort line of code.
#
################################################################################

################################################################################
# Prepares the polishing using FMLRC.
#
# 2025-05-16: add the memory and CPU tracking
#
# Arguments:
#   -a s1.fq
#   -b s2.fq
# Inputs:
#   s1.fq
#   s2.fq
# Outputs:
#   $${_arg_outdir}/msbwt
################################################################################

summarize_polish_memlog_file() {
  local _memlog_file="$1"
  local _summary_base="$2" # e.g., outdir/summary-polish.txt
  local _start_ts="${3:-}"
  local _end_ts="${4:-}"
  local _pipeline_log="${5:-}"

  local _timestamp
  _timestamp=$(date +"%Y%m%d_%H%M")
  local _summary_file="${_summary_base}.${_timestamp}.txt"

  local _memlog_basename
  _memlog_basename="$(basename "$_memlog_file")"

  awk -F',' -v start="$_start_ts" -v end="$_end_ts" -v file="$_memlog_basename" '
    function format_duration(seconds) {
      h = int(seconds / 3600)
      m = int((seconds % 3600) / 60)
      s = seconds % 60
      return sprintf("%d:%02d:%02d (%.2f h)", h, m, s, seconds / 3600)
    }

    NR == 1 { next }
    $2 > max {
      max = $2
      ts = $1
      load = $4
      disk = $5
      cmd = $6
    }
    END {
      cmd_str = "date -d @" ts " +\"%Y-%m-%d %H:%M:%S\""
      cmd_str | getline ts_human
      close(cmd_str)

      print "📊 Peak total PSS during polishing:"
      print "File:             " file
      print "Timestamp:        " ts " (" ts_human ")"
      printf "Peak PSS:         %d KB (%.2f GB)\n", max, max / 1048576

      if (start && end) {
        elapsed = end - start
        elapsed_str = format_duration(elapsed)
        print "Elapsed time:     " elapsed_str
      } else {
        print "Elapsed time:     (not available)"
      }

      print "CPU load (1min):  " load
      print "Disk free (GB):   " disk
      print "Max PSS command:  " cmd
    }
  ' "$_memlog_file" | tee "$_summary_file"

  if [[ -n "$_pipeline_log" && -r "$_pipeline_log" ]]; then
    echo -e "\n📋 Full pipeline log:" >>"$_summary_file"
    cat "$_pipeline_log" >>"$_summary_file"
  fi

  # Create/update symlink to point to the latest summary
  ln -sf "$(basename "$_summary_file")" "$_summary_base"

  echo "[INFO] Summary written to: $_summary_file"
  echo "[INFO] Symlink created:    $_summary_base -> $(basename "$_summary_file")"
}

# The version without memory-usage tracking
################################################################################
# Prepares the polishing using FMLRC.
#
# 2025-05-16: add the memory and CPU tracking
#
# Arguments:
#   -a s1.fq
#   -b s2.fq
# Inputs:
#   s1.fq
#   s2.fq
# Outputs:
#   $${_arg_outdir}/msbwt
################################################################################
function _run_polap_prepare-polishing { # prepare the polishing using FMLRC
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  local filesize
  _polap_set-variables-short-read
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  help_message=$(
    cat <<HEREDOC
# Prepares the polishing using FMLRC.
# Arguments:
#   -a ${_arg_short_read1}
#   -b ${_arg_short_read2}
#   or
#   --bioproject ${_arg_bioproject}
# Inputs:
#   ${_arg_short_read1}
#   ${_arg_short_read2}
# Outputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
# Precondition:
#   get-bioproject --bioproject ${_arg_bioproject}
Example: $(basename "$0") ${_arg_menu[0]} -a ${_arg_short_read1} [-b ${_arg_short_read2}]
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then
    ls -l "${_arg_outdir}/msbwt" >&2
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
    exit $EXIT_SUCCESS
  fi

  if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
    # Get the file size in bytes
    filesize=$(stat --format=%s "${_polap_var_outdir_msbwt}")
  else
    filesize=0
  fi

  if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
    _polap_log1_file "${_polap_var_outdir_msbwt_tar_gz}"
    if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
      _polap_log1_file "${_polap_var_outdir_msbwt}"
      _polap_log1 "  skipping the short-read polishing preparation."
    else
      tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
    fi
  elif ((filesize > 1024)); then
    # Check if the file size is greater than 100 KB (100 * 1024 bytes)
    _polap_log1_file "${_polap_var_outdir_msbwt}"
    _polap_log1 "  skipping the short-read polishing preparation."
  else

    # Initialize Conda
    # _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
    # source $HOME/miniconda3/etc/profile.d/conda.sh
    # conda activate polap-fmlrc
    #
    # source $HOME/miniconda3/bin/activate polap-fmlrc

    _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
    if ! run_check2; then
      echoerr "ERROR: change your conda environment to polap-fmlrc."
      echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
      echoerr "INFO: (base) $ conda activate polap-fmlrc"
      exit $EXIT_ERROR
    fi

    # check_file_existence "${_arg_short_read1}"

    _polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"

    # one or two files
    # fq or gzipped fq
    # exist or not
    #!/bin/bash

    _polap_log1 "  decompressing short-read data if necessary..."
    if [[ -s "${_arg_short_read1}" ]]; then
      _arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
    else
      _arg_short_read1=""
    fi
    if [[ -s "${_arg_short_read2}" ]]; then
      _arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
    else
      _arg_short_read2=""
    fi
    if [[ -n "${_arg_short_read1}" ]]; then
      _polap_log1 "    short-read1: ${_arg_short_read1}"
    else
      _polap_log1 "    short-read1: no such data file"
    fi
    if [[ -n "${_arg_short_read2}" ]]; then
      _polap_log1 "    short-read2: ${_arg_short_read2}"
    else
      _polap_log1 "    short-read2: no such data file"
    fi

    if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
      _polap_log1 "    short-read1 file: ${_arg_short_read1}"
    else
      _polap_log1 "    short-read1: no fastq or fq file: ${_arg_short_read1}"
    fi

    if [[ ${_arg_short_read2} = *.fastq || ${_arg_short_read2} = *.fq ]]; then
      _polap_log1 "    short-read2 file: ${_arg_short_read2}"
    else
      _polap_log1 "    short-read2: no fastq or fq file: ${_arg_short_read2}"
    fi

    rm -rf "${_arg_outdir}/msbwt"

    # 2025-05-27
    # Version 1
    # mkdir -p "${_polap_var_outdir_msbwt_dir}"
    if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
      _polap_log1 "  creating ${_arg_outdir}/msbwt ..."
      cat "${_arg_short_read1}" "${_arg_short_read2:-/dev/null}" |
        awk 'NR % 4 == 2' | sort | tr NT TN |
        ropebwt2 -LR 2>"${_polap_output_dest}" |
        tr NT TN |
        msbwt convert "${_arg_outdir}"/msbwt \
          >/dev/null 2>&1
    elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
      _polap_log1 "  creating ${_arg_outdir}/msbwt ..."
      zcat "${_arg_short_read1}" "${_arg_short_read2}" |
        awk 'NR % 4 == 2' | sort | tr NT TN |
        ropebwt2 -LR 2>"${_polap_output_dest}" |
        tr NT TN |
        msbwt convert "${_arg_outdir}"/msbwt \
          >/dev/null 2>&1
    else
      _polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
      _polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
    fi

    # Version 2
    # https://github.com/HudsonAlpha/rust-msbwt
    # mkdir -p "${_polap_var_outdir_msbwt_dir}"
    # if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
    # 	_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
    #
    # 	msbwt2-build \
    # 		-o "${_polap_var_outdir_msbwt}" \
    # 		"${_arg_short_read1}" "${_arg_short_read2:-}"
    #
    # elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
    # 	_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
    #
    # 	msbwt2-build \
    # 		-o "${_polap_var_outdir_msbwt}" \
    # 		"${_arg_short_read1}" "${_arg_short_read2:-}"
    #
    # else
    # 	_polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
    # 	_polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
    # fi

    conda deactivate
  fi

  _polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

function _run_polap_prepare-polishing_v2 { # prepare the polishing using FMLRC
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  local filesize
  _polap_set-variables-short-read
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  help_message=$(
    cat <<HEREDOC
# Prepares the polishing using FMLRC.
# Arguments:
#   -a ${_arg_short_read1}
#   -b ${_arg_short_read2}
#   or
#   --bioproject ${_arg_bioproject}
# Inputs:
#   ${_arg_short_read1}
#   ${_arg_short_read2}
# Outputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
# Precondition:
#   get-bioproject --bioproject ${_arg_bioproject}
Example: $(basename "$0") ${_arg_menu[0]} -a ${_arg_short_read1} [-b ${_arg_short_read2}]
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then
    ls -l "${_arg_outdir}/msbwt" >&2
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
    exit $EXIT_SUCCESS
  fi

  if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
    # Get the file size in bytes
    filesize=$(stat --format=%s "${_polap_var_outdir_msbwt}")
  else
    filesize=0
  fi

  if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
    _polap_log1_file "${_polap_var_outdir_msbwt_tar_gz}"
    if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
      _polap_log1_file "${_polap_var_outdir_msbwt}"
      _polap_log1 "  skipping the short-read polishing preparation."
    else
      tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
    fi
  elif ((filesize > 1024)); then
    # Check if the file size is greater than 100 KB (100 * 1024 bytes)
    _polap_log1_file "${_polap_var_outdir_msbwt}"
    _polap_log1 "  skipping the short-read polishing preparation."
  else

    # File preparation
    # check_file_existence "${_arg_short_read1}"

    # one or two files
    # fq or gzipped fq
    # exist or not
    #!/bin/bash

    _polap_log1 "  decompressing short-read data if necessary..."
    if [[ -s "${_arg_short_read1}" ]]; then
      _arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
    else
      _arg_short_read1=""
    fi
    if [[ -s "${_arg_short_read2}" ]]; then
      _arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
    else
      _arg_short_read2=""
    fi
    if [[ -n "${_arg_short_read1}" ]]; then
      _polap_log1 "    short-read1: ${_arg_short_read1}"
    else
      _polap_log1 "    short-read1: no such data file"
    fi
    if [[ -n "${_arg_short_read2}" ]]; then
      _polap_log1 "    short-read2: ${_arg_short_read2}"
    else
      _polap_log1 "    short-read2: no such data file"
    fi

    if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
      _polap_log1 "    short-read1 file: ${_arg_short_read1}"
    else
      _polap_log1 "    short-read1: no fastq or fq file: ${_arg_short_read1}"
    fi

    if [[ ${_arg_short_read2} = *.fastq || ${_arg_short_read2} = *.fq ]]; then
      _polap_log1 "    short-read2 file: ${_arg_short_read2}"
    else
      _polap_log1 "    short-read2: no fastq or fq file: ${_arg_short_read2}"
    fi

    rm -rf "${_arg_outdir}/msbwt"

    # Initialize Conda
    _polap_lib_conda-ensure_conda_env polap-fmlrc2 || exit 1
    if ! run_check2_v2; then
      echoerr "ERROR: change your conda environment to polap-fmlrc."
      echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
      echoerr "INFO: (base) $ conda activate polap-fmlrc"
      exit $EXIT_ERROR
    fi

    _polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"

    # Version 2
    # https://github.com/HudsonAlpha/rust-msbwt
    mkdir -p "${_polap_var_outdir_msbwt_dir}"
    if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
      _polap_log1 "  creating ${_arg_outdir}/msbwt ..."

      msbwt2-build \
        -o "${_polap_var_outdir_msbwt}" \
        "${_arg_short_read1}" "${_arg_short_read2:-}"

    elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
      _polap_log1 "  creating ${_arg_outdir}/msbwt ..."

      msbwt2-build \
        -o "${_polap_var_outdir_msbwt}" \
        "${_arg_short_read1}" "${_arg_short_read2:-}"

    else
      _polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
      _polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
    fi

    conda deactivate
  fi

  _polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Polishes using FMLRC.
#
# Arguments:
#   -p mt.0.fasta
#   -f mt.1.fa
# Inputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
################################################################################
function _run_polap_polish { # polish organelle genome sequences using FMLRC
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  help_message=$(
    cat <<HEREDOC
# Polish a draft sequence using FMLRC.
#
# Arguments:
#   -p ${_arg_unpolished_fasta}: a long-read draft genome assembly
#   -f ${_arg_final_assembly}: a final genome assembly sequence name
# Inputs:
#   ${_polap_var_outdir_msbwt}
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
Example: $(basename "$0") ${_arg_menu[0]} -p ${_arg_unpolished_fasta} -f ${_arg_final_assembly}
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
  if ! run_check2; then
    _polap_log0 "ERROR: change your conda environment to polap-fmlrc."
    _polap_log0 "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
    _polap_log0 "INFO: (base) $ conda activate polap-fmlrc"
    exit $EXIT_ERROR
  fi

  if [[ ! -s "${_polap_var_outdir_msbwt}" ]]; then
    _polap_log0 "ERROR: no msbwt at ${_polap_var_outdir_msbwt}"
    _polap_log0 "HINT: $(basename "$0") prepare-polishing [-a s1.fq] [-b s2.fq]"
    exit $EXIT_ERROR
  fi

  _polap_log1 "INFO: executing fmlrc on the draft sequence ${_arg_unpolished_fasta} ... be patient!"
  if [[ -s "${_arg_unpolished_fasta}" ]]; then
    # fmlrc -p "${_arg_threads}" "${_polap_var_outdir_msbwt}" "${_arg_unpolished_fasta}" "${_arg_final_assembly}" >/dev/null 2>&1
    # fmlrc -p "${_arg_threads}" "${_polap_var_outdir_msbwt}" "${_arg_unpolished_fasta}" "${_arg_final_assembly}" >/dev/null 2>&1

    # 2025-05-27
    # Version 1
    _polap_log3_pipe "fmlrc -p ${_arg_threads_fmlrc} ${_polap_var_outdir_msbwt} ${_arg_unpolished_fasta} ${_arg_final_assembly} >/dev/null 2>&1"

  else
    _polap_log0 "ERROR: no unpolished fasta file: [${_arg_unpolished_fasta}]"
    exit $EXIT_ERROR
  fi

  conda deactivate

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Polishes using FMLRC2.
#
# Arguments:
#   -p mt.0.fasta
#   -f mt.1.fa
# Inputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
################################################################################
function _run_polap_polish_v2 { # polish organelle genome sequences using FMLRC
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  help_message=$(
    cat <<HEREDOC
# Polish a draft sequence using FMLRC2.
#
# Arguments:
#   -p ${_arg_unpolished_fasta}: a long-read draft genome assembly
#   -f ${_arg_final_assembly}: a final genome assembly sequence name
# Inputs:
#   ${_polap_var_outdir_msbwt}
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
Example: $(basename "$0") ${_arg_menu[0]} -p ${_arg_unpolished_fasta} -f ${_arg_final_assembly}
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap-fmlrc2 || exit 1
  if ! run_check2_v2; then
    _polap_log0 "ERROR: change your conda environment to polap-fmlrc."
    _polap_log0 "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
    _polap_log0 "INFO: (base) $ conda activate polap-fmlrc"
    exit $EXIT_ERROR
  fi

  if [[ ! -s "${_polap_var_outdir_msbwt}" ]]; then
    _polap_log0 "ERROR: no msbwt at ${_polap_var_outdir_msbwt}"
    _polap_log0 "HINT: $(basename "$0") prepare-polishing [-a s1.fq] [-b s2.fq]"
    exit $EXIT_ERROR
  fi

  _polap_log1 "INFO: executing fmlrc2 on the draft sequence ${_arg_unpolished_fasta} ... be patient!"
  if [[ -s "${_arg_unpolished_fasta}" ]]; then

    # Version 2
    # https://github.com/HudsonAlpha/fmlrc2
    fmlrc2 \
      -t "${_arg_threads_fmlrc}" \
      ${_polap_var_outdir_msbwt} \
      "${_arg_unpolished_fasta}" \
      "${_arg_final_assembly}"

  else
    _polap_log0 "ERROR: no unpolished fasta file: [${_arg_unpolished_fasta}]"
    exit $EXIT_ERROR
  fi

  conda deactivate

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

# 2025-05-27
# This used to be original version 1 with sort changed in the main script.
# The other versions are modified so that we could keep two versions.
# sort --buffer-size=1G --temporary-directory=/tmp |
function _run_polap_prepare-polishing_sort_efficient { # prepare the polishing using FMLRC
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  local filesize
  _polap_set-variables-short-read
  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  help_message=$(
    cat <<HEREDOC
# Prepares the polishing using FMLRC.
# Arguments:
#   -a ${_arg_short_read1}
#   -b ${_arg_short_read2}
#   or
#   --bioproject ${_arg_bioproject}
# Inputs:
#   ${_arg_short_read1}
#   ${_arg_short_read2}
# Outputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
# Precondition:
#   get-bioproject --bioproject ${_arg_bioproject}
Example: $(basename "$0") ${_arg_menu[0]} -a ${_arg_short_read1} [-b ${_arg_short_read2}]
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then
    ls -l "${_arg_outdir}/msbwt" >&2
    # Disable debugging if previously enabled
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
    exit $EXIT_SUCCESS
  fi

  if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
    # Get the file size in bytes
    filesize=$(stat --format=%s "${_polap_var_outdir_msbwt}")
  else
    filesize=0
  fi

  if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
    _polap_log1_file "${_polap_var_outdir_msbwt_tar_gz}"
    if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
      _polap_log1_file "${_polap_var_outdir_msbwt}"
      _polap_log1 "  skipping the short-read polishing preparation."
    else
      tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
    fi
  elif ((filesize > 1024)); then
    # Check if the file size is greater than 100 KB (100 * 1024 bytes)
    _polap_log1_file "${_polap_var_outdir_msbwt}"
    _polap_log1 "  skipping the short-read polishing preparation."
  else

    # Initialize Conda
    source $HOME/miniconda3/etc/profile.d/conda.sh
    conda activate polap-fmlrc
    # source $HOME/miniconda3/bin/activate polap-fmlrc

    if ! run_check2; then
      echoerr "ERROR: change your conda environment to polap-fmlrc."
      echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
      echoerr "INFO: (base) $ conda activate polap-fmlrc"
      exit $EXIT_ERROR
    fi

    # check_file_existence "${_arg_short_read1}"

    _polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"

    # one or two files
    # fq or gzipped fq
    # exist or not
    #!/bin/bash

    _polap_log1 "  decompressing short-read data if necessary..."
    if [[ -s "${_arg_short_read1}" ]]; then
      _arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
    else
      _arg_short_read1=""
    fi
    if [[ -s "${_arg_short_read2}" ]]; then
      _arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
    else
      _arg_short_read2=""
    fi
    if [[ -n "${_arg_short_read1}" ]]; then
      _polap_log1 "    short-read1: ${_arg_short_read1}"
    else
      _polap_log1 "    short-read1: no such data file"
    fi
    if [[ -n "${_arg_short_read2}" ]]; then
      _polap_log1 "    short-read2: ${_arg_short_read2}"
    else
      _polap_log1 "    short-read2: no such data file"
    fi

    if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
      _polap_log1 "    short-read1 file: ${_arg_short_read1}"
    else
      _polap_log1 "    short-read1: no fastq or fq file: ${_arg_short_read1}"
    fi

    if [[ ${_arg_short_read2} = *.fastq || ${_arg_short_read2} = *.fq ]]; then
      _polap_log1 "    short-read2 file: ${_arg_short_read2}"
    else
      _polap_log1 "    short-read2: no fastq or fq file: ${_arg_short_read2}"
    fi

    rm -rf "${_arg_outdir}/msbwt"
    if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
      _polap_log1 "  creating ${_arg_outdir}/msbwt ..."
      cat "${_arg_short_read1}" "${_arg_short_read2:-/dev/null}" |
        awk 'NR % 4 == 2' |
        sort --buffer-size=1G --temporary-directory=/tmp |
        tr NT TN |
        ropebwt2 -LR 2>"${_polap_output_dest}" |
        tr NT TN |
        msbwt convert "${_arg_outdir}"/msbwt \
          >/dev/null 2>&1
    elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
      _polap_log1 "  creating ${_arg_outdir}/msbwt ..."
      zcat "${_arg_short_read1}" "${_arg_short_read2}" |
        awk 'NR % 4 == 2' |
        sort --buffer-size=1G --temporary-directory=/tmp |
        tr NT TN |
        ropebwt2 -LR 2>"${_polap_output_dest}" |
        tr NT TN |
        msbwt convert "${_arg_outdir}"/msbwt \
          >/dev/null 2>&1
    else
      _polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
      _polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
    fi

    conda deactivate
  fi

  _polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}
