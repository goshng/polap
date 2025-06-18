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
# Polap subcommand, disassemble, assembles ptDNAs using a sampling approach.
#
# See Also:
# polap-function-disassemble-seeds.sh
# polap-r-disassemble-man-benchmark-boxplots.R
# run-polap-function-disassemble.sh
# run-polap-r-disassemble.R
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

# Function to extract the number from the filename
_disassemble_extract_number_circular_path() {
  local filename="$1"
  # Extract the base name (remove leading paths)
  local basename=$(basename "$filename")
  # Extract the number using sed
  local number=$(echo "$basename" | sed -n 's/.*circular_path_\([0-9]*\)_concatenated.*/\1/p')
  echo "$number"
}

# Function to adjust a value of _alpha
# Example usage
# _alpha=5.5
# _delta=2.0
#
# echo "Original _alpha: $_alpha"
#
# # Increase _alpha
# _alpha=$(_disassemble_adjust_alpha_by_delta "$_alpha" "$_delta" "increase")
# echo "Updated _alpha (after increase): $_alpha"
#
# # Decrease _alpha
# _alpha=$(_disassemble_adjust_alpha_by_delta "$_alpha" "$_delta" "decrease")
# echo "Updated _alpha (after decrease): $_alpha"
#
# # Decrease _alpha to become negative and test revert
# _alpha=1.0
# _delta=2.0
# _alpha=$(_disassemble_adjust_alpha_by_delta "$_alpha" "$_delta" "decrease")
# echo "Updated _alpha (after revert): $_alpha"
# Function to adjust a value of _alpha
_disassemble_adjust_alpha_by_delta() {
  local _alpha=$1     # Current value of _alpha (float)
  local _delta=$2     # Delta value (float)
  local _direction=$3 # Direction: "increase" or "decrease"

  if [[ "$_direction" == "increase" ]]; then
    # Increase _alpha by _delta
    _alpha=$(echo "scale=5; $_alpha + $_delta" | bc)
  elif [[ "$_direction" == "decrease" ]]; then
    # Decrease _alpha by _delta
    _alpha=$(echo "scale=5; $_alpha - $_delta" | bc)
    # Check if _alpha is negative
    if (($(echo "$_alpha < 0" | bc -l))); then
      _alpha=$(echo "scale=5; $_alpha + $_delta" | bc)
    fi
  else
    die "Invalid direction. Use 'increase' or 'decrease'."
    return 1
  fi

  echo "$_alpha" # Output the updated value
}

# unzip a gzipped file leaving the input as is.
# do not unzip if the input is not a gzipped file.
#
# unzipped_file=$(_polap_gunzip_file "${_short_read1}")
# _rstatus="$?"
# if [[ "$_rstatus" -eq 0 ]]; then
#   _polap_log2 "  unzipped file: $unzipped_file"
#   _short_read1="$unzipped_file"
# fi
_polap_gunzip_file() {
  local input_file="$1"

  # Check if the input file exists
  if [[ ! -f "$input_file" ]]; then
    die "Error: File '$input_file' not found."
  fi

  # Check if the file is gzipped
  if file "$input_file" | grep -q "gzip compressed data"; then
    # Extract the file name without the .gz extension
    local output_file="${input_file%.gz}"

    # Unzip the file and keep the original
    if gunzip -c "$input_file" >"$output_file"; then
      echo "$output_file"
      return 0
    else
      die "ERROR: failed to unzip '$input_file'."
    fi
  else
    _polap_log3 "    file '$input_file' is not gzipped."
    return 1
  fi
}

# Estimate the genome size from a short-read dataset
#
# input1: a short-read fastq data
# input2: an output folder _outdir
# input3: redo default is on
# output: the genome size
# _outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
_polap_find-genome-size() {
  local _short_read1="$1"
  local _outdir="$2"
  local _redo="${3:-on}"
  local _rstatus
  local _outdir_genome_size
  local _outdir_jellyfish_out
  local _outdir_jellyfish_out_histo
  local unzipped_file
  local _EXPECTED_GENOME_SIZE
  local _expected_genome_size_bp

  # Set paths for bioproject data
  # source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
  # source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

  _polap_log2 "    estimate the genome size using the short-read data ${_short_read1} ..."

  if [[ -d "${_outdir}" ]]; then
    _polap_log3 "    output folder: ${_outdir}"
  else
    _polap_log3 mkdir -p "${_outdir}"
    mkdir -p "${_outdir}"
  fi
  check_file_existence "${_short_read1}"

  _polap_log3 "    input1: ${_short_read1}"
  _outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
  _outdir_jellyfish_out="${_outdir}/jellyfish_out"
  _outdir_jellyfish_out_histo="${_outdir}/jellyfish_out.histo"

  # See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
  if [ -s "${_outdir_genome_size}" ] && [ "${_redo}" = "off" ]; then
    _polap_log3 "    found: ${_outdir_genome_size}, so skipping the genome size estimation ..."
    _polap_log3_file "${_outdir_genome_size}"
  else
    if [ -s "${_outdir_jellyfish_out}" ] && [ "${_redo}" = "off" ]; then
      _polap_log3 "    found: ${_outdir_jellyfish_out}, so skipping the JellyFish counting ..."
      _polap_log3_file "${_outdir_jellyfish_out}"
    else
      unzipped_file=$(_polap_gunzip_file "${_short_read1}")
      _rstatus="$?"
      if [[ "$_rstatus" -eq 0 ]]; then
        _polap_log3 "    unzipped file: $unzipped_file"
        _short_read1="$unzipped_file"
      fi

      if [[ -s "${_short_read1}" ]]; then
        _polap_log3_pipe "command time -v jellyfish count \
					-t ${_arg_threads} -C -m 19 \
					-s ${_arg_jellyfish_s} \
					-o ${_outdir_jellyfish_out} \
					--min-qual-char=? \
          ${_short_read1} 2>${_outdir}/timing-jellyfish.txt"
      else
        die "ASSERT: we must have at least one short-read fastq file."
      fi
      check_file_existence "${_outdir_jellyfish_out}"
    fi

    check_file_existence "${_outdir_jellyfish_out}"
    _polap_log3_pipe "jellyfish histo \
			-o ${_outdir_jellyfish_out_histo} \
      ${_outdir_jellyfish_out} 2>${_outdir}/timing-jellyfish-histo.txt"
    _polap_log3_file "${_outdir_jellyfish_out_histo}"

    _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-jellyfish.R \
			${_outdir_jellyfish_out_histo} \
			${_outdir_genome_size}"
    # Check the exit status
    if [ $? -ne 0 ]; then
      # Take action if needed, e.g., logging, sending a notification, etc.
      _polap_log0 "ERROR: disk space might be full"
      die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
    fi
  fi
  _polap_log3 "    output: ${_outdir_genome_size}"
  _polap_log3_cat "${_outdir_genome_size}"

  local _EXPECTED_GENOME_SIZE=$(<"${_outdir_genome_size}")
  local _EXPECTED_GENOME_SIZE=${_EXPECTED_GENOME_SIZE%.*}
  local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
  _polap_log3 "    expected genome size using short-read data (bases): ${_expected_genome_size_bp}"

  _polap_log3_cmd rm -f "${_outdir_jellyfish_out}"
  _polap_log3_cmd rm -f "${_outdir_jellyfish_out_histo}"

  return 0
}

# _run_polap_step-disassemble-cflye <out> <long.fq> <genomesize> <alpha> <resume:on>
#
function _polap_disassemble-cflye {
  local out="${1}"
  local long_read="${2}"
  local expected_genome_size="${3}"
  local alpha="${4}"
  local resume="${5:-off}"
  local _command1

  if [[ "${_POLAP_RELEASE}" -eq 0 ]]; then
    _command1="command time -v $HOME/all/polap/Flye/bin/cflye"
  else
    _command1="command time -v cflye"
  fi
  _command1="command time -v cflye"

  _command1+=" \
          ${_arg_flye_data_type} \
          ${long_read} \
			    --out-dir ${out} \
			    --disjointig-min-coverage ${alpha} \
			    --threads ${_arg_threads}"

  _command1+=" \
			  --asm-coverage ${_arg_flye_asm_coverage} \
			  --genome-size ${expected_genome_size}"
  if [[ "${_arg_flye_asm_coverage}" -le 0 ]]; then
    die "ERROR: must be positive --flye-asm-coverage ${_arg_flye_asm_coverage}"
  fi

  if [[ "${_arg_contigger}" == "on" ]]; then
    _command1+=" \
		    --stop-after contigger"
  fi

  if [[ "${resume}" != "off" ]]; then
    _command1+=" \
		    --resume"
  fi

  if [[ "${_arg_debug}" == "on" ]]; then
    _command1+=" \
		      --debug"
  fi

  # _command1+=" \
  # 	    2>${_polap_output_dest}"
  _command1+=" \
		    2>${out}/timing-cflye.txt"

  # Initialize Conda and activate polap-cflye environment
  _polap_lib_conda-ensure_conda_env polap-cflye || exit 1

  _polap_log3_pipe "${_command1}"

  conda deactivate

  return 0
}

function _run_polap_step-disassemble-archive-cfile {
  local cfile="${1}"

  if [[ -s "${cfile}" ]]; then
    _polap_log3_cmd mkdir -p ${_arg_archive}/$(dirname "${cfile#*/}")
    _polap_log3_pipe "cp -p \
      ${cfile} \
      ${_arg_archive}/${cfile#*/}"
  fi
  return 0
}

function _disassemble_report1 {
  local _case_infer="${1:-0}"
  local s
  local s2
  local d
  local d1
  local d2

  s=${_disassemble_dir}/${_arg_disassemble_i}/1/summary1.txt
  s2=${_disassemble_dir}/${_arg_disassemble_i}/1/summary2.txt
  d=${_disassemble_dir}/${_arg_disassemble_i}/1/summary1.md
  # csvtk -t cut -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length \
  csvtk -t cut -f index,long_rate_sample,alpha,peak_ram_size_gb,draft_assembly_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length \
    "${s}" |
    csvtk -t round -n 4 -f long_rate_sample,coverage_ref,coverage_target |
    csvtk -t round -n 2 -f alpha |
    csvtk rename -t -f index,long_rate_sample,alpha,peak_ram_size_gb,draft_assembly_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length \
      -n I,Rate,Alpha,'Pmem','D',Time,'N','L','C','C_ref','C_target',Length \
      >"${s2}"
  # csvtk -t cut I,Rate,Alpha,'Pmem','G',Time,'N','L','C','C_ref','C_target',Length \

  # _polap_log0 "report1: case: ${_case_infer}"

  if [[ "${_case_infer}" == "0" ]]; then
    csvtk -t cut -f I,Rate,Alpha,'Pmem','D',Time,'N','L','C','C_ref','C_target',Length \
      "${s2}" |
      csvtk -t csv2md -a right - \
        >"${d}"
  else
    csvtk -t cut -f I,Rate,Alpha,'Pmem','D',Time,'N','L','C',Length \
      "${s2}" |
      csvtk -t csv2md -a right - \
        >"${d}"
  fi
  _polap_log1_head "${d}"
  d1="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-ordered.txt"
  d2="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-ordered.pdf"
  _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
        --table ${s} \
        --out ${d1} \
        --plot ${d2}"
  _rstatus="$?"
  if [[ "$_rstatus" -ne 0 ]]; then
    _polap_log0 "ERROR: run-polap-r-disassemble.R on ${s}"
  fi
}

function _disassemble_report2 {
  local _case_infer="${1:-0}"
  local s
  local s2
  local d
  local d1
  local d2

  # concatenate all summary1.txt to the summary table.
  s=${_disassemble_dir}/${_arg_disassemble_i}/2/summary1.txt
  s2=${_disassemble_dir}/${_arg_disassemble_i}/2/summary2.txt
  d=${_disassemble_dir}/${_arg_disassemble_i}/2/summary1.md
  csvtk -t cut -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length,pident \
    "${s}" |
    csvtk -t round -n 4 -f long_rate_sample,coverage_ref,coverage_target |
    csvtk -t round -n 2 -f alpha |
    csvtk rename -t -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length,pident \
      -n I,Rate,Alpha,'Pmem','G',Time,'N','L','C','C_ref','C_target',Length,'Identity' \
      >"${s2}"

  # _polap_log0 "report2: case: ${_case_infer}"

  if [[ "${_case_infer}" == "0" ]]; then

    csvtk -t cut -f I,Rate,Alpha,'Pmem','G',Time,'N','L','C','C_ref','C_target',Length,'Identity' \
      "${s2}" |
      csvtk -t csv2md -a right - \
        >"${d}"
  else
    csvtk -t cut -f I,Rate,Alpha,'Pmem','G',Time,'N','L','C',Length \
      "${s2}" |
      csvtk -t csv2md -a right - \
        >"${d}"
  fi
  _polap_log1_head "${d}"
  d1="${_disassemble_dir}/${_arg_disassemble_i}/2/summary1-ordered.txt"
  d2="${_disassemble_dir}/${_arg_disassemble_i}/2/summary1-ordered.pdf"
  _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
        --table ${s} \
        --out ${d1} \
        --plot ${d2}"
  _rstatus="$?"
  if [[ "$_rstatus" -ne 0 ]]; then
    _polap_log0 "ERROR: run-polap-r-disassemble.R on ${s}"
  fi
}

function _disassemble_report3y {
  s="${_disassemble_dir}/${_arg_disassemble_i}/3/summary1.txt"
  d="${_disassemble_dir}/${_arg_disassemble_i}/3/summary1y.md"
  if [[ -s "${s}" ]]; then
    csvtk -t cut -f index,sampling_datasize,short_rate_sample,time_prepare,memory_prepare,time_polishing,memory_polishing,length,pident \
      "${s}" |
      csvtk -t round -n 4 -f short_rate_sample |
      csvtk rename -t -f index,sampling_datasize,short_rate_sample,time_prepare,memory_prepare,time_polishing,memory_polishing,length,pident \
        -n I,Size,Rate,'Tp','Mp','Ts','Ms',Length,Pident |
      csvtk -t csv2md -a right - \
        >"${d}"
    _polap_log1_head "${d}"
    d1="${_disassemble_dir}/${_arg_disassemble_i}/3/summary1y-ordered.txt"
    d2="${_disassemble_dir}/${_arg_disassemble_i}/3/summary1y-ordered.pdf"
    _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
        --table ${s} \
        --out ${d1} \
        --plot ${d2} 2>/dev/null"
  fi
  return 0
}

# rewrite the stage 3 results
function _disassemble_report3 {
  local _brg_stage="${1:-3}"
  s="${_disassemble_dir}/${_arg_disassemble_i}/${_brg_stage}/summary1.txt"
  d="${_disassemble_dir}/${_arg_disassemble_i}/${_brg_stage}/summary1.md"
  if [[ -s "${s}" ]]; then
    csvtk -t cut -f index,short_rate_sample,sampling_datasize,short_sample_seed,time_prepare,memory_prepare,time_polishing,memory_polishing,pident,length \
      "${s}" |
      csvtk -t round -n 4 -f short_rate_sample |
      csvtk rename -t -f index,short_rate_sample,sampling_datasize,short_sample_seed,time_prepare,memory_prepare,time_polishing,memory_polishing,pident,length \
        -n I,Rate,Size,Seed,Tp,Mp,Ts,Ms,Pident,Length |
      csvtk -t csv2md -a right - \
        >"${d}"
    _polap_log1_head "${d}"
    d1="${_disassemble_dir}/${_arg_disassemble_i}/${_brg_stage}/summary1-ordered.txt"
    d2="${_disassemble_dir}/${_arg_disassemble_i}/${_brg_stage}/summary1-ordered.pdf"
    _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
		      --table ${s} \
		      --out ${d1} \
		      --plot ${d2} 2>/dev/null"
  fi
  return 0
}

_disassemble_find_folder_with_max_coverage() {
  local _var_mtdna="${1}"
  local max_value=-1
  local max_folder=""
  local j
  local value
  local folder
  local _ptdna

  shopt -s nullglob # Enable nullglob so that the loop gets no elements if no match
  for _ptdna in "${_var_mtdna}"/circular_path_*_concatenated.fa; do
    j=$(_disassemble_extract_number_circular_path "${_ptdna}")
    if [[ -f "${_var_mtdna}/${j}/coverage1.txt" ]]; then
      local value=$(<"${_var_mtdna}/${j}/coverage1.txt")
      if (($(echo "$value > $max_value" | bc -l))); then
        max_value=$value
        max_folder=$j
      fi
    fi
  done
  shopt -u nullglob # Restore default behavior

  if [[ -n $max_folder ]]; then
    echo "$max_folder" # Return the folder number with the largest value
  else
    _polap_log2 "  no valid coverage1.txt files found."
    return 1 # Indicate an error
  fi
}

# print the header of the polap-cflye main summary
function _polap_lib_table-cflye-main-summary-header {
  local _summary_file="${1}"

  # Summary
  #
  # long_total
  # long_rate_sample
  # long_sample_seed
  #
  # short_total
  # short_rate_sample
  # short_sample_seed
  #
  # _expected_genome_size
  #
  # peak_ram_size_gb
  #
  # alpha
  #
  # gfa_number_segments
  # gfa_number_links
  # gfa_total_segment_length
  #
  # note: annotation result e.g., number of MT/PT genes
  # note: seeds e.g., number of seeds
  #
  # num_circular_paths
  # num_circular_nodes
  #
  # coverage
  #
  # pident
  #
  # length
  # gc

  # Add header to the output file
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "index" \
    "sampling_datasize" \
    "long_total" "long_rate_sample" "long_sample_seed" \
    "short_total" "short_rate_sample" "short_sample_seed" \
    "expected_genome_size" \
    "peak_ram_size_gb" \
    "alpha" \
    "gfa_number_segments" "gfa_number_links" "gfa_total_segment_length" \
    "num_circular_paths" "num_circular_nodes" \
    "time" \
    "coverage_ref" \
    "coverage_target" \
    "pident" \
    "length" \
    "gc" \
    "j_candidate" \
    "draft_assembly_size" \
    >"${_summary_file}"
}

function _polap_lib_table-cflye-main-summary-content {
  local _summary_file="${1}"

  # Append the result to the output file
  # printf "index\tlong\trate\tsize\talpha\tmemory\tnsegments\tnlink\ttotalsegmentlength\tlength\tgc\tcoverage\n" >"${_summary1_file}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${_summary_i}" \
    "${_summary_sampling_datasize}" \
    "${_summary_long_total}" \
    "${_summary_long_rate_sample}" \
    "${_summary_long_sample_seed}" \
    "${_summary_short_total}" \
    "${_summary_short_rate_sample}" \
    "${_summary_short_sample_seed}" \
    "${_summary_genome_size}" \
    "${_summary_peak_ram_size_gb}" \
    "${_summary_pre_alpha}" \
    "${_summary_gfa_number_segments}" \
    "${_summary_gfa_number_links}" \
    "${_summary_gfa_total_segment_length}" \
    "${_summary_num_circular_paths}" \
    "${_summary_num_circular_nodes}" \
    "${_summary_elapsed_time}" \
    "${_summary_coverage_ref}" \
    "${_summary_coverage_target}" \
    "${_summary_pident}" \
    "${_summary_length}" \
    "${_summary_gc}" \
    "${_summary_j_candidate}" \
    "${_summary_draft_assembly_size}" \
    >>"${_summary_file}"
}

# print the header of the polap-cflye polish summary
function _polap_lib_table-cflye-polish-summary-header {
  local _summary_file="${1}"

  # add header to the output file
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "index" \
    "sampling_datasize" \
    "short_total" \
    "short_rate_sample" \
    "short_sample_seed" \
    "time_prepare" \
    "memory_prepare" \
    "time_polishing" \
    "memory_polishing" \
    "length" \
    "pident" \
    >"${_summary_file}"
}

function _polap_lib_table-cflye-polish-summary-content {
  local _summary_file="${1}"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${_summary_i}" \
    "${_summary_sampling_datasize}" \
    "${_summary_short_total}" \
    "${_summary_short_rate_sample}" \
    "${_summary_short_sample_seed}" \
    "${_summary_total_hours_msbwt}" \
    "${_summary_memory_gb_msbwt}" \
    "${_summary_total_hours_fmlrc}" \
    "${_summary_memory_gb_fmlrc}" \
    "${_summary_assembly_length}" \
    "${_summary_assembly_pident}" \
    >>"${_summary_file}"
}

# count the number of bases of a fastq and save the number in a text file
# arg1: long-read fastq
# arg2: long_total_length.txt
function _disassemble-step1 {
  local _input="$1"
  local _output="$2"

  if [[ -s "${_output}" &&
    "${_arg_redo}" == "off" ]]; then
    _polap_log2 "    found: ${_output}, skipping ..."
  else
    _polap_lib_fastq-total-length-of \
      "${_input}" \
      "${_output}"

  fi
  local _total_length_long=$(<"${_output}")
  local _total_length_long_bp=$(_polap_utility_convert_bp ${_total_length_long})
  _polap_log2 "    the size of the long-read dataset: ${_total_length_long_bp}"
}

# concatenate two or one short-read fastq files
# arg1: output concatenated fastq
# arg2: input short-read fastq 1
# arg3: input short-read fastq 2
function _disassemble-step2 {
  local _input_short_reads="${1}"
  local _input1_short_read="$2"
  local _input2_short_read="$3"

  if [[ -s "${_input_short_reads}" && "${_arg_redo}" = "off" ]]; then
    _polap_log2 "    found: ${_input_short_reads}, skipping ..."
  else
    if [[ -s "${_input1_short_read}" &&
      -s "${_input2_short_read}" ]]; then
      _polap_lib_fastq-concatenate-fastq-files "${_input_short_reads}" \
        "${_input1_short_read}" \
        "${_input2_short_read}"
    elif [[ -s "${_input2_short_read}" ]]; then
      _polap_lib_fastq-concatenate-fastq-files "${_input_short_reads}" \
        "${_input2_short_read}"
    else
      _polap_lib_fastq-concatenate-fastq-files "${_input_short_reads}" \
        "${_input1_short_read}"
    fi
  fi
}

# count bases of a short-read sequencing fastq file
# arg1: a short-read fastq file
function _disassemble-step3 {
  local _input="$1"
  local _output="$2"

  if [[ -s "${_output}" &&
    "${_arg_redo}" = "off" ]]; then
    _polap_log3 "    found: ${_output}, skipping ..."
  else
    _polap_lib_fastq-total-length-of \
      "${_input}" \
      "${_output}"
  fi
  local _total_length_short=$(<"${_output}")
  local _total_length_short_bp=$(_polap_utility_convert_bp ${_total_length_short})
  if [[ -n "${_total_length_short_bp}" ]]; then
    _polap_log3 "    the size of the short-read dataset: ${_total_length_short_bp}"
  else
    die "ERROR: short-read data size is not computed."
  fi
}

function _disassemble-step8 {
  local _outdir="${1}"
  local _long_read="${2}"
  local _summary_genome_size="${3}"
  local _alpha="${4}"
  local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
  local _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
  local _is_stop="off"

  local _summary_peak_ram_size_gb=-1
  local _summary_peak_ram=-1

  # no resume feature
  if [[ -d "${_outdir}/00-assembly" ]]; then
    _polap_log2 "  found: ${_outdir}/00-assembly, and cflye must have been tried so skipping it."
  else
    # cflye assemble stage
    _polap_disassemble-cflye \
      "${_outdir}" "${_long_read}" "${_summary_genome_size}" "${_alpha}"
  fi

  # INFO: I wanted to have a bash function for wrapping this part,
  # but it seems better what it is.
  local _memory_gb=-1
  local _total_hours=-1
  read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_outdir}/timing-cflye.txt")
  _summary_peak_ram_size_gb=$(_polap_lib_number-round-number "${_memory_gb}")

  # check the peak memory of cflye
  if [[ -d "${_outdir}/30-contigger" ]] &&
    [[ -s "${_contigger_edges_fasta}" ]]; then

    # _summary_peak_ram_size_gb=$(grep "Peak RAM usage" "${_outdir}/cflye.log" | awk 'NR==2 {print $7}')
    _summary_peak_ram=$(grep "Peak RAM usage" "${_outdir}/cflye.log" | awk 'NR==2 {print $7, $8}')

    if [[ -z "${_summary_peak_ram}" ]]; then
      die "ERROR: cflye failed to run."
    fi

    if ((_summary_peak_ram_size_gb < _arg_disassemble_memory)); then
      _polap_log2 "  cFlye: used memory is less than ${_arg_disassemble_memory} Gb ($i): $_summary_peak_ram_size_gb"
    else
      _polap_log2 "  cFlye: used memory is not less than ${_arg_disassemble_memory} Gb ($i): $_summary_peak_ram_size_gb"
      _polap_log2 "  exit the disassemble menu."
      _is_stop="on"
    fi

    # copy assembly_graph.gfa to graph_final.gfa if --no-contigger
    #
    if [[ "${_arg_contigger}" == "off" ]]; then
      if [[ -s "${_outdir}/assembly_graph.gfa" ]]; then
        if [[ ! -s "${_contigger_edges_gfa}.backup" ]]; then
          _polap_log2 "  backup ${_contigger_edges_gfa}"
          _polap_log3_cmd cp -p "${_contigger_edges_gfa}" "${_contigger_edges_gfa}.backup"
        fi
        _polap_log2 "  use ${_outdir}/assembly_graph.gfa by copying graph_final.gfa"
        _polap_log3_cmd cp -p "${_outdir}/assembly_graph.gfa" "${_contigger_edges_gfa}"
      else
        _polap_log0 "ERROR: cFlye did not finish with the polishing step."
        die "ERROR: no such file: ${_outdir}/assembly_graph.gfa"
      fi
    fi
  else
    _polap_log2 "  cFlye: no genome assembly"
  fi
  local arr=(
    "${_is_stop}"
    "${_summary_peak_ram_size_gb}"
  )
  echo "${arr[@]}"
}

function _disassemble-step9 {
  local _outdir="${1}"
  local _alpha="${2}"
  local _outdir_draft_assembly_size="${_outdir}/draft_assembly_size.txt"
  local _flye_draft_assembly_fasta="${_outdir}/00-assembly/draft_assembly.fasta"
  local _draft_assembly_size=0

  # get the total number of bases in the flye draft assembly sequences
  _polap_log3_cmd rm -f "${_outdir_draft_assembly_size}"
  if [[ -s "${_flye_draft_assembly_fasta}" ]]; then
    _polap_lib_fastq-total-length-of \
      "${_flye_draft_assembly_fasta}" \
      "${_outdir_draft_assembly_size}"

    _draft_assembly_size=$(<"${_outdir_draft_assembly_size}")
  fi

  # adjust alpha
  if ((_arg_disassemble_m < _draft_assembly_size)); then
    _alpha=$(_disassemble_adjust_alpha_by_delta "${_alpha}" \
      "${_arg_disassemble_delta}" \
      "increase")
    _polap_log2 "      alpha increased: ${_alpha}"
  else
    _alpha=$(_disassemble_adjust_alpha_by_delta "${_alpha}" \
      "${_arg_disassemble_delta}" \
      "decrease")
    _polap_log2 "      alpha decrease (if not negative): ${_alpha}"
  fi
  local arr=(
    "${_alpha}"
    "${_draft_assembly_size}"
  )
  echo "${arr[@]}"
}

function _disassemble-step10 {
  local _outdir="${1}"
  local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
  local _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"

  local _summary_gfa_number_segments=-1
  local _summary_gfa_number_links=-1
  local _summary_gfa_total_segment_length=-1

  if [[ -s "${_contigger_edges_gfa}" ]] &&
    [[ -s "${_contigger_edges_fasta}" ]]; then
    _polap_log3_pipe "gfatools stat -l ${_contigger_edges_gfa} \
          >${_outdir}/graph_final.txt 2>$_polap_output_dest"
    _summary_gfa_number_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')
    _summary_gfa_number_links=$(grep "Number of links:" "${_outdir}/graph_final.txt" | awk '{print $4}')
    _summary_gfa_total_segment_length=$(grep "Total segment length:" "${_outdir}/graph_final.txt" | awk '{print $4}')
  else
    if [[ ! -s "${_contigger_edges_gfa}" ]]; then
      _polap_log2 "    no input: ${_contigger_edges_gfa}"
    fi
    if [[ ! -s "${_contigger_edges_fasta}" ]]; then
      _polap_log2 "    no input: ${_contigger_edges_fasta}"
    fi
  fi

  local arr=(
    "${_summary_gfa_number_segments}"
    "${_summary_gfa_number_links}"
    "${_summary_gfa_total_segment_length}"
  )
  echo "${arr[@]}"
}

function _disassemble-step11 {
  local _outdir="${1}"

  local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
  local _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
  local _mtcontigname="${_outdir}/mt.contig.name"
  local _ga_annotation_all="${_outdir}/assembly_info_organelle_annotation_count-all.txt"

  if [[ -s "${_contigger_edges_gfa}" ]] &&
    [[ -s "${_contigger_edges_fasta}" ]]; then
    if [[ -s "${_ga_annotation_all}" ]]; then
      _polap_log2 "    found: ${_ga_annotation_all}, so skipping mt/pt annotation ..."
    else
      polap_annotate "${_contigger_edges_gfa}" "${_ga_annotation_all}"
    fi
  else
    # Check missing files in a loop
    for file in \
      "${_contigger_edges_gfa}" \
      "${_contigger_edges_fasta}"; do
      [[ ! -s "$file" ]] && _polap_log2 "    no input: $file"
    done
  fi
}

function _disassemble-step12 {
  local _outdir="${1}"

  local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
  local _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
  local _mtcontigname="${_outdir}/mt.contig.name"
  local _ga_annotation_all="${_outdir}/assembly_info_organelle_annotation_count-all.txt"

  if [[ -s "${_contigger_edges_gfa}" ]] &&
    [[ -s "${_contigger_edges_fasta}" ]] &&
    [[ -s "${_ga_annotation_all}" ]]; then
    # depth-filtering
    if [[ -s "${_mtcontigname}" ]]; then
      _polap_log2 "    found: ${_mtcontigname}, so skipping seeds ..."
    else
      polap_disassemble-seeds "${_contigger_edges_gfa}" \
        "${_ga_annotation_all}" \
        "${_mtcontigname}"
    fi
  else
    # Check missing files in a loop
    for file in \
      "${_contigger_edges_gfa}" \
      "${_contigger_edges_fasta}" \
      "${_ga_annotation_all}"; do
      [[ ! -s "$file" ]] && _polap_log2 "    no input: $file"
    done
  fi
}

# _summary_num_circular_paths
# -1 if no gfa
# -2 if too many segments
# some numbers if properly
function _disassemble-step13 {
  local _outdir="${1}"
  local _summary_num_circular_paths=-1
  local _summary_num_circular_nodes=-1
  local _fixed_upper_bound_number_segments=30 # the fixed number

  local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
  local _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
  local _mtcontigname="${_outdir}/mt.contig.name"
  local _var_mtdna="${_outdir}/52-mtdna"

  local _summary_gfa_number_segments=${_fixed_upper_bound_number_segments}
  if [[ -s "${_contigger_edges_gfa}" ]]; then
    gfatools stat -l "${_contigger_edges_gfa}" \
      >"${_outdir}/30-contigger/graph_final.txt" \
      2>"$_polap_output_dest"
    _summary_gfa_number_segments=$(grep "Number of segments:" "${_outdir}/30-contigger/graph_final.txt" | awk '{print $4}')
    if ((_summary_gfa_number_segments > _fixed_upper_bound_number_segments)); then
      _summary_num_circular_paths=-2
      _summary_num_circular_nodes=-2
    fi
  fi

  if ((_summary_gfa_number_segments <= _fixed_upper_bound_number_segments)) &&
    [[ -s "${_contigger_edges_gfa}" ]] &&
    [[ -s "${_contigger_edges_fasta}" ]] &&
    [[ -s "${_mtcontigname}" ]]; then
    _polap_log3_cmd rm -rf "${_var_mtdna}"
    _polap_log3_cmd mkdir -p "${_var_mtdna}"

    #

    _polap_log3_pipe "command time -v python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_var_mtdna} \
		        2>${_outdir}/timing-find-plastid.txt"
    if [[ -s "${_var_mtdna}/circular_path_count.txt" ]]; then
      _summary_num_circular_paths=$(<"${_var_mtdna}/circular_path_count.txt")
      _summary_num_circular_nodes=$(<"${_var_mtdna}/circular_path_nodes.txt")
      _polap_log3_cat "${_var_mtdna}/circular_path.txt"
      if [[ "${_summary_num_circular_paths}" -eq 4 ]] ||
        [[ "${_summary_num_circular_paths}" -eq 2 ]]; then
        _polap_log2 "    circular_path count is 2 or 4."
      else
        _polap_log2 "    circular_path count is neither 2 nor 4."
        _polap_log3_cmd rm -f "${_var_mtdna}"/circular_path_*_concatenated.fa
      fi
    fi
  else
    # Check missing files in a loop
    for file in \
      "${_contigger_edges_gfa}" \
      "${_contigger_edges_fasta}" \
      "${_mtcontigname}"; do
      [[ ! -s "$file" ]] && _polap_log2 "    no input: $file"
    done
  fi
  local arr=(
    "${_summary_num_circular_paths}"
    "${_summary_num_circular_nodes}"
  )
  echo "${arr[@]}"
}

function _disassemble-step14 {
  local _var_mtdna="${1}"
  local _circular_path_fasta="${_var_mtdna}/circular_path_1_concatenated.fa"
  local _arg_unpolished_fasta="${_var_mtdna}/ptdna.0.fa"
  local _arg_final_assembly="${_var_mtdna}/ptdna.1.fa"
  local _summary_coverage_ref=-1
  local _summary_coverage_target=-1
  local _summary_j_candidate=-1
  local num_circular_paths_file_count
  local _ptdir
  local j
  local _folder_with_max_coverage
  local _summary_coverage_ref
  local _summary_coverage_target
  local diff
  local _restarted_fasta

  num_circular_paths_file_count=$(find "${_var_mtdna}" -maxdepth 1 -name "circular_path_*_concatenated.fa" 2>/dev/null | wc -l)
  if [[ "${num_circular_paths_file_count}" -eq 4 ]] ||
    [[ "${num_circular_paths_file_count}" -eq 2 ]]; then
    _polap_log2 "    circular_path count is 2 or 4."
    _polap_log2 "    output1: unpolished: ${_arg_unpolished_fasta}"
    _polap_log2 "    output2: polished: ${_arg_final_assembly}"

    if [[ -n "${_arg_disassemble_c}" ]]; then
      _polap_log1 "  step 14-1: use alignment coverage to select the best ptDNA using reference: ${_arg_disassemble_c}"
      _polap_log2 "    input1: ${_arg_disassemble_c}"
      _polap_log2 "    input2: e.g., ${_circular_path_fasta}"
      _polap_log2 "    outdir: ${_var_mtdna}"

      # 14-1. compute the coverage for each candidate
      shopt -s nullglob # Enable nullglob so that the loop gets no elements if no match
      for _ptdna in "${_var_mtdna}"/circular_path_*_concatenated.fa; do
        j=$(_disassemble_extract_number_circular_path "${_ptdna}")
        _ptdir="${_var_mtdna}/${j}"
        _polap_log3_pipe "python ${_POLAPLIB_DIR}/run-polap-py-compare2ptdna.py \
    		      --seq1 ${_arg_disassemble_c} \
	    	      --seq2 ${_ptdna} \
		          --out ${_ptdir} \
		          2>${_polap_output_dest}"
      done
      shopt -u nullglob # Restore default behavior

      # 14-2. find the one with the max coverage
      _folder_with_max_coverage=$(_disassemble_find_folder_with_max_coverage "${_var_mtdna}")
      if [[ $? -eq 0 ]]; then
        # Calculate absolute difference using bc
        _summary_coverage_ref=$(<"${_var_mtdna}/${_folder_with_max_coverage}/coverage1.txt")
        _summary_coverage_target=$(<"${_var_mtdna}/${_folder_with_max_coverage}/coverage2.txt")
        local diff=$(echo "scale=5; a=(${_summary_coverage_ref} - ${_summary_coverage_target}); if (a < 0) a=-a; a" | bc)

        # Check if difference is less than 0.05
        if (($(echo "$diff < 0.05" | bc -l))); then
          _polap_log2 "  coverage_ref and coverage_target: difference ($diff) is less than 0.05"
          _summary_j_candidate="${_folder_with_max_coverage}"
        fi
      fi

      # 14-3. select one ptDNA
      if [[ "${_summary_j_candidate}" -ne -1 ]]; then

        cp -p "${_var_mtdna}/${_summary_j_candidate}/coverage1.txt" \
          "${_var_mtdna}/measure_ptdna-blast-alignment-coverage-rate1.txt"
        cp -p "${_var_mtdna}/${_summary_j_candidate}/coverage2.txt" \
          "${_var_mtdna}/measure_ptdna-blast-alignment-coverage-rate2.txt"
        _restarted_fasta="${_var_mtdna}/${_summary_j_candidate}/seq2_restarted.fasta"

        # polish the ptDNA if necessary
        #
        _polap_log3_cmd rm -f "${_arg_unpolished_fasta}"
        _polap_log3_cmd rm -f "${_arg_final_assembly}"

        if [[ -s "${_restarted_fasta}" ]]; then
          _polap_log2 "    ptDNA selected: ${_restarted_fasta}"
          if [[ "${_summary_num_circular_paths}" -eq 4 ]] || [[ "${_summary_num_circular_paths}" -eq 2 ]]; then
            _polap_log3_cmd cp -p "${_restarted_fasta}" "${_arg_unpolished_fasta}"
          else
            die "    but number of circular patths: ${_summary_num_circular_paths} neither 2 nor 4"
          fi
        else
          _polap_log2 "    no such ptDNA: ${_restarted_fasta}"
        fi

        if [[ "${_arg_polish}" == "on" ]]; then
          if [[ -s "${_arg_unpolished_fasta}" ]]; then
            _run_polap_polish
          else
            die "ERROR: no such file: ${_arg_unpolished_fasta}"
          fi
        else
          _polap_log3_cmd rm -f "${_arg_final_assembly}"
        fi
      else
        _polap_log2 "    (compare) no ptDNA selected"
      fi
    else
      _polap_log2 "    pick the first ptDNA: ${_circular_path_fasta}"
      _polap_log2 "      as the unpolished draft: ${_arg_unpolished_fasta}"
      _polap_log3_cmd cp -p "${_circular_path_fasta}" "${_arg_unpolished_fasta}"
      _polap_log3_cmd rm -f "${_arg_final_assembly}"
    fi
  else
    _polap_log2 "  no ptDNA: circular_path count is neither 2 nor 4."
    _polap_log3_cmd rm -f "${_arg_unpolished_fasta}"
    _polap_log3_cmd rm -f "${_arg_final_assembly}"
  fi
  # Code 1
  local arr=(
    "${_summary_coverage_ref}"
    "${_summary_coverage_target}"
    "${_summary_j_candidate}"
  )
  echo "${arr[@]}"
  #
  # Code 2
  # _a=(
  # 	"${_summary_coverage_ref}"
  # 	"${_summary_coverage_target}"
  # 	"${_summary_j_candidate}"
  # )
}

function _disassemble-step15 {
  local var_mtdna="${1}"
  local _summary_pident
  local _ptdna="${_var_mtdna}/ptdna.1.fa"
  local _ptdir="${_var_mtdna}/b"

  _summary_pident=-1
  if [[ -s "${_var_mtdna}/b/pident.txt" ]] &&
    [[ "${_arg_test}" == "on" ]]; then
    _polap_log2 "    found: ${_var_mtdna}/b/pident.txt, so skipping percent identity ..."
    _summary_pident=$(<"${_var_mtdna}/b/pident.txt")
  else
    _summary_pident=0
    if [[ -s "${_ptdna}" ]]; then
      _polap_log3_pipe "python ${_POLAPLIB_DIR}/run-polap-py-compare2ptdna.py \
    		    --seq1 ${_arg_disassemble_c} \
	    	    --seq2 ${_ptdna} \
		        --out ${_ptdir} \
		        2>$_polap_output_dest"
      if [[ -s "${_var_mtdna}/b/pident.txt" ]]; then
        _summary_pident=$(<"${_var_mtdna}/b/pident.txt")
      fi
    else
      _polap_log2 "    no such polished ptDNA assembled: ${_ptdna}"
    fi
  fi
  echo "${_summary_pident}"
}

# local _a=($(_disassemble-step16 "${_var_mtdna}"))
# _nseq="${_a[0]}"
# _summary_gc="${_a[1]}"
# _summary_length="${_a[2]}"
function _disassemble-step16 {
  local var_mtdna="${1}"
  local _nseq
  local _summary_gc
  local _summary_length
  local ptdna_file

  ptdna_file="${_var_mtdna}/ptdna.1.fa"
  if [[ ! -s "${ptdna_file}" ]]; then
    ptdna_file="${_var_mtdna}/ptdna.0.fa"
  fi
  _polap_log2 "    input1: ${ptdna_file}"

  _nseq=-1
  _summary_gc=-1
  _summary_length=-1
  # cflye other steps
  if [[ -s "${ptdna_file}" ]]; then
    # NOTE: subset gfa with the seeds
    # TODO: a circular path

    _nseq=$(seqkit stats -Ta "${ptdna_file}" |
      csvtk cut -t -f num_seqs |
      csvtk del-header |
      head -n 1)
    [[ "${_nseq}" == 1 ]] || die "ERROR: ${_nseq} only one sequence in ${ptdna_file}"
    _summary_gc=$(seqkit stats -Ta "${ptdna_file}" |
      csvtk cut -t -f "GC(%)" |
      csvtk del-header |
      head -n 1)
    _summary_length=$(seqkit stats -Ta "${ptdna_file}" |
      csvtk cut -t -f sum_len |
      csvtk del-header |
      head -n 1)
  else
    _polap_log2 "  no ptDNA assembled: ${ptdna_file}"
  fi
  local arr=("${_nseq}" "${_summary_gc}" "${_summary_length}")
  echo "${arr[@]}"
}

function _disassemble-redownsample {

  local input_file date input1 output1 long_read_bp genome_size
  local long_read_coverage target_coverage sampling_rate random_seed

  input_file="${_ga_outdir}/lx.txt" # change if your input comes from a different file

  local -A redown # this array is local to this function
  _disassemble_parse_redownsample "$input_file" redown

  # Use the array locally
  _polap_log0 "Input1: ${redown[input1]}"
  _polap_log0 "Output1: ${redown[output1]}"
  _polap_log0 "random_seed: ${redown[random_seed]}"
  _polap_log0 "sampling_rate: ${redown[sampling_rate]}"

  # cmd
  seqkit sample \
    -p "${redown[sampling_rate]}" \
    -s "${redown[random_seed]}" \
    "${redown[input1]}" \
    -o "${redown[output1]}" 2>${_polap_output_dest}

  input_file="${_ga_outdir}/sx.txt" # change if your input comes from a different file
  _disassemble_parse_redownsample "$input_file" redown
  _polap_log0 "input1: ${redown[input1]}"
  _polap_log0 "input2: ${redown[input2]}"
  _polap_log0 "output1: ${redown[output1]}"
  _polap_log0 "output2: ${redown[output2]}"
  _polap_log0 "random_seed: ${redown[random_seed]}"
  _polap_log0 "sampling_rate: ${redown[sampling_rate]}"

  seqtk sample \
    -s "${redown[random_seed]}" \
    "${redown[input1]}" \
    "${redown[sampling_rate]}" \
    >"${redown[output1]}" 2>${_polap_output_dest}

  seqtk sample \
    -s "${redown[random_seed]}" \
    "${redown[input2]}" \
    "${redown[sampling_rate]}" \
    >"${redown[output2]}" 2>${_polap_output_dest}

  cat "${redown[input1]}" \
    "${redown[input2]}" \
    >"${_ga_outdir}/s.fq"
}

_disassemble_parse_redownsample() {
  local input_file="$1"
  local -n result="$2" # name-ref to caller's associative array

  while IFS= read -r line; do
    line="$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

    case "$line" in
    "coverage compute date:"*)
      result[date]="${line#*: }"
      ;;
    "input1:"*)
      result[input1]="${line#*: }"
      ;;
    "input2:"*)
      result[input2]="${line#*: }"
      ;;
    "output1:"*)
      result[output1]="${line#*: }"
      ;;
    "output2:"*)
      result[output2]="${line#*: }"
      ;;
    "short-read1:"*)
      val="${line#*: }"
      result[short_read1]="${val% (bp)}"
      ;;
    "short-read2:"*)
      val="${line#*: }"
      result[short_read2]="${val% (bp)}"
      ;;
    "short-read:"*)
      val="${line#*: }"
      result[short_read_total]="${val% (bp)}"
      ;;
    "short-read coverage:"*)
      val="${line#*: }"
      result[short_read_coverage]="${val%x}"
      ;;
    "long-read coverage:"*)
      val="${line#*: }"
      result[long_read_coverage]="${val%x}"
      ;;
    "genome size:"*)
      val="${line#*: }"
      result[genome_size]="${val% (bp)}"
      ;;
    "target coverage:"*)
      val="${line#*: }"
      result[target_coverage]="${val%x}"
      ;;
    "sampling rate:"*)
      result[sampling_rate]="${line#*: }"
      ;;
    "random seed:"*)
      result[random_seed]="${line#*: }"
      ;;
    esac
  done <"$input_file"
}

# stage 0
# input: input sequencing data files
# output:
# 0/l.fq
# 0/s_1.fq
# 0/s_2.fq
#
# output:
# o/short_expected_genome_size.txt
# o/l.fq.txt
# o/s1.fq.txt
# o/s2.fq.txt
# o/0/lx.txt
# o/0/sx.txt
# o/0/target-coverage ... txt
#
# created and deleted:
# o/s.fq
function _disassemble-stage0 {
  _polap_log1 "  subsample the input total read data using a given target coverage: ${_arg_downsample}x"

  local _polap_output_dest="/dev/null"
  local _infile="${_arg_long_reads}"
  local _infile1="${_arg_short_read1}"
  local _infile2="${_arg_short_read2}"
  local _outfile="${_arg_outdir}/${_arg_inum}/l.fq"
  local _outfile1="${_arg_outdir}/${_arg_inum}/s_1.fq"
  local _outfile2="${_arg_outdir}/${_arg_inum}/s_2.fq"
  local _ga_input_short_reads="${_arg_outdir}/${_arg_inum}/s.fq"

  local _ga_outdir="${_arg_outdir}/${_arg_inum}"
  local _total_input_short_reads="${_arg_outdir}/s.fq"

  # long-read downsampling
  # short-read downsampling
  #
  # concatenate the short-read data
  # estimate the genome size using the input short-read data
  # compute the total size of the input long-read data
  # use target coverage to determine the sampling rate
  # if the rate is less than 1, then subsample it with the rate.
  # Otherwise, link the input data with subsampled data file name.

  local _outdir_genome_size="${_arg_outdir}/short_expected_genome_size.txt"

  # concatenate the short-read data
  # if [[ -s "${_outdir_genome_size}" ]]; then
  # 	_polap_log2 "    found: genome size: ${_outdir_genome_size}"
  # else
  # 	_disassemble-step2 "${_total_input_short_reads}" "${_infile1}" "${_infile2}"
  # fi

  # estimate the genome size using the input short-read data
  local _v
  if [[ -z "${_arg_genomesize}" ]]; then
    _polap_log2 "    no genome size option"
    if [[ -s "${_outdir_genome_size}" ]]; then
      _polap_log2 "    found: genome size: ${_outdir_genome_size}"
    else
      _disassemble-step2 "${_total_input_short_reads}" "${_infile1}" "${_infile2}"
      _polap_find-genome-size \
        "${_total_input_short_reads}" \
        "${_arg_outdir}"
    fi
    _v=$(<"${_outdir_genome_size}")
  else
    _polap_log2 "    genome size option: ${_arg_genomesize}"
    _v="${_arg_genomesize}"
  fi

  # count the bases in the input long-read data file
  if [[ -s "${_arg_outdir}/l.fq.txt" ]]; then
    _polap_log2 "    found: ${_arg_outdir}/l.fq.txt"
  else
    _polap_lib_fastq-total-length-of "${_arg_long_reads}" "${_arg_outdir}/l.fq.txt"
  fi
  local _l=$(<"${_arg_outdir}/l.fq.txt")

  # compute the long-read coverage
  local _coverage_long=$(echo "scale=5; ${_l} / ${_v}" | bc)
  # compute the long-read data down-sampling rate accordingly
  local _rate
  if [[ "${_arg_downsample}" == "0" ]]; then
    _rate=1
  else
    _rate=$(echo "scale=5; ${_arg_downsample} / ${_coverage_long}" | bc)
  fi

  local _msg1="coverage compute date: $(date)"
  _msg1+="
    input1: ${_infile}
    output1: ${_outfile}
    long-read: ${_l} (bp)
    genome size: ${_v} (bp)
    long-read coverage: ${_coverage_long}x
    target coverage: ${_arg_downsample}x
    sampling rate: ${_rate}"

  _polap_log2 "${_msg1}"

  echo "${_msg1}" >>"${_arg_outdir}/${_arg_inum}/lx.txt"
  echo "${_arg_downsample}" >"${_arg_outdir}/${_arg_inum}/target-coverage-${_arg_downsample}x.txt"

  # create a downsampled long-read data
  # if there is no one.
  if [[ -s "${_outfile}" ]]; then
    _polap_log0 "  found1: ${_outfile}"
  else

    local result=$(echo "$_rate < 1" | bc)

    if [ "$result" -eq 1 ]; then
      # echo "The rate value is less than 1"
      if [[ "${_arg_dry}" == "off" ]]; then
        _polap_lib_random-get
        _seed=${_polap_var_random_number}
        echo "random seed: ${_seed}" >>"${_arg_outdir}/${_arg_inum}/lx.txt"

        seqkit sample \
          -p "${_rate}" \
          -s "${_seed}" \
          "${_infile}" \
          -o ${_outfile} 2>${_polap_output_dest}
      fi
    else
      # echo "The value is not less than 1"
      _polap_log2 "    sampling rate is not less than 1: ${_rate}"
      _polap_log2 "    no subsampling of input: ${_infile}"
      _polap_log1 "    no subsampling: ${_outfile}"
      ln -s "$PWD/${_infile}" "${_outfile}"
    fi
  fi

  # compute the total size of the input short-read data and the sampling rate
  if [[ -s "${_arg_outdir}/s1.fq.txt" ]]; then
    _polap_log2 "    found: ${_arg_outdir}/s1.fq.txt"
  else
    _polap_lib_fastq-total-length-of "${_arg_short_read1}" "${_arg_outdir}/s1.fq.txt"
  fi
  if [[ -s "${_arg_outdir}/s2.fq.txt" ]]; then
    _polap_log2 "    found: ${_arg_outdir}/s2.fq.txt"
  else
    _polap_lib_fastq-total-length-of "${_arg_short_read2}" "${_arg_outdir}/s2.fq.txt"
  fi
  local _s1=$(<"${_arg_outdir}/s1.fq.txt")
  local _s2=$(<"${_arg_outdir}/s2.fq.txt")
  local _s=$((_s1 + _s2))
  if [[ "${_s1}" == "${_s2}" ]]; then
    _polap_log1 "Two of the pair are the same in the number of reads."
  else
    _polap_log0 "  short-read1: ${_s1} (bp)"
    _polap_log0 "  short-read2: ${_s2} (bp)"
    _polap_log0 "WARNING: two of the pair are different in the number of reads."
  fi

  local _coverage_short=$(echo "scale=5; ${_s} / ${_v}" | bc)
  local _rate
  if [[ "${_arg_downsample}" == "0" ]]; then
    _rate=1
  else
    _rate=$(echo "scale=5; ${_arg_downsample} / ${_coverage_short}" | bc)
  fi

  local _msg1="coverage compute date: $(date)"
  _msg1+="
    input1: ${_infile1}
    input2: ${_infile2}
    output1: ${_outfile1}
    output2: ${_outfile2}
    short-read1: ${_s1} (bp)
    short-read2: ${_s2} (bp)
    short-read: ${_s} (bp)
    short-read coverage: ${_coverage_short}x
    genome size: ${_v} (bp)
    target coverage: ${_arg_downsample}x
    sampling rate: ${_rate}"

  _polap_log2 "${_msg1}"

  echo "${_msg1}" >>"${_arg_outdir}/${_arg_inum}/sx.txt"

  # create a downsampled short-read data file
  # if there is no one.
  if [[ -s "${_ga_input_short_reads}" ]]; then
    _polap_log0 "  found2: ${_ga_input_short_reads}"
  else

    local result=$(echo "$_rate < 1" | bc)

    if [ "$result" -eq 1 ]; then
      # echo "The rate value is less than 1"

      # Example:
      # seqtk sample -s100 read1.fq 0.1 >sub1.fq
      # seqtk sample -s100 read2.fq 0.1 >sub2.fq

      if [[ "${_arg_dry}" == "off" ]]; then
        _polap_lib_random-get
        _seed=${_polap_var_random_number}
        echo "random seed: ${_seed}" >>"${_arg_outdir}/${_arg_inum}/sx.txt"

        seqtk sample \
          -s"${_seed}" \
          "${_infile1}" \
          "${_rate}" \
          >"${_outfile1}"

        seqtk sample \
          -s"${_seed}" \
          "${_infile2}" \
          "${_rate}" \
          >"${_outfile2}"

      fi
    else
      # echo "The value is not less than 1"
      _polap_log2 "  sampling rate is not less than 1: ${_rate}"
      _polap_log2 "  no subsampling of input: ${_infile1}"
      _polap_log2 "  no subsampling of input: ${_infile2}"
      _polap_log1 "  no subsampling: ${_outfile1}"
      _polap_log1 "  no subsampling: ${_outfile2}"
      ln -s "$PWD/${_infile1}" "${_outfile1}"
      ln -s "$PWD/${_infile2}" "${_outfile2}"
    fi
    _polap_log1 "   concatenate the two downsampled short-read data files"
    _disassemble-step2 "${_ga_input_short_reads}" "${_outfile1}" "${_outfile2}"
    check_file_existence "${_ga_input_short_reads}"
    # _polap_log2 "     clean-up the two downsampled short-read data files"
    # _polap_log3_cmd rm -f "${_outfile1}" "${_outfile2}"
  fi

  # delete the s.fq
  # We need to save the disk space because s.fq at the outdir can be huge.
  # I do not find any other use of this s.fq except for the countig the bases.
  # _total_input_short_reads="${_arg_outdir}/s.fq"
  if [[ -s "${_total_input_short_reads}" ]]; then
    _polap_log3_cmd rm -f "${_total_input_short_reads}"
  fi
}

# stage 1
function _disassemble-stage1 {
  _disassemble_i_stage="${_disassemble_i}/1"

  _polap_log1 "    determine the case ..."
  if [[ "${_arg_disassemble_c_is}" == "off" ]]; then
    _polap_log2 "    case: infer"
    _arg_contigger="on"
    _arg_polish="off"
  else
    if [[ "${_arg_disassemble_align_reference}" == "off" ]]; then
      _polap_log2 "    case: compare"
      _arg_contigger="off"
      _arg_polish="on"
    else
      _polap_log2 "    case: check"
      _arg_contigger="on"
      _arg_polish="off"
    fi
  fi
  _polap_log2 "    contigger: ${_arg_contigger}"
  _polap_log2 "    polish: ${_arg_polish}"

  if [[ -z "${_arg_disassemble_s}" ]] && [[ -z "${_arg_disassemble_beta}" ]]; then
    _polap_log1 "    no option: --disassemble-s or --disassemble-beta"
    _polap_log2 "      -> use a range of rates for subsampling long-read sequencing data"

    # short-read data
    if [[ -s "${_arg_short_read1}" ]] ||
      [[ -s "${_arg_short_read2}" ]]; then
      _polap_log1 "    use short-read data for genome size estimates"
      if [[ -z "${_arg_steps_include}" ]]; then
        _arg_steps_include="1-16"
      fi
    else
      _polap_log0 "    no short-read data specified, so use a range of genome sizes"
      die "ERROR: not implemented yet!"
    fi
    _polap_log3_cmd mkdir -p "${_disassemble_i}/1"
    _disassemble_make_params_stage1 "${_disassemble_i}/1/params.txt"
    _run_polap_step-disassemble 1
    _rstatus="$?"
    [[ "$_rstatus" -ne 0 ]] && return "$_rstatus"
  else
    # --disassemble-s -> only one sample or single s
    _polap_log1 "    --disassemble-s -> no iteration or no stage 1"
    _arg_steps_include="1,3"
    _run_polap_step-disassemble 1
    local _ga_outdir="${_arg_outdir}/${_arg_inum}"
    local _ga_long_total_length=$(<"${_ga_outdir}/long_total_length.txt")
    if [[ -z "${_arg_disassemble_s}" ]]; then
      if [[ -z "${_arg_disassemble_beta}" ]]; then
        die "BUG: not possible case"
      else
        _arg_disassemble_s=$(echo "scale=0; ${_ga_long_total_length} * ${_arg_disassemble_beta} / 1" | bc)
      fi
    else
      if [[ -z "${_arg_disassemble_beta}" ]]; then
        _arg_disassemble_beta=$(echo "scale=3; ${_arg_disassemble_s} / ${_ga_long_total_length}" | bc)
      fi
    fi
    _polap_log1 "  subsample size: ${_arg_disassemble_s}"
    _polap_log1 "  subsample rate: ${_arg_disassemble_beta}"
  fi
}

# stage 2
function _disassemble-stage2 {
  # select the size and alpha using stage 1
  _disassemble_i_stage="${_disassemble_i}/1"
  _selected_index_stage1="${_disassemble_i}/1/selected-index.txt"
  local _j_best_stage1=-1

  if [[ -z "${_arg_disassemble_s}" ]]; then
    _summary1_ordered="${_disassemble_i_stage}/summary1-ordered.txt"
    _summary1_txt="${_disassemble_i_stage}/summary1.txt"
    if [[ -s "${_summary1_txt}" ]]; then

      if [[ "${_arg_disassemble_align_reference}" == "off" ]] &&
        [[ "${_arg_disassemble_c_is}" == "on" ]]; then
        # case: compare
        _disassemble_report1
      else
        # case: check or infer
        _disassemble_report1 true
      fi
      #
      # _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
      #      --table ${_disassemble_i_stage}/summary1.txt \
      #      --out ${_disassemble_i_stage}/summary1-ordered.txt \
      #      --plot ${_disassemble_i_stage}/summary1-ordered.pdf \
      #      2>${_polap_output_dest}"
      # _rstatus="$?"
      # if [[ "$_rstatus" -ne 0 ]]; then
      # 	_polap_log0 "ERROR: run-polap-r-disassemble.R on ${_disassemble_i_stage}/summary1.txt"
      # fi
    else
      die "ERROR: no such file stage1 summary: ${_summary1_txt}"
    fi

    if [[ -s "${_summary1_ordered}" ]]; then

      # Use awk to extract the desired columns and save the output to a Bash variable
      output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $5, $9, $11}' "${_summary1_ordered}")
      # n=$(awk 'END {print NR - 1}' "${_summary1_ordered}")
      n=$(awk 'NR > 1 && $0 !~ /^#/ {count++} END {print count}' "${_summary1_ordered}")
      if ((n < 2)); then
        _polap_log0 "  the effective sample size: ${n}: ${_summary1_ordered}"
        _polap_log0 "ERROR: the number of potential ptDNA assemblies is less than 2"
        _polap_log0 "  suggestion: increase the maximum memory requirement --disassemble-memory"
        return "${_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES}"
      elif ((n < 5)); then
        _polap_log0 "  the effective sample size: ${n}: ${_summary1_ordered}"
        _polap_log0 "  warning: the number of potential ptDNA assemblies is too small to select one: $n < 5"
        _polap_log0 "  suggestion: increase the step size --disassemble-n"
        _polap_log0 "  suggestion: increase the maximum memory requirement --disassemble-memory"
        #     _polap_log0 "  suggestion: decrease the max subsample size rate --disassemble-p"
        # _polap_log0 "  suggestion: decrease the subsample size --disassemble-b"
        # _polap_log0 "  suggestion: decrease the maximum draft genome size for update alpha --disassemble-m"
      fi
      # Use read to split the output into individual variables
      read -r index size rate randomseed genomesize alpha <<<"$output"
      # Print the extracted variables
      _polap_log2 "    best Index: $index"
      _polap_log2 "    rate: $rate"
      _polap_log2 "    size: $size"
      _polap_log2 "    alpha: $alpha"
      _polap_log2 "    genome size: $genomesize"
      _j_best_stage1="${index}"
      _arg_disassemble_s="$size"
      _arg_disassemble_alpha="$alpha"
      _arg_disassemble_n_is="on"
    else
      _polap_log0 "ERROR: no such file: ${_summary1_ordered}"
      _polap_log0 "ERROR: subsampling and assembly in stage 1 did not produce enough assemblies."
      _polap_log0 "  suggestion: increase the max subsample size rate --disassemble-p"
      _polap_log0 "  suggestion: increase the subsample size --disassemble-b"
      _polap_log0 "  suggestion: increase the downsample size --downsample"
      die "SUGGESTION: increase the sampling size --disassemble-n"
    fi
  else
    _polap_log1 "  --disassemble-s -> skip stage 1"
    _polap_log1 "  --disassemble-s -> no iteration or no stage 1"
  fi
  _polap_log2 "  --disassemble-alpha must have been set as well: ${_arg_disassemble_alpha}"

  echo "${_j_best_stage1}" >"${_selected_index_stage1}"

  # runs the replicates or previously stage 2
  _arg_steps_include=""
  _disassemble_i_stage="${_disassemble_i}/2"

  _polap_log1 "    flye runs through the end of the polishing in the stage 2 of replicates"
  _arg_contigger="off"
  if [[ "${_arg_disassemble_c_is}" == "off" ]]; then
    _polap_log2 "    case: infer"
    _arg_polish="off"
  else
    if [[ "${_arg_disassemble_align_reference}" == "off" ]]; then
      _polap_log2 "    case: compare"
      _arg_polish="on"
    else
      _polap_log2 "    case: check"
      _arg_polish="off"
    fi
  fi
  _polap_log1 "    contigger: ${_arg_contigger}"
  _polap_log1 "    polish: ${_arg_polish}"

  # stage 2:
  # _polap_log0 "  stage 2: assemble with a single set of parameters"
  if [[ -z "${_arg_disassemble_s}" ]]; then
    _polap_log0 "ERROR: subsample size must be given either by the stage 1 or option --disassemble-s"
    die "ERROR: no --disassemble-s"
  fi
  if [[ -z "${_arg_disassemble_alpha}" ]]; then
    _polap_log0 "ERROR: subsampling minimum coverage alpha must be given either by the stage 1 or option --disassemble-s"
    die "ERROR: no --disassemble-alpha"
  fi
  _arg_disassemble_a="${_arg_disassemble_s}"
  _arg_disassemble_b="${_arg_disassemble_s}"
  _arg_disassemble_b_is="on"
  _arg_disassemble_n="${_arg_disassemble_r}"
  _arg_disassemble_n_is="on"
  _polap_log1 "  input1: --disassemble-s=${_arg_disassemble_s}"

  if [[ -s "${_arg_short_read1}" ]] ||
    [[ -s "${_arg_short_read2}" ]]; then
    _polap_log2 "  genome size estimated by the short-read data"
    if [[ -z "${_arg_steps_include}" ]]; then
      _arg_steps_include="1-16"
      _arg_steps_exclude="9" # no changes to alpha
      check_file_existence "${_ga_input_short_reads}"
    fi
  else
    _polap_log1 "    no short-read data, so use a range of genome sizes"
    die "ERROR: not implemented yet!"
  fi

  _polap_log3_cmd mkdir -p "${_disassemble_i}/2"
  _disassemble_make_params_stage2 "${_disassemble_i}/2/params.txt"
  _run_polap_step-disassemble 2
  _rstatus="$?"
  [[ "$_rstatus" -ne 0 ]] && return "$_rstatus"
}

# stage 3
function _disassemble-stage3 {
  # stage 3

  _disassemble_i_stage="${_disassemble_i}/2"
  _summary1_ordered="${_disassemble_i_stage}/summary1-ordered.txt"
  _summary1_txt="${_disassemble_i_stage}/summary1.txt"
  _final_assembly="${_disassemble_i}/pt.1.fa"
  _unpolished_final_assembly="${_disassemble_i}/pt.0.fa"
  _index_file="${_disassemble_i}/stage2-selected-index.txt"
  _selected_index_stage2="${_disassemble_i}/2/selected-index.txt"
  _selected_index_stage3="${_disassemble_i}/3/selected-index.txt"

  if [[ -s "${_summary1_txt}" ]]; then
    if [[ "${_arg_disassemble_align_reference}" == "off" ]] &&
      [[ "${_arg_disassemble_c_is}" == "on" ]]; then
      # case: compare
      _disassemble_report2
    else
      # case: check or infer
      _disassemble_report2 true
    fi
    #
    # _polap_log1 "  draw the length distribution of the ptDNA sequences in stage 2"
    # _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
    #       --table ${_disassemble_i_stage}/summary1.txt \
    #       --out ${_disassemble_i_stage}/summary1-ordered.txt \
    #       --plot ${_disassemble_i_stage}/summary1-ordered.pdf \
    #       2>${_polap_output_dest}"
  else
    die "ERROR: no such file stage2 summary: ${_summary1_txt}"
  fi

  local _j_best_stage2=-1
  local _j_best_stage3=-1
  _var_mtdna=""
  _arg_unpolished_fasta=""
  if [[ "${_arg_disassemble_align_reference}" == "off" ]]; then
    if [[ -s "${_summary1_ordered}" ]]; then
      _polap_log1 "  determine the best j in stage 2 using the length distribution"

      # Use awk to extract the desired columns and save the output to a Bash variable
      output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $5, $9, $11}' "${_summary1_ordered}")
      # n=$( awk 'END {print NR - 1}' "${_summary1_ordered}")
      n=$(awk 'NR > 1 && $0 !~ /^#/ {count++} END {print count}' "${_summary1_ordered}")
      if ((n < 2)); then
        _polap_log0 "  the effective sample size: ${n}: ${_summary1_ordered}"
        _polap_log0 "ERROR: the number of potential ptDNA assemblies is less than 2"
        _polap_log0 "  suggestion: delete ${_disassemble_dir}/infer-1/2"
        _polap_log0 "  suggestion: increase the maximum memory requirement --disassemble-memory"
        return "${_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES}"
      elif ((n < 5)); then
        _polap_log0 "  the effective sample size: ${n}: ${_summary1_ordered}"
        _polap_log0 "  warning: the number of potential ptDNA assemblies is too small to select one"
        _polap_log0 "  suggestion: increase the replicate size --disassemble_r"
        _polap_log0 "  suggestion: increase the maximum memory requirement --disassemble-memory"
      fi
      # Use read to split the output into individual variables
      read -r index size rate randomseed genomesize alpha <<<"$output"
      # Print the extracted variables
      _polap_log1 "  best index: $index"
      _j_best_stage2="${index}"
    else
      die "ERROR: no such file: ${_summary1_ordered}"
    fi
    echo "${_j_best_stage2}" >"${_selected_index_stage2}"

    # no check with align-reference
    # j: should be determined
    # the first sequence is used as ptdna.0.fa
    _var_mtdna="${_disassemble_i_stage}/${_j_best_stage2}/52-mtdna"
    _circular_path_fasta="${_var_mtdna}/circular_path_1_concatenated.fa"
    _arg_unpolished_fasta="${_var_mtdna}/ptdna.0.fa"
    _arg_final_assembly="${_var_mtdna}/ptdna.1.fa"
    cp -p "${_circular_path_fasta}" "${_arg_unpolished_fasta}"
    _polap_log1 "  the selected unpolished ptDNA: ${_circular_path_fasta}"
  else
    # Stage: check
    #
    # Copy of "check" menu
    #
    # use the reference to select one of the four ptDNA candidates

    # do not place the pt.1.fa and pt.0.fa as the final yet
    # we just determine the candidate in 52-mtdna folder.

    _polap_log1 "  we use the ordered table to select one with 4 candidates"
    _polap_log1 "  stage check: use the reference to select one of the four ptDNA candidates"
    rm -f "${_index_file}"
    if [[ -s "${_summary1_ordered}" ]]; then
      j=$(awk -F'\t' '$15 == 4 {print $1; exit}' "${_summary1_ordered}")
      _j_best_stage2="${j}"
      if [[ -z "${j}" ]]; then
        _polap_log1 "WARNING: no such index with 4 paths in file: ${_summary1_ordered}"
        _polap_log1 "  suggestion: increase the replicate size --disassemble-r"

        # fall back to 2-fragment cases
        j=$(awk -F'\t' '$15 == 2 {print $1; exit}' "${_summary1_ordered}")
        _j_best_stage2="${j}"
        if [[ -z "${j}" ]]; then
          _polap_log0 "ERROR: no such index with 2 paths in file: ${_summary1_ordered}"
          _polap_log0 "  suggestion: no way of comparing it with ptGAUL ptDNA"
          return 1
        fi

      else
        _polap_log1 "    index with 4 paths: ${j}"
      fi

      # TODO: check this part.
      # check if the chosen index is the best one
      # the top best: index of summary1-ordered.txt
      # selected-index.txt: the chosen index
      local _output_line=$(awk -F'\t' 'NR==2 {print $1}' "${_summary1_ordered}")
      read -r _top_index <<<"$_output_line"
      if [[ "${_top_index}" != "${_j_best_stage2}" ]]; then
        _polap_log2 "    the top best: index of summary1-ordered.txt: ${_top_index}"
        _polap_log2 "    selected-index.txt: the chosen index: ${_j_best_stage2}"
        _polap_log0 "Warning: the selected index of assembly length nearest the mode of assembly length distribution is different from the chosen index."
        _polap_log0 "Warning: the percent identity of the assembly might not be near 100% because the selected plastid genome structure is different from the reference although they are actually very similar."
      fi
    else
      _polap_log0 "  you have no ptDNA assembly from the stage 2"
      _polap_log0 "  suggestion: increase the replicate size --disassemble-r"
      _polap_log0 "ERROR: no such file: ${_summary1_ordered}"
      return 1
    fi
    echo "${_j_best_stage2}" >"${_selected_index_stage2}"

    _var_mtdna="${_disassemble_i}/2/${_j_best_stage2}/52-mtdna"
    _polap_log1 "    best choice of mtdna dir with 4 paths: ${_var_mtdna}"

    # make sure that there are 4 paths in the selected ptDNA assembly
    if [[ -s "${_var_mtdna}/circular_path_count.txt" ]]; then
      _summary_num_circular_paths=$(<"${_var_mtdna}/circular_path_count.txt")
      _summary_num_circular_nodes=$(<"${_var_mtdna}/circular_path_nodes.txt")
      _polap_log3_cat "${_var_mtdna}/circular_path.txt"
      if [[ "${_summary_num_circular_paths}" -eq 4 ]]; then
        _polap_log2 "    check: circular_path count is 4."
      elif [[ "${_summary_num_circular_paths}" -eq 2 ]]; then
        _polap_log2 "    check: circular_path count is 2."
      else
        _polap_log0 "INFO: Increase the sample size in stage 2"
        die "ERROR: circular_path count is neither 2 nor 4: see ${_var_mtdna}"
      fi
    fi

    _circular_path_fasta="${_var_mtdna}/circular_path_1_concatenated.fa"
    _arg_unpolished_fasta="${_var_mtdna}/ptdna.0.fa"
    _arg_final_assembly="${_var_mtdna}/ptdna.1.fa"
    # Copy of:
    # 14. Pick one sequence based on the reference, if any.
    # Use the first sequence if there is no such reference.
    #
    # for each ptDNA of the potential ptDNA sequences,
    #   select one candidate ptDNA
    #   rearrange it so that we could do a pairwise sequence alignment
    _summary_coverage_ref=-1
    _summary_coverage_target=-1
    _summary_j_candidate=-1

    # if _polap_contains_step 14 "${_step_array[@]}"; then
    _polap_log1 "  checking stage or step 14 of step-disassemble"
    _polap_log1 "    : choose one of the multiple candidate ptDNA sequences"
    _polap_log2 "    input1: ${_circular_path_fasta}"

    # FIXME: disassemble example cannot call itself because of this
    # running in a subshell. Output from _disassemble-step14 gets out to
    # the three return values.
    # Use other way to get the 3 return values not using echo at the end of
    # the function.
    #
    # Code 1:
    # local _a=($(_disassemble-step14 "${_var_mtdna}"))
    #
    # Code 2:
    # local _a # return of function: _disassemble-step14
    # _disassemble-step14 "${_var_mtdna}"
    local _a=($(_disassemble-step14 "${_var_mtdna}" | tail -n 1))
    _summary_coverage_ref="${_a[0]}"
    _summary_coverage_target="${_a[1]}"
    _summary_j_candidate="${_a[2]}"

    _polap_log2 "    output1: ${_summary_coverage_ref}"
    _polap_log2 "    output2: ${_summary_coverage_target}"
    _polap_log2 "    output3: one of the 4 candidates: ${_summary_j_candidate}"

    if [[ -s "${_arg_unpolished_fasta}" ]]; then
      _polap_log2 "    output4: ${_arg_unpolished_fasta}"
    fi
    if [[ -s "${_arg_final_assembly}" ]]; then
      _polap_log2 "    output5: ${_arg_final_assembly}"
    fi

    _polap_log1 "  the selected unpolished ptDNA: ${_var_mtdna}/${_summary_j_candidate}/seq2_restarted.fasta"
  fi

  if [[ -s "${_arg_unpolished_fasta}" ]]; then
    _polap_log1 "  the unpolished ptDNA: ${_arg_unpolished_fasta}"
  else
    die "ERROR: no such file: ${_arg_unpolished_fasta}"
  fi

  # revert to the user provided options
  _arg_disassemble_a="${_b_disassemble_a}"
  _arg_disassemble_b="${_b_disassemble_b}"
  _arg_disassemble_s="${_b_disassemble_s}"
  _arg_disassemble_p="${_b_disassemble_p}"
  _arg_disassemble_n="${_b_disassemble_n}"
  _arg_disassemble_b_is="${_b_disassemble_b_is}"
  _arg_disassemble_n_is="${_b_disassemble_n_is}"
  # if _polap_contains_step 6 "${_stage_array[@]}"; then
  _disassemble_i_stage="${_disassemble_i}/3"

  # choose the source ptDNA from stage 2
  # best index:
  # j_candidate:
  if [[ "${_j_best_stage2}" != "-1" ]]; then
    local _target_ptdna="${_disassemble_i}/2/${_j_best_stage2}/52-mtdna/ptdna.0.fa"
  else
    die "ERROR: no such best j from the stage 2"
  fi

  _polap_log1 "  the selected unpolished ptDNA: ${_target_ptdna}"
  _arg_unpolished_fasta="${_target_ptdna}"

  if [[ "${_arg_disassemble_simple_polishing}" == "off" ]]; then

    _polap_log1 "  simple-polishing: ${_arg_disassemble_simple_polishing}"

    if [[ -s "${_arg_unpolished_fasta}" ]]; then
      _run_polap_polish-disassemble
    else
      _polap_log0 "WARNING: no such file: ${_disassemble_i}/pt.0.fa"
    fi

    # step 6
    _disassemble_i_stage="${_disassemble_i}/3"

    _disassemble_report3

    _summary1_ordered="${_disassemble_dir}/${_arg_disassemble_i}/3/summary1-ordered.txt"

    if [[ -s "${_summary1_ordered}" ]]; then
      # Use awk to extract the desired columns and save the output to a Bash variable
      output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $5, $6, $7, $8, $9, $10}' "${_summary1_ordered}")
      n=$(awk 'END {print NR - 1}' "${_summary1_ordered}")
      if ((n < 5)); then
        _polap_log0 "  warning: the number of potential ptDNA assemblies is too small to select one"
        _polap_log0 "  suggestion: increase the replicate size --disassemble-r"
      elif ((n < 2)); then
        _polap_log0 "ERROR: the number of potential ptDNA assemblies is less than 2"
        return "${_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES}"
      fi
      # Use read to split the output into individual variables
      read -r index size rate randomseed memory1 time1 memory2 time2 length <<<"$output"
      # Print the extracted variables
      _polap_log1 "  best index: $index"
      _j_best_stage3="${index}"
      _var_mtdna="${_disassemble_i_stage}/${_j_best_stage3}"
      if [[ -s "${_var_mtdna}/ptdna.1.fa" ]]; then
        if [[ "${_arg_disassemble_align_reference}" == "off" ]]; then
          cp -p "${_var_mtdna}/ptdna.1.fa" "${_disassemble_i}/pt.subsample-polishing.1.fa"
        else
          cp -p "${_var_mtdna}/ptdna.1.fa" "${_disassemble_i}/pt.subsample-polishing.reference.aligned.1.fa"
        fi
      else
        die "ERROR: polished ptDNA: no such file: ${_var_mtdna}/ptdna.1.fa"
      fi
    else
      die "ERROR: no such file: ${_summary1_ordered}"
    fi
    echo "${_j_best_stage3}" >"${_selected_index_stage3}"

  else
    # Copy of "polishing" menu
    _arg_final_assembly="${_disassemble_dir}/${_arg_disassemble_i}/pt.1.fa"
    if [[ "${_arg_disassemble_align_reference}" == "off" ]]; then
      _arg_final_assembly="${_disassemble_dir}/${_arg_disassemble_i}/pt.simple-polishing.1.fa"
    else
      _arg_final_assembly="${_disassemble_dir}/${_arg_disassemble_i}/pt.simple-polishing.reference.aligned.1.fa"
    fi
    _polap_log1 "polishing ptDNA: ${_arg_unpolished_fasta}"
    _run_polap_polish
    _polap_log1 "full short-read data polished assembly: ${_arg_final_assembly}"
  fi
}

_disassemble_make_params_txt() {
  local _disassemble_params_file="${1}"

  rm -f "${_disassemble_params_file}"
  printf "%s: %s\n" "I" "${_arg_disassemble_i}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "P" "${_arg_disassemble_p}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "N" "${_arg_disassemble_n}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "R" "${_arg_disassemble_r}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "A" "${_arg_disassemble_a}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "B" "${_arg_disassemble_b}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "M" "${_arg_disassemble_m}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "D" "${_arg_disassemble_delta}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "Alpha" "${_arg_disassemble_alpha}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "Memory" "${_arg_disassemble_memory}" >>"${_disassemble_params_file}"
}

_disassemble_make_params_stage1() {
  local _disassemble_params_file="${1}"

  rm -f "${_disassemble_params_file}"
  printf "%s: %s\n" "--contigger" "${_arg_contigger}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "--polish" "${_arg_polish}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "I" "${_arg_disassemble_i}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "P" "${_arg_disassemble_p}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "N" "${_arg_disassemble_n}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "R" "${_arg_disassemble_r}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "A" "${_arg_disassemble_a}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "B" "${_arg_disassemble_b}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "M" "${_arg_disassemble_m}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "D" "${_arg_disassemble_delta}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "Alpha" "${_arg_disassemble_alpha}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "Memory" "${_arg_disassemble_memory}" >>"${_disassemble_params_file}"
}

_disassemble_make_params_stage2() {
  local _disassemble_params_file="${1}"

  rm -f "${_disassemble_params_file}"
  printf "%s: %s\n" "S" "${_arg_disassemble_s}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "Alpha" "${_arg_disassemble_alpha}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "--contigger" "${_arg_contigger}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "--polish" "${_arg_polish}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "I" "${_arg_disassemble_i}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "P" "${_arg_disassemble_p}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "N" "${_arg_disassemble_n}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "R" "${_arg_disassemble_r}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "A" "${_arg_disassemble_a}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "B" "${_arg_disassemble_b}" >>"${_disassemble_params_file}"
}

_disassemble_make_params_stage3() {
  local _disassemble_params_file="${1}"

  rm -f "${_disassemble_params_file}"
  printf "%s: %s\n" "I" "${_arg_disassemble_i}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "P" "${_arg_disassemble_p}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "N" "${_arg_disassemble_n}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "R" "${_arg_disassemble_r}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "A" "${_arg_disassemble_a}" >>"${_disassemble_params_file}"
  printf "%s: %s\n" "B" "${_arg_disassemble_b}" >>"${_disassemble_params_file}"
}

################################################################################
# It is not as straightforward as it might seem. The simplest approach
# is to test plastid genome assembly by subsampling long-read data. A
# more complex approach involves using a slightly modified version of
# Flye, called cflye, which provides coverage options for generating
# disjointigs. Short-read data is still used to get a rough estimate of
# the overall genome size, but once subsampling is applied, the size of
# the genome assembled from long reads can change. The sizes of both
# the long-read and short-read data are defined by their total number
# of bases. When subsampling, we reduce the amount of short-read data
# proportionally to match the size of the subsampled long-read data, then
# use JellyFish (a k-mer counter) to predict the total size of the genome
# assembled from these long reads.
#
# If the total long-read data is smaller than the short-read data,
# we simply subsample the short reads to match the long-read dataset
# size. However, if the long-read data is larger than the short-read
# data, we use all the short reads. Because our goal is to estimate the
# approximate genome size that can be assembled from subsampled long reads,
# some level of error is to be expected.
#
# In the Flye genome assembly program, disjointigs are generated, and
# those whose read coverage falls below a certain cutoff are filtered
# out. We set this initial cutoff using the size of the long-read data (A
# = 30 Mb, where Mb stands for megabases) and alpha (α = 0.1). We then
# compare the total length (S) of the disjointigs against a fixed value
# M. If S is greater than M, alpha is increased by an amount delta (δ =
# 0.1); if S is less than M, alpha is decreased by the same delta. Along
# with this process, after each cycle we also increase the long-read data
# size A by d (for example, 10 Mb) and repeat these cycles until A reaches
# B, for a total of n cycles. The values of A, B, and n are predetermined,
# and d is derived from these values. For instance, if A = 30 Mb, B =
# 120 Mb, and n = 10, then d = 10 Mb.
#
# At the end of each cycle, we extract and evaluate potential plastid
# genome candidates from the newly assembled genome. To qualify as a plastid
# candidate, the assembly must consist of a single circular contig or three
# contigs that form a circular structure. We then assess the size of this
# circular contig, the number of expected genes, and so on. These values
# can be plotted as Y-coordinates against the cycle number on the X-axis
# in a scatter plot, allowing us to visualize progress across cycles.
#
#### Step 1: Estimating Long-Read Dataset Size
#
# The process starts by determining the total length of the long-read
# dataset by counting the total number of bases present in the sequencing
# reads. This information is crucial as it defines the upper limit for
# subsampling and subsequent analyses.
#
#### Step 2 & 3: Preparing Short-Read Data or presteps 1 and 2
#
# When short-read data is available, the sequencing reads from different
# files are first combined into a single dataset to ensure uniformity. The
# total number of bases in this combined short-read dataset is then
# calculated. This base count is used to balance the amount of short-read
# data with the subsampled long-read data in later steps.
#
#### Step 4: Defining Subsampling Range
#
# A range of data sizes is defined to ensure a comprehensive assembly
# across varying input sizes. The minimum dataset size is set by comparing
# the total base counts of the long-read and short-read datasets, selecting
# the smaller of the two. If a maximum subsample size is not provided, it
# defaults to the size of the smaller dataset. This range is then divided
# into evenly spaced intervals, creating a systematic plan for subsampling
# and subsequent assembly.
#
#### Step 5: Subsampling Long-Read Data
#
# To explore the effect of varying coverage levels on assembly quality,
# smaller subsets of the long-read dataset are created. The proportion
# of reads to be included in each subset is determined by dividing the
# target size by the total dataset size. Random selection of reads ensures
# unbiased sampling, and the process is repeated across the predefined
# range of data sizes.
#
#### Step 6: Subsampling Short-Read Data
#
# Short-read data is similarly subsampled to match the size of the
# corresponding long-read subset. This ensures balanced input for assembly,
# which is important for accurate genome size estimation and assembly
# quality.
#
#### Step 7: Genome Size Estimation
#
# If a reference genome size is not available, the genome size is
# estimated using the subsampled short-read data. A k-mer counting tool
# is employed to generate a histogram of k-mer frequencies, which is then
# analyzed to derive an estimated genome size. This value serves as a
# reference for determining coverage thresholds in subsequent assembly
# steps.
#
#### Step 8: Assembly Using Modified Genome Assembler
#
# The assembly process uses a specialized version of a genome
# assembler tailored for handling plastid genomes. Parameters such as
# the estimated genome size and minimum coverage thresholds are supplied
# to guide the assembly. The assembler produces intermediate outputs
# known as disjointigs, representing contiguous sequences before final
# polishing. Memory usage is monitored to prevent excessive resource
# consumption.
#
#### Step 9: Adjusting Coverage Threshold
#
# Post-assembly, the total length of the assembled genome is compared
# to a predefined threshold. If the length exceeds this threshold, the
# coverage requirement is increased to filter out low-coverage regions in
# the next iteration. Conversely, if the length is below the threshold,
# the coverage requirement is reduced. This iterative adjustment helps
# refine the assembly.
#
#### Step 10: Generating Assembly Graph Statistics
#
# Key statistics of the assembly graph, such as the number of contigs and
# total segment length, are calculated. These metrics provide an overview
# of the structural complexity of the assembly and help in evaluating
# its completeness.
#
#### Step 11: Annotating Contigs
#
# The assembled contigs are examined for the presence of plastid-specific
# genes. Contigs with a high density of plastid genes are prioritized for
# further analysis, while those containing predominantly non-plastid genes
# are excluded.
#
#### Step 12: Selecting Seed Contigs
#
# Seed contigs are selected based on specific criteria to identify
# candidates most likely to represent the plastid genome. This process
# involves several substeps:
#
# - **Step 12-1: Preselection of Contigs**
#   Contigs are initially filtered based on gene density. Contigs with
# a sufficient number of plastid genes and minimal contamination from
# non-plastid genes are retained.
#
# - **Step 12-2: Determining Depth Range**
#   The sequencing depth of selected contigs is analyzed, and an
# appropriate depth range is defined. Contigs outside this range are
# excluded to ensure that only those with typical plastid genome coverage
# are retained.
#
# - **Step 12-3: Filtering the Assembly Graph**
#   The assembly graph is filtered by removing contigs and connections
# that fall outside the selected depth range. This step simplifies the
# graph and focuses on high-confidence regions.
#
# - **Step 12-4: Identifying Connected Components**
#   The filtered graph is analyzed to find clusters of interconnected
# contigs. These clusters represent potential plastid genome sequences. The
# contigs forming these clusters are selected as seed contigs for further
# analysis.
#
#### Step 13: Extracting Circular Genome Sequences
#
# Circular genome sequences are identified by tracing paths through the
# assembly graph that form closed loops. Only complete circular structures
# are retained, as plastid genomes are typically circular. The extracted
# sequences are saved for further evaluation.
#
#### Step 14: Selecting the Best Genome Candidate
#
# When a reference plastid genome is available, the extracted candidates
# are compared against it. The candidate with the highest sequence
# similarity and coverage is selected as the final assembly. In cases where
# no reference is provided, the first complete circular sequence is chosen.
#
#### Step 15: Finalizing the Assembly
#
# The final plastid genome assembly is validated by calculating key
# metrics such as total length and GC content. These metrics are logged,
# and the assembled sequence is saved as the final result. This completes
# the iterative workflow, yielding a high-quality plastid genome assembly.
#
################################################################################
function _run_polap_disassemble {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  local exit_code
  local i
  local j
  local k
  local c
  local n
  local s
  local d
  local col1
  local cfile
  local _rstatus
  local index
  local rate
  local size
  local alpha
  local genomesize
  local output
  local long_read
  local randomseed
  local _disassemble_dir
  local _summary_table
  local _ga_input_short_reads
  local _ptdna
  local _ptdir
  local _ref_ptdna
  local summary_j_ordered
  local _disassemble_i
  local _index_table
  local _mtcontigname
  local _include
  local _exclude
  local _stage_array
  local _disassemble_i
  local _disassemble_i_stage
  local _summary1_ordered
  local _disassemble_params_file
  local _run_type

  if [[ "${_arg_jellyfish_s_is}" == "off" ]]; then
    _arg_jellyfish_s="2G"
  fi
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  local _ga_outdir="${_arg_outdir}/${_arg_inum}"
  local _disassemble_dir="${_ga_outdir}/disassemble"
  local _ga_input_short_reads="${_ga_outdir}/s.fq"

  # if ! run_check_flye; then
  # 	_polap_log0 "Suggestion: (polap) $ conda install goshng::cflye"
  # 	exit $EXIT_ERROR
  # fi

  # plastid genome is assembled; not for mitochondrial genome assembly
  _arg_plastid="on"

  # convert --disassemble-a and --disassemble-b

  # base unit conversion to a regular positive number: 1k -> 1000
  _arg_disassemble_s=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_s}")
  _arg_disassemble_a=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_a}")
  _arg_disassemble_b=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_b}")
  _arg_disassemble_m=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_m}")
  _arg_disassemble_s_max=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_s_max}")
  _arg_genomesize_a=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize_a}")
  _arg_genomesize_b=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize_b}")
  if [[ -n "${_arg_genomesize}" ]]; then
    _arg_genomesize=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize}")
  fi

  if [[ "${_arg_disassemble_q_is}" == "off" ]]; then
    _arg_disassemble_q=${_arg_disassemble_p}
  fi

  help_message=$(
    cat <<HEREDOC
Plastid genome assembly by subsampling long-read data without references

Inputs
------

- long-read data: ${_arg_long_reads} (default: l.fq)
- short-read data 1: ${_arg_short_read1} (default: s1.fa)
- short-read data 2: ${_arg_short_read2} (default: s2.fa)

Main arguments
--------------

-l ${_arg_long_reads}: long-read fastq file
-a ${_arg_short_read1}: short-read fastq file 1
-b ${_arg_short_read2}: short-read fastq file 2

--downsample ${_arg_downsample}: maximum genome coverage to downsample
--disassemble-n ${_arg_disassemble_n}: the number of steps in stage 1
--disassemble-p ${_arg_disassemble_p}: the maximum percent of long-read data
--disassemble-r ${_arg_disassemble_r}: the number of replicates in stages 2/3

Two use cases:
1. case inference: default
2. case check: compare one selected from stage 2 (no alignment in stages 1 and 2)
  --disassemble-c and --disassemble-align-reference

For the cases inference or check we use use either short-read 
subsampling-based polishing or short-read no-subsampling (simple) polishing
  --simple-polishing on or off

Outputs
-------

- plastid genome assembly: ${_arg_outdir}/ptdna.${_arg_inum}.fa

Arguments
---------

-o ${_arg_outdir}: output folder
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file 1
-b ${_arg_short_read2}: a short-read fastq data file 2
-i ${_arg_inum}: output number creating folder: ${_arg_outdir}/${_arg_inum}
-t ${_arg_threads}: the number of CPU cores
--disassemble-i ${_arg_disassemble_i}: the index in disassemble any string
--downsample ${_arg_downsample}: maximum genome coverage to downsample
--disassemble-a ${_arg_disassemble_a}: int for base pairs or .float for rate the smallest base pairs for a subsampling range
--disassemble-b ${_arg_disassemble_b}: int for base pairs or .float for rate the largest base pairs for a subsampling range
--disassemble-p ${_arg_disassemble_p}: the percentile of the largest long read, -a/-b or -p
--disassemble-n ${_arg_disassemble_n}: the number of steps
--disassemble-m ${_arg_disassemble_m}: the upper bound for a Flye assembly
--disassemble-memory ${_arg_disassemble_memory}: the maximum memory in Gb
--disassemble-alpha ${_arg_disassemble_alpha}: the starting Flye's disjointig coverage
--disassemble-delta ${_arg_disassemble_delta}: the move size of alpha (0.1 - 1.0)
--disassemble-c <FASTA>: a single reference sequence in FASTA
--random-seed <arg>: 5-digit number or 0 for random seed
--disassemble-r ${_arg_disassemble_r}: the number of replicates

Menus
-----

- help: display this help message
- view: show some results
- downsample: delete output files 
- archive: archive output files leaving out too large files
- ptgaul: extract ptDNA from ${_arg_outdir}/ptgaul

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1} -b ${_arg_short_read2}

Examples
--------
$(basename "$0") ${_arg_menu[0]} example <NUMBER>

NUMBER=1
--------
$(basename "$0") x-ncbi-fetch-sra --sra SRR7153095
$(basename "$0") x-ncbi-fetch-sra --sra SRR7161123

NUMBER=2
--------
cp -s SRR7153095.fastq l.fq
cp -s SRR7161123_1.fastq s1.fq
cp -s SRR7161123_2.fastq s2.fq
$(basename "$0") disassemble

NUMBER=3
--------
$(basename "$0") get-mtdna --plastid --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta o/ptdna-reference.fa

NUMBER=4
--------
$(basename "$0") disassemble --disassemble-i 1 --stages-include 3 \
  -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq \
  --disassemble-align-reference --disassemble-c o/ptdna-reference.fa

NUMBER=5
--------
mkdir -p o/0/mafft
$(basename "$0") mafft-mtdna -a o/ptdna-reference.fa \
  -b o/0/disassemble/2/pt.subsample-polishing.reference.aligned.1.fa \
  -o o/0/mafft >o/0/mafft/log.txt
cat o/0/mafft/pident.txt

NUMBER=6
--------
$(basename "$0") disassemble --disassemble-i 2 \
  -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq \
  --disassemble-align-reference --disassemble-c o/ptdna-reference.fa

NUMBER=7
--------
$(basename "$0") disassemble --disassemble-i 3 \
  -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq \
  --disassemble-c o/ptdna-reference.fa

HEREDOC
  )

  dev_message=$(
    cat <<HEREDOC
Plastid genome assembly by subsampling long-read data without references

Internals for development version
---------------------------------
--start-index <INDEX>: from <INDEX>
--end-index <INDEX>: to <INDEX> - 1
--steps-include <STPES>: STEPS can be 1,2,3 or 1-15
--steps-exclude <STPES>: STEPS can be 1,2,3 or 1-15
--stages-include <STAGES> for internal use
--stages-exclude <STAGES> for internal use
--disassemble-stop-after stage1: stop after stage 1
--disassemble-stop-after stage2: stop after stage 2

Examples
--------
$(basename "$0") ${_arg_menu[0]} example <NUMBER>

NUMBER=8
--------
bash ptgaul/ptGAUL.sh -o o-ptgaul -r ptdna-Eucalyptus_pauciflora-known.fa -g 180000 -l long.fastq
cp -pr o-ptgaul/result_3000 o/ptgaul

More menu examples
------------------
$(basename "$0") disassemble view
$(basename "$0") disassemble view 0
$(basename "$0") disassemble view x
$(basename "$0") disassemble view best 0
$(basename "$0") disassemble view best x
$(basename "$0") disassemble report
$(basename "$0") disassemble report coverage <- with known ptDNA
$(basename "$0") disassemble best 0 46
$(basename "$0") disassemble best x 2
$(basename "$0") disassemble archive
$(basename "$0") disassemble archive polishing <- backup msbwt as well
$(basename "$0") disassemble polishing
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
  [[ ${_arg_menu[1]} == "dev" ]] && _polap_echo0 "${dev_message}" && return 0
  [[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then

    _disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"
    if [[ "${_arg_menu[2]}" == "1" ]]; then
      _polap_log0_cat "${_disassemble_i}/1/summary1.md"
    fi

    if [[ "${_arg_menu[2]}" == "2" ]]; then
      _polap_log0_cat "${_disassemble_i}/2/summary1.md"
    fi

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  # menu: polishing
  # no need of a known ptDNA; use one of the paths
  # 1. copy the selected j's ptdna.0.fa
  # disassemble/4/2/j/52-mtdna/ptdna.0.fa -> disassemble/4/pt.0.fa
  # 2. polish the pt.0.fa -> pt.1.fa
  if [[ "${_arg_menu[1]}" == "polishing" ]]; then
    _disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"
    _index_file="${_disassemble_i}/stage2-selected-index.txt"

    rm -f "${_index_file}"
    # best j
    local _summary1_ordered="${_disassemble_i}/2/summary1-ordered.txt"
    if [[ -s "${_summary1_ordered}" ]]; then
      j=$(awk 'NR==2 {print $1}' "${_summary1_ordered}")
      echo "${j}" >"${_index_file}"
    else
      _polap_log0 "ERROR: no such file: ${_summary1_ordered}"
      return 1
    fi

    _var_mtdna="${_disassemble_i}/2/${j}/52-mtdna"
    _polap_log0 "Best: ${_var_mtdna}"

    _selected_fasta="${_var_mtdna}/ptdna.0.fa"
    _arg_unpolished_fasta="${_disassemble_dir}/${_arg_disassemble_i}/pt.0.fa"

    cp -p "${_selected_fasta}" "${_arg_unpolished_fasta}"

    _arg_final_assembly="${_disassemble_dir}/${_arg_disassemble_i}/pt.1.fa"
    _polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
    _run_polap_polish
    _polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  # menu: report
  if [[ "${_arg_menu[1]}" == "report" ]]; then

    # report 1 --disassemble-i 3
    if [[ "${_arg_menu[2]}" == "1" ]]; then
      if [[ "${_arg_menu[3]}" == "thirdfile" ]]; then
        _disassemble_report1
      else
        _disassemble_report1 true
      fi
    fi

    # report 2 --disassemble-i 3
    if [[ "${_arg_menu[2]}" == "2" ]]; then
      if [[ "${_arg_menu[3]}" == "thirdfile" ]]; then
        _disassemble_report2
      else
        _disassemble_report2 true
      fi
    fi

    if [[ "${_arg_menu[2]}" == "3" ]]; then
      _disassemble_report3
    fi

    if [[ "${_arg_menu[2]}" == "3x" ]]; then
      _disassemble_report3x "${_arg_menu[3]}"
    fi

    # report 3 --disassemble-i 3
    if [[ "${_arg_menu[2]}" == "4" ]]; then
      # concatenate all summary1.txt to the summary table.
      s0="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1.txt"
      d1="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-scatter.txt"
      d2="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-scatter.pdf"
      _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
        --table ${s0} \
        --out ${d1} \
        --plot ${d2} \
			  --scatter"
    fi

    if [[ "${_arg_menu[2]}" == "9" ]]; then
      # concatenate all summary1.txt to the summary table.
      for i in "${_disassemble_dir}"/*/; do
        if [ -d "$i" ]; then
          i=$(basename "${i%/}")
          if [[ -s "${_disassemble_dir}/${i}/summary1.txt" ]]; then
            if [[ "${_arg_menu[2]}" == "coverage" ]]; then
              _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
        --table ${_disassemble_dir}/${i}/summary1.txt \
        --out ${_disassemble_dir}/${i}/summary1-ordered.txt \
        --plot ${_disassemble_dir}/${i}/summary1-ordered.pdf \
			  --coverage"
            else
              _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
        --table ${_disassemble_dir}/${i}/summary1.txt \
        --out ${_disassemble_dir}/${i}/summary1-ordered.txt \
        --plot ${_disassemble_dir}/${i}/summary1-ordered.pdf"
            fi
            _polap_log0_column ${_disassemble_dir}/${i}/summary1-ordered.txt
          fi
        fi
      done
    fi

    return 0
  fi

  # menu: best
  if [[ "${_arg_menu[1]}" == "best" ]]; then
    _disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"
    _disassemble_i_stage="${_disassemble_i}/2"
    _selected_index_stage2="${_disassemble_i}/2/selected-index.txt"

    local _j_best_stage2=$(<"${_selected_index_stage2}")

    _polap_log0 "best ptDNA graph: ${_disassemble_i_stage}/${_j_best_stage2}/assembly_graph.gfa"

    return 0
  fi

  # menu: bandage
  if [[ "${_arg_menu[1]}" == "bandage" ]]; then
    _disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"
    _disassemble_i_stage="${_disassemble_i}/2"
    _selected_index_stage2="${_disassemble_i}/2/selected-index.txt"

    local _j_best_stage2=$(<"${_selected_index_stage2}")
    local _gfa_best_stage2="${_disassemble_i_stage}/${_j_best_stage2}/assembly_graph.gfa"
    local _bandage_best_stage2="${_disassemble_i}/pt.1.gfa"
    local _png_best_stage2="${_disassemble_i}/pt.1.png"

    cp -p "${_gfa_best_stage2}" "${_bandage_best_stage2}"

    Bandage image "${_gfa_best_stage2}" \
      "${_png_best_stage2}" \
      --colour uniform \
      --unicolpos \#EEEEEE \
      --singlearr \
      --toutline 1 \
      --names \
      --lengths \
      --depth \
      --fontsize 3

    _polap_log0 "Polap ptDNA gfa: ${_bandage_best_stage2}"
    _polap_log1 "Polap ptDNA png: ${_png_best_stage2}"

    return 0
  fi

  # menu: archive
  # archive the disassemble analysis
  if [[ "${_arg_menu[1]}" == "archive" ]]; then
    _polap_log1 "archiving ${_arg_outdir} to ${_arg_archive} ... upto ${_arg_max_filesize}"

    _arg_menu[1]="cflye"
    _run_polap_archive

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  if [[ "${_arg_menu[1]}" == "downsample" ]]; then

    _disassemble-stage0

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  if [[ "${_arg_menu[1]}" == "redownsample" ]]; then

    _disassemble-redownsample

    # read lx.txt
    # read sx.txt

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  # menu: ptgaul - 2025-06-04
  # ptgaul menu should be in a different module: polap-cmd-ptgau.sh.

  # menu: ptgaul
  # ptgaul
  if [[ "${_arg_menu[1]}" == "ptgaul" ]]; then

    _contigger_edges_gfa="${_arg_outdir}/ptgaul/flye_cpONT/assembly_graph.gfa"
    _outdir="${_arg_outdir}/ptgaul/flye_cpONT/ptdna"

    _mtcontigname="${_outdir}/mt.contig.name"
    _arg_unpolished_fasta="${_outdir}/circular_path_1_concatenated.fa"
    _arg_final_assembly="${_outdir}/pt.1.fa"

    if [[ "${_arg_menu[2]}" == "outfile" ]]; then
      _polap_log3_cmd mkdir -p "${_outdir}"
      gfatools stat -l "${_contigger_edges_gfa}" \
        >"${_outdir}/graph_final.txt" \
        2>"$_polap_output_dest"
      _summary_gfa_number_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')

      if [[ "${_summary_gfa_number_segments}" -eq 3 ]]; then
        _polap_log1 "ptGAUL assembly's segments: ${_summary_gfa_number_segments}"
        _polap_log1 "you may use ptgaul menu:"
        _polap_log1 "  ptgaul 1"
        _polap_log1 "  ptgaul 1 1"
        _polap_log1 "  ptgaul 1 3"
        _polap_log1 "  ptgaul 2"
        _polap_log1 "  ptgaul 3"
        # mt.contig.name
        _polap_log1 "create file: ${_mtcontigname}"
        local string="edge_1,edge_2,edge_3"
        echo "$string" | tr ',' '\n' >"${_mtcontigname}"
        # extract ptDNA from ptGAUL's result
        _polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"
        # polish one of the ptDNA sequences
        _polap_log1 "  polishing ptDNA: ${_arg_unpolished_fasta}"
        _run_polap_polish
        exit_code=$?
        # If a function or command is killed due to an out-of-memory (OOM) event,
        # it usually exits with signal 137 (SIGKILL).
        if [ $exit_code -eq 137 ]; then
          die "polishing was killed due to OOM: an out-of-memory event"
        fi
        if [[ -s "${_arg_final_assembly}" ]]; then
          _polap_log1 "  ptGAUL polished assembly: ${_arg_final_assembly}"
        else
          die "ERROR: no such ptGAUL polished assembly: ${_arg_final_assembly}"
        fi
      elif [[ "${_summary_gfa_number_segments}" -eq 1 ]]; then
        _polap_log1 "  ptGAUL assembly's segments: ${_summary_gfa_number_segments}"
        _arg_unpolished_fasta="${_arg_outdir}/ptgaul/ptGAUL_final_assembly/final_assembly.fasta"
        _polap_log1 "  ptGAUL assembly: ${_arg_unpolished_fasta}"
        _polap_log1 "  polishing ptDNA: ${_arg_unpolished_fasta}"
        _run_polap_polish
        exit_code=$?
        # If a function or command is killed due to an out-of-memory (OOM) event,
        # it usually exits with signal 137 (SIGKILL).
        if [ $exit_code -eq 137 ]; then
          die "polishing was killed due to OOM: an out-of-memory event"
        fi
        if [[ -s "${_arg_final_assembly}" ]]; then
          _polap_log1 "  ptGAUL polished assembly: ${_arg_final_assembly}"
        else
          die "ERROR: no such ptGAUL polished assembly: ${_arg_final_assembly}"
        fi
      else
        # extract three unique edges from the gfa for a ptDNA candidate
        #
        # L	edge_3	+	edge_4	+	0M	RC:i:41
        # L	edge_3	+	edge_4	-	0M	RC:i:40
        # L	edge_3	-	edge_5	-	0M	RC:i:42
        # L	edge_3	-	edge_5	+	0M	RC:i:41
        # P	contig_5	edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3-	*
        # P	contig_1	edge_1+	*
        # P	contig_2	edge_2+	*
        #
        bash ${_POLAPLIB_DIR}/run-polap-sh-extract-three-edges-of-ptdna.sh \
          "${_contigger_edges_gfa}" \
          "${_mtcontigname}"
        if [[ -s "${_mtcontigname}" ]]; then
          # extract ptDNA from ptGAUL's result
          _polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"
          # polish one of the ptDNA sequences
          _polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
          _run_polap_polish
          exit_code=$?
          # If a function or command is killed due to an out-of-memory (OOM) event,
          # it usually exits with signal 137 (SIGKILL).
          if [ $exit_code -eq 137 ]; then
            die "polishing was killed due to OOM: an out-of-memory event"
          fi
          if [[ -s "${_arg_final_assembly}" ]]; then
            _polap_log1 "  ptGAUL polished assembly: ${_arg_final_assembly}"
          else
            die "ERROR: no such ptGAUL polished assembly: ${_arg_final_assembly}"
          fi
        else
          _polap_log0 "ptGAUL assembly's segment counts: ${_summary_gfa_number_segments}"
          _polap_log0 "you may not use ptgaul menu at the moment."
          _polap_log0 "check: ${_contigger_edges_gfa}"
          _polap_log0 "create: ${_mtcontigname} with e.g., edge_1 that is originated from the ptDNA."
          _polap_log0 "then execute:disassemble ptgaul 2"
          _polap_log0 "polap disassemble ptgaul 2 -o ${_arg_outdir}"
          _polap_log0 "polap disassemble ptgaul 3 -o ${_arg_outdir}"
        fi
      fi
    fi

    # _mtcontigname
    if [[ "${_arg_menu[2]}" == "1" ]]; then
      if [[ "${_arg_menu[3]}" == "1" ]]; then
        _polap_log0 "create file: ${_mtcontigname}"
        local string="edge_1"
        echo "$string" | tr ',' '\n' >"${_mtcontigname}"
      elif [[ "${_arg_menu[3]}" == "3" ]]; then
        _polap_log0 "create file: ${_mtcontigname}"
        local string="edge_1,edge_2,edge_3"
        echo "$string" | tr ',' '\n' >"${_mtcontigname}"
      fi
    fi

    # extract
    if [[ "${_arg_menu[2]}" == "2" ]]; then
      _polap_log0 "extract ptDNA: ${_arg_unpolished_fasta}"
      _polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"
    fi

    # polishing
    if [[ "${_arg_menu[2]}" == "3" ]]; then
      _polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
      _run_polap_polish
      _polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"
    fi

    if [[ "${_arg_menu[2]}" == "4" ]]; then
      for ((i = 1; i <= 4; i++)); do
        _arg_unpolished_fasta="${_outdir}/circular_path_${i}_concatenated.fa"
        _arg_final_assembly="${_outdir}/pt.${i}.fa"
        _polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
        _run_polap_polish
        _polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"
      done
    fi

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  # menu: polish
  if [[ "${_arg_menu[1]}" == "polish" ]]; then

    _run_polap_polish-disassemble

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  # menu: example
  if [[ "${_arg_menu[1]}" == "example" ]]; then
    if [[ "${_arg_menu[2]}" == "1" ]]; then
      _polap_log0 "Downloading sequencing data from NCBI SRA database: SRR7153095 and SRR7161123 ... takes time depending on your network speed ... "
      _polap_log0 "  it takes time depending on your network speed ... be patient!"
      _polap_log0 "  long-read SRA: SRR7153095 ..."
      "$0" x-ncbi-fetch-sra --sra SRR7153095
      _polap_log0 "  short-read SRA: SRR7161123 ..."
      "$0" x-ncbi-fetch-sra --sra SRR7161123
    fi

    if [[ "${_arg_menu[2]}" == "2" ]]; then
      _polap_log0 "Assembling the plastid genome by subsampling NCBI SRA database: SRR7153095 and SRR7161123 ..."
      cp -s SRR7153095.fastq l.fq
      cp -s SRR7161123_1.fastq s1.fq
      cp -s SRR7161123_2.fastq s2.fq
      "$0" disassemble \
        -l SRR7153095.fastq \
        -a SRR7161123_1.fastq \
        -b SRR7161123_2.fastq
      _polap_log0 "  plastid genome assembly: o/ptdna.0.fa"
    fi

    if [[ "${_arg_menu[2]}" == "3" ]]; then
      _polap_log0 "Downloading plastid genome sequences for Eucalyptus pauciflora ... -> o/ptdna-reference.fa"
      "$0" get-mtdna --plastid --species "Eucalyptus pauciflora"
      cp o/00-bioproject/2-mtdna.fasta o/ptdna-reference.fa
    fi

    if [[ "${_arg_menu[2]}" == "4" ]]; then
      # subshell problem
      # local _a=($(_disassemble-step14 "${_var_mtdna}"))
      _polap_log0 "Assembling the plastid genome and check it with the reference ..."
      _polap_log0 "  execute the following:"
      "$0" disassemble \
        --stages-include 3 \
        --disassemble-i 1 \
        -l SRR7153095.fastq \
        -a SRR7161123_1.fastq \
        -b SRR7161123_2.fastq \
        --disassemble-align-reference \
        --disassemble-c o/ptdna-reference.fa
    fi

    if [[ "${_arg_menu[2]}" == "5" ]]; then
      _polap_log0_cmd mkdir -p o/0/mafft
      "$0" mafft-mtdna \
        -a o/ptdna-reference.fa \
        -b o/0/disassemble/2/pt.subsample-polishing.reference.aligned.1.fa \
        -o o/0/mafft \
        >o/0/mafft/log.txt
      _polap_log0 "see o/0/mafft/pident.txt"
      _polap_log0_cat "o/0/mafft/pident.txt"
    fi

    if [[ "${_arg_menu[2]}" == "6" ]]; then
      # subshell problem
      # local _a=($(_disassemble-step14 "${_var_mtdna}"))
      _polap_log0 "Assembling the plastid genome and check it with the reference ..."
      _polap_log0 "  execute the following:"
      "$0" disassemble \
        --disassemble-i 2 \
        -l SRR7153095.fastq \
        -a SRR7161123_1.fastq \
        -b SRR7161123_2.fastq \
        --disassemble-align-reference \
        --disassemble-c o/ptdna-reference.fa
    fi

    if [[ "${_arg_menu[2]}" == "7" ]]; then
      _polap_log0 "Assembling the plastid genome and comparing it with the ptGAUL assembly ..."
      _polap_log0 "  execute the following:"
      "$0" disassemble \
        --disassemble-i 3 \
        -l SRR7153095.fastq \
        -a SRR7161123_1.fastq \
        -b SRR7161123_2.fastq \
        --disassemble-c o/ptdna-reference.fa
    fi

    if [[ "${_arg_menu[2]}" == "8" ]]; then
      _polap_log0 "Assembling the plastid genome using ptGAUL ..."
      bash ${_POLAPLIB_DIR}/polap-ptGAUL1.sh \
        -o ptgaul \
        -r ptdna-reference.fa \
        -g 160000 \
        -l SRR7153095.fastq
      mv ptgaul/result_3000 o/ptgaul
      rm ptgaul

      _polap_log0 "Preparing the short-read polishing: This polishing preparation takes long ..."
      $0 prepare-polishing \
        -a SRR7161123_1.fastq \
        -b SRR7161123_2.fastq

      _polap_log0 "Polishing the ptGAUL draft sequence ..."
      "$0" disassemble ptgaul
    fi

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  _polap_lib_conda-ensure_conda_env polap-cflye || exit 1
  if ! run_check_flye; then
    _polap_log0 "Suggestion: (polap) $ conda install goshng::cflye"
    exit $EXIT_ERROR
  fi
  conda deactivate

  # main

  # save the user provided options
  local _b_disassemble_a="${_arg_disassemble_a}"
  local _b_disassemble_b="${_arg_disassemble_b}"
  local _b_disassemble_s="${_arg_disassemble_s}"
  local _b_disassemble_p="${_arg_disassemble_p}"
  local _b_disassemble_n="${_arg_disassemble_n}"
  local _b_disassemble_b_is="${_arg_disassemble_b_is}"
  local _b_disassemble_n_is="${_arg_disassemble_n_is}"

  # select stages
  # FIXME: sometimes, this is set to different thing
  # even though --stages-include was set correcly.
  # local _stages="--stages-include 1-3"
  # command time -v ${_polap_cmd} disassemble \
  # CORRECT:
  # 	${_stages} \
  # WRONG:
  # 	"${_stages}" \
  # _polap_log3 "  BUG1: --stages-include: ${_arg_stages_include}"
  # _polap_log3 "  BUG2: --stages-include is: ${_arg_stages_is}"
  if [[ -z "${_arg_stages_include}" ]]; then
    _arg_stages_include="0-3"
    _arg_stages_exclude=""
  fi
  # _polap_log3 "  BUG3: --stages-include: ${_arg_stages_include}"

  if [[ -z "${_arg_disassemble_c}" ]]; then
    # case: infer
    _run_type="inference"
    _polap_log2 "case: infer"
  else
    # case: check or compare
    if [[ "${_arg_disassemble_align_reference}" == "off" ]]; then
      # compare
      _run_type="compare"
      _polap_log2 "case: compare"
      _polap_log2 "  turn on --disassemble-simple-polishing"
      _arg_disassemble_simple_polishing="on"
    else
      # check
      _run_type="check" # check the inference
      _polap_log2 "case: check"
    fi
  fi

  _polap_log0 "Assemble plastid genomes by subsampling long-read data (i: ${_arg_inum}, disassemble index: ${_arg_disassemble_i}, run type: ${_run_type}, full short-read polishing: ${_arg_disassemble_simple_polishing})"

  # if [[ "${_arg_menu[1]}" == "stage1" ]]; then
  # 	_arg_stages_include="0-1"
  # 	_arg_stages_exclude=""
  # elif [[ "${_arg_menu[1]}" == "stage2" ]]; then
  # 	_arg_stages_include="2-4"
  # 	_arg_stages_exclude=""
  # elif [[ "${_arg_menu[1]}" == "stage3" ]]; then
  # 	_arg_stages_include="5-6"
  # 	_arg_stages_exclude=""
  # fi
  _include="${_arg_stages_include}"
  _exclude="${_arg_stages_exclude}" # Optional range or list of steps to exclude
  _stage_array=()
  _stage_array=($(_polap_parse_steps "${_include}" "${_exclude}"))
  _polap_log2 "  executing stages: ${_stage_array[*]}"

  # FIXME: input data - l.fq, s1.fq, s2.fq might worth including them.
  # currently we strictly use the options or not specified then we won't use
  # them.
  #
  # BioProject
  # -l -a -b options
  # default -l -a -b files
  _polap_log1 "  input1: long-read: ${_arg_long_reads}"
  _polap_log1 "  input2: short-read1: ${_arg_short_read1}"
  _polap_log1 "  input3: short-read2: ${_arg_short_read2}"
  [[ -s "${_arg_long_reads}" ]] || return ${_POLAP_ERR_CMD_OPTION_LONGREAD}
  [[ -s "${_arg_short_read1}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
  [[ -s "${_arg_short_read2}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
  if [[ "${_arg_long_reads_is}" == "off" ]]; then
    _polap_log0 "WARNING: we use the default long-read: ${_arg_long_reads}"
    _polap_log0 "  you do not use -l option for a long-read FASTQ file"
  fi
  if [[ "${_arg_short_read1_is}" == "off" ]]; then
    _polap_log0 "WARNING: we use the default short-read1: ${_arg_short_read1}"
    _polap_log0 "  you do not use -a option for a short-read FASTQ file"
  fi
  if [[ "${_arg_short_read2_is}" == "off" ]]; then
    _polap_log0 "WARNING: we use the default short-read2: ${_arg_short_read2}"
    _polap_log0 "  you do not use -b option for a short-read FASTQ file"
  fi

  # if [[ "${_arg_short_read1_is}" == "off" ]]; then
  # 	_polap_log1 "  input2: no short-read1"
  # else
  # 	_polap_log1 "  input2: short-read1: ${_arg_short_read1}"
  # 	[[ -s "${_arg_short_read1}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
  # fi
  # if [[ "${_arg_short_read2_is}" == "off" ]]; then
  # 	_polap_log1 "  input3: no short-read2"
  # else
  # 	_polap_log1 "  input3: short-read2: ${_arg_short_read2}"
  # 	[[ -s "${_arg_short_read2}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
  # fi

  # Iteration over p x n for each i in --disassemble-i
  #
  # output: o/disassemble + --disassemble-i
  # output: o/disassemble/i
  _disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"

  # Downsample the input data with a target coverage
  # input: long- and short-read data
  # output1: tmp/l50x.fq
  # output2: tmp/s50x_1.fq
  # output3: tmp/s50x_2.fq
  _polap_log3_cmd mkdir -p "${_disassemble_i}"
  if _polap_contains_step 0 "${_stage_array[@]}"; then
    _polap_log1 "  stage 0: downsample the input data, if possible, to ${_arg_downsample}x"

    if [[ -d "${_disassemble_i}" ]]; then
      _polap_log1 "  delete the previous disassemble: ${_disassemble_i}"
      _polap_log3_cmd rm -rf "${_disassemble_i}"
    fi
    _polap_log3_cmd mkdir -p "${_disassemble_i}"
    _polap_log1 "  output: ${_disassemble_i}"

    _disassemble-stage0
  fi

  # short-read polishing preparation
  # output: o/msbwt
  if [[ "${_arg_disassemble_simple_polishing}" == "on" ]]; then
    _run_polap_prepare-polishing
  fi

  # Stage 1
  # output: o/disassemble/i/1
  if _polap_contains_step 1 "${_stage_array[@]}"; then
    _polap_log1 "  stage 1: determine the sample size and the minimum read coverage: run type: ${_run_type}"

    # output: o/disassemble/i/params.txt
    _disassemble_params_file="${_disassemble_i}/params.txt"
    _disassemble_make_params_txt "${_disassemble_params_file}"

    _disassemble-stage1
  fi

  if [[ -n "${_arg_disassemble_stop_after}" ]]; then
    if [[ "${_arg_disassemble_stop_after}" == "stage1" ]]; then
      return 0
    fi
  fi

  # use summary1 to select one
  #
  if _polap_contains_step 2 "${_stage_array[@]}"; then
    _polap_log1 "  stage 2: assemble with subsample-replicate: run type: ${_run_type}"
    _disassemble-stage2
  fi

  if [[ -n "${_arg_disassemble_stop_after}" ]]; then
    if [[ "${_arg_disassemble_stop_after}" == "stage2" ]]; then
      return 0
    fi
  fi

  if _polap_contains_step 3 "${_stage_array[@]}"; then
    if [[ "${_arg_disassemble_simple_polishing}" == "off" ]]; then
      _polap_log1 "  stage 3: assemble with subsample-polishing: run type: ${_run_type}"
    else
      _polap_log1 "  stage 3: assemble with simple-polishing: run type: ${_run_type}"
    fi
    _disassemble-stage3

    # link the output
    local _s="${_disassemble_i}/pt.subsample-polishing.1.fa"
    local _d="${_arg_outdir}/ptdna.${_arg_inum}.fa"
    if [[ -s "${_s}" ]]; then
      ln -sfn $(realpath "${_s}") "${_d}"
      _polap_log0 "  Final plastid genome assembly (subsampling-based): ${_d}"
    fi

    local _s="${_disassemble_i}/pt.subsample-polishing.reference.aligned.1.fa"
    local _d="${_arg_outdir}/ptdna.ref.${_arg_inum}.fa"
    if [[ -s "${_s}" ]]; then
      ln -sfn $(realpath "${_s}") "${_d}"
      _polap_log0 "  Final plastid genome assembly (subsampling-based and reference-aligned): ${_d}"
    fi

    local _s="${_disassemble_i}/pt.simple-polishing.1.fa"
    local _d="${_arg_outdir}/ptdna.simple.${_arg_inum}.fa"
    if [[ -s "${_s}" ]]; then
      ln -sfn $(realpath "${_s}") "${_d}"
      _polap_log0 "  Final plastid genome assembly (full short-read polishing): ${_d}"
    fi

    local _s="${_disassemble_i}/pt.simple-polishing.reference.aligned.1.fa"
    local _d=${_arg_outdir}/ptdna.simple.ref.${_arg_inum}.fa
    if [[ -s "${_s}" ]]; then
      ln -sfn $(realpath "${_s}") "${_d}"
      _polap_log0 "  Final plastid genome assembly (full short-read polishing and reference-aligned): ${_d}"
    fi
  fi

  # Disable debugging if previously enabled
  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

function _run_polap_step-disassemble {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  _arg_plastid="on"
  local _msg_level=2
  local i
  local j
  local col1 col2 col3
  local _command1
  local _disassemble_dir
  local _ga_input_short_reads
  local _include
  local _exclude
  local _datasize_array
  local _datasize_array_length
  local _i_last
  local _index_table
  local sampling_datasize
  local genomesize
  local _size_index
  local _alpha
  local _summary_long_total
  local _summary_short_total
  local _summary_long_rate_sample
  local _summary_short_rate_sample
  local _summary_long_sample_seed
  local _summary_pre_alpha
  local _summary_elapsed_time
  local _outdir
  local _sampling_datasize_bp
  local _total_length
  local _total_short_read
  local _short_read1
  local _long_read
  local _is_stop
  local _summary1_file
  local resume=""
  local _summary_num_circular_paths
  local _summary_num_circular_nodes
  local _summary_genome_size
  local _contigger_edges_gfa
  local _contigger_edges_fasta
  local _summary_peak_ram_size_gb
  local _summary_peak_ram
  local _outdir_draft_assembly_size
  local _draft_assembly_size
  local _draft_assembly_size_bp
  local _sampling_datasize_bp
  local _summary_gfa_number_segments
  local _summary_gfa_number_links
  local _summary_gfa_total_segment_length
  local _summary_coverage_ref
  local _summary_coverage_target
  local _summary_pident
  local _summary_gc
  local _summary_length
  local _summary_j_candidate
  local _summary_draft_assembly_size=-1
  local coverage
  local _ga_annotation_all
  local _mtcontigname
  local _ptdna
  local _ptdir
  local ptdna_file
  local nseq
  local _total_iterations
  local time_per_iteration
  local remaining_iterations
  local remaining_time
  local _terminal_width
  local status
  local _var_mtdna
  local _restarted_fasta
  local _circular_path_fasta
  local _folder_with_max_coverage
  local _step_array

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  source "${_POLAPLIB_DIR}/polap-function-disassemble-seeds.sh"

  help_message=$(
    cat <<HEREDOC
Plastid genome assembly step-by-step subsampling long-read data

Inputs
------

- long-read data: ${_arg_long_reads}
- short-read data 1: ${_arg_short_read1}
- short-read data 2 (optional): ${_arg_short_read2}

Outputs
-------

- plastid genome assembly
- trace plots for the features of plastid genome assemblies
- ${_ga_outdir}/disassemble/<number>

Arguments
---------

-o ${_arg_outdir}
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file 1
-b ${_arg_short_read2} (optional): a short-read fastq data file 2

Menus
-----

Usages
------

Examples
--------
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
  [[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  if [[ "${_arg_menu[1]}" == "reset" ]] ||
    [[ "${_arg_menu[1]}" == "clean" ]]; then

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  if [[ "${_arg_menu[1]}" == "example" ]]; then

    # Disable debugging if previously enabled
    _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
    [ "$_POLAP_DEBUG" -eq 1 ] && set +x
    return 0
  fi

  local _start_time=$(date +%s)
  # _polap_lib_timing-step-reset

  i=0
  _summary_gfa_number_segments=0
  _summary_gfa_number_links=0
  _summary_gfa_total_segment_length=0
  _summary_num_circular_paths=0
  _summary_num_circular_nodes=0
  length=0
  _summary_gc=0
  coverage=0

  local _ga_outdir="${_arg_outdir}/${_arg_inum}"
  local _arg1="${1:-$_arg_disassemble_j}"
  _disassemble_dir="${_ga_outdir}/disassemble/${_arg_disassemble_i}/${_arg1}"

  local ga_long_reads="${_ga_outdir}/l.fq"
  local ga_short_read1="${_ga_outdir}/s_1.fq"
  local ga_short_read2="${_ga_outdir}/s_2.fq"

  _ga_input_short_reads="${_ga_outdir}/s.fq"
  _index_table="${_disassemble_dir}/index.txt"
  _short_read1=${_disassemble_dir}/s.fq
  _long_read=${_disassemble_dir}/l.fq.gz
  _summary1_file="${_disassemble_dir}/summary1.txt"

  _polap_log1 "  assemble a plastid genome without reference"

  if [[ "${_arg1}" -eq 1 ]] || [[ "${_arg1}" -eq 2 ]]; then
    _polap_log2 "    disassemble stage ${_arg1}."
  else
    _polap_log0 "ERROR: disassemble stage must be 1 or 2."
  fi

  if [[ -s "${ga_long_reads}" ]]; then
    _polap_log2 "    input1: long-read: ${ga_long_reads}"
  else
    _polap_log2 "    input1: no long-read"
  fi

  if [[ -s "${_ga_input_short_reads}" ]]; then
    _polap_log2 "    found: short-read1 + short-read2: ${_ga_input_short_reads}"
  else
    if [[ -s "${ga_short_read1}" ]]; then
      die "ERROR: input2: ${ga_short_read1} must have been deleted."
    fi

    set -o
    if [[ -s "${ga_short_read2}" ]]; then
      die "ERROR: input3: ${ga_short_read2} must have been deleted."
    fi
  fi

  _polap_log2 "    input4: option --plastid: ${_arg_plastid}"

  if [[ -n "${_arg_species}" ]]; then
    _polap_log2 "    input6: species: ${_arg_species}"
    _run_polap_get-mtdna
    _arg_disassemble_c="${_arg_outdir}/00-bioproject/2-mtdna.fasta"
  elif [[ -n "${_arg_disassemble_c}" ]]; then
    _polap_log2 "    input6: compare to fasta: ${_arg_disassemble_c}"
    check_file_existence "${_arg_disassemble_c}"
  else
    _polap_log2 "    input6: neither species nor fasta file is given"
    _polap_log2 "    inference of ptDNA"
    _polap_log2 "    one of the sequences of an assembly graph is chosen"
  fi
  _polap_log2 "  output: ${_disassemble_dir}"

  # select steps
  local _include="${_arg_steps_include}"
  local _exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
  _step_array=()

  _step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))
  _rstatus="$?"
  if [[ "${_rstatus}" -ne 0 ]]; then
    _polap_log0 "ERROR: Error parsing steps."
    return "${_POLAP_ERR_CMD_OPTION_STEPS}"
  fi

  _polap_log2 "  executing steps: ${_step_array[*]}"

  # steps
  # 1. long-read count
  # 2. short-read concatenate
  # 3. short-read count
  # 4. loop index
  # 5. long-read sample
  # 6. short-read sample
  # 7. genome size estimate
  # 8. flye run
  # 9. alpha adjustment
  # 10. assembly summary
  # 11. annotate
  # 12. seeds
  # 13. ptDNA extraction
  # 14. best ptDNA
  # 15. ptDNA evaluation
  # 16.

  # 1. long-read count
  if _polap_contains_step 1 "${_step_array[@]}"; then
    _polap_log1 "  step 1: count the total number of bases in the input long-read data"
    _polap_log2 "    input1: downsampled long-read fastq: ${ga_long_reads}"
    _polap_log2 "    output1: long-read length: ${_ga_outdir}/long_total_length.txt"
    if [[ -s "${_ga_outdir}/long_total_length.txt" ]]; then
      _polap_log2 "    found: ${_ga_outdir}/long_total_length.txt"
    else
      _disassemble-step1 "${ga_long_reads}" "${_ga_outdir}/long_total_length.txt"
    fi
  fi
  check_file_existence "${_ga_outdir}/long_total_length.txt"
  local _ga_long_total_length=$(<"${_ga_outdir}/long_total_length.txt")

  # concatenation done during the downsampling at -stage0 function
  #
  # 2. short-read concatenate
  # combine short-read dataset files
  # create o/s.fq to deal with a single or two short-read files
  # overwrite the s.fq
  # if _polap_contains_step 2 "${_step_array[@]}"; then
  # 	_polap_log1 "  step 2: concatenate the short-read data if there are two"
  # 	_polap_log2 "    input1: command option short-read1: ${ga_short_read1}"
  # 	_polap_log2 "    input2: command option short-read2: ${ga_short_read2}"
  # 	_polap_log2 "    output: concatenated short-reads: ${_ga_input_short_reads}"
  # 	if [[ -s "${_ga_input_short_reads}" ]]; then
  # 		_polap_log2 "    found: concatenated short-reads: ${_ga_input_short_reads}"
  # 		if [[ -s "${ga_short_read1}" ]] || [[ -s "${ga_short_read2}" ]]; then
  # 			_polap_log0 "WARNING: we should not have such file: ${ga_short_read1}"
  # 			_polap_log0 "WARNING: we should not have such file: ${ga_short_read2}"
  # 			die "ERROR: we must have deleted each short-read data files."
  # 		fi
  # 	else
  # 		die "ERROR: we must have ${_ga_input_short_reads} in stage 0."
  # 		_disassemble-step2 "${_ga_input_short_reads}" "${ga_short_read1}" "${ga_short_read2}"
  # 		# delete each short reads to save the disk space
  # 		_polap_log3_cmd rm -f "${ga_short_read1}" "${ga_short_read2}"
  # 	fi
  # fi

  # 3. short-read count
  if _polap_contains_step 3 "${_step_array[@]}"; then
    _polap_log1 "  step 3: count the total number of bases in the concatenated short-read data"
    _polap_log2 "    input1: concatenated short-reads: ${_ga_input_short_reads}"
    _polap_log2 "    output1: ${_ga_outdir}/short_total_length.txt"
    check_file_existence "${_ga_input_short_reads}"
    if [[ -s "${_ga_outdir}/short_total_length.txt" ]]; then
      _polap_log2 "    found: ${_ga_outdir}/short_total_length.txt"
    else
      _disassemble-step3 "${_ga_input_short_reads}" "${_ga_outdir}/short_total_length.txt"
    fi
  fi
  check_file_existence "${_ga_outdir}/short_total_length.txt"
  local _ga_short_total_length=$(<"${_ga_outdir}/short_total_length.txt")

  # init
  _alpha="${_arg_disassemble_alpha}"
  if [ -d "${_disassemble_dir}" ] && [ "${_arg_redo}" = "on" ]; then
    _polap_log3_cmd rm -rf "${_disassemble_dir}"
  fi
  _polap_log3_cmd mkdir -p "${_disassemble_dir}"

  # 4. loop index
  if _polap_contains_step 4 "${_step_array[@]}"; then
    _polap_log1 "  step 4: prepare for the loop iteration"

    if [[ "${_arg_disassemble_b_is}" == "on" ]]; then
      _polap_lib_array-make-index \
        "${_index_table}" \
        "${_arg_disassemble_n}" \
        "${_arg_disassemble_a}" \
        "${_arg_disassemble_b}"
    else
      # local _total_length_long=$(<"${_ga_outdir}/long_total_length.txt")
      local _disassemble_b=$(_polap_utility_compute_percentage \
        "${_arg_disassemble_p}" \
        "${_ga_long_total_length}")
      _polap_lib_array-make-index \
        "${_index_table}" \
        "${_arg_disassemble_n}" \
        "${_arg_disassemble_a}" \
        "${_disassemble_b}"
    fi

    # read the 2nd column or sampling sizes into an array
    # local _array_index_table
    # mapfile -t _array_index_table < <(cut -f2 "${_index_table}")
    # local _total_iterations="${#_array_index_table[@]}"
    # _total_iterations="${_arg_disassemble_n}"
    # _polap_lib_timing-reset "${_total_iterations}"

    _polap_lib_table-cflye-main-summary-header "${_summary1_file}"
    _is_stop="off"
  else
    # no iteration
    rm -f "${_index_table}"
    touch "${_index_table}"
    _total_iterations=0
    _is_stop="on"
  fi

  _i_last=0
  _terminal_width=$(tput cols) # Get terminal width
  cat <"${_index_table}" | while IFS=$'\t' read -r i sampling_datasize genomesize; do
    local _summary_sampling_datasize="${sampling_datasize}"
    local _summary_i="${i}"
    if [[ "${_is_stop}" == "on" ]]; then
      break
    fi

    if ((i < _arg_start_index)); then
      continue
    fi

    if [[ -n "${_arg_disassemble_s_max}" ]]; then
      if [[ "${sampling_datasize}" -gt "${_arg_disassemble_s_max}" ]]; then
        break
      fi
    fi

    _polap_lib_timing-set
    _i_last=$((_i_last + 1))
    _summary_pre_alpha="${_alpha}"

    # Init variables in the loop
    _outdir=${_disassemble_dir}/$i
    _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
    _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
    _mtcontigname="${_outdir}/mt.contig.name"
    _ga_annotation_all="${_outdir}/assembly_info_organelle_annotation_count-all.txt"
    _var_mtdna="${_outdir}/52-mtdna"
    _ptdna="${_var_mtdna}/circular_path_1_concatenated.fa"

    _sampling_datasize_bp=$(_polap_utility_convert_bp "${sampling_datasize}")
    _genomesize_bp=$(_polap_utility_convert_bp "${genomesize}")

    _polap_log2 "  begin: i=${i}:S=${_sampling_datasize_bp}:alpha=${_summary_pre_alpha}"

    _polap_log3_cmd mkdir -p "${_outdir}"

    # 5. Sample long-read data
    #
    _summary_long_total=-1
    _summary_long_rate_sample=-1
    _summary_long_sample_seed=-1
    if _polap_contains_step 5 "${_step_array[@]}"; then
      _polap_log1 "  step 5: subsample long-read data"
      _summary_long_total=$(<"${_ga_outdir}/long_total_length.txt")
      _summary_long_rate_sample=$(echo "scale=10; ${sampling_datasize}/${_summary_long_total}" | bc)
      _polap_log2 "    input1: long read data: ${ga_long_reads} (${_summary_long_total} bp)"
      _polap_lib_random-get
      _summary_long_sample_seed=${_polap_var_random_number}
      _polap_log3 "    input2: sampling rate: ${_summary_long_rate_sample}"
      _polap_log3 "    input3: random seed: ${_summary_long_sample_seed}"
      _polap_log2 "    output: ${_long_read}"
      _polap_log3_pipe "seqkit sample \
            -p ${_summary_long_rate_sample} \
            -s ${_summary_long_sample_seed} \
            ${ga_long_reads} \
            -o ${_long_read} 2>${_polap_output_dest}"

      # WARNING: the actual sample size is not equal to the sampling datasize.
      # We might remove the sampling datasize or note that it is intended to be
      # such a sample size by seqkit. But, the value is not actual sampling
      # size but rough estimate. We could count the sampled long-read data.
      # However, it would take computing time without much meaning.
      # So, we do not count the the number of bases in the sampled long-read
      # data.
      # _summary_sampling_datasize="${sampling_datasize}"
    fi

    # 6. Sample short-read data
    # We utilize the identical sampling rate employed
    # by the long-read sampling method, which is set to default.
    #
    _summary_short_total=-1
    _summary_short_rate_sample=-1
    _summary_short_sample_seed=-1
    if _polap_contains_step 6 "${_step_array[@]}"; then
      _polap_log1 "  step 6: subsample the short-read data"
      _polap_log2 "    input1: short read data: ${_ga_input_short_reads}"
      _polap_log2 "    output: ${_short_read1}"
      check_file_existence "${_ga_input_short_reads}"
      # _summary_short_total=$(<"${_ga_outdir}/short_total_length.txt")
      _summary_short_rate_sample="${_summary_long_rate_sample}"
      _polap_lib_random-get
      _summary_short_sample_seed=${_polap_var_random_number}
      _polap_log3 "    input1: sampling rate: ${_summary_short_rate_sample}"
      _polap_log3 "    input2: random seed: ${_summary_short_sample_seed}"
      _polap_log3_pipe "seqkit sample \
            -p ${_summary_short_rate_sample} \
            -s ${_summary_short_sample_seed} \
            ${_ga_input_short_reads} \
            -o ${_short_read1} 2>${_polap_output_dest}"
    fi

    # 7. Estimate the genome assembly size
    #
    _summary_genome_size=-1
    if _polap_contains_step 7 "${_step_array[@]}"; then
      _polap_log1 "  step 7: genome size using the subsample of the short-read data"
      check_file_existence "${_short_read1}"
      _polap_log2 "    input1: ${_short_read1}"
      _polap_log2 "    output: ${_outdir}/short_expected_genome_size.txt"
      _polap_find-genome-size \
        "${_short_read1}" \
        "${_outdir}"
      _summary_genome_size=$(<"${_outdir}/short_expected_genome_size.txt")
      _polap_log2 "    output: estimated genome size: ${_summary_genome_size} bp"
    fi

    # 8. flye
    #
    _summary_peak_ram_size_gb=-1
    if _polap_contains_step 8 "${_step_array[@]}"; then
      _polap_log1 "  step 8: cFlye on the subsample of the long-read data"
      _polap_log2 "    input1: long-read: ${_long_read}"
      _polap_log2 "    input2: genome size: ${_summary_genome_size}"
      _polap_log2 "    input3: alpha: ${_alpha}"
      _polap_log2 "    outdir: ${_outdir}"
      local _a=($(_disassemble-step8 \
        "${_outdir}" \
        "${_long_read}" \
        "${_summary_genome_size}" \
        "${_alpha}" | tail -n 1))
      _is_stop="${_a[0]}"
      _summary_peak_ram_size_gb="${_a[1]}"
    fi

    # 9. Adjust the alpha
    #
    _outdir_draft_assembly_size="${_outdir}/draft_assembly_size.txt"
    if _polap_contains_step 9 "${_step_array[@]}"; then
      _polap_log1 "  step 9: adjust the alpha"
      _polap_log2 "    input1: ${_outdir}/00-assembly/draft_assembly.fasta"
      local _a=($(_disassemble-step9 "${_outdir}" "${_alpha}" | tail -n 1))
      _alpha="${_a[0]}"
      _summary_draft_assembly_size="${_a[1]}"
      _polap_log2 "    output1: ${_summary_pre_alpha} -> ${_alpha}"
      _polap_log2 "    output2: draft_assembly_size (disjointigs): ${_summary_draft_assembly_size}"
    else
      _summary_draft_assembly_size="-1"
    fi

    # 10. Summary of assembly
    #
    _summary_gfa_number_segments=-1
    _summary_gfa_number_links=-1
    _summary_gfa_total_segment_length=-1
    if _polap_contains_step 10 "${_step_array[@]}"; then
      _polap_log1 "  step 10: summary statistics of gfa assembly graph"
      _polap_log2 "    input1: ${_contigger_edges_gfa}"
      _polap_log2 "    output: ${_outdir}/graph_final.txt"
      local _a=($(_disassemble-step10 "${_outdir}" | tail -n 1))
      _summary_gfa_number_segments="${_a[0]}"
      _summary_gfa_number_links="${_a[1]}"
      _summary_gfa_total_segment_length="${_a[2]}"
      _polap_log2 "    output1: ${_summary_gfa_number_segments}"
      _polap_log2 "    output2: ${_summary_gfa_number_links}"
      _polap_log2 "    output3: ${_summary_gfa_total_segment_length}"
    fi

    # 11. Annotate
    #
    if _polap_contains_step 11 "${_step_array[@]}"; then
      _polap_log1 "  step 11: annotate to find the contig with more plastid genes than mitochondrial"
      _polap_log2 "    input1: ${_contigger_edges_gfa}"
      _polap_log2 "    output: ${_ga_annotation_all}"
      _disassemble-step11 "${_outdir}"
    fi

    # 12. seeds
    # create mt.contig.name.txt
    if _polap_contains_step 12 "${_step_array[@]}"; then
      _polap_log1 "  step 12: seeds"
      _polap_log2 "    input1: ${_contigger_edges_gfa}"
      _polap_log2 "    output: ${_mtcontigname}"
      _disassemble-step12 "${_outdir}"
    fi

    # 13. Extract ptDNA sequences
    #
    # locate the circular paths
    # extract candidate circular ptDNA sequences
    _summary_num_circular_paths=-1
    _summary_num_circular_nodes=-1
    if _polap_contains_step 13 "${_step_array[@]}"; then
      _polap_log1 "  step 13: extract ptDNA sequences from the gfa assembly graph"
      _polap_log2 "    input1: ${_contigger_edges_gfa}"
      _polap_log2 "    input2: ${_mtcontigname}"
      _polap_log2 "    outdir: ${_var_mtdna}"
      local _a=($(_disassemble-step13 "${_outdir}" | tail -n 1))
      _summary_num_circular_paths="${_a[0]}"
      _summary_num_circular_nodes="${_a[1]}"
      _polap_log2 "    output1: ${_summary_num_circular_paths}"
      _polap_log2 "    output2: ${_summary_num_circular_nodes}"
    fi

    # 14. Pick one sequence based on the reference, if any.
    # Use the first sequence if there is no such reference.
    #
    # for each ptDNA of the potential ptDNA sequences,
    #   select one candidate ptDNA
    #   rearrange it so that we could do a pairwise sequence alignment
    _summary_coverage_ref=-1
    _summary_coverage_target=-1
    _summary_j_candidate=-1
    _circular_path_fasta="${_var_mtdna}/circular_path_1_concatenated.fa"
    _arg_unpolished_fasta="${_var_mtdna}/ptdna.0.fa"
    _arg_final_assembly="${_var_mtdna}/ptdna.1.fa"
    if _polap_contains_step 14 "${_step_array[@]}"; then
      _polap_log1 "  step 14: choose one of the multiple candidate ptDNA sequences"
      _polap_log2 "    input1: ${_circular_path_fasta}"

      # put in a bash function
      # INFO: tail -n 1 is necessary because we need to take the last
      # return value not something else that might have been printed out
      # in the function.
      local _a=($(_disassemble-step14 "${_var_mtdna}" | tail -n 1))
      _summary_coverage_ref="${_a[0]}"
      _summary_coverage_target="${_a[1]}"
      _summary_j_candidate="${_a[2]}"
      _polap_log2 "    output1: ${_summary_coverage_ref}"
      _polap_log2 "    output2: ${_summary_coverage_target}"
      _polap_log2 "    output3: ${_summary_j_candidate}"
      _polap_log2 "    output4: ${_arg_unpolished_fasta}"
      _polap_log2 "    output5: ${_arg_final_assembly}"

    fi

    # ptDNA BLAST global alignment
    # pident
    _summary_pident=-1
    if _polap_contains_step 15 "${_step_array[@]}"; then
      _polap_log1 "  step 15: percent identity comparison"
      _polap_log2 "    output: ${_var_mtdna}/b/pident.txt"
      _polap_log2 "    output: ${_var_mtdna}/b/blast_results_2.txt"
      _summary_pident=$(_disassemble-step15 "${_var_mtdna}" | tail -n 1)
    fi

    # ptDNA summary
    #
    _nseq=-1
    _summary_length=-1
    _summary_gc=-1
    if _polap_contains_step 16 "${_step_array[@]}"; then
      _polap_log1 "  step 16: ptDNA summary"
      _polap_log2 "    input1: ${_var_mtdna}/ptdna.0.fa"
      _polap_log2 "    input2: ${_var_mtdna}/ptdna.1.fa"

      # put in a bash function
      local _a=($(_disassemble-step16 "${_var_mtdna}" | tail -n 1))
      _nseq="${_a[0]}"
      _summary_gc="${_a[1]}"
      _summary_length="${_a[2]}"
      _polap_log2 "    output1: ${_nseq}"
      _polap_log2 "    output2: ${_summary_gc}"
      _polap_log2 "    output3: ${_summary_length}"
    fi

    _summary_elapsed_time="$(_polap_get_elapsed_time)"

    _polap_lib_table-cflye-main-summary-content "${_summary1_file}"

    # Progress
    # Calculate elapsed time and remaining time
    #
    # stage 1: weighting the total iterations
    # stage 2: not weight on the total iterations
    local status=$(_polap_lib_timing-step "${i}" "${_arg_disassemble_n}" "${_arg1}" "stage $_arg1")

    # Display the progress and remaining time on the same line if no verbose
    if [ "${_arg_verbose}" -eq "1" ]; then
      printf "\r%-${_terminal_width}s" " " >&3
      _polap_log0_ne "\r$status, elapsed time: $(_polap_get_elapsed_time ${_start_time})"
    fi

    # determine the sample size and the Flye's alpha
    _polap_log1 "  end: index $i: sample size: $_sampling_datasize_bp alpha: $_summary_pre_alpha cflye peak memory (GB): ${_summary_peak_ram_size_gb}"
  done

  _polap_log0 ""

  if [[ "${_arg_test}" == "off" ]]; then
    rm -f "${_short_read1}"
    rm -f "${_long_read}"
  fi

  return 0
}

# input1: 1.fa
# input2: 2.fa
# output: output/pident.txt
function _polap_mafft-compare-two {
  local _fa1="${_final_assembly}"
  local _fa2="${_prev_final_assembly}"
  local _out1="${_outdir}"

  cat "${_fa1}" >"${_out1}/in.fa"
  cat "${_fa2}" >>"${_out1}/in.fa"

  _polap_log3_pipe "mafft \
    --auto \
    ${_out1}/in.fa \
    >${_out1}/out.mafft \
    2>${_polap_output_dest}"

  _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-mafft.R \
    --input ${_out1}/out.mafft \
    --out ${_out1}/out.txt"

  local pident_mafft=$(grep -oP '(?<=Percent Identity: )\S+' "${_out1}/out.txt")
  local pident_mafft=${pident_mafft%\%}

  echo "$pident_mafft" >"${_out1}/pident.txt"
}

function _run_polap_polish-disassemble {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  local _index_table

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  source "${_POLAPLIB_DIR}/polap-function-disassemble-seeds.sh"

  help_message=$(
    cat <<HEREDOC
plastid genome assembly polishing step-by-step via subsampling long-read data

HEREDOC
  )

  # display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
  [[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

  local _start_time=$(date +%s)
  # _polap_lib_timing-step-reset
  local i=0
  local _summary_gfa_number_segments=0
  local _summary_gfa_number_links=0
  local _summary_gfa_total_segment_length=0
  local _summary_num_circular_paths=0
  local _summary_num_circular_nodes=0
  length=0
  _summary_gc=0
  coverage=0

  local _unpolished_fasta
  local _final_assembly
  local _disassemble_dir
  local _ga_input_short_reads
  local _short_read1
  local _long_read
  local _summary1_file

  local _ga_outdir="${_arg_outdir}/${_arg_inum}"
  _unpolished_fasta="${_arg_unpolished_fasta}"
  _ga_input_short_reads="${_ga_outdir}/s.fq"
  _disassemble_dir="${_ga_outdir}/disassemble/${_arg_disassemble_i}/3"
  _index_table="${_disassemble_dir}/index.short.txt"
  _short_read1=${_disassemble_dir}/s.fq
  _long_read=${_disassemble_dir}/l.fq.gz
  _summary1_file="${_disassemble_dir}/summary1.txt"
  check_file_existence "${_ga_outdir}/long_total_length.txt"
  check_file_existence "${_ga_outdir}/short_total_length.txt"
  local _ga_long_total_length=$(<"${_ga_outdir}/long_total_length.txt")
  local _ga_short_total_length=$(<"${_ga_outdir}/short_total_length.txt")

  _polap_log1 "  polish a plastid genome without reference"

  # select steps
  if [[ -z "${_arg_steps_include}" ]]; then
    _arg_steps_include="4,6,7,8"
    _arg_steps_exclude=""
  fi
  local _include="${_arg_steps_include}"
  local _exclude="${_arg_steps_exclude}" # optional range or list of steps to exclude
  _step_array=()

  _step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))
  _rstatus="$?"
  if [[ "${_rstatus}" -ne 0 ]]; then
    _polap_log0 "ERROR: error parsing steps."
    return "${_POLAP_ERR_CMD_OPTION_STEPS}"
  fi

  # 4. loop index
  local _disassemble_b
  local total_size_data
  if _polap_contains_step 4 "${_step_array[@]}"; then
    _polap_log1 "  step 4: prepare for the loop iteration for subsample polishing"
    _polap_log3_cmd mkdir -p "${_disassemble_dir}"

    # total_size_data=$(<"${_polap_var_outdir_short_total_length}")
    _disassemble_b=$(_polap_utility_compute_percentage \
      "${_arg_disassemble_q}" \
      "${_ga_short_total_length}")

    _polap_lib_array-make-index \
      "${_index_table}" \
      "${_arg_disassemble_r}" \
      "${_arg_disassemble_a}" \
      "${_disassemble_b}"

    # read the 2nd column or sampling sizes into an array
    # local _array_index_table
    # mapfile -t _array_index_table < <(cut -f2 "${_index_table}")
    # local _total_iterations="${#_array_index_table[@]}"
    # _polap_lib_timing-reset "${_total_iterations}"

    _polap_lib_table-cflye-polish-summary-header "${_summary1_file}"
    _is_stop="off"
  fi

  _i_last=0
  _terminal_width=$(tput cols) # Get terminal width
  cat <"${_index_table}" | while IFS=$'\t' read -r i sampling_datasize genomesize; do
    local _summary_sampling_datasize="${sampling_datasize}"
    local _summary_i="${i}"
    if [[ "${_is_stop}" == "on" ]]; then
      break
    fi

    if ((i < _arg_start_index)); then
      continue
    fi

    if [[ -n "${_arg_disassemble_s_max}" ]]; then
      if [[ "${sampling_datasize}" -gt "${_arg_disassemble_s_max}" ]]; then
        break
      fi
    fi

    _polap_lib_timing-set
    _i_last=$((_i_last + 1))
    _outdir=${_disassemble_dir}/$i
    _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
    _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
    _mtcontigname="${_outdir}/mt.contig.name"
    _ga_annotation_all="${_outdir}/assembly_info_organelle_annotation_count-all.txt"
    _var_mtdna="${_outdir}/52-mtdna"
    _ptdna="${_var_mtdna}/circular_path_1_concatenated.fa"

    _sampling_datasize_bp=$(_polap_utility_convert_bp "${sampling_datasize}")
    _genomesize_bp=$(_polap_utility_convert_bp "${genomesize}")
    _polap_log2 "  i=${i}:S=${_sampling_datasize_bp}"
    _polap_log3_cmd mkdir -p "${_outdir}"

    # begin: steps

    # 6. Sample short-read data
    # We utilize the identical sampling rate employed
    # by the long-read sampling method, which is set to default.
    #
    _summary_short_total=-1
    _summary_short_rate_sample=-1
    _summary_short_sample_seed=-1
    if _polap_contains_step 6 "${_step_array[@]}"; then
      check_file_existence "${_ga_input_short_reads}"
      _polap_log1 "  step 6: subsample the short-read data"
      # _summary_short_total=$(<"${_polap_var_outdir_short_total_length}")
      _polap_log1 "    use sample size: ${sampling_datasize} bp"
      _summary_short_rate_sample=$(echo "scale=10; ${sampling_datasize}/${_ga_short_total_length}" | bc)
      _polap_lib_random-get
      _summary_short_sample_seed=${_polap_var_random_number}
      _polap_log2 "    input1: sampling rate: ${_summary_short_rate_sample}"
      _polap_log2 "    input2: random seed: ${_summary_short_sample_seed}"
      _polap_log2 "    input3: short read data: ${_ga_input_short_reads}"
      _polap_log2 "    output: ${_short_read1}"
      _polap_log3_pipe "seqkit sample \
          -p ${_summary_short_rate_sample} \
          -s ${_summary_short_sample_seed} \
          ${_ga_input_short_reads} \
          -o ${_short_read1} 2>${_polap_output_dest}"
    fi

    local _summary_memory_gb_msbwt
    local _summary_total_hours_msbwt
    local _msbwt_dir="${_disassemble_dir}/msbwt"
    local _msbwt="${_disassemble_dir}/msbwt/comp_msbwt.npy"
    if _polap_contains_step 7 "${_step_array[@]}"; then
      _polap_log1 "  step 7: creating ${_msbwt_dir} ..."

      rm -rf "${_msbwt_dir}"
      command time -v bash "${_POLAPLIB_DIR}"/polap-build-msbwt.sh \
        "${_short_read1}" \
        "${_msbwt_dir}" \
        2>"${_outdir}/timing-msbwt.txt"

      # Time and memory
      read -r _summary_memory_gb_msbwt _summary_total_hours_msbwt < <(_polap_lib_timing-parse-timing "${_outdir}/timing-msbwt.txt")
      _polap_log2 "    msbtw memory: ${_summary_memory_gb_msbwt}"
      _polap_log2 "    msbtw time: ${_summary_total_hours_msbwt}"
    fi

    local _summary_memory_gb_fmlrc
    local _summary_total_hours_fmlrc
    local _summary_assembly_length
    if _polap_contains_step 8 "${_step_array[@]}"; then
      _polap_log1 "  step 8: polishing ${_unpolished_fasta}"
      # Initialize Conda
      source $HOME/miniconda3/etc/profile.d/conda.sh
      conda activate polap-fmlrc

      _final_assembly="${_outdir}/ptdna.1.fa"
      _polap_log2 "    input1: ${_unpolished_fasta}"
      _polap_log2 "    input2: ${_msbwt_dir}"
      _polap_log2 "    output1: ${_final_assembly}"
      command time -v fmlrc -p "${_arg_threads}" \
        "${_msbwt}" \
        "${_unpolished_fasta}" \
        "${_final_assembly}" \
        2>"${_outdir}/timing-fmlrc.txt"
      # 2>&1

      conda deactivate

      # Time and memory
      read -r _summary_memory_gb_fmlrc _summary_total_hours_fmlrc < <(_polap_lib_timing-parse-timing "${_outdir}/timing-fmlrc.txt")
      _summary_assembly_length=$(grep -v "^>" "${_final_assembly}" | tr -d '\n' | wc -c)
      _polap_log2 "    fmlrc memory: ${_summary_memory_gb_fmlrc}"
      _polap_log2 "    fmlrc time: ${_summary_total_hours_fmlrc}"
      _polap_log2 "    length: ${_summary_assembly_length}"
    fi

    # stop if the _summary_memory_gb_msbwt is greater than the max memory
    # stop if the _summary_memory_gb_fmlrc is greater than the max memory
    if (($(echo "$_summary_memory_gb_msbwt  < $_arg_disassemble_memory" | bc -l))); then
      _polap_log2 "  cFlye: msbwt used memory is less than ${_arg_disassemble_memory} Gb ($i): ${_summary_memory_gb_msbwt}"
    else
      _polap_log2 "  cFlye: msbwt used memory is not less than ${_arg_disassemble_memory} Gb ($i): ${_summary_memory_gb_msbwt}"
      _polap_log2 "  exit the polish disassemble menu."
      _is_stop="on"
    fi

    if (($(echo "$_summary_memory_gb_fmlrc  < $_arg_disassemble_memory" | bc -l))); then
      _polap_log2 "  cFlye: fmlrc used memory is less than ${_arg_disassemble_memory} Gb ($i): ${_summary_memory_gb_fmlrc}"
    else
      _polap_log2 "  cFlye: fmlrc used memory is not less than ${_arg_disassemble_memory} Gb ($i): ${_summary_memory_gb_fmlrc}"
      _polap_log2 "  exit the polish disassemble menu."
      _is_stop="on"
    fi

    rm -rf "${_msbwt_dir}"

    # compare one and its previous sequences
    if ((_summary_i > 0)); then
      local _prev_summary_i=$((_summary_i - 1))
      local _prev_outdir=${_disassemble_dir}/${_prev_summary_i}
      local _prev_final_assembly="${_prev_outdir}/ptdna.1.fa"
      # _polap_log0 "${_final_assembly}"
      # _polap_log0 "${_prev_final_assembly}"

      _polap_mafft-compare-two \
        "${_final_assembly}" \
        "${_prev_final_assembly}" \
        "${_outdir}"

      _summary_assembly_pident=$(<"${_outdir}/pident.txt")
    else
      _summary_assembly_pident="-1"
    fi
    # _polap_log0 "${_summary_assembly_pident}"

    _polap_lib_table-cflye-polish-summary-content "${_summary1_file}"

    # end: steps

    # Progress
    # Calculate elapsed time and remaining time
    #
    local status=$(_polap_lib_timing-step "${i}" "${_arg_disassemble_r}" 2 "stage 3")

    # Display the progress and remaining time on the same line
    if [ "${_arg_verbose}" -eq "1" ]; then
      printf "\r%-${_terminal_width}s" " " >&3
      _polap_log0_ne "\r$status, elapsed time: $(_polap_get_elapsed_time ${_start_time})"
    fi

    # determine the sample size and the Flye's alpha
    _polap_log1 "  end: index $i: sample size: $_sampling_datasize_bp memory:"
  done

  _polap_log0 ""

  if [[ "${_arg_test}" == "off" ]]; then
    rm -f "${_short_read1}"
    rm -f "${_long_read}"
  fi

  return 0
}
