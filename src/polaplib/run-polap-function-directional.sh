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

################################################################################
# The idea is as follows. One intriguing yet somewhat challenging aspect is the use of directional reads in genome assembly. This approach could be useful when the results of plant mitochondrial genome assemblies are difficult to interpret.
#
# Let us consider repetitive sequences, specifically direct repeats and inverted repeats. In an assembly graph, multiple sequences are linked by edges. Following a path along these edges corresponds to interpreting that path as the assembled genome. If the assembly graph is simple, identifying this path is straightforward. A plastid genome (chloroplast genome), for instance, may form a simplest circular structure if a single sequence is connected end to end. More commonly, however, there are three sequences connected by four edges, composed of the large subunit, small subunit, and an inverted repeat. In the case of an inverted repeat, if we regard it as two strands (similar to a DNA double helix), it can essentially be treated as a single sequence. In reality, chloroplast genomes often consist of four sequences arranged in a circular form.
#
# However, if an assembled mitochondrial genome is made up of multiple sequences linked in a complex arrangement, tracing a path through the graph—like one would do with a simpler chloroplast genome—may not be easy. One possible approach is to search for the simplest path in the assembly graph, but confirming whether such a path actually exists is not trivial. Instead, we aim to investigate whether we can simplify the genome assembly graph by using that plausible path as seed contigs.
#
# The aforementioned direct and inverted repeats are included in the assembled genome if the assembly path traverses them more than once. However, direct repeats lie in the same orientation along the path, while inverted repeats are in the opposite orientation. If we only follow the completed path in one direction, direct repeats remain repetitive sequences distinguishable by their flanking regions, whereas inverted repeats, despite being repetitive, differ entirely in orientation. This implies that in a unidirectional assembly process, direct repeats may still be indistinguishable, while inverted repeats might be resolved.
#
# The issue is that raw nucleotide sequences are typically generated without preserving directional information, making it rare to find genome assembly tools that handle directional reads. One exception is directional RNA sequencing, which does preserve orientation and thus factors it into transcript assembly. However, because standard genome assembly does not generate directionally specific sequences, such tools have generally not been necessary. For nuclear genomes, incorporating directional information into assembly likely offers little benefit. But for smaller genomes, such as plant mitochondria or possibly bacterial genomes—where assembly can be complicated by repetitive sequences—taking directionality into account could be advantageous in certain cases. In essence, direct repeats would still remain in the graph, but repeats like inverted repeats might be resolved, simplifying the assembly graph overall.
#
# There are two specific challenges: (1) generating directional read data, and (2) developing or adapting a tool capable of performing genome assembly using these directional reads. We generate directional read data by selecting sequences that map in the same orientation to seed contigs. Among available genome assembly tools, we modified Flye to utilize this directional data in the assembly process.
#
# As a proof of concept, we tested our method on a plant mitochondrial genome that was otherwise difficult to interpret, and found that it produced a somewhat simpler assembly graph. Going forward, it seems promising to apply this method to bacterial genome assembly as well.
################################################################################
function _run_polap_directional {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

  help_message=$(
    cat <<HEREDOC
Plant mitochondrial genome assembly using directional long-read data

1. Reference generating or whole-genome assembly or subsampling and WGA
2. Select directional seed contigs: e.g., edge_1+,edge_2- not just edge_1,edge_2
3. Map reads with only one directional
4. Assemble the directional reads using dflye.
5. Test the data

Inputs
------

- long-read data: ${_arg_long_reads}
- short-read data: ${_arg_short_read1}

Outputs
-------

whole-genome assembly or subsampled whole-genome assembly
let's first select reads

Arguments
---------
-o ${_arg_outdir}
--directional-i <name>: name the directional
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file

-m ${_arg_min_read_length}: the long-read sequence length threshold
-t ${_arg_threads}: the number of CPU cores
-c ${_arg_coverage}: the coverage option
--random-seed <arg>: 5-digit number

view
----

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1}
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
  [[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

  # Display the content of output files
  if [[ "${_arg_menu[1]}" == "view" ]]; then

    return 0
  fi

  _polap_log0 "starting the directional genome assembly on ${_arg_outdir} ..."

  #

  # funcition: like assemble2
  #
  # seeds: prepare seeds
  # map: map reads
  # reads: select reads and test reads
  # flye3 or dflye

  # _run_polap_directional-prepare-seeds
  # _run_polap_directional-map-reads
  # _run_polap_directional-select-reads

  # INFO:
  # _polap_directional-dflye -> call this dga
  # disassemble combined wga/oga and assemble modules
  # the v0.2 separtes wga/oga and assemble modules
  # _run_polap_directional -> like _run_polap_assemble2
  # so, we could exectue dflye in _run_polap_dflye in dga module
  # run_polap_dflye -> like run_polap_flye2

  # _run_polap_dflye

  # local out="${1}"
  # local long_read="${2}"
  # local expected_genome_size="${3}"

  local step_array

  # select steps
  if [[ -z "${_arg_steps_include}" ]]; then
    _arg_steps_include="1-6"
    _arg_steps_exclude=""
  fi

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
  local _directional_dir="${_arg_outdir}/${_arg_jnum}/directional/${_arg_directional_i}"

  _polap_log3_cmd mkdir -p "${_directional_dir}"

  if _polap_contains_step 1 "${_step_array[@]}"; then
    _polap_directional-prepare-seeds
  fi

  if _polap_contains_step 2 "${_step_array[@]}"; then
    _polap_directional-map-reads
  fi

  # create index
  if _polap_contains_step 3 "${_step_array[@]}"; then
    _polap_directional-prepare-loop-iteration
  fi

  local _index_table="${_directional_dir}/index.txt"
  local i
  local omega
  while IFS=$'\t' read -r i omega; do
    _polap_log0 "index: ${i} - ${omega}"

    local _directional_dir_i="${_directional_dir}/${i}"
    _polap_log3_cmd mkdir -p "${_directional_dir_i}"

    if _polap_contains_step 4 "${_step_array[@]}"; then
      _polap_directional-reads "${i}" "${omega}"
    fi

    if _polap_contains_step 5 "${_step_array[@]}"; then
      _polap_directional-subsample "${i}" "${omega}"
    fi

    if _polap_contains_step 6 "${_step_array[@]}"; then
      _polap_directional-flye "${i}" "${omega}"
    fi

    if _polap_contains_step 7 "${_step_array[@]}"; then
      _polap_log0 "step 7: dflye summary"
    fi

    if _polap_contains_step 8 "${_step_array[@]}"; then
      _polap_log0 "step 8: dflye plot"
    fi

  done <"${_index_table}"

  _polap_log1 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}
