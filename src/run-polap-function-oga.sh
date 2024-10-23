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

################################################################################
# Selects reads mapped on a genome assembly.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -r ${_arg_pair_min}: minimum minimap2 alignment length for a pair of contigs
#   -x ${_arg_bridge_min}: minimum long-read length for connecting the pair of contigs
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
# Inputs:
#   ${MTCONTIGNAME}
#   ${_polap_var_contigger_edges_fasta}
# Outputs:
#   ${MTSEEDSDIR}
#   ${MTDIR}/contig.fa
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/contig.paf
#   ${MTDIR}/contig.tab
#   ${MTSEEDSDIR}/1.names
#   ${MTSEEDSDIR}/1.fq.gz
#   ${MTSEEDSDIR}/2.fq.gz
################################################################################
function _run_polap_select-reads() { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-base.sh"     # '.' means 'source'
	source "$script_dir/polap-variables-mtcontig.sh" # '.' means 'source'
	source "$script_dir/polap-variables-oga.sh"      # '.' means 'source'
	source "$script_dir/polap-variables-ga.sh"       # '.' means 'source'

	local MTDIR="${_polap_var_oga}"                              # target: ${_polap_var_oga}
	local MTSEEDSDIR="${_polap_var_oga}/seeds"                   # ${_polap_var_seeds} for oga-class
	local MTCONTIGNAME="${_polap_var_ga}/mt.contig.name-${JNUM}" # ${_polap_var_mtcontigname}

	# for contigs
	#	_polap_var_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Select long reads mapped on a genome assembly.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number (or 0)
#   -j $JNUM: destination Flye organelle assembly number
#   --rwx ${_arg_rwx}: set the option values of -r, -x, -w to the same one
# Inputs:
#   ${MTCONTIGNAME}
#   ${_polap_var_contigger_edges_fasta}
# Outputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/1.fq.gz <- single-mapped reads
#   ${MTSEEDSDIR}/2.fq.gz <- reduced single-mapped reads, if we have many such reads
#   ${MTSEEDSDIR}/3.fq.gz <- reduced single-mapped reads plus pair-mapped reads
Example: $0 ${_arg_menu[0]} -i ${INUM} -j ${JNUM} --rwx ${_arg_rwx}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "selecting long-reads mapped on the seed contigs in file: ${MTCONTIGNAME} and ${_polap_var_contigger_edges_fasta}"
	_polap_log1 "  input1: ${MTCONTIGNAME}"
	_polap_log1 "  input2: ${_polap_var_contigger_edges_fasta}"

	if [ ! -s "${MTCONTIGNAME}" ]; then
		_polap_log0 "ERROR: no such mt.contig.name file: ${MTCONTIGNAME}"
		exit $EXIT_ERROR
	fi

	if [ ! -s "${_polap_var_contigger_edges_fasta}" ]; then
		_polap_log0 "ERROR: no assembly fasta file: ${_polap_var_contigger_edges_fasta}"
		exit $EXIT_ERROR
	fi

	# if [ -d "${MTDIR}" ] && [ "${_arg_yes}" = "off" ]; then
	if [ -d "${MTDIR}" ] && [ "${_arg_redo}" = "off" ]; then
		while true; do
			read -r -p "Folder [${MTDIR}] already exists. Do you want to replace it? [y/n] " yn
			case $yn in
			[Yy]*)
				_polap_log3_cmd rm -rf "${MTDIR}"
				_polap_log0 "  folder ${MTDIR} is deleted."
				break
				;;
			[Nn]*)
				_polap_log0 "  folder ${MTDIR} is not replaced."
				_polap_log0 "  you might want a new mt.contig.name file for flye2 step."
				exit $EXIT_FAIL
				;;
			*) _polap_log0 "Please answer yes or no." ;;
			esac
		done
	else
		_polap_log3_cmd rm -rf "${MTDIR}"
		_polap_log1 "  ${MTDIR} is deleted if there is one."
	fi

	_polap_log1 "  creating ${ODIR}/organelle-assembly_${INUM}-${JNUM}"
	echo "$CMD" >"${ODIR}/organelle-assembly_${INUM}-${JNUM}"

	_polap_log1 "  creates ${MTSEEDSDIR}"
	_polap_log3_cmd mkdir -p "${MTSEEDSDIR}"
	_polap_log3_cmd ln -s "$PWD/${_polap_var_base_nk_fq_gz}" -t "${MTDIR}"
	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input2: ${MTCONTIGNAME}"
	_polap_log2 "    output: ${MTDIR}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    --threads ${_arg_threads} \
    -f ${MTCONTIGNAME} \
		${_polap_var_contigger_edges_fasta} \
		-o ${MTDIR}/contig.fa \
		2>${_polap_output_dest}"

	# we could circularize a single contig.
	local contig_count=$(wc -l <"${MTCONTIGNAME}")
	_polap_log1 "  number of seed contigs: ${contig_count}"
	if [[ ${_arg_circularize} == "on" ]]; then
		if [ "$contig_count" -eq 1 ]; then

			_polap_log1 "  circularizes the single seed contig ..."
			_polap_log2 "    input1: ${MTDIR}"
			_polap_log2 "    input2: ${MTDIR}/contig.fa"
			_polap_log2 "    input3: ${MTCONTIGNAME}"
			_polap_log2 "    output1: ${MTCONTIGNAME}-backup"
			_polap_log2 "    output2: (new) ${MTCONTIGNAME}"
			_polap_log2 "    output3: ${MTDIR}/contig.fa-backup"
			_polap_log2 "    output4: (new) ${MTDIR}/contig.fa"
			_polap_log3_cmd bash "$script_dir/run-polap-sh-half-cut.sh" "${MTDIR}" "${MTDIR}/contig.fa" "${MTCONTIGNAME}"

			# cp "${MTDIR}/contig.fa" "${MTDIR}/contig.fa-backup"
			# seqkit fx2tab --length --name "${MTDIR}"/contig.fa -o "${MTDIR}"/contig.fa.len >/dev/null 2>&1
			# A=$(cut -f2 "${MTDIR}"/contig.fa.len)
			# B=$(echo "scale=0; $A/2" | bc)
			# C=$((B + 1))
			# seqkit subseq -r 1:"$B" "${MTDIR}"/contig.fa -o "${MTDIR}"/c1.fa >/dev/null 2>&1
			# seqkit subseq -r "$C":"$A" "${MTDIR}"/contig.fa -o "${MTDIR}"/c2.fa >/dev/null 2>&1
			# cat "${MTDIR}"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "${MTDIR}"/contig.fa >/dev/null 2>&1
			# cp "${MTCONTIGNAME}" "${MTCONTIGNAME}"-backup
			# echo -e "edge_1\nedge_2" >"${MTCONTIGNAME}"

			_polap_log1 "    creating new ${MTDIR}/contig.fa and ${MTCONTIGNAME}"

		else
			_polap_log0 "Not implemented yet!"
			exit $EXIT_ERROR
			# "$script_dir"/run-polap-single.R "${MTSEEDSDIR}"/contig.tab "${MTSEEDSDIR}" "${_arg_single_min}" >/dev/null 2>&1
			# cat "${MTSEEDSDIR}"/single.names | sort | uniq >"${MTSEEDSDIR}"/1.names
			# echo "INFO: creates long read single name in ${MTSEEDSDIR}/1.names"
		fi
	fi

	_polap_log1 "  please, wait for a long-read data selection ..."
	_polap_log1 "    $INUM (source) -> $JNUM (target) ..."
	_polap_log1 "    bridge=${_arg_bridge_min}"
	_polap_log1 "    p_mapping=${_arg_pair_min}"
	_polap_log1 "    s_mapping=${_arg_single_min}"
	_polap_log1 "    min_len_read=${_arg_min_read_length}"

	CONTIG_LENGTH=$(seqkit stats -Ta "${MTDIR}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	echo "$CONTIG_LENGTH" >"${MTDIR}"/contig_total_length.txt
	_polap_log1 "  organelle genome size based on the seed contig selection: $CONTIG_LENGTH"

	_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
	_polap_log2 "    input1: ${MTDIR}/contig.fa"
	_polap_log2 "    input2: ${_polap_var_base_nk_fq_gz}"
	_polap_log2 "    output: ${MTDIR}/contig.paf"
	if [[ -s "${MTDIR}"/contig.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
		_polap_log1 "  found: ${MTDIR}/contig.paf, skipping the minimap2 mapping step ..."
	else
		_polap_log3_pipe "minimap2 -cx map-ont \
      ${MTDIR}/contig.fa \
      ${_polap_var_base_nk_fq_gz} \
      -t ${_arg_threads} \
      -o ${MTDIR}/contig.paf \
      >${_polap_output_dest} 2>&1"
	fi

	_polap_log1 "  filtering out reads less than ${_arg_min_read_length} bp and converting PAF to TAB ..."
	_polap_log2 "    input1: ${MTDIR}/contig.paf"
	_polap_log2 "    output: ${MTDIR}/contig.tab"
	# cut -f1-11 "${MTDIR}"/contig.paf | awk -v minlength="${_arg_min_read_length}" '{if ($2>=minlength) {print}}' >"${MTDIR}"/contig.tab
	_polap_log3_cmd bash "$script_dir/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${MTDIR}/contig.paf" "${MTDIR}/contig.tab"

	# MT: MPAIR=3000 MBRIDGE=3000 MSINGLE=3000
	# PT: MPAIR=1000 MBRIDGE=5000 MSINGLE=0
	_polap_log1 "  selecting seeds that meet the selection criteria ..."
	_polap_log2 "    bridge (POLAP pairs bridge minimum)=${_arg_bridge_min}"
	_polap_log2 "    p_mapping (POLAP pairs alignment minimum)=${_arg_pair_min}"
	_polap_log2 "    s_mapping (POLAP single alignment minimum)=${_arg_single_min}"
	_polap_log2 "    min_len_read=${_arg_min_read_length}"
	_polap_log2 "    input1: ${MTCONTIGNAME}"
	_polap_log2 "    input2: ${MTDIR}/contig.tab"
	_polap_log2 "    output: ${MTSEEDSDIR}/single.names for reads mapped on a single contig"
	_polap_log2 "    output: ${MTSEEDSDIR}/<edge1>-<edge2>.name for reads mapped on contig pairs"
	_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-pairs.R \
		${MTCONTIGNAME} \
		${MTDIR}/contig.tab \
		${MTSEEDSDIR} \
		${_arg_pair_min} \
    ${_arg_bridge_min} \
    ${_arg_single_min} \
    >${_polap_output_dest} 2>&1"
	# "$script_dir"/run-polap-pairs.R "${MTCONTIGNAME}" ${MTDIR}/contig.tab ${MTSEEDSDIR} ${_arg_pair_min} ${_arg_bridge_min} ${_arg_single_min} >/dev/null 2>&1

	# cat "${MTSEEDSDIR}/"*".name" "${MTSEEDSDIR}"/single.names | sort | uniq >"${MTSEEDSDIR}"/1.names
	_polap_log1 "  creates names of reads mapped on a single name in ${MTSEEDSDIR}/1.names"
	_polap_log2 "    input1: ${MTSEEDSDIR}/single.names"
	_polap_log2 "    output: ${MTSEEDSDIR}/1.names"
	_polap_log3_pipe "cat ${MTSEEDSDIR}/single.names | sort | uniq >${MTSEEDSDIR}/1.names"

	# seqkit grep --threads ${_arg_threads} -f "${MTSEEDSDIR}"/1.names ${_polap_var_base_nk_fq_gz} -o "${MTSEEDSDIR}"/1.fq.gz >/dev/null 2>&1
	_polap_log1 "  creating reads mapped on a single contig in ${MTSEEDSDIR}/1.fq.gz"
	_polap_log2 "    input1: ${_polap_var_base_nk_fq_gz}"
	_polap_log2 "    input2: ${MTSEEDSDIR}/1.names"
	_polap_log2 "    output: ${MTSEEDSDIR}/1.fq.gz"
	seqtk subseq "${_polap_var_base_nk_fq_gz}" "${MTSEEDSDIR}"/1.names | gzip >"${MTSEEDSDIR}"/1.fq.gz

	_polap_log1 "  computing the total size of the single-mapped reads ..."
	_polap_log2 "    input1: ${MTSEEDSDIR}/1.fq.gz"
	local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
	_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
	local TOTAL_LENGTH=$(seqkit stats -Ta "${MTSEEDSDIR}"/1.fq.gz | csvtk cut -t -f "sum_len" | csvtk del-header)
	local _total_length_bp=$(_polap_utility_convert_bp ${TOTAL_LENGTH})
	_polap_log2 "    result1 (total size of reads mapped on a single contig): ${_total_length_bp}"
	local EXPECTED_ORGANELLE_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
	_polap_log2 "    result2 (expected organelle coverage): ${EXPECTED_ORGANELLE_COVERAGE}x"

	_polap_log1 "  reducing the single-mapped reads upto the coverage of ${_arg_coverage}x"
	_polap_log2 "    input1: ${MTSEEDSDIR}/1.fq.gz"
	_polap_log2 "    output: ${MTSEEDSDIR}/2.fq.gz"
	if [[ "${_arg_test}" == "on" ]]; then
		_polap_log0 "    OPTION: --test : No reduction of the test long-read data"
		_polap_log3_cmd ln -s "$(realpath "${MTSEEDSDIR}"/1.fq.gz)" "${MTSEEDSDIR}"/2.fq.gz
	elif [[ "${_arg_coverage_check}" == "off" ]]; then
		_polap_log0 "    OPTION: --no-coverage-check : No reduction of the long-read data"
		_polap_log3_cmd ln -s "$(realpath "${MTSEEDSDIR}"/1.fq.gz)" "${MTSEEDSDIR}"/2.fq.gz
	else
		if [ "$EXPECTED_ORGANELLE_COVERAGE" -lt "${_arg_coverage}" ]; then
			_polap_log1 "    no reduction of the long-read data because $EXPECTED_ORGANELLE_COVERAGE < ${_arg_coverage}"
			_polap_log3_cmd ln -s "$(realpath "${MTSEEDSDIR}"/1.fq.gz)" "${MTSEEDSDIR}"/2.fq.gz
		else
			_polap_log1 "    SUGGESTION: you might want to increase the minimum read lengths (use --rwx or -m) because you have enough long-read data."
			RATE=$(echo "scale=10; ${_arg_coverage}/$EXPECTED_ORGANELLE_COVERAGE" | bc)
			_polap_log0 "  long-read data reduction by rate of $RATE <= COV[${_arg_coverage}] / long-read organelle coverage[$EXPECTED_ORGANELLE_COVERAGE]"
			_polap_log1 "    sampling long-read data by $RATE ... wait ..."
			# seqkit sample -p "$RATE" "${MTSEEDSDIR}/1.fq.gz" -o "${MTSEEDSDIR}/2.fq.gz" >/dev/null 2>&1
			local seed=${_arg_random_seed:-$RANDOM}
			_polap_log1 "    random seed for reducing the single-mapped long-read data: ${seed}"
			_polap_log3_pipe "seqkit sample \
        -p ${RATE} \
        -s ${seed} \
        ${MTSEEDSDIR}/1.fq.gz \
        -o ${MTSEEDSDIR}/2.fq.gz 2>${_polap_output_dest}"
		fi
	fi

	# 1.fq.gz: reads mapped on a single read
	# 2.fq.gz: a reduced data of 1.fq.gz
	#          or a reduced data of 1.fq.gz plus the inter-mapped reads, if any
	# single.names.2: read names from the reduced data
	_polap_log1 "  adding pair-mapped bridging long reads to the single-mapped data ..."
	local C=$(ls -1 "${MTSEEDSDIR}/"*".name" 2>/dev/null | wc -l)
	if [ "$C" != 0 ]; then
		_polap_log2 "    input1: combinations of $C bridging reads"
		_polap_log1 "  creating a list of reduced reads mapped on a single contig"
		_polap_log2 "    input1: ${MTSEEDSDIR}/2.fq.gz"
		_polap_log2 "    output: ${MTSEEDSDIR}/single.names.2"
		_polap_log3_cmd seqkit seq -n -i "${MTSEEDSDIR}"/2.fq.gz -o "${MTSEEDSDIR}"/single.names.2

		_polap_log3 "  FIXME: _polap_log3_cmd or _polap_log3_pipe for cat *"
		cat "${MTSEEDSDIR}/"*".name" "${MTSEEDSDIR}"/single.names.2 | sort | uniq >"${MTSEEDSDIR}/1.names.2"
		# seqkit grep --threads ${_arg_threads} -f "${MTSEEDSDIR}"/1.names.2 ${_polap_var_base_nk_fq_gz} -o "${MTSEEDSDIR}"/2.fq.gz >/dev/null 2>&1

		_polap_log1 "  extracting single- and pair-mapped reads ..."
		_polap_log2 "    input1: ${_polap_var_base_nk_fq_gz}"
		_polap_log2 "    input1: ${MTSEEDSDIR}/1.names.2"
		_polap_log2 "    output: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log3_pipe "seqtk subseq ${_polap_var_base_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz"
	else
		_polap_log1 "    no pair-mapped reads: 2.fq.gz -> 3.fq.gz"
		_polap_log3_cmd ln -s "$(realpath "${MTSEEDSDIR}"/2.fq.gz)" "${MTSEEDSDIR}"/3.fq.gz
	fi
	_polap_log1 "  organelle reads in ${MTSEEDSDIR}/3.fq.gz"

	# put the backup to the original
	if [[ ${_arg_circularize} == "on" ]]; then
		_polap_log1 "  circularize option: putting seed contigs back to the original mt.contig.name and ${MTDIR}/contig.fa"
		if [[ -s "${MTCONTIGNAME}"-backup ]]; then
			mv "${MTCONTIGNAME}"-backup "${MTCONTIGNAME}"
			mv "${MTDIR}/contig.fa-backup" "${MTDIR}/contig.fa"
		else
			echo "DEV: not implemented yet"
			exit $EXIT_ERROR
		fi
	fi

	_polap_log1 NEXT: $0 flye2 -o "$ODIR" -j "$JNUM" -t "${_arg_threads}" -c "${_arg_coverage}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Executes Flye for an organelle-genome assembly
#
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
# Outputs:
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/30-contigger/contigs.fasta
#   ${MTDIR}R}/30-contigger/contigs_stats.txt
#   ${MTDIR}R}/30-contigger/graph_final.fasta
#   ${MTDIR}/30-contigger/graph_final.gfa
################################################################################
function _run_polap_flye2() { # executes Flye for an organelle-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-oga.sh"

	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="${MTDIR}"/seeds

	help_message=$(
		cat <<HEREDOC
# Executes Flye for an organelle-genome assembly
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
# Outputs:
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/30-contigger/contigs.fasta
#   ${MTDIR}/30-contigger/contigs_stats.txt
#   ${MTDIR}/30-contigger/graph_final.fasta
#   ${MTDIR}/30-contigger/graph_final.gfa
Example: $0 ${_arg_menu[0]} -j <arg>
Example: $0 ${_arg_menu[0]} -j <arg> polishing
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_polap_log0 "flye2 assembly with polishing on the reads mapped on seed contigs ..."
	else
		_polap_log0 "flye2 assembly on the reads mapped on seed contigs ..."
	fi
	_polap_log1 "  input1: ${MTDIR}/contig.fa"
	_polap_log1 "  input2: ${MTSEEDSDIR}/3.fq.gz"

	if [ ! -s "${MTDIR}/contig.fa" ]; then
		echoall "ERROR: no selected-contig file [${MTDIR}/contig.fa]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "${MTSEEDSDIR}/3.fq.gz" ]; then
		echoall "ERROR: no long-read file: ${MTSEEDSDIR}/3.fq.gz"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	CONTIG_LENGTH=$(seqkit stats -Ta "${MTDIR}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	echo "$CONTIG_LENGTH" >"${MTDIR}"/contig_total_length.txt
	_polap_log1 "  organelle genome size based on the seed contig selection: $CONTIG_LENGTH"

	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_polap_log1 "  executing the organelle-genome long-read flye polishing assembly on $JNUM ..."
		_polap_log1 "    input1: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log1 "    output: ${MTDIR}/assembly_graph.gfa"
	else
		_polap_log1 "  executing the organelle-genome assembly using flye on $JNUM ..."
		_polap_log1 "    input1: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log1 "    output: ${MTDIR}/30-contigger/graph_final.gfa"
	fi
	local _command1="flye \
    --nano-raw ${MTSEEDSDIR}/3.fq.gz \
		--out-dir ${MTDIR} \
		--threads ${_arg_threads} \
		--asm-coverage ${_arg_flye_asm_coverage} \
		--genome-size $CONTIG_LENGTH"
	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_command1+=" \
		  --resume"
	else
		_command1+=" \
		  --stop-after contigger"
	fi
	_command1+=" \
		2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_polap_log0 "  output: the long-read polished assembly graph: $PWD/${_polap_var_oga_assembly_graph_gfa}"
		_polap_log0 "  DO: extract a draft organelle genome sequence (mt.0.fasta)"
		_polap_log0 "      from the polished assembly graph before short-read polishing"
		_polap_log1 "NEXT: $0 prepare-polishing -a ${_arg_short_read1} -b ${_arg_short_read2}"
	else
		_polap_log0 "  output: the assembly fasta: ${_polap_var_oga_contigger_edges_fasta}"
		_polap_log0 "  output: the assembly graph: ${_polap_var_oga_contigger_edges_gfa}"
		jnum_next=$((JNUM + 1))
		_polap_log1 "  create and edit $ODIR/$JNUM/mt.contig.name-${jnum_next}"
		_polap_log1 "NEXT: $0 assemble2 -o ${ODIR} -i ${JNUM} -j ${jnum_next}"
		_polap_log1 or you could finish with Flye organelle-genome assembly with its polishing stage.
		_polap_log1 "NEXT: $0 flye-polishing -o ${ODIR} -j $${JNUM}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Finishes a Flye organelle-genome assembly upto flye-polishing step
#
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
#   ${MTDIR}/30-contigger
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
################################################################################
function _run_polap_flye-polishing() { # finish a Flye organelle-genome assembly upto flye-polishing step
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="${MTDIR}"/seeds

	help_message=$(
		cat <<HEREDOC
# Finishes the Flye organelle-genome assembly.
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/3.fq.gz
#   ${MTDIR}/30-contigger
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
Example: $0 ${_arg_menu[0]} -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && return

	if [ ! -s "${MTDIR}/contig.fa" ]; then
		echoall "ERROR: no selected-contig file [${MTDIR}/contig.fa]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "${MTSEEDSDIR}/3.fq.gz" ]; then
		echoall "ERROR: no long-read file [${MTSEEDSDIR}/3.fq.gz]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	_arg_menu[1]="polishing"
	_run_polap_flye2

	# CONTIG_LENGTH=$(seqkit stats -Ta "${MTDIR}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	# _polap_log1 "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"
	#
	# _polap_log1 "INFO: polishing the organelle-genome assembly using flye ... be patient!"
	# _polap_log0 "please, wait for Flye long-read polishing of the organelle-genome assembly on $JNUM ..."
	# _polap_log3_pipe "flye --nano-raw ${MTSEEDSDIR}/3.fq.gz \
	# 	--out-dir ${MTDIR} \
	# 	--threads ${_arg_threads} \
	# 	--asm-coverage ${_arg_flye_asm_coverage} \
	# 	--genome-size $CONTIG_LENGTH \
	# 	--resume \
	# 	2>$_polap_output_dest"
	#
	# _polap_log0 "  output: the long-read polished assembly graph: $PWD/${MTDIR}/assembly_graph.gfa"
	# _polap_log0 "  DO: extract a draft organelle genome sequence (mt.0.fasta) from the polished assembly graph"

	# echo "column -t $ODIR/assembly_info_organelle_annotation_count.txt"
	# echoall NEXT: $(basename $0) check-coverage [-p ${_arg_unpolished_fasta}]
	# _polap_log1 NEXT: "$(basename "$0")" prepare-polishing -a "${_arg_short_read1}" -b "${_arg_short_read2}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Reports the organelle-genome assembly results.
################################################################################
function _run_polap_report-assembly() { # report an organelle-genome assembly result
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-wga.sh" # '.' means 'source'

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Reports the organelle-genome assembly results.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject
# Inputs:
#   ${ODIR}: output folder for BioProject
# Outputs:
Example: $(basename $0) ${_arg_menu[0]} [-o ${ODIR}]
Example: report-assembly -o PRJDB10540a 2>&1 | tr '\n' '\t' | sed 's/\t$/\n/'
HEREDOC
	)

	# Set variables for file paths
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'
	source "$script_dir/polap-variables-base.sh"       # '.' means 'source'

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_log0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log1_log "reporting the organelle-genome assembly at ${ODIR} ..."

	if [ -d "${ODIR}" ]; then
		_polap_log2_file "the main output folder: ${ODIR}"
	else
		_polap_log2 "ERROR: no such output folder; use -o option"
		exit $EXIT_SUCCESS
	fi

	_polap_log0 $(cut -f1 "${_polap_var_bioproject_txt}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_sra_long_read}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_sra_short_read}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_species}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")

	# wc -l "${_polap_var_wga}/mt.contig.name-"? | awk '$2 != "total" {print $1}' | head -5 >&2

	# for i in "${_arg_select_contig_numbers[@]}"; do
	# 	# Call the function corresponding to the current number (index is i-1)
	# 	INUM="${i}"
	# 	source "$script_dir/polap-variables-oga.sh" # '.' means 'source'
	# 	_polap_log0 $(cat "${_polap_var_wga}/mt.contig.name-${i}" | wc -l)
	# 	_polap_log0_cat "${_polap_var_mtdna_compare}"
	# done

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
		local FDIR="${_polap_var_wga}"
		JNUM="${i}"
		file="$FDIR"/mt.contig.name-$JNUM

		unique_file=$(is_unique "$file")
		if [ $? -eq 0 ]; then
			# If unique, add it to the unique_files array
			unique_files+=("$file")
			echo "$file is unique."

			MTCONTIGNAME="$file"
			INUM="${i}"
		else
			_polap_log1 "$file is the same as $unique_file."
			INUM="${unique_file##*-}"
		fi
		source "$script_dir/polap-variables-oga.sh" # '.' means 'source'
		_polap_log0 $(cat "${_polap_var_wga}/mt.contig.name-${i}" | wc -l)
		_polap_log0_cat "${_polap_var_mtdna_compare}"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
