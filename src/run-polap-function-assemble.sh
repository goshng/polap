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
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Runs the whole-genome assembly.
#
# Defaults:
#   l.fq
#   s1.fq
#   s2.fq
#   number of threads: ${_arg_threads}
#   assembly coverage: ${_arg_coverage}
#
# Outputs:
#   ${_arg_outdir}/long_total_length.txt
#   ${_arg_outdir}/jellyfish_out.histo
#   ${_arg_outdir}/short_expected_genome_size.txt
#   ${_polap_var_outdir_nk_fq_gz}
#   $FDIR
################################################################################
function _run_polap_assemble1 { # whole-genome genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Flye whole-genome assembly using the long-read data.
#
# Arguments:
#   -o ${_arg_outdir}
#   -l ${_arg_long_reads}: a long-read fastq data file
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file
#   -m ${_arg_min_read_length}: the long-read sequence length threshold
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the coverage option
#   --no-reduction-reads: no data reduction in a whole-genome assembly
#   -g <arg>: computed by find-genome-size menu or given by users
#   --flye-asm-coverage ${_arg_flye_asm_coverage}: Flye option --asm-coverage
#   --random-seed <arg>: 5-digit number
#   --redo: do not use previously generated intermediate results
# Inputs:
#   ${_arg_long_reads}: a long-read fastq 
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file
# Outputs:
#   ${_polap_var_outdir_s1_fq_stats}
#   ${_polap_var_outdir_s2_fq_stats}
#   ${_polap_var_outdir_long_total_length}
#   ${_polap_var_outdir_genome_size}
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_polap_var_outdir_lk_fq_gz}
#   ${_polap_var_wga_contigger_edges_gfa}
#   ${_polap_var_ga_contigger_edges_stats}
# View:
#   Check file: ${_polap_var_wga_contigger_edges_gfa}
Example: $(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir} -l ${_arg_long_reads} -a ${_arg_short_read1} -b ${_arg_short_read2}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -s "${_polap_var_wga_contigger_edges_gfa}" ]]; then
			_polap_log0_file "${_polap_var_wga_contigger_edges_gfa}"
		else
			_polap_log0 "No such file: ${_polap_var_wga_contigger_edges_gfa}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "starting the whole-genome assembly on ${_arg_outdir} ..."
	_polap_log1 "  output1: ${_polap_var_outdir_s1_fq_stats}"
	_polap_log1 "  output1: ${_polap_var_outdir_s2_fq_stats}"
	_polap_log1 "  output2: ${_polap_var_outdir_long_total_length}"
	_polap_log1 "  output3: ${_polap_var_outdir_genome_size}"
	_polap_log1 "  output4: ${_polap_var_outdir_lk_fq_gz}"
	_polap_log1 "  output5: ${_polap_var_outdir_nk_fq_gz}"
	_polap_log1 "  output6: ${_polap_var_wga_contigger_edges_gfa}"
	_polap_log1 "  output7: ${_polap_var_wga_contigger_edges_stats}"

	# Skip flye1 if you want
	if [[ "${_arg_flye}" == "on" ]]; then
		if [[ -d "${_polap_var_wga}" ]]; then
			if [[ "${_arg_redo}" = "on" ]]; then
				_polap_log3_cmd rm -rf "${_polap_var_wga}"
				_polap_log3_cmd mkdir -p "${_polap_var_wga}"
			else
				if confirm "Do you want to do the whole-genome assembly, which will delete ${_polap_var_wga}?"; then
					_polap_log0 "  deleting and creating ${_polap_var_wga} ..."
					_polap_log3_cmd rm -rf "${_polap_var_wga}"
					_polap_log3_cmd mkdir -p "${_polap_var_wga}"
				else
					_polap_log0 "You have cancelled the whole-genome assembly."
					return
				fi
			fi
		else
			_polap_log0 "  creating ${_polap_var_wga} ..."
			_polap_log3_cmd mkdir -p "${_polap_var_wga}"
		fi
	fi

	if [ -s "${_polap_var_outdir_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping total-length-long ..."
	else
		_run_polap_total-length-long
	fi

	if [ -s "${_polap_var_outdir_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping find-genome-size ..."
	else
		_run_polap_find-genome-size
	fi

	# prepare-polishing early on to delete short-read data files
	if [ -s "${_polap_var_outdir_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_outdir_msbwt_tar_gz} ... later when we polish it with the short-read data."
			# tar zxf "${_polap_var_outdir_msbwt_tar_gz}"
		else
			_polap_log0 "  Do the preparation of short-read polishing ... early"
			check_file_existence "${_arg_short_read1}"
			_run_polap_prepare-polishing
		fi
	fi

	if [ -s "${_polap_var_outdir_nk_fq_gz}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping reduce-data ..."
	else
		_run_polap_reduce-data
	fi

	if [ -s "${_polap_var_outdir_lk_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping summary-reads ..."
	else
		_run_polap_summary-reads
	fi

	if [[ "${_arg_stopafter}" == "data" ]]; then
		return 0
	fi

	if [ -s "${_polap_var_wga_contigger_edges_gfa}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping flye1 ..."
	else
		if [[ "${_arg_flye}" == "on" ]]; then
			_run_polap_flye1
		else
			_polap_log2 "  skipping flye1 ..."
		fi
	fi

	if [ -s "${_polap_var_ga_contigger_edges_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping edges-stats ..."
	else
		_run_polap_edges-stats
		if [ $? -eq $RETURN_FAIL ]; then
			_polap_log0 "ERROR: edges-stats step failed."
			return $RETURN_FAIL
		fi
	fi

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping annotation ..."
	else
		_run_polap_annotate
		if [ $? -eq $RETURN_FAIL ]; then
			_polap_log0 "ERROR: annotation step failed."
			return $RETURN_FAIL
		fi
	fi

	_polap_log1 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Selects and assembles long-read data.
# Arguments:
#   -i 0
#   -j 1
#   -o o
#   -t ${_arg_threads}
#   -m ${_arg_min_read_length}
#   MPAIR
#   MBRIDGE
#   COV
#   ${_arg_circularize}
# Inputs:
# Outputs:
################################################################################
function _run_polap_assemble2 { # organelle-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	FDIR="${_arg_outdir}"/${_arg_inum}

	_polap_var_mtcontigname="$FDIR"/mt.contig.name-"${_arg_jnum}"

	help_message=$(
		cat <<HEREDOC
# Organelle-genome assembly using seed contig sequences.
#
# Arguments:
#   -i ${_arg_inum}: index of the source of an organelle-genome assembly
#   -j ${_arg_jnum}: index of the target organelle-genome assembly
#   -w ${_arg_single_min}: minimum mapping length for read selection
#   -m ${_arg_min_read_length}: minimum read length
#   --polap-reads: uses the POLAP read selection not ptGAUL's
#   -c ${_arg_coverage}: maximum coverage of reads 
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by seed contig size or given by users
#   --no-coverage-check: no data reduction in an organelle-genome assembly
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
#   ${_polap_var_ga_contigger_edges_gfa} 
#     if no such file: ${_polap_var_ga_contigger_edges_fasta}
# Outputs:
#   ${_polap_var_oga_assembly_graph_gfa}
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum} -j ${_arg_jnum} -w ${_arg_single_min}
Example: $(basename "$0") ${_arg_menu[0]} --polap-reads
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
	fi

	check_file_existence "${_polap_var_mtcontigname}"
	check_file_existence "${_polap_var_ga_contigger_edges_fasta}"
	check_file_existence "${_polap_var_outdir_lk_fq_gz}"

	_run_polap_map-reads
	if [[ "${_arg_polap_reads}" == "on" ]]; then
		_arg_menu[1]="polap-reads"

		local _n=$(wc -l <"${_polap_var_mtcontigname}")
		if [[ "${_n}" -eq 1 ]]; then
			_polap_log0 "  single seed contig; read-selection type change: polap-reads -> intra-reads"
			_arg_menu[1]="intra-reads"
		fi
	else
		_arg_menu[1]="infile"
	fi
	_run_polap_select-reads
	_run_polap_flye2

	_polap_log1 "NEXT: $(basename "$0") prepare-polishing -o ${_arg_outdir} -i ${_arg_inum} -j ${_arg_jnum}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# NOT IMPLEMENTED YET!
# because we need a manual long-read selection step.
# You could execute this menu with option --test.
#
# Runs the organelle-genome assembly.
# Arguments:
#   -o ${_arg_outdir}
#   -l ${_arg_long_reads}: a long-read fastq data file
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file
# Inputs:
#   ${_arg_long_reads}: a long-read fastq
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file
# Outputs:
################################################################################
function _run_polap_assemble { # whole-genome and then organelle-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Run the POLAP organelle-genome assembly with sequencing data.
#
# Steps:
#   1. assemble1
#   2. seeds -> multiple mt.contig.name files
#   3. assemble2 -> assemble as many organelle genomes as mt.contig.name files
# 
# Arguments:
#   -o ${_arg_outdir}: output folder
#   -l ${_arg_long_reads}: a long-read fastq data file
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file
# Inputs:
#   ${_arg_long_reads}: a long-read fastq 
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file
# Outputs:
#   ${_polap_var_oga_assembly_graph_gfa}
Example: $(basename $0) -l ${_arg_long_reads} -a ${_arg_short_read1} -b ${_arg_short_read2}
Example: $(basename $0) --test
Example: $(basename $0) ${_arg_menu[0]} -o o2
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "assembling the organelle genome ..."
	_polap_log1 "  input1: ${_arg_long_reads}"
	_polap_log1 "  input2: ${_arg_short_read1}"
	_polap_log1 "  input3: ${_arg_short_read2}"
	_polap_log1 "  output1: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log1 "  output2: ${_polap_var_ga_annotation_all}"
	_polap_log1 "  output3: ${_polap_var_mtcontigname}"
	_polap_log1 "  output4: ${_polap_var_oga_contigger_edges_gfa}"
	_polap_log1 "  output5: ${_polap_var_oga_assembly_graph_gfa}"

	# Delete all except the polap.log
	if [[ "${_arg_redo}" = "on" ]]; then
		rm -rf "${_arg_outdir}/${_arg_inum}" "${_arg_outdir}/${_arg_jnum}"
		rm -f "${_arg_outdir}"/jellyfish_out* \
			"${_arg_outdir}"/*.txt \
			"${_arg_outdir}"/*.stats \
			"${_arg_outdir}"/*.gz
	fi

	# Run assembly, annotation, and contig selection steps
	if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]] &&
		[[ "${_arg_redo}" == "off" ]]; then
		_polap_log1 "  found: ${_polap_var_ga_contigger_edges_gfa}, skipping the whole-genome assembly"
	else
		check_file_existence "${_arg_long_reads}"
		check_file_existence "${_arg_short_read1}"
		_run_polap_assemble1
	fi

	if [[ -s "${_polap_var_mtcontigname}" ]] &&
		[[ "${_arg_redo}" == "off" ]]; then
		_polap_log0 "  found: ${_polap_var_mtcontigname}, skipping seed contig selection"
	else
		_run_polap_seeds
	fi

	# Loop over all files that match the pattern "mt.contig.name-<number>"
	for file in "${_polap_var_ga}"/mt.contig.name-*; do
		# Extract the <number> part using parameter expansion
		file=$(basename $file)
		local number="${file#mt.contig.name-}"

		_arg_jnum="${number}"
		source "$script_dir/polap-variables-common.sh"
		if [[ -s "${_polap_var_oga_contigger_edges_gfa}" ]] &&
			[[ "${_arg_redo}" == "off" ]]; then
			_polap_log0 "  found: ${_polap_var_oga_contigger_edges_gfa}, skipping organelle-genome assembly"
		else
			_run_polap_assemble2
		fi
	done

	# We need a way of extracting mtDNA sequences.
	# _run_polap_select-mtdna
	# _run_polap_preepare-polishing
	# _run_polap_polish

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
