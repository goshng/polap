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
set +u; [[ -n "${!_POLAP_INCLUDE_}" ]] && return 0; set -u
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
#   $ODIR/long_total_length.txt
#   $ODIR/jellyfish_out.histo
#   $ODIR/short_expected_genome_size.txt
#   ${_polap_var_base_nk_fq_gz}
#   $FDIR
################################################################################
function _run_polap_assemble1() { # whole-genome genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-base.sh" # '.' means 'source'
	source "$script_dir/polap-variables-wga.sh"  # '.' means 'source'
	LRNK="${_polap_var_base_nk_fq_gz}"

	help_message=$(
		cat <<HEREDOC
# Run the whole-genome assembly.
#
# Arguments:
#   -o ${ODIR}
#   -l ${_arg_long_reads}: a long-read fastq data file
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file
#   -m ${_arg_min_read_length}: the long-read sequence length threshold
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${_arg_long_reads}: a long-read fastq 
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file
# Outputs:
#   ${_polap_var_base_fq_stats}
#   ${_polap_var_base_long_total_length}
#   ${_polap_var_base_genome_size}
#   ${_polap_var_base_nk_fq_gz}
#   ${_polap_var_wga_contigger_gfa}
Example: $0 ${_arg_menu[0]} [-o ${ODIR}] [-l ${_arg_long_reads}] [-a ${_arg_short_read1}] [-b ${_arg_short_read2}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -s "${_polap_var_wga_contigger_gfa}" ]]; then
			_polap_log0_file "${_polap_var_wga_contigger_gfa}"
		else
			_polap_log0 "No such file: ${_polap_var_wga_contigger_gfa}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x; return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "starting the whole-genome assembly on ${ODIR} ..."

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

	if [ -s "${_polap_var_base_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping summary-reads ..."
	else
		if [[ "${_arg_flye}" == "on" ]]; then
			_run_polap_summary-reads
		else
			_polap_log2 "  skipping summary-reads ..."
		fi
	fi

	if [ -s "${_polap_var_base_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping total-length-long ..."
	else
		_run_polap_total-length-long
	fi

	if [ -s "${_polap_var_base_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping find-genome-size ..."
	else
		_run_polap_find-genome-size
	fi

	if [ -s "${_polap_var_base_nk_fq_gz}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping reduce-data ..."
	else
		_run_polap_reduce-data
	fi

	check_file_existence "${_polap_var_base_nk_fq_gz}"

	if [ -s "${_polap_var_wga_contigger_gfa}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping flye1 ..."
	else
		if [[ "${_arg_flye}" == "on" ]]; then
			_run_polap_flye1
		else
			_polap_log2 "  skipping flye1 ..."
		fi
	fi

	_polap_log1 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x; return 0
}

################################################################################
# FIXME: separate command parser
################################################################################
function ncbi_command_parse() {
	BIOPRJ=""
	SPECIES=""
	SRA=""
	while getopts "b:s:r:o:" option; do
		case $option in
		b) BIOPRJ=$OPTARG ;;
		s) SPECIES=$OPTARG ;;
		r) SRA=$OPTARG ;;
		o) ODIR=$OPTARG ;;
		# h) usage_of ;;
		\?) # incorrect option
			echo "Error: Invalid option; try -h option"
			exit
			;;
		esac
	done
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
function _run_polap_assemble2() { # organelle-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-base.sh" # '.' means 'source'
	source "$script_dir/polap-variables-ga.sh"   # '.' means 'source'
	LRNK="$ODIR/nk.fq.gz"
	FDIR="$ODIR"/$INUM
	local MTDIR="$ODIR"/$JNUM
	local MTSEEDSDIR="${MTDIR}"/seeds

	MTCONTIGNAME="$FDIR"/mt.contig.name-"$JNUM"

	# for contigs
	#	assembly_graph_final_fasta=o/30-contigger/contigs.fasta
	#	for edges
	# _polap_var_contigger_edges_fasta="$FDIR"/30-contigger/graph_final.fasta

	help_message=$(
		cat <<HEREDOC
# Selects reads mapped on a genome assembly and assembles an organelle genome.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -r ${_arg_pair_min}: minimum minimap2 alignment length for a pair of contigs
#   -x ${_arg_bridge_min}: minimum long-read length for connecting the pair of contigs
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTCONTIGNAME}
#   ${_polap_var_contigger_edges_fasta}
# Outputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/1.names
#   ${MTSEEDSDIR}/2.fq.gz
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/30-contigger/contigs.fasta
#   ${MTDIR}/30-contigger/contigs_stats.txt
#   ${MTDIR}/30-contigger/graph_final.fasta
#   ${MTDIR}/30-contigger/graph_final.gfa
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [-r|--pair-min <arg>] [-x|--bridge-min <arg>] [-w|--single-min <arg>] [-t|--threads <arg>] [-c|--coverage <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	check_file_existence "${MTCONTIGNAME}"
	check_file_existence "${_polap_var_contigger_edges_fasta}"
	check_file_existence "${_polap_var_base_nk_fq_gz}"

	_polap_log1 "NEXT: $(basename $0) select-reads -o $ODIR -i $INUM -j $JNUM"
	_run_polap_select-reads
	if [[ "${_arg_flye}" == "on" ]]; then
		_run_polap_flye2
	else
		_polap_log2 "  skipping flye2 ..."
	fi

	[ "$DEBUG" -eq 1 ] && set +x; return 0
}

################################################################################
# NOT IMPLEMENTED YET!
# because we need a manual long-read selection step.
# You could execute this menu with option --test.
#
# Runs the organelle-genome assembly.
# Arguments:
#   -o $ODIR
#   -l ${_arg_long_reads}: a long-read fastq data file
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file
# Inputs:
#   ${_arg_long_reads}: a long-read fastq
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
################################################################################
function _run_polap_assemble() { # whole-genome and then organelle-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-base.sh" # '.' means 'source'
	source "$script_dir/polap-variables-wga.sh"  # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Runs the POLAP organelle-genome assembly with sequencing data.
# 
# Arguments:
#   -o $ODIR: output folder (default: o)
#   -l ${_arg_long_reads}: a long-read fastq data file
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file
#   --rwx <arg>: mapping length for read selection
#   -m <arg>: minimum read length
#   -c <arg>: maximum coverage of reads 
#   -t <arg>: number of CPUs
#   -i <arg>: index of the source of an organelle-genome assembly
#   -j <arg>: index of the target organelle-genome assembly
#   --no-reduction-reads: no data reduction in a whole-genome assembly
#   --no-coverage-check: no data reduction in an organelle-genome assembly
#   --random-seed <arg>: 5-digit number
#   --flye-asm-coverage <arg>: Flye --asm-coverage
#   --redo: do not use previously created intermediate results
# Inputs:
#   ${_arg_long_reads}: a long-read fastq 
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
Example: $(basename $0) ${_arg_menu[0]} --test
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Not delete the output directory.
	if [[ "${_arg_redo}" = "on" ]]; then
		_polap_log3_cmd rm -rf "${ODIR}"
	fi
	mkdir -p "$ODIR"

	# Run assembly, annotation, and contig selection steps
	if [ -s "${_polap_var_wga_contigger_gfa}" ]; then
		_polap_log1 "  skipping the whole-genome assembly"
	else
		check_file_existence "${_arg_long_reads}"
		_run_polap_assemble1
	fi

	if [ -s "${_polap_var_wga_annotation}" ]; then
		_polap_log1 "  skipping the organelle annotation on the whole-genome"
	else
		_run_polap_edges-stats
		_run_polap_annotate
	fi

	# Select seed contigs
	if [[ "${_arg_test}" = "on" ]]; then
		_arg_plastid="on"
	fi
	_arg_menu[1]="auto"
	_run_polap_seeds

	_run_polap_assemble2
	_run_polap_flye-polishing

	return

	INUM="${i}" _run_polap_edges-stats
	INUM="${i}" _run_polap_annotate
	JNUM="${i}" _run_polap_flye-polishing
	INUM="${i}" _run_polap_select-mtdna

	if [ -s "${_polap_var_base_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_base_msbwt_tar_gz} ... later when we polish it with the short-read data."
			tar -zxf "${_polap_var_base_msbwt_tar_gz}" -C "${ODIR}"
		else
			_polap_log1 "  Do the preparation of short-read polishing ... early"
			check_file_existence "${_arg_short_read1}"
			check_file_existence "${_arg_short_read2}"
			_run_polap_prepare-polishing
		fi
	fi

	# Run the polishing step
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Define the paths for mtDNA sequences to be polished
		PA="${ODIR}/${i}/mt.0.fasta"
		FA="${ODIR}/${i}/mt.1.fa"

		if [ -s "${PA}" ] && [ -s "${_polap_var_base_msbwt}" ]; then
			_run_polap_polish
		fi
	done

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x; return 0
}
