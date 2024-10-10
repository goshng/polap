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
# Runs the whole-genome assembly.
# Defaults:
#   l.fq
#   s1.fq
#   s2.fq
#   number of threads: $NT
#   assembly coverage: $COV
# Outputs:
#   $ODIR/long_total_length.txt
#   $ODIR/jellyfish_out.histo
#   $ODIR/short_expected_genome_size.txt
#   $LRNK
#   $FDIR
################################################################################
function _run_polap_assemble1() {
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
# Runs the whole-genome assembly.
# Arguments:
#   -o $ODIR
#   -l $LR: a long-read fastq data file
#   -a $SR1: a short-read fastq data file
#   -b $SR2: another short-read fastq data file
#   -m $MR: the long-read sequence length threshold
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $LR: a long-read fastq 
#   $SR1: a short-read fastq data file
#   $SR2: another short-read fastq data file
ta file
# Outputs:
#   $ODIR/long_total_length.txt
#   $ODIR/short_expected_genome_size.txt
#   $LRNK
#   $FDIR/30-contigger/contigs.fasta
#   $FDIR/30-contigger/contigs_stats.txt
#   $FDIR/30-contigger/graph_final.fasta
#   $FDIR/30-contigger/graph_final.gfa
Example: $(basename $0) ${_arg_menu[0]} [-o|--outdir <arg>] [-l|--long-reads <arg>] [-a|--short-read1 <arg>] [-b|--short-read2 <arg>] [-t|--threads <arg>] [-c|--coverage <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		if [[ -s "${_polap_var_wga_contigger_gfa}" ]]; then
			_polap_log0_file "${_polap_var_wga_contigger_gfa}"
		fi
		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	# check_file_existence "${LR}"
	# check_file_existence "${SR1}"
	# check_file_existence "${SR2}"

	source "$script_dir/polap-variables-base.sh" # '.' means 'source'

	if [ -s "${_polap_var_base_fq_stats}" ]; then
		_polap_log2 "  skipping summary-reads ..."
	else
		_run_polap_summary-reads
	fi

	if [ -s "${_polap_var_base_long_total_length}" ]; then
		_polap_log2 "  skipping total-length-long ..."
	else
		_run_polap_total-length-long
	fi

	if [ -s "${_polap_var_base_genome_size}" ]; then
		_polap_log2 "  skipping find-genome-size ..."
	else
		_run_polap_find-genome-size
	fi

	if [ -s "${_polap_var_base_nk_fq_gz}" ]; then
		_polap_log2 "  skipping reduce-data ..."
	else
		_run_polap_reduce-data
	fi

	check_file_existence "${_polap_var_base_nk_fq_gz}"

	_run_polap_flye1

	_polap_log1 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
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
#   -t $NT
#   -m $MR
#   MPAIR
#   MBRIDGE
#   COV
#   CIRCULARIZE
# Inputs:
# Outputs:
################################################################################
function _run_polap_assemble2() {
	[ "$DEBUG" -eq 1 ] && set -x

	LRNK="$ODIR/nk.fq.gz"
	MR=$_arg_min_read_length
	FDIR="$ODIR"/$INUM
	ADIR="$FDIR"/50-annotation
	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="$MTDIR"/seeds

	MTCONTIGNAME="$FDIR"/mt.contig.name-"$JNUM"

	# for contigs
	#	assembly_graph_final_fasta=o/30-contigger/contigs.fasta
	#	for edges
	assembly_graph_final_fasta="$FDIR"/30-contigger/graph_final.fasta

	help_message=$(
		cat <<HEREDOC
# Selects reads mapped on a genome assembly and assembles an organelle genome.
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -r $MPAIR: minimum minimap2 alignment length for a pair of contigs
#   -x $MBRIDGE: minimum long-read length for connecting the pair of contigs
#   -w $MSINGLE: minimum minimap2 alignment length for a single contig
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $MTCONTIGNAME
#   ${assembly_graph_final_fasta}
# Outputs:
#   $MTDIR/contig.fa
#   $MTSEEDSDIR/1.names
#   $MTSEEDSDIR/2.fq.gz
#   $MTDIR/contig_total_length.txt
#   $MTDIR/30-contigger/contigs.fasta
#   $MTDIR/30-contigger/contigs_stats.txt
#   $MTDIR/30-contigger/graph_final.fasta
#   $MTDIR/30-contigger/graph_final.gfa
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [-r|--pair-min <arg>] [-x|--bridge-min <arg>] [-w|--single-min <arg>] [-t|--threads <arg>] [-c|--coverage <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	check_file_existence "${MTCONTIGNAME}"
	check_file_existence "${assembly_graph_final_fasta}"
	check_file_existence "${LRNK}"

	echoerr "NEXT: $(basename $0) select-reads -o $ODIR -i $INUM -j $JNUM"
	_run_polap_select-reads
	_run_polap_flye2

	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# NOT IMPLEMENTED YET!
# because we need a manual long-read selection step.
# You could execute this menu with option --test.
#
# Runs the organelle-genome assembly.
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
################################################################################
function _run_polap_assemble() {
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
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Not delete the output directory.
	mkdir -p "$ODIR"

	# Run assembly, annotation, and contig selection steps
	if [ -s "${_polap_var_wga_contigger_gfa}" ]; then
		_polap_log1 "  skipping the whole-genome assembly"
	else
		check_file_existence "${LR}"
		_run_polap_assemble1
	fi

	if [ -s "${_polap_var_wga_annotation}" ]; then
		_polap_log1 "  skipping the organelle annotation on the whole-genome"
	else
		_run_polap_annotate
	fi

	# Select seed contigs
	_run_polap_select-contigs

	# Loop over numbers from 1 to 5
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Call the function corresponding to the current number (index is i-1)
		INUM=0
		FDIR="${ODIR}/${INUM}"
		JNUM="${i}"

		MTCONTIGNAME="$FDIR"/mt.contig.name-$JNUM
		# check the mt.contig.name-1
		if [ -s "$MTCONTIGNAME" ]; then
			# Run secondary assembly, polishing, and mtDNA selection steps
			_polap_log1_file "${MTCONTIGNAME}"
			_run_polap_assemble2
			INUM="${i}" _run_polap_annotate
			INUM="${i}" _run_polap_flye-polishing
			INUM="${i}" _run_polap_select-mtdna
		else
			_polap_log1 "LOG: $MTCONTIGNAME is empty for select-contig type $i ..."
		fi
	done

	if [ -s "${_polap_var_base_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_base_msbwt_tar_gz} ... later when we polish it with the short-read data."
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

		if [ -s "${PA}" ] && [ -s "${_polap_var_base_msbwt}" ]; then
			_run_polap_polish
		fi
	done

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
