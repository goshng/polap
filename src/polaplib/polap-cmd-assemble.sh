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
# Assemble subcommands is the key subcommand of the polap.
# It assembles the whole-genome and organelle-genome sequences.
# We tested assemble1 and assemble2 functions only.
# assemble-wrange is a copy of assemble2 with _run_polap_test-reads, which
# is used in polap-data-dflye.
# assemble is a wrapper for assemble1 and assemble2. It needs automatic seed
# contig selection, which is not tested yet.
#
# Functions:
# function _run_polap_assemble1 { # whole-genome genome assembly
# function _run_polap_assemble2 { # organelle-genome assembly
# function _run_polap_assemble-wrange { # organelle-genome assembly
# function _run_polap_assemble { # whole-genome and then organelle-genome assembly
#
# TODO:
#   - implement _run_polap_assemble-wrange
#   - implement _run_polap_assemble
#
# See Also:
#   polaplib/run-polap-function-directional.sh
#   polaplib/polap-cmd-bioproject.sh
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
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<EOF
Name:
  polap assemble1 - assemble whole-genome sequences

Synopsis:
  polap assemble1 [options]

Description:
  polap assemble1 runs the whole-genome assembly using Flye. 

Options:
  -l FASTQ
    reads data file [default: ${_arg_long_reads}]
  
  -a FASTQ
    short-read data file 1 [default: ${_arg_short_read1}]

  -b FASTQ
    short-read data file 2 [default: ${_arg_short_read2}]

  -o OUTDIR
    output folder [default: ${_arg_outdir}]

  -m INT
    minimum read length [default: ${_arg_min_read_length}]

  -t INT
    the number of CPU cores [default: ${_arg_threads}]

  -c INT
    the coverage option [default: ${_arg_coverage}]

  -g <arg>
    computed by find-genome-size polap command or given by users

  --reduction-reads [default: ${_arg_reduction_reads}]
    data reduction in a whole-genome assembly

  --redo [default: ${_arg_redo}]
    do not use previously generated intermediate results

  --stopafter {data,flye1} [default: ${_arg_stopafter}]
    stop after data or flye1 step

Outputs:
  ${_polap_var_outdir_s1_fq_stats}

  ${_polap_var_outdir_s2_fq_stats}

  ${_polap_var_outdir_long_total_length}

  ${_polap_var_outdir_genome_size}

  ${_polap_var_outdir_nk_fq_gz}

  ${_polap_var_outdir_lk_fq_gz}

  ${_polap_var_wga_contigger_edges_gfa}

  ${_polap_var_ga_contigger_edges_stats}

Examples:
  Get organelle genome sequences:
    polap assemble -l l.fq -a s1.fq -b s2.fq -o outdir

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show3 "assemble1" help_message || return 0

	# Display help message
	# if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
	# 	local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
	# 	man "$manfile" >&3
	# 	rm -f "$manfile"
	# 	return
	# fi

	# Display help message
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
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1

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
		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_run_polap_find-genome-size-for-pacbio
		else
			_run_polap_find-genome-size
		fi
	fi

	# TODO: we leave the polishing step later after the organelle-genome assembly
	# It is because polishing preparation and polishing itself can use
	# much more memory than genome assembly.
	# prepare-polishing early on to delete short-read data files
	# if [ -s "${_polap_var_outdir_msbwt}" ]; then
	# 	_polap_log1 "  skipping the preparation of short-read polishing ..."
	# else
	# 	if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
	# 		_polap_log1 "  decompressing ${_polap_var_outdir_msbwt_tar_gz} ... later when we polish it with the short-read data."
	# 		# tar zxf "${_polap_var_outdir_msbwt_tar_gz}"
	# 	else
	# 		_polap_log0 "  Do the preparation of short-read polishing ... early"
	# 		check_file_existence "${_arg_short_read1}"
	# 		_run_polap_prepare-polishing
	# 	fi
	# fi

	if [ -s "${_polap_var_outdir_nk_fq_gz}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping reduce-data ..."
	else
		_run_polap_reduce-data
	fi

	if [ -s "${_polap_var_outdir_l_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
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
			# Allow to use nk.fq.gz not -l arguments
			_arg_long_reads_is="off"
			_run_polap_flye1
			if [ $? -eq $RETURN_FAIL ]; then
				_polap_log0 "ERROR: flye1 step failed."
				return $RETURN_FAIL
			fi
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

	conda deactivate

	_polap_log1 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
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
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	FDIR="${_arg_outdir}"/${_arg_inum}

	_polap_var_mtcontigname="$FDIR"/mt.contig.name-"${_arg_jnum}"
	local i=0

	help_message=$(
		cat <<EOF
Name:
  polap assemble2 - assemble organelle genome sequences

Synopsis:
  polap assemble2 [options]

Description:
  polap assemble2 runs the organelle-genome assembly using Flye. 

Options:
  -i INT
    index of the source of an organelle-genome assembly [default: ${_arg_inum}]
  
  -j INT
    index of the target organelle-genome assembly [default: ${_arg_jnum}]
  
  -w INT
    minimum mapping length for read selection [default: ${_arg_single_min}]

  -c INT
    the coverage option [default: ${_arg_coverage}]

  -t INT
    the number of CPU cores [default: ${_arg_threads}]

  --polap-reads [default: ${_arg_polap_reads}]
    uses the POLAP read selection not ptGAUL's

  --coverage-check [default: ${_arg_coverage_check}]

  --no-coverage-check: no data reduction in an organelle-genome assembly

  -l FASTQ
    reads data file [default: ${_arg_long_reads}]
  
  -o OUTDIR
    output folder [default: ${_arg_outdir}]

  -m INT
    minimum read length [default: ${_arg_min_read_length}]

Inputs:
  ${_polap_var_mtcontigname}

  ${_polap_var_ga_contigger_edges_fasta}

  ${_polap_var_ga_contigger_edges_gfa} 

  if no such file: ${_polap_var_ga_contigger_edges_fasta}

Outputs:
  ${_polap_var_oga_assembly_graph_gfa}

Examples:
  Get organelle genome sequences using minimum mapping length 6000:
    polap assemble2 -l l.fq -o outdir -w 6000

  Get organelle genome sequences using seed contig index 0 and target organelle-genome assembly index 1:
    polap assemble2 -i 0 -j 1 -l l.fq -o outdir

  Get organelle genome sequences using polap read selection:
    polap assemble2 --polap-reads -l l.fq -o outdir

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
	fi

	check_file_existence "${_polap_var_mtcontigname}"
	check_file_existence "${_polap_var_ga_contigger_edges_fasta}"
	# 2025-07-15
	# allow -l option for assemble2
	#
	# check_file_existence "${_polap_var_outdir_lk_fq_gz}"

	_polap_lib_conda-ensure_conda_env polap || exit 1
	_polap_log3 "i: $i"
	_run_polap_map-reads
	_polap_log3 "i: $i"
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
	_polap_log3 "i: $i"
	_run_polap_select-reads
	_polap_log3 "i: $i"
	_run_polap_flye2
	_polap_log3 "i: $i"

	_polap_log1 "NEXT: $(basename "$0") prepare-polishing -o ${_arg_outdir} -i ${_arg_inum} -j ${_arg_jnum}"

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# A copy of _run_polap_assemble2
# replaces:
# _run_polap_select-reads
# _run_polap_flye2
# with:
# _run_polap_test-reads
#
# we need one just like assemble2 but with _run_polap_test-reads not _run_polap_select-reads.
function _run_polap_assemble-wrange { # organelle-genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	FDIR="${_arg_outdir}"/${_arg_inum}

	_polap_var_mtcontigname="$FDIR"/mt.contig.name-"${_arg_jnum}"

	help_message=$(
		cat <<HEREDOC
# Organelle-genome assembly using seed contig sequences with a range -s.
# no assembly but test-reads only, which includes flye assembly.
#
# Arguments:
#   -i ${_arg_inum}: index of the source of an organelle-genome assembly
#   -j ${_arg_jnum}: index of the target organelle-genome assembly
#   -s, --select-read-range <start,end,count>
#   --start-index <index>: to start at somewhere not start
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
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum} -j ${_arg_jnum} --select-read-range 3000,27000,5
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

	# 2025-06-13: done for polap-data-dflye
	# test-reads only
	_run_polap_test-reads
	# _run_polap_flye2

	_polap_log1 "NEXT: $(basename "$0") prepare-polishing -o ${_arg_outdir} -i ${_arg_inum} -j ${_arg_jnum}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# Version: v0.1.0
# Function: recreate_dir_preserve_log
# Usage: recreate_dir_preserve_log <dir_path> [log_filename]
# Example: recreate_dir_preserve_log "$_arg_outdir" polap.log

recreate_dir_preserve_log() {
	local dir="$1"
	local logname="${2:-polap.log}"
	local logfile="${dir}/${logname}"
	local tmp_log=""

	# If log exists, move it temporarily
	if [[ -f "$logfile" ]]; then
		tmp_log="$(mktemp)"
		mv "$logfile" "$tmp_log"
	fi

	# Recreate directory
	rm -rf "$dir"
	mkdir -p "$dir"

	# Restore log
	if [[ -n "$tmp_log" ]]; then
		mv "$tmp_log" "$logfile"
	fi
}

################################################################################
# 2025-10-29
# We may be ready to implement this.
# We will use miniassemble and Oatk's extraction subcommands to implement this.
# We will also use a hybrid approach of polishing.
# Many different functions are intergrated into this main subcommand.
# We will take ONT and short-read datasets to assemble plastid and mitochondrial
# genome sequences and polish them as a final form.
#
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
#   ${_arg_long_reads}.mt.1.fasta
#   ${_arg_long_reads}.pt.1.fasta
################################################################################
function _run_polap_assemble { # whole-genome and then organelle-genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local polap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  polap assemble - assemble plant organelle-genome sequences using ONT long-read data

Synopsis:
  polap assemble [options]

Description:
  polap assemble subcommand uses minimap2 and miniasm to generate seed contigs,
which are fed into Flye to finalize the assembly. It also uses Oatk's pathfinder
to extract organelle genome sequences, which are polished using Racon, fmlrc2, and
polypolish.

Options:
  -l FASTQ
    reads data file
	
  -a FASTQ
    short-read data file 1

  -b FASTQ
    short-read data file 2

  -o OUTDIR
    output folder

  -o, --outdir: output folder name (default: ${_arg_outdir})
    The option '-o' or '--outdir' specifies the output folder name, 
    with a default value of 'o'. The output folder typically contains input 
    files that are long-read and short-read data files. Input data files can 
    be specified using the options provided by -l, -a, and -b.

  -l, --long-reads: long-reads data file in fastq format (default: ${_arg_long_reads})
    The option '-l' or '--long-reads' specifies the location of a long-reads 
    data file in fastq format, with a default filename of 'l.fq'.

  -a, --short-read1: short-read fastq file 1 (default: ${_arg_short_read1})
    The option '-a' or '--short-read1' specifies the first short-read fastq 
    file to be used, with a default value of "s1.fq".

  -b, --short-read2: short-read fastq file 2 (default: ${_arg_short_read2})
    The option '-b' or '--short-read2' specifies a short-read fastq file 2, 
    with a default value of 's2.fq'. The second short-read data file, 
    if provided, is considered optional.
    (Note: not tested yet; -a & -b are required.)

Outputs:
  ${_polap_var_ga_contigger_edges_stats}

Examples:
  Get organelle genome sequences:
    polap assemble -l l.fq -a s1.fq -b s2.fq -o outdir --prefix t1
    Check t1.pt.gfa, t1.mt.gfa, t1.pt.fa, and t1.mt.fa

See also:
  version

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi

EOF
	)

	# Display help message
	_polap_lib_help-maybe-show3 "$polap_cmd" help_message || return 0

	_polap_log0 "assembling the plant organelle genome ..."
	_polap_log0 "  long-read: ${_arg_long_reads}"
	_polap_log0 "  short-read1: ${_arg_short_read1}"
	_polap_log0 "  short-read2: ${_arg_short_read2}"

	# Delete all except the polap.log
	# Deleting output folder and files are not a good idea
	# because it accidentally deletes unintended folders.
	# recreate_dir_preserve_log "${_arg_outdir}" "polap.log"
	# find "${_arg_outdir}" -mindepth 1 ! -name 'polap.log' -exec rm -rf {} +

	_run_polap_miniassemble

	local pt_gfa="${_arg_prefix}.pt.gfa"
	local mt_gfa="${_arg_prefix}.mt.gfa"
	rm -f "${pt_gfa}"
	rm -f "${mt_gfa}"
	cp -p "${_arg_outdir}/pt.1.gfa" "${pt_gfa}"
	cp -p "${_arg_outdir}/mt.1.gfa" "${mt_gfa}"

	_arg_infile1="${_arg_outdir}/pt.1.gfa"
	_arg_infile2="${_arg_outdir}/mt.1.gfa"
	_run_polap_extract

	local pt_fa="${_arg_prefix}.pt.fa"
	local mt_fa="${_arg_prefix}.mt.fa"

	_arg_plastid="on"
	_arg_infile="${_arg_outdir}/extract/oatk.pltd.ctg.fasta"
	_arg_outfile="${_arg_outdir}/pt.1.fasta"
	if [[ -s "${_arg_infile}" ]]; then
		_run_polap_polish-longshort
		rm -rf "${pt_fa}"
		cp -p "${_arg_outfile}" "${pt_fa}"
	fi

	_arg_plastid="off"
	_arg_infile="${_arg_outdir}/extract/oatk.mito.ctg.fasta"
	_arg_outfile="${_arg_outdir}/mt.1.fasta"
	if [[ -s "${_arg_infile}" ]]; then
		_run_polap_polish-longshort
		rm -f "${mt_fa}"
		cp -p "${_arg_outfile}" "${mt_fa}"
	fi

	# Report lines
	_polap_log0_file_exist "ptDNA assembly format" "${pt_gfa}"
	_polap_log0_file_exist "mtDNA assembly format" "${mt_gfa}"
	_polap_log0_file_exist "ptDNA sequence in fasta format" "${pt_fa}"
	_polap_log0_file_exist "mtDNA sequence in fasta format" "${mt_fa}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
