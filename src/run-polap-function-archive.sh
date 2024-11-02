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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Archive the ${ODIR} folder to ${_arg_archive}
################################################################################
function _run_polap_archive() { # archive a POLAP output folder for later use
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	if [ "${_arg_short_read1_is}" = "on" ]; then
		_arg_archive="${_arg_short_read1}"
	fi

	local FDIR="$ODIR"/$INUM
	local MTCONTIGNAME="$FDIR"/mt.contig.name-"$JNUM"
	local _polap_var_source_0="${_polap_var_wga}"
	local _polap_var_source_0_30_contigger="$FDIR"/"$JNUM"/mtcontigs

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Archive the ${ODIR} folder to ${_arg_archive}
#
# Arguments:
#   -o ${ODIR}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
#
# Inputs:
#   ${ODIR}
#
# Outputs:
#   ${_arg_archive}
#
Example: $(basename $0) ${_arg_menu[0]} -o <folder> --archive <folder>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0_log "archiving ${ODIR} to ${_arg_archive} ..."

	local source_folders=(
		./0-bioproject
		./0/30-contigger
	)

	local source_files0=(
		./0-bioproject/2-mtdna.fasta
		./0-bioproject/1-species.txt
		./0-bioproject/2-mtdna.accession
		./0-bioproject/1-sra-short-read.tsv
		./0-bioproject/1-taxonomy.txt
		./0-bioproject/1-mtdna.fasta
		./0-bioproject/1-taxon-id.txt
		./0-bioproject/1-sra-long-read.tsv
		./0-bioproject/1-runinfo.all
		./0-bioproject/1-runinfo.tsv
		./0-bioproject/1-mtdna.fasta.stats
		./0-bioproject/1-passed.txt
		./short_expected_genome_size.txt
		./long_total_length.txt
		./fq.stats
		./nk.fq.stats
		./nk.fq.gz
		./0/30-contigger/graph_final.gfa
		./0/30-contigger/contigs_stats.txt
		./0/30-contigger/contigs.fasta
		./0/30-contigger/graph_final.fasta
		./0/assembly_info_organelle_annotation_count-all.txt
		./0/assembly_info_organelle_annotation_count.txt
		./0/contig-annotation-table.txt
	)

	local source_files1=(
		./30-contigger/graph_final.gfa
		./30-contigger/contigs_stats.txt
		./30-contigger/contigs.fasta
		./30-contigger/graph_final.fasta
		./mt.0.fasta
		./mt.1.fa
	)

	local source_files2=(
		./1/mt.0.fasta
		./1/contig-annotation-table.txt
		./1/30-contigger
		./1/30-contigger/graph_final.gfa
		./1/30-contigger/scaffolds_links.txt
		./1/30-contigger/contigs_stats.txt
		./1/30-contigger/contigs.fasta
		./1/30-contigger/contigs.fasta.fai
		./1/30-contigger/graph_final.fasta
		./1/30-contigger/graph_final.gv
		./1/assembly_graph.gv
		./1/contig_total_length.txt
		./1/assembly.fasta
		./1/mt.0.edges
		./1/mt.1.fa
		./1/assembly_info_organelle_annotation_count-all.txt
		./1/contig.fa
		./1/assembly_info_organelle_annotation_count.txt
		./1/params.json
		./1/assembly_graph.gfa
		./1/contig.tab
		./0/mt.contig.name-6
		./0/mt.contig.name-5
		./0/mt.contig.name-3
		./0/mt.contig.name-1
		./0/mt.contig.name-2
		./organelle-assembly_0-3
		./organelle-assembly_0-2
		./organelle-assembly_0-6
		./organelle-assembly_0-5
		./organelle-assembly_0-4
	)

	if [ -d "${_arg_archive}" ]; then
		if confirm "Do you want to redo? It will delete the folder ${_arg_archive}"; then
			_polap_log0 "deleting the folder: ${_arg_archive} ..."
			rm -rf "${_arg_archive}"
		else
			_polap_log0 "You cancelled the operation."
			_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
			[ "$DEBUG" -eq 1 ] && set +x
			exit $EXIT_SUCCESS
		fi
	fi

	_polap_log1 "creating folders ..."
	for f in "${source_folders[@]}"; do
		_polap_log2 "mkdir -p ${_arg_archive}/${f}"
		mkdir -p "${_arg_archive}/${f}"
	done

	_polap_log1 "compressing the long-read data file ..."
	if [ -s "${_polap_var_base_sra_long_fastq}" ]; then
		gzip -c "${_polap_var_base_sra_long_fastq}" \
			>"${_arg_archive}/l.fq.gz"
	fi
	_polap_log1 "compressing the short-read polishing data file ..."
	if [ -s "${ODIR}/msbwt/comp_msbwt.npy" ]; then
		cd "${ODIR}"
		tar zcf ../"${_arg_archive}/msbwt.tar.gz" msbwt
		cd -
	fi

	cp -p "${ODIR}"/log* "${_arg_archive}/"

	_polap_log1 "copying folders ..."
	for f in "${source_files0[@]}"; do
		_polap_log2 "cp -p" "${ODIR}/${f}" "${_arg_archive}/${f}"
		cp -p "${ODIR}/${f}" "${_arg_archive}/${f}"
	done

	for gfa in $(find "${ODIR}" -name graph_final.gfa); do
		local b="${gfa%/30-contigger/graph_final.gfa}"
		# Extract the final directory name from the path (e.g., "0" or "1")
		b=$(basename "${b}")
		# Skip if the directory path is "o/0"
		if [ "${b}" == "0" ]; then
			continue
		fi
		_polap_log2 "mkdir -p" "${_arg_archive}/${b}/30-contigger"
		mkdir -p "${_arg_archive}/${b}/30-contigger"
		for f in "${source_files1[@]}"; do
			_polap_log2 "cp -p" "${ODIR}/${b}/${f}" "${_arg_archive}/${b}/${f}"
			cp -p "${ODIR}/${b}/${f}" "${_arg_archive}/${b}/${f}"
		done
		cp -p "${_polap_var_wga}/mt.contig.name-${b}" "${_arg_archive}/0/"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Clean up the ${ODIR}.
################################################################################
function _run_polap_cleanup() { # cleanup an POLAP output folder
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Clean up the ${ODIR}.
#
# Arguments:
#   -o ${ODIR}: the output directory
#   -j ${JNUM}: the target assembly
# Inputs:
#   ${ODIR}
# Outputs:
#
Example: $0 ${_arg_menu[0]} -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log3_cmd rm -f "${_polap_var_base_lk_fq_gz}"
	_polap_log3_cmd rm -rf "${_polap_var_oga}"/{00-assembly,10-consensus,20-repeat,40-polishing}
	_polap_log3_cmd rm -rf "${_polap_var_oga_contig}"/{contig.paf,contig.tab}
	for _pread_sel in ptgaul-intra-base-length single-intra-base-length combined-intra-base-length; do
		_polap_log3_cmd rm -rf "${_polap_var_oga_reads}/${_pread_sel}"
		_polap_log3_cmd rm -rf "${_polap_var_oga_seeds}/${_pread_sel}"
		_polap_log3_cmd rm -rf "${_polap_var_oga_sample}/${_pread_sel}"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

source "$script_dir/run-polap-function-utilities.sh"
source "$script_dir/run-polap-function-menus.sh"

################################################################################
# Initializes polap analysis in a starting folder,
#
# creating an output folder.
# Arguments:
#   -o $ODIR
# Inputs: nothing
# Outputs:
#   $ODIR
################################################################################
function _run_polap_init() { # initialize an output folder
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
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
#   -o ${ODIR}: the output folder
# Inputs: none
# Outputs:
#   ${ODIR}
Example: $0 ${_arg_menu[0]} -o ${ODIR}
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
		if [[ -d "${ODIR}" ]]; then
			ls -l "${ODIR}" >&3
		else
			_polap_log0 "No such output folder: ${ODIR}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	mkdir -p "${ODIR}"
	_polap_log0 "creating output folder [$ODIR] if no such folder exists ..."
	if [ "$ODIR" != "o" ]; then
		_polap_log1 "  Use -o $ODIR option in all subsequent analysis"
		_polap_log1 "  because your output folder is not the default of 'o'."
	fi
	_run_polap_make-menus
	_log_command_versions

	_polap_log1 "NEXT: $0 summary-reads -o ${ODIR} -l ${_arg_long_reads}"
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# View the polap log file.
################################################################################
function _run_polap_log() { # display the polap log
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Display the polap log.
#
# Arguments:
#   -o ${ODIR}: the output folder
# Inputs: none
# Outputs:
#   ${ODIR}
Example: $0 ${_arg_menu[0]} -o ${ODIR}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -d "${ODIR}" ]]; then
			ls -l "${ODIR}" >&3
		else
			_polap_log0 "No such output folder: ${ODIR}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	less "${LOG_FILE}" >&3

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
