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
	local _polap_var_mtcontigname="$FDIR"/mt.contig.name-"$JNUM"
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
# Inputs:
#   ${ODIR}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
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
	if [ -s "${_polap_var_outdir_sra_long_fastq}" ]; then
		gzip -c "${_polap_var_outdir_sra_long_fastq}" \
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
# To filter a GFA file by a specific depth range
################################################################################
function _polap_archive_gfa-depth-filtered() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
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
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-depthfilter-gfa.R \
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
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_filter-gfa-with-edge() { # archive a POLAP output folder for later use
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-package-common.sh"

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
#   -o ${ODIR}: the source output folder to archive
# Inputs:
#   ${ODIR}
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0 "reducing the whole-genome assembly with the annotated and selected MT contigs ..."
	_polap_log1 "  input1: ${_polap_var_ga_annotation_table}"
	_polap_log1 "  input2: ${_polap_var_mtcontigname}"

	INUM="0"
	source "$script_dir/polap-variables-common.sh"

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
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-contig2edge.R \
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
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_package() { # archive a POLAP output folder for later use
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-package-common.sh"

	if [ "${_arg_short_read1_is}" = "on" ]; then
		_arg_archive="${_arg_short_read1}"
	fi

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Package the ${ODIR} folder to ${_arg_archive}
#
# Arguments:
#   -o ${ODIR}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive
# Inputs:
#   ${ODIR}
# Outputs:
#   ${_arg_archive}
Example: $(basename "$0") ${_arg_menu[0]} -o <folder1> -a <folder2>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0 "packaging ${ODIR} to ${_arg_archive} ..."

	# Create a directory structure that mirrors the existing hierarchy by duplicating all folders and subfolders.
	rsync -a -f"+ */" -f"- *" "${ODIR}"/ "${_arg_archive}"/

	# Copy all contents from the base directory into a new location.
	cp -fp "${_polap_var_outdir_genome_size}" \
		"${_ppack_var_base_genome_size}" 2>/dev/null
	cp -fp "${_polap_var_outdir_long_total_length}" \
		"${_ppack_var_base_long_total_length}" 2>/dev/null
	rsync -a "${_polap_var_project}"/ "${_ppack_var_project}"/ 2>/dev/null

	# copy mt.contig.name files
	for dir in $(find ${ODIR} -maxdepth 1 -type d -regex '.*/[0-9]+$' | sort); do
		local _assembly_number=$(basename "$dir")

		INUM="${_assembly_number}"
		JNUM="${_assembly_number}"
		source "$script_dir/polap-variables-common.sh"
		source "$script_dir/polap-package-common.sh"

		if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
			# Copy all contents from the "whole-genome" assembly folder into a new location.
			_polap_log1 "    step 3-1: creating GFA without sequence data"
			_polap_log1 "      input1: ${_polap_var_ga_contigger_edges_gfa}"
			_polap_log2 "      output: ${_ppack_var_ga_contigger_edges_gfa}"
			_polap_log3_pipe "gfatools view \
		    -S ${_polap_var_ga_contigger_edges_gfa} \
		    >${_ppack_var_assembly_graph_final_gfa}.without.sequences \
		    2>$_polap_output_dest"
		fi

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

				_polap_log3_pipe "gfatools view \
		      -l @${_file_selected_edges} \
    	  	${_polap_var_ga_contigger_edges_gfa} \
    	  	2>$_polap_output_dest \
    	  	>${_ppack_var_ga_contigger}/graph_final-${number}.gfa"

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
			INUM="${_assembly_number}"
			JNUM="${_assembly_number}"
			source "$script_dir/polap-variables-common.sh"
			source "$script_dir/polap-package-common.sh"
			# Your processing code here
			rsync -a "${_polap_var_oga}/01-contig"/*.txt "${_ppack_var_oga}/01-contig/" 2>/dev/null
			rsync -a "${_polap_var_oga}/06-summary/" "${_ppack_var_oga}/06-summary/" 2>/dev/null
			rsync -a "${_polap_var_oga}/07-plot/" "${_ppack_var_oga}/07-plot/" 2>/dev/null
			cp -fp "${_polap_var_oga_assembly_graph_gfa}" "${_ppack_var_oga}" 2>/dev/null
			cp -fp "${_polap_var_oga}"/*.png "${_ppack_var_oga}" 2>/dev/null
		fi
	done

	cp -pf "${LOG_FILE}" "${_arg_archive}"

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
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR}
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

	_polap_log1 "NEXT: $(basename "$0") summary-reads -o ${ODIR} -l ${_arg_long_reads}"
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
#   --log <FILE>
# Inputs:
#   ${LOG_FILE}
# Outputs:
#   a page view of the log file: ${LOG_FILE}
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR}
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
		return 0
	fi

	less "${LOG_FILE}" >&3

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
