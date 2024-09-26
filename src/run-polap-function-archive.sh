################################################################################
# Archive the ${ODIR} folder to ${_arg_archive}
################################################################################
function _run_polap_archive() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'

	if [ "${_arg_short_read1_is}" = "on" ]; then
		_arg_archive="${_arg_short_read1}"
	fi
	local FDIR="$ODIR"/$INUM
	local MTCONTIGNAME="$FDIR"/mt.contig.name-"$JNUM"
	local _polap_var_source_0="${ODIR}/0"
	local _polap_var_source_0_30_contigger="$FDIR"/"$JNUM"/mtcontigs

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Archive the ${ODIR} folder to ${_arg_archive}
#
# Arguments:
#   -o ${ODIR}: the source output folder to archive
#   -a ${_arg_archive}: the target output folder for the archive [default: $_arg_archive]
#
# Inputs:
#   ${ODIR}
#
# Outputs:
#   ${_arg_archive}
#
Example: $(basename $0) ${_arg_menu[0]} [-o|--outdir <folder>] [-a|--archive <folder>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

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
			_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
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
		cp -p "${ODIR}/0/mt.contig.name-${b}" "${_arg_archive}/0/"
	done

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
