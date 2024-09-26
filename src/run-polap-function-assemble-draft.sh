################################################################################
# Template for an external shell script
#
# You could use this function template to create a new menu.
# Rename template and delete x in the name. You could execute such menu
# but such menus are not created as empty files by make-menus menu.
################################################################################
function _run_polap_assemble-draft() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

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
		check_file_existence "${LRNK}"
		_run_polap_assemble1
	fi

	if [ -s "${_polap_var_wga_annotation}" ]; then
		_polap_log1 "  skipping the organelle annotation on the whole-genome"
	else
		_run_polap_annotate
	fi

	# Select seed contigs
	_run_polap_select-contigs

	# Array to store the names of the original files
	files=($(ls "${ODIR}/0/mt.contig.name-"?))

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
		INUM=0
		FDIR="${ODIR}/${INUM}"
		JNUM="${i}"
		file="$FDIR"/mt.contig.name-$JNUM
		mkdir -p "${ODIR}/${JNUM}"

		unique_file=$(is_unique "$file")
		if [ $? -eq 0 ]; then
			# If unique, add it to the unique_files array
			unique_files+=("$file")
			echo "$file is unique."

			MTCONTIGNAME="$file"
			# check the mt.contig.name-1
			if [ -s "$MTCONTIGNAME" ]; then
				# Run secondary assembly, polishing, and mtDNA selection steps
				_polap_log1_file "${MTCONTIGNAME}"

				if [ -s "${ODIR}/${JNUM}/30-contigger/graph_final.gfa" ]; then
					_polap_log2 "  skipping organelle-genome assembly -i 0 -j ${JNUM} ..."
				else
					_polap_log2 "  not skipping organelle-genome assembly -i 0 -j ${JNUM} ..."
					_arg_yes="on"
					_run_polap_assemble2
				fi

				if [ -s "${ODIR}/${JNUM}/contig-annotation-table.txt" ]; then
					_polap_log2 "  skipping organelle-genome annotation -i ${JNUM} ..."
				else
					INUM="${i}" _run_polap_annotate
				fi

				if [ -s "${ODIR}/${JNUM}/assembly_graph.gfa" ]; then
					_polap_log2 "  skipping organelle-genome long-read polishing -i ${JNUM} ..."
				else
					INUM="${i}" _run_polap_flye-polishing
				fi

				if [ -s "${ODIR}/${JNUM}/mt.0.fasta" ]; then
					_polap_log2 "  skipping organelle-genome extraction -i ${JNUM} ..."
				else
					INUM="${i}" _run_polap_select-mtdna
				fi
			else
				_polap_log1 "LOG: $MTCONTIGNAME is empty for select-contig type $i ..."
			fi
		else
			_polap_log1 "$file is the same as $unique_file."
		fi
	done

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
