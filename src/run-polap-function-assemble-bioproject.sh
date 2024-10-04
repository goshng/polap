################################################################################
# Assembles organelle genome from a BioProject using long and short reads.
# Compares the known and assembled mtDNA using BLAST.
################################################################################
function _run_polap_assemble-bioproject() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables for file paths
	source "$script_dir/polap-variables-base.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'

	# Set the output directory for the current job number
	if [[ -s "${_polap_var_bioproject_txt}" ]]; then
		_arg_bioproject=$(<"${_polap_var_bioproject_txt}")
		local BIOPRJ="${_arg_bioproject}"
	else
		_polap_log0 "NOTE: you have no ${_polap_var_bioproject_txt}"
		_polap_log0 "  save BioProject ID in ${_polap_var_bioproject_txt}"
		_polap_log0 "  if you already have the folder ${_polap_var_bioproject}"
		if [ "${_arg_short_read2_is}" = "on" ]; then
			_arg_bioproject="${SR2}"
		fi

		if [ -z "${_arg_bioproject}" ]; then
			_polap_log0 "ERROR: use --bioproject option"
			exit $EXIT_SUCCESS
		fi
	fi

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Runs the organelle-genome assembly for a given BioProject ID.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject
#   -b ${BIOPRJ}: BioProject ID for creating o/0-bioproject
# Inputs:
#   ${ODIR}
# Outputs:
#   mt.1.fa
Example: $(basename $0) ${_arg_menu[0]} -o ${ODIR} [-b ${BIOPRJ}]
HEREDOC
	)

	LRNK="${ODIR}/nk.fq.gz"

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_log0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "assembling organelle genomes of BioProject ${BIOPRJ} ..."

	if [ -d "${ODIR}" ]; then
		_polap_log2_file "the main output folder: ${ODIR}"
	else
		_polap_log2 "  no output folder; creating ${ODIR}"
		mkdir -p "${ODIR}"
	fi

	# Get the long-read information from NCBI
	if [ -s "${_polap_var_bioproject_sra_long_read}" ]; then
		_polap_log2 "  using SRA info at: ${_polap_var_bioproject_sra_long_read}"
		_polap_log2_cat "${_polap_var_bioproject_sra_long_read}"
	else
		_run_polap_get-bioproject
	fi

	# Fetch the long-read dataset
	if [ -s "${LRNK}" ]; then
		_polap_log1 "we use the long-read data: ${LRNK}"
	elif [[ -s "${_polap_var_base_l_fq_gz}" ]]; then
		_polap_log1 "we use the long-read data: ${_polap_var_base_l_fq_gz}"
	elif [ -s "${_polap_var_bioproject_sra_long_read}" ]; then
		local SRA=$(cut -f1 "${_polap_var_bioproject_sra_long_read}")
		LR="${ODIR}/${SRA}.fastq"
		if [ -s "${LR}" ]; then
			_polap_log2_file "${LR}"
		elif [ -s "${LR}.gz" ]; then
			_polap_log2_file "${LR}.gz is being extracted ..."
			gunzip "${LR}.gz"
		else
			"$script_dir"/run-polap-ncbitools fetch sra "$SRA"
			mv "$SRA.fastq" "${ODIR}/"
		fi
	else
		die "ERROR: no long-read dataset for the BioProject: $BIOPRJ"
	fi

	# Fetch the short-read dataset
	if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
		_polap_log1 "we use the short-read data in: ${_polap_var_base_msbwt_tar_gz}"
		# check_file_existence "${_polap_var_base_genome_size}"
	elif [[ -s "${_polap_var_base_msbwt}" ]]; then
		_polap_log1 "we use the short-read data in: ${_polap_var_base_msbwt}"
		# check_file_existence "${_polap_var_base_genome_size}"
	elif [ -s "${_polap_var_bioproject_sra_short_read}" ]; then
		local SRA=$(cut -f1 "${_polap_var_bioproject_sra_short_read}")
		SR1="${ODIR}/${SRA}_1.fastq"
		SR2="${ODIR}/${SRA}_2.fastq"
		if [ -s "${SR1}" ]; then
			_polap_log2_file "${SR1}"
		elif [ -s "${SR1}.gz" ]; then
			_polap_log2_file "${SR1}.gz is being extracted ..."
			gunzip "${SR1}.gz"
		fi

		if [ -s "${SR2}" ]; then
			_polap_log2_file "${SR2}"
		elif [ -s "${SR2}.gz" ]; then
			_polap_log2_file "${SR2}.gz is being extracted ..."
			gunzip "${SR2}.gz"
		fi
		# if [ -s "${SR1}" ]; then
		# 	seqkit stats -T "${SR1}" >"${_polap_var_base_fq_stats}"
		# fi
		# if [ -s "${SR2}" ]; then
		# 	seqkit stats -T "${SR2}" >>"${_polap_var_base_fq_stats}"
		# fi

		if [ ! -s "${SR1}" ] || [ ! -s "${SR2}" ]; then
			_polap_log1 "  downloading the paired-end short-read data: $SRA"
			"$script_dir"/run-polap-ncbitools fetch sra "$SRA"
			mv "${SRA}_1.fastq" "${ODIR}/"
			mv "${SRA}_2.fastq" "${ODIR}/"
		fi
	else
		die "ERROR: no read short-read dataset for the BioProject: $BIOPRJ"
	fi

	# check_file_existence "${LR}"
	# check_file_existence "${SR1}"
	# check_file_existence "${SR2}"
	#
	# if [ -s "${_polap_var_base_fq_stats}" ]; then
	# 	_polap_log2_file "${_polap_var_base_fq_stats}"
	# 	_polap_log3_cat "${_polap_var_base_fq_stats}"
	# else
	# 	_run_polap_summary-reads
	# fi

	# _run_polap_total-length-long
	# _run_polap_find-genome-size
	# _run_polap_reduce-data

	if [ -s "${_polap_var_base_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_base_msbwt_tar_gz} ... later when we polish it with the short-read data."
			# tar zxf "${_polap_var_base_msbwt_tar_gz}"
		else
			_polap_log1 "  Do the preparation of short-read polishing ... early"
			check_file_existence "${SR1}"
			check_file_existence "${SR2}"
			_run_polap_prepare-polishing
		fi
	fi

	_run_polap_assemble-draft

	if [ -s "${_polap_var_base_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_base_msbwt_tar_gz} ..."
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
		check_file_existence "${_polap_var_base_msbwt}"

		if [ -s "${PA}" ]; then
			if [ -s "${FA}" ]; then
				_polap_log1 "  skipping the short-read polishing ..."
			else
				_polap_log1 "  mtDNA is being polished for select-contig type ${i}..."
				_run_polap_polish
			fi
		else
			_polap_log1 "  no mtDNA candidate for select-contig type ${i}..."
		fi
	done

	# Download the known reference mtDNA sequence in fasta format if available
	# -o PRJNA914763
	# INUM=0
	# if [ -s "${_polap_var_bioproject_mtdna_fasta2}" ]; then
	# 	_polap_log2 "  skipping downloading mtDNA from NCBI"
	# 	_polap_log2 "  we use the mtDNA: ${_polap_var_bioproject_mtdna_fasta2}"
	# else
	# 	_run_polap_get-mtdna
	#
	# fi

	if [[ -s "${_polap_var_bioproject_mtdna_fasta2}" ]]; then

		# Compare the known mtDNA and the assembled one.
		for i in "${_arg_select_contig_numbers[@]}"; do
			# Define the paths for mtDNA sequences to be polished
			FA="${ODIR}/${i}/mt.1.fa"

			if [ -s "${FA}" ]; then
				INUM="${i}"
				_run_polap_compare-mtdna
			else
				_polap_log1 "  skipping the short-read polishing ..."
				INUM="${i}"
				source "$script_dir/polap-variables-oga.sh" # '.' means 'source'
				local n1=$(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")
				local l1="0"
				local l2="0"
				local c1="0"
				printf "%s\t%d\t%d\t%f\n" ${n1} ${l1} ${l2} ${c1} >"${_polap_var_mtdna_compare}"
			fi
		done
	else
		_polap_log0 "No known mtDNA for $(<${_polap_var_bioproject_species})"
	fi

	# Finalize the assemble-bioproject function.
	_polap_log1 "rsync back to thorne ..."
	touch "${ODIR}/log-assembled.txt"
	# touch "$BIOPRJ/log-copying-to-xxx.txt"
	# touch "$BIOPRJ/log-not-started-yet.txt"
	if [ "$(hostname)" = "thorne" ]; then
		if [[ "${PWD}" = "/home/goshng/run/polap" ]]; then
			cp -pr "${ODIR}/" "/home/goshng/all/polap/figshare/bioprojects/"
			touch "/home/goshng/all/polap/figshare/bioprojects/log-assembled-${ODIR}.txt"
		fi
	else
		touch "log-assembled-${ODIR}.txt"
		scp -pq "log-assembled-${ODIR}.txt" thorne:"$PWD/"
		rsync -a -e ssh "${ODIR}/" thorne:"$PWD/${ODIR}/"
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
