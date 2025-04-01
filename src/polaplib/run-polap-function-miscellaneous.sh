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

function _run_polap_report-table {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Report info.
#
Example: $(basename "$0") ${_arg_menu[0]} Brassica_rapa
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	local _project="${_arg_menu[1]}"
	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		local _project=$(basename $PWD)
	fi
	_polap_log0 "Project: ${_project}"

	local _l=$(basename "$(readlink -f l.fq)" .fastq)
	_polap_log0 "Long SRA: ${_l}"
	local _s=$(basename "$(readlink -f s1.fq)" .fastq)
	_s=${_s%_1}
	_polap_log0 "Short SRA: ${_s}"

	_polap_utility_get_contig_length "${_arg_outdir}/1/mt.0.fasta" "${_arg_outdir}/1/mt.0.fasta.len"
	local _l=$(<"${_arg_outdir}/1/mt.0.fasta.len")
	_polap_log0 "mtDNA total length: ${_l}"

	_arg_menu[1]=$(basename "${_project}" | sed 's/_/ /g')
	_run_polap_get-taxonomy-species

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_get-taxonomy-species {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Get taxonomy.
#
Example: $(basename "$0") ${_arg_menu[0]} "Brassica rapa"
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	local _species="${_arg_menu[1]}"
	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		local _species=$(basename "$PWD" | sed 's/_/ /g')
	fi

	MDATA="${_arg_outdir}/00-bioproject/taxonomy.txt"
	_polap_log0 "species: ${_species}"
	local _v=$(esearch -db taxonomy -query "${_species}" | efetch -format xml | xtract -pattern Taxon -block "LineageEx/Taxon" -sep ";" -element ScientificName)
	_polap_log0 "taxonomy: ${_v}"

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

################################################################################
# FIXME: something are complicated. Do we need it?
################################################################################
function _run_polap_x-bioproject2 {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	mkdir ${_arg_outdir}
	# BIOPRJ=$_arg_bioproject

	MDATA="${_arg_outdir}/accessions_demo.txt"
	esearch -db bioproject -query $BIOPRJ | elink -target biosample | efetch -format docsum | xtract.Linux -pattern DocumentSummary -block Ids -element Id -group SRA >${MDATA}
	SRSs=$(cat ${MDATA} | while read LINE; do
		NCOL=$(echo ${LINE} | wc -w)
		ACC=$(echo ${LINE} | cut -d ' ' -f ${NCOL})
		echo ${ACC}
	done | xargs)

	# SRS_str=$(join_by ' OR ' ${SRSs[@]})
	# echoall $SRS_str
	# esearch -db SRA -query "\"" $SRS_str "\"" | efetch -format runinfo | csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName | csvtk pretty 1>&2

	for SRS in ${SRSs[@]}; do
		esearch -db SRA -query ${SRS} | efetch -format runinfo
		# esearch -db SRA -query ${SRS} | efetch -format runinfo | tail -n +2
	done >${_arg_outdir}/accessions_demo.tab
	# csvtk cut -f Run,bases,LibraryName,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName | csvtk pretty 1>&2

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

###############################################################################
# Feteches SRA data file.
# Arguments:
#   --sra SRR10190639
# Outputs:
#   SRR10190639.fastq
###############################################################################
function _run_polap_x-ncbi-fetch-sra {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Feteches SRA data file.
# Arguments:
#   --sra SRR10190639
# Outputs:
#   SRR10190639.fastq
Example: $(basename $0) ${_arg_menu[0]} --sra <arg>
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if ! run_check1; then
		error_polap_conda
		exit $EXIT_ERROR
	fi

	if [ -z "$_arg_sra" ]; then
		echoerr "ERROR: no --sra option is used."
		exit $EXIT_SUCCESS
	fi

	SRA=$_arg_sra
	_polap_log0 fetching "$SRA".fastq from NCBI SRA database ... please wait ...
	"${_POLAPLIB_DIR}"/run-polap-ncbitools fetch sra "$SRA"

	_polap_log1 You have a file called "$SRA".fastq and a folder named "$SRA"
	_polap_log1 if your download try is successful. Then, you would want to delete
	_polap_log1 the folder because we need only the fastq file.

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

################################################################################
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --sra SRR10190639
################################################################################
function _run_polap_x-ncbi-fetch-sra-runinfo {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --sra SRR10190639
# Outputs:
#   bases
Example: $(basename $0) ${_arg_menu[0]} --sra <arg>
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "$_arg_sra" ]; then
		echoerr "ERROR: no --sra option is used."
		exit $EXIT_SUCCESS
	fi

	echoerr "counting bases in SRA database and its downloaded FASTQ files of SRA [${_arg_sra}] ..."
	bases=$(esearch -db sra -query "${_arg_sra}" |
		efetch -format runinfo |
		csvtk cut -f "bases" |
		csvtk del-header)
	echoall "SRA: ${_arg_sra}: ${bases} (bp)"

	bases2=$(seqkit stats -T "${_arg_sra}"*.fastq 2>/dev/null | csvtk cut -t -f "sum_len" | paste -s -d+ - | bc)
	echoall "FASTQ: ${_arg_sra}*.fastq: ${bases2} (bp)"

	if [ "${bases}" -eq "${bases2}" ]; then
		echoall "LOG: SRA ${_arg_sra} and its FASTQ files: total bases match."
	else
		echoall "ERROR: SRA ${_arg_sra} and its FASTQ files: total bases do not match."
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_x-get-sra-info {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
Get SRA info.

esearch -db sra -query <SRA> | efetch -format runinfo | csvtk cut -f Run | csvtk pretty

Arguments:
--sra SRR10190639

Outputs:
bases

Example: $(basename $0) ${_arg_menu[0]} --sra <arg>
HEREDOC
	)

	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0

	if [ -z "$_arg_sra" ]; then
		_polap_log0 "ERROR: no --sra option"
		return 0
	fi

	esearch -db sra -query "${_arg_sra}" |
		efetch -format runinfo |
		csvtk cut -f Run,bases,LibraryStrategy,LibrarySource,LibraryLayout,Platform,ScientificName |
		csvtk pretty >&3

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_x-check-disk-space {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	local _polap_return=${RETURN_SUCCESS}

	help_message=$(
		cat <<HEREDOC
# Roughly compute the space for a POLAP analysis.
# Arguments:
#   long-read sra ID
#   short-read sra ID
Example: $(basename $0) ${_arg_menu[0]} SRR10190639 SRR10250248
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if [[ "${_arg_menu[1]}" == "infile" ]] || [[ "${_arg_menu[2]}" == "outfile" ]]; then
		_polap_echo0 "${help_message}"
		return
	fi

	local _sra1="${_arg_menu[1]}"
	local _sra2="${_arg_menu[2]}"

	local _bases1=$(esearch -db sra -query "${_sra1}" |
		efetch -format runinfo |
		csvtk cut -f "bases" |
		csvtk del-header)

	local _bases2=$(esearch -db sra -query "${_sra2}" |
		efetch -format runinfo |
		csvtk cut -f "bases" |
		csvtk del-header)

	local _available_disk_space=$(df --block-size=1 --output=avail ${_arg_outdir} | tail -1)

	local _required_disk_space=$(((_bases1 + _bases2) * 5))

	local _required_disk_space_bp=$(_polap_utility_convert_bp ${_required_disk_space})
	_polap_log0 "required disk space: ${_required_disk_space_bp}"
	local _available_disk_space_bp=$(_polap_utility_convert_bp ${_available_disk_space})
	_polap_log0 "available disk space: ${_available_disk_space_bp}"

	if ((_required_disk_space > _available_disk_space)); then
		_polap_return=${_POLAP_ERR_NO_DISK_SPACE}
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	return ${_polap_return}
}

function _run_polap_x-link-sra {
	# Loop through all .fastq files in the current directory
	for file in *.fastq; do
		# Check if the file is a long-read (ending in just .fastq without _1 or _2)
		if [[ "$file" =~ ^[DES]RR[0-9]+\.fastq$ ]]; then
			ln -s "$file" l.fq
			_polap_log0 "Created soft link for long-read file: $file -> l.fq"
		# Check if the file is a short-read pair _1.fastq
		elif [[ "$file" =~ ^[DES]RR[0-9]+_1\.fastq$ ]]; then
			ln -s "$file" s1.fq
			_polap_log0 "Created soft link for short-read file 1: $file -> s1"
		# Check if the file is a short-read pair _2.fastq
		elif [[ "$file" =~ ^[DES]RR[0-9]+_2\.fastq$ ]]; then
			ln -s "$file" s2.fq
			_polap_log0 "Created soft link for short-read file 2: $file -> s2"
		fi
	done
}

################################################################################
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --species
################################################################################
function _run_polap_x-ncbi-fetch-mtdna-genbank {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --species species-name
# Inputs:
#   species-name
# Outputs:
#   species-name.mt.gb
Example: $(basename "$0") ${_arg_menu[0]} --species <arg>
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "$_arg_species" ]; then
		echoerr "ERROR: no --species option is used."
	else
		S="${_arg_species// /-}"
		esearch -db nuccore -query "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${_arg_species}[Organism]" |
			efetch -format gb >"${S}".mt.gb
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_x-ncbi-fetch-cpdna-genbank {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --species species-name
# Inputs:
#   species-name
# Outputs:
#   species-name.mt.gb
Example: $(basename "$0") ${_arg_menu[0]} --species <arg>
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "$_arg_species" ]; then
		echoerr "ERROR: no --species option is used."
	else
		S="${_arg_species// /-}"
		esearch -db nuccore -query "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${_arg_species}[Organism]" |
			efetch -format gb >"${S}".mt.gb
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}
################################################################################
# Fetches mtDNA genome sequence by accession.
# Arguments:
#   --accession
# Inputs:
#   accession ID
# Outputs:
#   <accession>.fa
################################################################################
function _run_polap_x-ncbi-fetch-mtdna-nucleotide {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --accession <arg>
# Inputs:
#   accession ID
# Outputs:
#   <accession>.fa
Example: $(basename "$0") ${_arg_menu[0]} --accession <arg>
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "$_arg_accession" ]; then
		echoerr "ERROR: no --accession option is used."
	else
		esearch -db nuccore -query "${_arg_accession}[ACCN]" </dev/null |
			efetch -format fasta >"${_arg_accession}".fa
		_polap_log1 "NEXT: $(basename "$0") align-two-dna-sequences --query mt.1.fa --subject ${_arg_accession}.fa"
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

################################################################################
# Uses BLAST to align two very similar DNA sequences.
# Before align the two sequences, rearrange the assembled query sequences
# so that they are alignable using multiple sequence alignment tools
# like ClustalW.
# Arguments:
#   --query ${_arg_query}
#   --subject ${_arg_subject}
# Inputs:
#   query: mt.1.fa or the assembled sequence
#   subject: a known mtDNA sequence
# Outputs:
#   pairwise-alignment.txt
################################################################################
function _run_polap_x-align-two-dna-sequences {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Uses BLAST to align two very similar DNA sequences.
# Before align the two sequences, rearrange the assembled query sequences
# so that they are alignable using multiple sequence alignment tools
# like ClustalW.
# Arguments:
#   --query ${_arg_query}
#   --subject ${_arg_subject}
# Inputs:
#   query: mt.1.fa or the assembled sequence
#   subject: a known mtDNA sequence
# Outputs:
#   pairwise-alignment.txt
Example: $(basename "$0") ${_arg_menu[0]} [--query <arg>] [--subject <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "$_arg_query" -a -z "$_arg_subject" ]; then
		echoerr "ERROR: no --query and --subject option are used."
		echoerr "INFO: --query mt.1.fa --subject NCBI.fa"
		echoerr "INFO: seqkit seq -p -r mt.1.fa -o mt.1r.fa"
		echoerr "INFO: seqkit restart -i <POS> mt.1.fa -o mt.2.fa"
	else
		blastn -query "$_arg_query" -subject "$_arg_subject" >pairwise-alignment.txt
		echoerr "INFO: seqkit seq -p -r mt.1.fa -o mt.1r.fa"
		echoerr "INFO: seqkit restart -i <POS> mt.1.fa -o mt.2.fa"
		echoerr "INFO: $(basename "$0") clustal --query mt.2.fa --subject NCBI-ACC.fa"
		echoerr see pairwise-alignment.txt
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

################################################################################
# Uses clustalw to align two DNA sequences.
# After rearranging parts of a query assembled sequence,
# we use clustalw to align the two to check how similar they are.
# Arguments:
#   --query ${_arg_query}
#   --subject ${_arg_subject}
# Inputs:
#   query: mt.1.fa or the assembled sequence
#   subject: a known mtDNA sequence
################################################################################
function _run_polap_x-clustal {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Uses clustalw to align two DNA sequences.
# After rearranging parts of a query assembled sequence,
# we use clustalw to align the two to check how similar they are.
# Arguments:
#   --query ${_arg_query}
#   --subject ${_arg_subject}
# Inputs:
#   query: mt.1.fa or the assembled sequence
#   subject: a known mtDNA sequence
Example: $(basename "$0") ${_arg_menu[0]} [--query <arg>] [--subject <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [[ -s x.aln ]]; then
		echoall "INFO: clustalw alignemnt is found."
	else
		local MT1=$_arg_query
		shift
		local MT2=$_arg_subject
		shift
		command -v clustalw2 >/dev/null 2>&1 || {
			echo >&2 "clustalw2: not installed"
			exit 1
		}
		cat "$MT1" "$MT2" >x.fa
		clustalw2 x.fa
	fi
	# TODO: works for only Linux.
	GAP1=$(sed -n '4~4p' x.aln | tr -s ' ' | cut -d' ' -f2 | tr -d '\n' | tr -dc '-' | wc -c)
	GAP2=$(sed -n '5~4p' x.aln | tr -s ' ' | cut -d' ' -f2 | tr -d '\n' | tr -dc '-' | wc -c)
	MATCH=$(sed -n '6~4p' x.aln | tr -dc '*' | wc -c)
	TOTAL=$(sed -n '4~4p' x.aln | tr -s ' ' | cut -d' ' -f2 | tr -d '\n' | wc -c)
	GAP=$((GAP1 + GAP2))
	MISMATCH=$((TOTAL - MATCH - GAP))
	echoall "INFO: pairwise alignment:   length: $TOTAL"
	echoall "INFO: pairwise alignment:    match: $MATCH"
	echoall "INFO: pairwise alignment: mismatch: $MISMATCH"
	echoall "INFO: pairwise alignment:     gaps: $GAP"
	LENGTH1=$(seqkit stats -Ta "$_arg_subject" | csvtk cut -t -f "sum_len" | csvtk del-header)
	LENGTH2=$(seqkit stats -Ta "$_arg_query" | csvtk cut -t -f "sum_len" | csvtk del-header)
	NAME1=$(seqkit seq -n -i "$_arg_subject")
	NAME2=$(seqkit seq -n -i "$_arg_query")
	echoall "INFO: Sequence1: $NAME1"
	echoall "INFO: SequenceLen1: $LENGTH1"
	echoall "INFO: Sequence2: $NAME2"
	echoall "INFO: SequenceLen2: $LENGTH2"

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

################################################################################
# Checks the coverage using samtools.
# Arguments:
#   -p mt.0.fasta
# Inputs:
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
################################################################################
function _run_polap_x-check-coverage {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	LRNK="${_arg_outdir}/nk.fq.gz"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# NOT IMPLEMENTED YET
# Checks the coverage on a draft genome using samtools.
# Arguments:
#   -p ${_arg_unpolished_fasta}: a draft genome
# Inputs:
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_arg_unpolished_fasta}
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-p|--unpolished-fasta <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "${_arg_unpolished_fasta}" ] && [ -z "${_polap_var_outdir_nk_fq_gz}" ]; then
		echoerr "ERROR: no -p option are used."
		echoerr "INFO: --p mt.0.fasta"
	else
		echo "INFO: executing minimap2 and samtools for checking the long reads coverage on the ${_arg_unpolished_fasta} ... be patient!"
		if [[ ! -s "${_arg_unpolished_fasta}" ]]; then
			echoall "ERROR: no such file ${_arg_unpolished_fasta}"
			exit $EXIT_ERROR
		fi

		if [[ ! -s "${_polap_var_outdir_nk_fq_gz}" ]]; then
			echoall "ERROR: no such file ${_polap_var_outdir_nk_fq_gz}"
			exit $EXIT_ERROR
		fi

		minimap2 -t "${_arg_threads}" -ax map-ont "${_arg_unpolished_fasta}" "${_polap_var_outdir_nk_fq_gz}" 2>/dev/null |
			samtools view -u 2>/dev/null |
			samtools sort -o "${_arg_outdir}"/1.bam \
				>/dev/null 2>&1
		samtools coverage -A -w 32 "${_arg_outdir}"/1.bam 1>&2

		_polap_log1 INFO: conda env create -f "${_POLAPLIB_DIR}"/environment-fmlrc.yaml
		_polap_log1 INFO: conda activate polap-fmlrc
		_polap_log1 "NEXT: $(basename "$0") prepare-polishing [-a s1.fq] [-b s2.fq]"
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_x-help {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	print_x-help 1>&2

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_x-link-fastq {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	# Loop through all .fastq files in the current directory
	for file in *.fastq *.fq; do
		# Check if the file is a long-read (ending in just .fastq without _1 or _2)
		if [[ "$file" =~ ^[DES]RR[0-9]+\.fastq$ ]]; then
			ln -s "$file" l.fq
			_polap_log0 "Created soft link for long-read file: $file -> l.fq"
		# Check if the file is a short-read pair _1.fastq
		elif [[ "$file" =~ ^[DES]RR[0-9]+_1\.fastq$ ]]; then
			ln -s "$file" s1.fq
			_polap_log0 "Created soft link for short-read file 1: $file -> s1"
		# Check if the file is a short-read pair _2.fastq
		elif [[ "$file" =~ ^[DES]RR[0-9]+_2\.fastq$ ]]; then
			ln -s "$file" s2.fq
			_polap_log0 "Created soft link for short-read file 2: $file -> s2"
		fi
	done

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	return 0
}

function _run_polap_x-prepend-gplv3 {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	help_message=$(
		cat <<HEREDOC
# Prepend GPLv3 copyright to a file.
Example: $(basename "$0") ${_arg_menu[0]} <file>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Define your files
	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "ERROR: $(basename "$0") ${_arg_menu[0]} <file>"
		return $RETURN_SUCCESS
	fi

	local file1="${_arg_menu[1]}" # File to prepend

	# Create a temporary file that will hold the combined content
	local temp_file=$(mktemp)

	# Combine the content of file1 and file2 into the temporary file
	cat <(head -16 "${_POLAPLIB_DIR}/run-polap-function-template.sh") "$file1" >"$temp_file"

	# Replace file2 with the new combined content
	mv "$temp_file" "$file1"

	_polap_log0 "$file1" now has a header of the POLAP GPLv3 copyright.

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	return 0
}

# copy seed test data
function _run_polap_x-copy-seed-test {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local _id=$(basename $PWD)
	local _base_name="${_polap_var_ga_annotation_all##*/}"
	local _base_name="${_base_name%.*}"
	_polap_log0_cmd cp "${_polap_var_ga_annotation_all}" "$HOME/all/polap/github/trash/${_base_name}_${_id}.txt"

	echo "input1 <- file.path(input_dir0, \"${_base_name}_${_id}.txt\")" >>"${_POLAPLIB_DIR}/run-polap-r-determine-depth-range.txt"

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	return 0
}

function _polap_gfatools-gfa2fasta {
	if [[ ! -s "${_polap_var_ga_contigger_edges_fasta}" ]]; then
		if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
			_polap_log3_pipe "gfatools gfa2fa \
		    ${_polap_var_ga_contigger_edges_gfa} \
		    >${_polap_var_ga_contigger_edges_fasta} \
        2>${_polap_output_dest}"
		else
			return ${_POLAP_ERR_NO_EDGES_GFA}
		fi
	fi
}

function _run_polap_get-revision1 {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	cp -p "${_POLAPLIB_DIR}/polap-revision1.sh" .

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
	return 0
}
