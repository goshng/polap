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
# FIXME: something are complicated. Do we need it?
################################################################################
function _run_polap_x-bioproject2() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	mkdir $ODIR
	BIOPRJ=$_arg_bioproject

	MDATA="$ODIR/accessions_demo.txt"
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
	done >$ODIR/accessions_demo.tab
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
function _run_polap_x-ncbi-fetch-sra() {
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
	"$script_dir"/run-polap-ncbitools fetch sra "$SRA"

	echoerr You have a file called "$SRA".fastq and a folder named "$SRA"
	echoerr if your download try is successful. Then, you would want to delete
	echoerr the folder because we need only the fastq file.

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

################################################################################
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --sra SRR10190639
################################################################################
function _run_polap_x-ncbi-fetch-sra-runinfo() {
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

################################################################################
# Fetches mtDNA genome sequence by species name.
# Arguments:
#   --species
################################################################################
function _run_polap_x-ncbi-fetch-mtdna-genbank() {
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

################################################################################
# Fetches mtDNA genome sequence by accession.
# Arguments:
#   --accession
# Inputs:
#   accession ID
# Outputs:
#   <accession>.fa
################################################################################
function _run_polap_x-ncbi-fetch-mtdna-nucleotide() {
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
		echoall "NEXT: $(basename "$0") align-two-dna-sequences --query mt.1.fa --subject ${_arg_accession}.fa"
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
function _run_polap_x-align-two-dna-sequences() {
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
function _run_polap_x-clustal() {
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
#   $LRNK
#   $PA
# Outputs:
#   $FA
################################################################################
function _run_polap_x-check-coverage() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	LRNK="$ODIR/nk.fq.gz"

	help_message=$(
		cat <<HEREDOC
# NOT IMPLEMENTED YET
# Checks the coverage on a draft genome using samtools.
# Arguments:
#   -p $PA: a draft genome
# Inputs:
#   $LRNK
#   $PA
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-p|--unpolished-fasta <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ -z "$PA" ] && [ -z "$LRNK" ]; then
		echoerr "ERROR: no -p option are used."
		echoerr "INFO: --p mt.0.fasta"
	else
		echo "INFO: executing minimap2 and samtools for checking the long reads coverage on the $PA ... be patient!"
		if [[ ! -s "$PA" ]]; then
			echoall "ERROR: no such file $PA"
			exit $EXIT_ERROR
		fi

		if [[ ! -s "$LRNK" ]]; then
			echoall "ERROR: no such file $LRNK"
			exit $EXIT_ERROR
		fi

		minimap2 -t "$NT" -ax map-ont "$PA" "$LRNK" 2>/dev/null |
			samtools view -u 2>/dev/null |
			samtools sort -o "$ODIR"/1.bam \
				>/dev/null 2>&1
		samtools coverage -A -w 32 "$ODIR"/1.bam 1>&2

		_polap_log1 INFO: conda env create -f "$WDIR"/environment-fmlrc.yaml
		_polap_log1 INFO: conda activate polap-fmlrc
		_polap_log1 "NEXT: $(basename "$0") prepare-polishing [-a s1.fq] [-b s2.fq]"
	fi

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}

function _run_polap_x-help() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	print_x-help 1>&2

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}