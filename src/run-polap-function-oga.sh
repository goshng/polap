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
# Selects reads mapped on a genome assembly.
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -r $MPAIR: minimum minimap2 alignment length for a pair of contigs
#   -x $MBRIDGE: minimum long-read length for connecting the pair of contigs
#   -w $MSINGLE: minimum minimap2 alignment length for a single contig
# Inputs:
#   $MTCONTIGNAME
#   ${assembly_graph_final_fasta}
# Outputs:
#   $MTSEEDSDIR
#   $MTDIR/contig.fa
#   $MTDIR/contig_total_length.txt
#   $MTDIR/contig.paf
#   $MTDIR/contig.tab
#   $MTSEEDSDIR/1.names
#   $MTSEEDSDIR/1.fq.gz
#   $MTSEEDSDIR/2.fq.gz
################################################################################
function _run_polap_select-reads() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	LRNK="$ODIR/nk.fq.gz"
	MR=$_arg_min_read_length
	FDIR="$ODIR"/$INUM
	ADIR="$FDIR"/50-annotation
	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="$MTDIR"/seeds

	MTCONTIGNAME="$FDIR"/mt.contig.name-$JNUM

	# for contigs
	#	assembly_graph_final_fasta=o/30-contigger/contigs.fasta
	#	for edges
	assembly_graph_final_fasta="$FDIR"/30-contigger/graph_final.fasta

	help_message=$(
		cat <<HEREDOC
# Selects reads mapped on a genome assembly.
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -r $MPAIR: minimum minimap2 alignment length for a pair of contigs
#   -x $MBRIDGE: minimum long-read length for connecting the pair of contigs
#   -w $MSINGLE: minimum minimap2 alignment length for a single contig
# Inputs:
#   $MTCONTIGNAME
#   ${assembly_graph_final_fasta}
# Outputs:
#   $MTDIR/contig.fa
#   $MTSEEDSDIR/1.names
#   $MTSEEDSDIR/2.fq.gz
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [-r|--pair-min <arg>] [-x|--bridge-min <arg>] [-w|--single-min <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "$MTCONTIGNAME" ]; then
		echoall "ERROR: no such mt.contig.name file: $MTCONTIGNAME"
		exit $EXIT_ERROR
	fi

	if [ ! -s "${assembly_graph_final_fasta}" ]; then
		echoall "ERROR: no assembly fasta file: ${assembly_graph_final_fasta}"
		exit $EXIT_ERROR
	fi

	if [ -d "$MTDIR" ] && [ "${_arg_yes}" = "off" ]; then
		while true; do
			read -r -p "Folder [$MTDIR] already exists. Do you want to replace it? [y/n] " yn
			case $yn in
			[Yy]*)
				rm -rf "$MTDIR"
				echo "INFO: $MTDIR is deleted."
				break
				;;
			[Nn]*)
				echoall "INFO: [$MTDIR] is not replaced."
				echoerr "you might want a new mt.contig.name file for flye2 step."
				exit $EXIT_FAIL
				;;
			*) echo "Please answer yes or no." ;;
			esac
		done
	else
		rm -rf "$MTDIR"
		echo "INFO: $MTDIR is deleted if there is one."
	fi

	echo "$CMD" >"$ODIR"/organelle-assembly_"${INUM}"-"${JNUM}"

	mkdir -p "$MTSEEDSDIR"
	ln -s "$PWD"/"$ODIR"/nk.fq.gz -t "$MTDIR"
	echoall "INFO: extracts contigs from the assembly ${assembly_graph_final_fasta}"
	echo "INFO: uses mt.contig.name at $MTCONTIGNAME"
	echoerr "FILE: reading contigs in $MTCONTIGNAME"
	cat "$MTCONTIGNAME" >&2
	echoerr "please, wait for a long-read data selection ... $INUM -> $JNUM ... bridge=$MBRIDGE p_mapping=$MPAIR s_mapping=$MSINGLE min_len_read=$MR"
	seqkit grep --threads "$NT" -f "$MTCONTIGNAME" "${assembly_graph_final_fasta}" -o "$MTDIR"/contig.fa >/dev/null 2>&1

	contig_count=$(wc -l <"$MTCONTIGNAME")
	if [[ $CIRCULARIZE == "on" ]]; then
		if [ "$contig_count" -eq 1 ]; then
			seqkit fx2tab --length --name "$MTDIR"/contig.fa -o "$MTDIR"/contig.fa.len >/dev/null 2>&1
			A=$(cut -f2 "$MTDIR"/contig.fa.len)
			B=$(echo "scale=0; $A/2" | bc)
			C=$((B + 1))
			seqkit subseq -r 1:"$B" "$MTDIR"/contig.fa -o "$MTDIR"/c1.fa >/dev/null 2>&1
			seqkit subseq -r "$C":"$A" "$MTDIR"/contig.fa -o "$MTDIR"/c2.fa >/dev/null 2>&1
			cat "$MTDIR"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "$MTDIR"/contig.fa >/dev/null 2>&1
			cp "$MTCONTIGNAME" "$MTCONTIGNAME"-backup
			echo -e "edge_1\nedge_2" >"$MTCONTIGNAME"
			echo "INFO: creates new $MTDIR/contig.fa and $MTCONTIGNAME"
		else
			echo "DEV: not implemented yet"
			exit $EXIT_ERROR
			# "$WDIR"/run-polap-single.R "$MTSEEDSDIR"/contig.tab "$MTSEEDSDIR" "$MSINGLE" >/dev/null 2>&1
			# cat "$MTSEEDSDIR"/single.names | sort | uniq >"$MTSEEDSDIR"/1.names
			# echo "INFO: creates long read single name in $MTSEEDSDIR/1.names"
		fi
	fi
	echo "DATA: $MTDIR/contig.fa is created."

	CONTIG_LENGTH=$(seqkit stats -Ta "$MTDIR"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	echo "$CONTIG_LENGTH" >"$MTDIR"/contig_total_length.txt
	echo "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"

	if [[ -s "$MTDIR"/contig.paf ]]; then
		echo "DATA: previously created $MTDIR/contig.paf is used without executing minimap2."
	else
		minimap2 -cx map-ont "$MTDIR"/contig.fa "$LRNK" -t "$NT" -o "$MTDIR"/contig.paf >/dev/null 2>&1
		echo "DATA: $LRNK is used to select reads."
		echo "DATA: $MTDIR/contig.paf is created."
	fi

	cut -f1-11 "$MTDIR"/contig.paf | awk -v minlength="$MR" '{if ($2>=minlength) {print}}' >"$MTDIR"/contig.tab
	echo "DATA: minimum length of long reads in the read selection: $MR"
	echo "DATA: $MTDIR/contig.tab is created."

	# MT: MPAIR=3000 MBRIDGE=3000 MSINGLE=3000
	# PT: MPAIR=1000 MBRIDGE=5000 MSINGLE=0
	Rscript --vanilla "$WDIR"/run-polap-pairs.R "$MTCONTIGNAME" "$MTDIR"/contig.tab "$MTSEEDSDIR" "$MPAIR" "$MBRIDGE" "$MSINGLE" >/dev/null 2>&1
	# "$WDIR"/run-polap-pairs.R "$MTCONTIGNAME" $MTDIR/contig.tab $MTSEEDSDIR $MPAIR $MBRIDGE $MSINGLE >/dev/null 2>&1
	echo "OPTION polap pairs alignment minimum: $MPAIR"
	echo "OPTION polap pairs bridge minimum: $MBRIDGE"
	echo "DATA: pair contig names in $MTSEEDSDIR are created."
	echo "DATA: single contig name in $MTSEEDSDIR is created."

	# cat "$MTSEEDSDIR/"*".name" "$MTSEEDSDIR"/single.names | sort | uniq >"$MTSEEDSDIR"/1.names
	cat "$MTSEEDSDIR"/single.names | sort | uniq >"$MTSEEDSDIR"/1.names
	echo "INFO: creates long read names and the single name in $MTSEEDSDIR/1.names"

	# seqkit grep --threads $NT -f "$MTSEEDSDIR"/1.names $LRNK -o "$MTSEEDSDIR"/1.fq.gz >/dev/null 2>&1
	seqtk subseq "$LRNK" "$MTSEEDSDIR"/1.names | gzip >"$MTSEEDSDIR"/1.fq.gz
	echo "DATA: organelle reads in $MTSEEDSDIR/1.fq.gz"

	TOTAL_LENGTH=$(seqkit stats -Ta "$MTSEEDSDIR"/1.fq.gz | csvtk cut -t -f "sum_len" | csvtk del-header)
	EXPECTED_ORGANELLE_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
	echo "INFO: expected coverage: ${EXPECTED_ORGANELLE_COVERAGE}x"

	if [[ "${_arg_test}" == "on" ]]; then
		echoall "OPTION: --test : No reduction of the test long-read data"
		ln -s "$(realpath "$MTSEEDSDIR"/1.fq.gz)" "$MTSEEDSDIR"/2.fq.gz
	elif [[ "${_arg_coverage_check}" == "off" ]]; then
		echoall "OPTION: --no-coverage-check : No reduction of the long-read data"
		ln -s "$(realpath "$MTSEEDSDIR"/1.fq.gz)" "$MTSEEDSDIR"/2.fq.gz
	else
		if [ "$EXPECTED_ORGANELLE_COVERAGE" -lt "$COV" ]; then
			echoall "LOG: No reduction of the long-read data because $EXPECTED_ORGANELLE_COVERAGE < $COV"
			ln -s "$(realpath "$MTSEEDSDIR"/1.fq.gz)" "$MTSEEDSDIR"/2.fq.gz
		else
			echoall "SUGGESTION: you might want to increase the minimum read lengths because you have enough long-read data."
			RATE=$(echo "scale=10; $COV/$EXPECTED_ORGANELLE_COVERAGE" | bc)
			echoall "LOG: long-read data reduction by rate of $RATE <= COV[$COV] / long-read organelle coverage[$EXPECTED_ORGANELLE_COVERAGE]"
			echoall "sampling long-read data by $RATE ... wait ..."
			# seqkit sample -p "$RATE" "$MTSEEDSDIR/1.fq.gz" -o "$MTSEEDSDIR/2.fq.gz" >/dev/null 2>&1
			local seed=${_arg_seed:-$RANDOM}
			_polap_log3_cmd "seqkit sample -p ${RATE} -s ${seed} ${MTSEEDSDIR}/1.fq.gz -o ${MTSEEDSDIR}/2.fq.gz 2>${_polap_output_dest}"
			echoall "DATA: a reduced long-read data $MTSEEDSDIR/2.fq.gz is created"
		fi
	fi

	C=$(ls -1 "$MTSEEDSDIR/"*".name" 2>/dev/null | wc -l)
	if [ "$C" != 0 ]; then
		echo "INFO: bridging reads exist: combinations of $C."
		seqkit seq -n -i "$MTSEEDSDIR"/2.fq.gz >"$MTSEEDSDIR"/single.names.2
		cat "$MTSEEDSDIR/"*".name" "$MTSEEDSDIR"/single.names.2 | sort | uniq >"$MTSEEDSDIR"/1.names.2
		# seqkit grep --threads $NT -f "$MTSEEDSDIR"/1.names.2 $LRNK -o "$MTSEEDSDIR"/2.fq.gz >/dev/null 2>&1
		seqtk subseq "$LRNK" "$MTSEEDSDIR"/1.names.2 | gzip >"$MTSEEDSDIR"/2.fq.gz
	fi
	echoall "DATA: organelle reads in $MTSEEDSDIR/2.fq.gz"

	# put the backup to the original
	if [[ $CIRCULARIZE == "on" ]]; then
		if [[ -s "$MTCONTIGNAME"-backup ]]; then
			mv "$MTCONTIGNAME"-backup "$MTCONTIGNAME"
		else
			echo "DEV: not implemented yet"
			exit $EXIT_ERROR
		fi
	fi

	echoerr NEXT: "$(basename "$0")" flye2 -o "$ODIR" -j "$JNUM" -t "$NT" -c "$COV"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Executes Flye for an organelle-genome assembly
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $MTDIR/contig.fa
#   $MTSEEDSDIR/2.fq.gz
# Outputs:
#   $MTDIR/contig_total_length.txt
#   $MTDIR/30-contigger/contigs.fasta
#   $MTDIR/30-contigger/contigs_stats.txt
#   $MTDIR/30-contigger/graph_final.fasta
#   $MTDIR/30-contigger/graph_final.gfa
################################################################################
function _run_polap_flye2() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	echo "INFO: organelle-genome assembly on $JNUM"

	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="$MTDIR"/seeds

	help_message=$(
		cat <<HEREDOC
# Executes Flye for an organelle-genome assembly
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $MTDIR/contig.fa
#   $MTSEEDSDIR/2.fq.gz
# Outputs:
#   $MTDIR/contig_total_length.txt
#   $MTDIR/30-contigger/contigs.fasta
#   $MTDIR/30-contigger/contigs_stats.txt
#   $MTDIR/30-contigger/graph_final.fasta
#   $MTDIR/30-contigger/graph_final.gfa
Example: $(basename $0) ${_arg_menu[0]} [-j|--jnum <arg>] [-t|--threads <arg>] [-c|--coverage <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "$MTDIR/contig.fa" ]; then
		echoall "ERROR: no selected-contig file [$MTDIR/contig.fa]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "$MTSEEDSDIR/2.fq.gz" ]; then
		echoall "ERROR: no long-read file [$MTSEEDSDIR/2.fq.gz]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	CONTIG_LENGTH=$(seqkit stats -Ta "$MTDIR"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	echo "$CONTIG_LENGTH" >"$MTDIR"/contig_total_length.txt
	echo "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"

	echo "INFO: executing the organelle-genome assembly using flye ... be patient!"
	echoerr "please, wait for an organelle-genome assembly on $JNUM ..."
	flye --nano-raw "$MTSEEDSDIR"/2.fq.gz \
		--out-dir "$MTDIR" \
		--threads "$NT" \
		--asm-coverage "$COV" \
		--genome-size "$CONTIG_LENGTH" \
		--stop-after contigger \
		>/dev/null 2>&1
	echoall "CHECK: assembly graph "$PWD/$MTDIR"/30-contigger/graph_final.gfa"
	# echo "column -t $ODIR/assembly_info_organelle_annotation_count.txt"

	jnum_next=$((JNUM + 1))
	echoall Create and edit $ODIR/$JNUM/mt.contig.name-${jnum_next}
	echoall NEXT: "$(basename "$0")" assemble2 -o "$ODIR" -j ${jnum_next}
	echoall or you could finish with Flye organelle-genome assembly with its polishing stage.
	echoall NEXT: "$(basename "$0")" flye-polishing -o "$ODIR" -j "$JNUM"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Finishes the Flye organelle-genome assembly.
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $MTDIR/contig.fa
#   $MTSEEDSDIR/2.fq.gz
#   $MTDIR/30-contigger
# Outputs:
#   $MTDIR/assembly_graph.gfa
################################################################################
function _run_polap_flye-polishing() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	echo "INFO: polishing organelle-genome assembly on $JNUM"

	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="$MTDIR"/seeds

	help_message=$(
		cat <<HEREDOC
# Finishes the Flye organelle-genome assembly.
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $MTDIR/contig.fa
#   $MTSEEDSDIR/2.fq.gz
#   $MTDIR/30-contigger
# Outputs:
#   $MTDIR/assembly_graph.gfa
Example: $(basename "$0") ${_arg_menu[0]} [-j|--jnum <arg>] [-t|--threads <arg>] [-c|--coverage <arg>]
HEREDOC
	)

	if [[ ${_arg_menu[1]} == "help" ]]; then
		echoerr "${help_message}"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "$MTDIR/contig.fa" ]; then
		echoall "ERROR: no selected-contig file [$MTDIR/contig.fa]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "$MTSEEDSDIR/2.fq.gz" ]; then
		echoall "ERROR: no long-read file [$MTSEEDSDIR/2.fq.gz]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	CONTIG_LENGTH=$(seqkit stats -Ta "$MTDIR"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	echo "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"

	echo "INFO: polishing the organelle-genome assembly using flye ... be patient!"
	echoerr "please, wait for Flye long-read polishing of the organelle-genome assembly on $JNUM ..."
	flye --nano-raw "$MTSEEDSDIR"/2.fq.gz \
		--out-dir "$MTDIR" \
		--threads "$NT" \
		--asm-coverage "$COV" \
		--genome-size "$CONTIG_LENGTH" \
		--resume \
		>/dev/null 2>&1
	echoall "CHECK: the long-read polished assembly graph $PWD/$MTDIR/assembly_graph.gfa"
	echoerr "DO: extract a draft organelle genome sequence (mt.0.fasta) from the polished assembly graph"
	# echo "column -t $ODIR/assembly_info_organelle_annotation_count.txt"
	# echoall NEXT: $(basename $0) check-coverage [-p $PA]
	echoall NEXT: "$(basename "$0")" prepare-polishing -a "$SR1" -b "$SR2"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
