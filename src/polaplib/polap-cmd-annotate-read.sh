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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_polap_annotate-read-pt {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap annotate-read - annotate reads with organelle genes

Synopsis:
  polap annotate-read [options]

Description:
  polap annotate-read uses plastid and organelle genes to annotate reads
  using minimap2.

Options:
  -l FASTQ
    reads data file

  --annotate-read-min-mapq INT
    minimum mapping quality for reads (0 ~ 255)

  --annotate-read-min-identity FLOAT
    minimum identity for reads (0 ~ 1.00)

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

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

	local annotatedir="${_arg_outdir}"/annotate-read-pt
	mkdir -p "${annotatedir}"
	local pt_table="${annotatedir}"/pt-contig-annotation-depth-table.txt
	local mt_table="${annotatedir}"/contig-annotation-depth-table.txt
	local all_table="${annotatedir}"/assembly_info_organelle_annotation_count-all.txt

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${all_table}"
			;;
		mt)
			_polap_log0_column "${mt_table}"
			;;
		pt-table | pt)
			_polap_log0_column "${pt_table}"
			;;
		table | *)
			_polap_log0_column "${mt_table}"
			;;
		esac
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# if [[ -s "${pt_table}" && -s "${mt_table}" && -s "${at_table}" ]]; then
	# 	_polap_log1 "Annotation tables were already prepared."
	# 	return 0
	# fi

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.fna

	_polap_log2 "map reads on mitochondrial genes using minimap2"
	_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${MTAA} \
      ${_arg_long_reads} \
      -t ${_arg_threads} \
      -o ${annotatedir}/mt.paf \
      >${_polap_output_dest} 2>&1"

	_polap_log2 "map reads on plastid genes using minimap2"
	_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${PTAA} \
      ${_arg_long_reads} \
      -t ${_arg_threads} \
      -o ${annotatedir}/pt.paf \
      >${_polap_output_dest} 2>&1"

	_polap_log2 "create the annotation table for mt and pt genes"
	if [[ -s "${annotatedir}"/mt.paf ]] &&
		[[ -s "${annotatedir}"/pt.paf ]]; then
		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-reads.R \
    --mt ${annotatedir}/mt.paf \
    --pt ${annotatedir}/pt.paf \
    --output ${annotatedir} \
    --min-mapq ${_arg_annotate_read_min_mapq} \
    --min-identity ${_arg_annotate_read_min_identity} \
	  >${_polap_output_dest} 2>&1"
	else
		_polap_log0 "No minimap2 results!"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_annotate-read-mt {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap annotate-read-mt - annotate reads with organelle genes for mitochondrial genome assembly

Synopsis:
  polap annotate-read-mt [options]

Description:
  polap annotate-read-mt uses organelle genes to annotate reads using minimap2.

Options:
  -l FASTQ
    reads data file

  --annotate-read-min-mapq INT
    minimum mapping quality for reads (0 ~ 255)

  --annotate-read-min-identity FLOAT
    minimum identity for reads (0 ~ 1.00)

Examples:
  Get organelle genome sequences:
    polap annotate-read-mt -l l.fq

TODO:
  Dev.

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

	local annotatedir="${_arg_outdir}"/annotate-read-mt
	mkdir -p "${annotatedir}/at"
	local pt_table="${annotatedir}"/pt-contig-annotation-depth-table.txt
	local mt_table="${annotatedir}"/contig-annotation-depth-table.txt
	local at_table="${annotatedir}"/at/contig-annotation-depth-table.txt
	local all_table="${annotatedir}"/assembly_info_organelle_annotation_count-all.txt

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${all_table}"
			;;
		mt)
			_polap_log0_column "${mt_table}"
			;;
		pt-table | pt)
			_polap_log0_column "${pt_table}"
			;;
		table | *)
			_polap_log0_column "${mt_table}"
			;;
		esac
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# if [[ -s "${pt_table}" && -s "${mt_table}" && -s "${at_table}" ]]; then
	# 	_polap_log1 "Annotation tables were already prepared."
	# 	return 0
	# fi

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.fna
	# local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.3k.c80.fna
	local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.c80.fna

	_polap_log1 "  step 1: map reads on mitochondrial genes using minimap2"
	_polap_log2 "    input1: ${_arg_minimap2_data_type}"
	_polap_log2 "    input2: ${MTAA}"
	_polap_log2 "    input3: ${_arg_long_reads}"
	_polap_log2 "    output: ${annotatedir}/mt.paf"
	_polap_log3_cmdout minimap2 -cx \
		${_arg_minimap2_data_type} \
		${MTAA} \
		${_arg_long_reads} \
		-t ${_arg_threads} \
		-o ${annotatedir}/mt.paf

	_polap_log1 "  step 2: map reads on plastid genes using minimap2"
	_polap_log2 "    input1: ${_arg_minimap2_data_type}"
	_polap_log2 "    input2: ${PTAA}"
	_polap_log2 "    input3: ${_arg_long_reads}"
	_polap_log2 "    output: ${annotatedir}/pt.paf"
	_polap_log3_cmdout minimap2 -cx \
		${_arg_minimap2_data_type} \
		${PTAA} \
		${_arg_long_reads} \
		-t ${_arg_threads} \
		-o ${annotatedir}/pt.paf

	_polap_log1 "  step3: create the annotation table for mt and pt genes"
	_polap_log2 "    input1: ${annotatedir}/mt.paf"
	_polap_log2 "    input2: ${annotatedir}/pt.paf"
	_polap_log2 "    input3: --min-mapq ${_arg_annotate_read_min_mapq}"
	_polap_log2 "    input4: --min-identity ${_arg_annotate_read_min_identity}"
	_polap_log2 "    input5: --min-pt 1"
	_polap_log2 "    output: ${annotatedir}"
	if [[ -s "${annotatedir}"/mt.paf ]] &&
		[[ -s "${annotatedir}"/pt.paf ]]; then
		_polap_log3_cmdout Rscript --vanilla "${_POLAPLIB_DIR}/polap-r-reads.R" \
			--mt "${annotatedir}"/mt.paf \
			--pt "${annotatedir}"/pt.paf \
			--output "${annotatedir}" \
			--min-mapq "${_arg_annotate_read_min_mapq}" \
			--min-identity "${_arg_annotate_read_min_identity}"
	else
		_polap_log0 "No minimap2 results!"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_annotate-read-nt {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polapannotate-read-nt - annotate reads with organelle genes

Synopsis:
  polapannotate-read-nt [options]

Description:
  polapannotate-read-nt uses plastid and organelle genes to annotate reads
  using minimap2.

Options:
  -l FASTQ
    reads data file

  --annotate-read-min-mapq INT
    minimum mapping quality for reads (0 ~ 255)

  --annotate-read-min-identity FLOAT
    minimum identity for reads (0 ~ 1.00)

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

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

	local annotatedir="${_arg_outdir}"/annotate-read-nt
	mkdir -p "${annotatedir}/at"
	local pt_table="${annotatedir}"/pt-contig-annotation-depth-table.txt
	local mt_table="${annotatedir}"/contig-annotation-depth-table.txt
	local at_table="${annotatedir}"/at/contig-annotation-depth-table.txt
	local all_table="${annotatedir}"/assembly_info_organelle_annotation_count-all.txt

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${all_table}"
			;;
		mt)
			_polap_log0_column "${mt_table}"
			;;
		pt-table | pt)
			_polap_log0_column "${pt_table}"
			;;
		table | *)
			_polap_log0_column "${mt_table}"
			;;
		esac
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# if [[ -s "${pt_table}" && -s "${mt_table}" && -s "${at_table}" ]]; then
	# 	_polap_log1 "Annotation tables were already prepared."
	# 	return 0
	# fi

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.fna
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.fna
	# local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.3k.c80.fna
	local NTAA="${_POLAPLIB_DIR}"/polap-mt.noncds.c80.fna

	_polap_log2 "map reads on mitochondrial genes using minimap2"
	_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${MTAA} \
      ${_arg_long_reads} \
      -t ${_arg_threads} \
      -o ${annotatedir}/mt.paf \
      >${_polap_output_dest} 2>&1"

	_polap_log2 "map reads on plastid genes using minimap2"
	_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${PTAA} \
      ${_arg_long_reads} \
      -t ${_arg_threads} \
      -o ${annotatedir}/pt.paf \
      >${_polap_output_dest} 2>&1"

	_polap_log2 "create the annotation table for mt and pt genes"
	if [[ -s "${annotatedir}"/mt.paf ]] &&
		[[ -s "${annotatedir}"/pt.paf ]]; then
		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-reads.R \
    --mt ${annotatedir}/mt.paf \
    --pt ${annotatedir}/pt.paf \
    --output ${annotatedir} \
    --min-mapq ${_arg_annotate_read_min_mapq} \
    --min-identity ${_arg_annotate_read_min_identity} \
	  >${_polap_output_dest} 2>&1"
	else
		_polap_log0 "No minimap2 results!"
	fi

	_polap_log2 "map reads on mitochondrial noncoding regions using minimap2"
	# v0.5.2.3
	# -k13 -w5 -m20 -p0.6 -N20 \
	_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${NTAA} \
      ${_arg_long_reads} \
      -t ${_arg_threads} \
      -o ${annotatedir}/nt.paf \
      >${_polap_output_dest} 2>&1"

	_polap_log2 "create the annotation table for mito noncoding and pt genes"
	if [[ -s "${annotatedir}"/nt.paf ]] &&
		[[ -s "${annotatedir}"/pt.paf ]]; then
		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-reads.R \
    --mt ${annotatedir}/nt.paf \
    --pt ${annotatedir}/pt.paf \
    --output ${annotatedir}/at \
    --min-mapq ${_arg_annotate_read_min_mapq} \
    --min-identity ${_arg_annotate_read_min_identity} \
    --min-pt 1 \
	  >${_polap_output_dest} 2>&1"
	else
		_polap_log0 "No minimap2 results!"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_assemble-annotated-read {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap assemble-annotated-read - annotate reads with organelle genes

Synopsis:
  polap assemble-annotated-read [options]

Description:
  polap assemble-annotated-read ...

Options:
  --plastid

  --animal

  -l FASTQ
    reads data file

  --annotate-read-min-mapq INT
    minimum mapping quality for reads (0 ~ 255)

  --annotate-read-min-identity FLOAT
    minimum identity for reads (0 ~ 1.00)

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

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

	local annotatedir="${_arg_outdir}"/annotate-read-mt
	if [[ "${_arg_plastid}" == "on" ]]; then
		annotatedir="${_arg_outdir}"/annotate-read-pt
	fi
	mkdir -p "${annotatedir}/at"
	local pt_table="${annotatedir}"/pt-contig-annotation-depth-table.txt
	local mt_table="${annotatedir}"/contig-annotation-depth-table.txt
	local at_table="${annotatedir}"/at/contig-annotation-depth-table.txt
	local all_table="${annotatedir}"/assembly_info_organelle_annotation_count-all.txt

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${all_table}"
			;;
		mt)
			_polap_log0_column "${mt_table}"
			;;
		pt-table | pt)
			_polap_log0_column "${pt_table}"
			;;
		table | *)
			_polap_log0_column "${mt_table}"
			;;
		esac
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# seed reads
	if [[ "${_arg_plastid}" == "on" ]]; then
		# for plastid genome assembly
		if [[ ! -s "${pt_table}" ]]; then
			_polap_log0 "ERROR: no pt table: ${pt_table}"
			return 1
		fi

		Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
			--table "${pt_table}" \
			--length 3e+7 \
			--output "${annotatedir}"/pt.id.txt
		# tail -n +2 "${pt_table}" | cut -f1 >"${annotatedir}"/pt.id.txt

		seqtk subseq "${_arg_long_reads}" "${annotatedir}"/pt.id.txt >"${annotatedir}"/pt.fq

		# subsample the data so that pt.fq is less than 30Mb
		rm -f "${annotatedir}"/pt0.fq
		_polap_lib_fastq-sample-to \
			"${annotatedir}"/pt.fq "${annotatedir}"/pt0.fq "30m"

		# flye v2.9.6
		if _polap_lib_version-check_flye_version; then
			rm -rf "${annotatedir}"/pt
			_polap_log0 "flye 2.9.6 assembly on ptDNA"
			flye "${_arg_flye_data_type}" "${annotatedir}"/pt0.fq -t "${_arg_threads}" \
				--out-dir "${annotatedir}"/pt 2>"${_polap_output_dest}"
			if [[ -s "${annotatedir}"/pt/assembly_graph.gfa ]]; then
				_polap_log1 "output: PT assembly: ${annotatedir}/pt/assembly_graph.gfa"
			else
				_polap_log0 "output: no PT assembly"
			fi
		else
			echo "Flye 2.9.6 is required. Aborting."
			exit 1
		fi
		# Go to: readassemble-ont-pt_genus_species() { at polap-lib-data.sh
		# Go to: readassemble-ont-pt-iterate_genus_species() { at polap-lib-data.sh
	else
		# for mitochondrial genome assembly

		if [[ ! -s "${mt_table}" ]]; then
			_polap_log0 "ERROR: no mt table: ${mt_table}"
			return 1
		fi

		if [[ "${_arg_animal}" == "off" ]]; then

			# Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
			# 	--table "${pt_table}" \
			# 	--length 3e+8 \
			# 	--output "${annotatedir}"/pt0mt.id.txt
			# _polap_log1 "pt0mt.id: $(cat "${annotatedir}"/pt0mt.id.txt | wc -l)"

			Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
				--table "${mt_table}" \
				--length 1e+8 \
				--output "${annotatedir}"/mt0.id.txt
			_polap_log1 "mt0.id: $(cat "${annotatedir}"/mt0.id.txt | wc -l)"

			# Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
			# 	--table "${at_table}" \
			# 	--length 1e+8 \
			# 	--output "${annotatedir}"/at0.id.txt
			# _polap_log1 "at0.id: $(cat "${annotatedir}"/at0.id.txt | wc -l)"

			# subtract pt.id.txt from at0.id.txt
			# grep -vFf "${annotatedir}"/pt0mt.id.txt "${annotatedir}"/at0.id.txt \
			# 	>"${annotatedir}"/at.id.txt
			# _polap_log1 "at.id (at0 - pt0mt): $(cat "${annotatedir}"/at.id.txt | wc -l)"
			#
			# cat "${annotatedir}"/mt0.id.txt "${annotatedir}"/at.id.txt |
			# 	sort | uniq >"${annotatedir}"/mt.id.txt
			# _polap_log1 "mt.id (mt0 + at): $(cat "${annotatedir}"/mt.id.txt | wc -l)"

			# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
			rm -f "${annotatedir}"/mt0.fq
			seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt0.fq

			# subsample the data so that mt.fq is less than 100Mb
			rm -f "${annotatedir}"/mt.fq
			_polap_lib_fastq-sample-to \
				"${annotatedir}"/mt0.fq "${annotatedir}"/mt.fq "100m"

			# flye v2.9.6
			if _polap_lib_version-check_flye_version; then
				rm -rf "${annotatedir}"/mt
				_polap_log0 "flye 2.9.6 assembly on mtDNA"
				flye "${_arg_flye_data_type}" "${annotatedir}"/mt.fq -t "${_arg_threads}" \
					--out-dir "${annotatedir}"/mt 2>"${_polap_output_dest}"
				if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
					_polap_log1 "output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
				else
					_polap_log0 "output: no MT assembly"
				fi
			else
				echo "Flye 2.9.6 is required. Aborting."
				exit 1
			fi
		# Go to: readassemble-ont-mt_genus_species() { at polap-lib-data.sh
		# Go to: readassemble-ont-mt-iterate_genus_species() { at polap-lib-data.sh
		else
			# for animal mitochondrial genome assembly
			Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
				--table "${mt_table}" \
				--length 1e+6 \
				--output "${annotatedir}"/mt0.id.txt

			# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
			seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt.fq

			# subsample the data so that mt.fq is less than 100Mb
			rm -f "${annotatedir}"/mt0.fq
			_polap_lib_fastq-sample-to \
				"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "1m"

			# flye v2.9.6
			if _polap_lib_version-check_flye_version; then
				rm -rf "${annotatedir}"/mt
				_polap_log0 "flye 2.9.6 assembly on mtDNA"
				flye "${_arg_flye_data_type}" "${annotatedir}"/mt0.fq -t "${_arg_threads}" \
					--out-dir "${annotatedir}"/mt 2>"${_polap_output_dest}"
				if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
					_polap_log1 "output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
				else
					_polap_log0 "output: no MT assembly"
				fi
			else
				echo "Flye 2.9.6 is required. Aborting."
				exit 1
			fi
		fi
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_assemble-annotated-read-nt {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap assemble-annotated-read - annotate reads with organelle genes

Synopsis:
  polap assemble-annotated-read [options]

Description:
  polap assemble-annotated-read ...

Options:
  --plastid

  --animal

  -l FASTQ
    reads data file

  --annotate-read-min-mapq INT
    minimum mapping quality for reads (0 ~ 255)

  --annotate-read-min-identity FLOAT
    minimum identity for reads (0 ~ 1.00)

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

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

	local annotatedir="${_arg_outdir}"/annotate-read-nt
	if [[ "${_arg_plastid}" == "on" ]]; then
		annotatedir="${_arg_outdir}"/annotate-read-pt
	fi
	mkdir -p "${annotatedir}/at"
	local pt_table="${annotatedir}"/pt-contig-annotation-depth-table.txt
	local mt_table="${annotatedir}"/contig-annotation-depth-table.txt
	local at_table="${annotatedir}"/at/contig-annotation-depth-table.txt
	local all_table="${annotatedir}"/assembly_info_organelle_annotation_count-all.txt

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${all_table}"
			;;
		mt)
			_polap_log0_column "${mt_table}"
			;;
		pt-table | pt)
			_polap_log0_column "${pt_table}"
			;;
		table | *)
			_polap_log0_column "${mt_table}"
			;;
		esac
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# seed reads
	if [[ "${_arg_plastid}" == "on" ]]; then
		# for plastid genome assembly
		if [[ ! -s "${pt_table}" ]]; then
			_polap_log0 "ERROR: no pt table: ${pt_table}"
			return 1
		fi

		Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
			--table "${pt_table}" \
			--length 3e+7 \
			--output "${annotatedir}"/pt.id.txt
		# tail -n +2 "${pt_table}" | cut -f1 >"${annotatedir}"/pt.id.txt

		seqtk subseq "${_arg_long_reads}" "${annotatedir}"/pt.id.txt >"${annotatedir}"/pt.fq

		# subsample the data so that pt.fq is less than 30Mb
		rm -f "${annotatedir}"/pt0.fq
		_polap_lib_fastq-sample-to \
			"${annotatedir}"/pt.fq "${annotatedir}"/pt0.fq "30m"

		# flye v2.9.6
		if _polap_lib_version-check_flye_version; then
			rm -rf "${annotatedir}"/pt
			_polap_log0 "flye 2.9.6 assembly on ptDNA"
			flye "${_arg_flye_data_type}" "${annotatedir}"/pt0.fq -t "${_arg_threads}" \
				--out-dir "${annotatedir}"/pt 2>"${_polap_output_dest}"
			if [[ -s "${annotatedir}"/pt/assembly_graph.gfa ]]; then
				_polap_log1 "output: PT assembly: ${annotatedir}/pt/assembly_graph.gfa"
			else
				_polap_log0 "output: no PT assembly"
			fi
		else
			echo "Flye 2.9.6 is required. Aborting."
			exit 1
		fi
	else
		# for mitochondrial genome assembly

		if [[ ! -s "${mt_table}" ]]; then
			_polap_log0 "ERROR: no mt table: ${mt_table}"
			return 1
		fi

		if [[ "${_arg_animal}" == "off" ]]; then

			Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
				--table "${pt_table}" \
				--length 3e+9 \
				--output "${annotatedir}"/pt0mt.id.txt
			_polap_log1 "pt0mt.id: $(cat "${annotatedir}"/pt0mt.id.txt | wc -l)"

			Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
				--table "${mt_table}" \
				--length 1e+9 \
				--output "${annotatedir}"/mt0.id.txt
			_polap_log1 "mt0.id: $(cat "${annotatedir}"/mt0.id.txt | wc -l)"

			Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
				--table "${at_table}" \
				--length 1e+9 \
				--output "${annotatedir}"/at0.id.txt
			_polap_log1 "at0.id: $(cat "${annotatedir}"/at0.id.txt | wc -l)"

			# subtract pt.id.txt from at0.id.txt
			grep -vFf "${annotatedir}"/pt0mt.id.txt "${annotatedir}"/at0.id.txt \
				>"${annotatedir}"/at.id.txt
			_polap_log1 "at.id (at0 - pt0mt): $(cat "${annotatedir}"/at.id.txt | wc -l)"

			cat "${annotatedir}"/mt0.id.txt "${annotatedir}"/at.id.txt |
				sort | uniq >"${annotatedir}"/mt.id.txt
			_polap_log1 "mt.id (mt0 + at): $(cat "${annotatedir}"/mt.id.txt | wc -l)"

			# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
			rm -f "${annotatedir}"/mt.fq
			seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt.id.txt >"${annotatedir}"/mt.fq

			# subsample the data so that mt.fq is less than 100Mb
			rm -f "${annotatedir}"/mt0.fq
			_polap_lib_fastq-sample-to \
				"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "100m"

			# flye v2.9.6
			if _polap_lib_version-check_flye_version; then
				rm -rf "${annotatedir}"/mt
				_polap_log0 "flye 2.9.6 assembly on mtDNA"
				flye "${_arg_flye_data_type}" "${annotatedir}"/mt0.fq -t "${_arg_threads}" \
					--out-dir "${annotatedir}"/mt 2>"${_polap_output_dest}"
				if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
					_polap_log1 "output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
				else
					_polap_log0 "output: no MT assembly"
				fi
			else
				echo "Flye 2.9.6 is required. Aborting."
				exit 1
			fi
		# Go to: readassemble-ont-nt_genus_species() { at polap-lib-data.sh
		# Go to: readassemble-ont-nt-iterate_genus_species() { at polap-lib-data.sh
		else
			# for animal mitochondrial genome assembly
			Rscript --vanilla "${_POLAPLIB_DIR}"/polap-r-filter-organelle-reads.R \
				--table "${mt_table}" \
				--length 1e+6 \
				--output "${annotatedir}"/mt0.id.txt

			# tail -n +2 "${mt_table}" | cut -f1 >"${annotatedir}"/mt.id.txt
			seqtk subseq "${_arg_long_reads}" "${annotatedir}"/mt0.id.txt >"${annotatedir}"/mt.fq

			# subsample the data so that mt.fq is less than 100Mb
			rm -f "${annotatedir}"/mt0.fq
			_polap_lib_fastq-sample-to \
				"${annotatedir}"/mt.fq "${annotatedir}"/mt0.fq "1m"

			# flye v2.9.6
			if _polap_lib_version-check_flye_version; then
				rm -rf "${annotatedir}"/mt
				_polap_log0 "flye 2.9.6 assembly on mtDNA"
				flye "${_arg_flye_data_type}" "${annotatedir}"/mt0.fq -t "${_arg_threads}" \
					--out-dir "${annotatedir}"/mt 2>"${_polap_output_dest}"
				if [[ -s "${annotatedir}"/mt/assembly_graph.gfa ]]; then
					_polap_log1 "output: MT assembly: ${annotatedir}/mt/assembly_graph.gfa"
				else
					_polap_log0 "output: no MT assembly"
				fi
			else
				echo "Flye 2.9.6 is required. Aborting."
				exit 1
			fi
		fi
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
