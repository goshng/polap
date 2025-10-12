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

function _run_polap_mtpt {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local polap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - annotate rougly reads with organelle genes

Synopsis:
  polap ${polap_cmd} [options]

Description:
  polap ${polap_cmd} uses plastid and organelle genes to annotate reads
  using minimap2.

Options:
  -l FASTQ
    reads data file

Examples:
  Get organelle genome sequences:
    polap ${polap_cmd} -l l.fq

TODO:
  Dev.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

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

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" == "stats" ]]; then

		_polap_lib_conda-ensure_conda_env polap || exit 1

		# Minimal: recursively find files named 'mtpt.tsv' under the root folder
		Rscript --vanilla "${_POLAPLIB_DIR}/scripts/mtpt_turnover_fig1.R" \
			--root "$PWD" \
			--pattern mtpt.tsv \
			--out out_mtpt_summary
		conda deactivate

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_lib_conda-ensure_conda_env polap-evo || exit 1

	# _polap_log0 bash "${_POLAPLIB_DIR}/polap-bash-mt-stoich-fractions.sh" \
	# 	-f "${_arg_mt_ref}" \
	# 	-l "${_arg_long_reads}" \
	# 	-o "${_arg_outdir}/stoich/" \
	# 	--label S1 \
	# 	--threads 12 \
	# 	--min-repeat-len 1000 \
	# 	--min-repeat-pid 90 \
	# 	--flank 800 \
	# 	--min-mapq 20 \
	# 	--min-span 250 \
	# 	--prefer mummer \
	# 	--step all
	#
	# _polap_log0 bash "${_POLAPLIB_DIR}/polap-bash-mt-isomers-mtpt.sh" \
	# 	-f "${_arg_mt_ref}" \
	# 	-c "${_arg_pt_ref}" \
	# 	-l "${_arg_long_reads}" \
	# 	-o "${_arg_outdir}/mtpt/" \
	# 	--label S1 --threads 12 \
	# 	--min-repeat-len 1000 --min-repeat-pid 90 \
	# 	--flank 800 --min-mapq 20 --min-span 250 --prefer mummer \
	# 	--mtpt-min-len 150 --mtpt-min-pid 85 --mtpt-recent 97 --mtpt-intermediate 90 \
	# 	--enrich-window 1000 --enrich-perm 2000 --step all
	# edit these three paths

	local MT="${_arg_mt_ref}"
	local CP="${_arg_pt_ref}"
	local OUT="${_arg_outdir}/mtpt"

	mkdir -p "$OUT/mtpt"

	# makeblastdb -in "$MT" -dbtype nucl -out "$OUT/mtpt/mt"
	# blastn -task megablast \
	# 	-db "$OUT/mtpt/mt" \
	# 	-query "$CP" \
	# 	-evalue 1e-5 -dust no -soft_masking false \
	# 	-perc_identity 75 -word_size 11 \
	# 	-outfmt "6 qseqid sseqid pident length qstart qend sstart send qcovs" \
	# 	>"$OUT/mtpt/raw.blast6.tsv"
	#
	# python3 "$_POLAPLIB_DIR/scripts/polap_py_mtpt_scan.py" \
	# 	--blast6 "$OUT/mtpt/raw.blast6.tsv" \
	# 	--mt-fasta "$MT" \
	# 	--min-len 150 \
	# 	--recent 97 --intermediate 90 \
	# 	--out-tsv "$OUT/mtpt/mtpt.tsv" \
	# 	--out-bed "$OUT/mtpt/mtpt.bed"

	python3 "$_POLAPLIB_DIR/polap-py-mtpt-verify-prepare.py" \
		--fasta "$MT" \
		--reads "${_arg_long_reads}" \
		--mtpt-tsv "$OUT/mtpt/mtpt.tsv" \
		--flank 800 \
		--threads 12 \
		--preset map-ont \
		--out "$OUT/mtpt/verify"

	python3 "$_POLAPLIB_DIR/polap-py-mtpt-verify-batch.py" \
		--fasta "$MT" \
		--reads "${_arg_long_reads}" \
		--mtpt-tsv "$OUT/mtpt/mtpt.tsv" \
		--out "$OUT/mtpt" \
		--mode both

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
