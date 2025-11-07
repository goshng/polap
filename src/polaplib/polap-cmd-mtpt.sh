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

	local MT="${_arg_mt_ref}"
	local PT="${_arg_pt_ref}"
	local OUT="${_arg_outdir}/mtpt"

	mkdir -p "$OUT/mtpt"

	# python3 "$_POLAPLIB_DIR/polap-py-mtpt-verify-prepare.py" \
	# 	--fasta "$MT" \
	# 	--reads "${_arg_long_reads}" \
	# 	--mtpt-tsv "$OUT/mtpt/mtpt.tsv" \
	# 	--flank 800 \
	# 	--threads 12 \
	# 	--preset map-ont \
	# 	--out "$OUT/mtpt/verify"
	#
	# python3 "$_POLAPLIB_DIR/polap-py-mtpt-verify-batch.py" \
	# 	--fasta "$MT" \
	# 	--reads "${_arg_long_reads}" \
	# 	--mtpt-tsv "$OUT/mtpt/mtpt.tsv" \
	# 	--out "$OUT/mtpt" \
	# 	--mode both

	# ----------------------------------------------------------------------------
	# Run per-species
	# ----------------------------------------------------------------------------
	# 1) init one species (repeat per species)
	bash "$_POLAPLIB_DIR/polap-bash-init-analysis.sh" \
		-s "${_arg_outdir}" \
		--copy \
		--mt "$MT" \
		--pt "$PT"

	# 2) Step 1 per species (MTPT detector)
	bash "$_POLAPLIB_DIR/polap-bash-step1-mtpt-detector.sh" \
		-s "${_arg_outdir}" \
		-H OatkDB/v20230921/embryophyta_pltd.fam \
		-t 8 \
		--min-len 150 \
		--min-pid 0.85

	# 3) Step 2 (pool all species’ MTPTs)
	# bash "$_POLAPLIB_DIR/polap-bash-step2-mtpt-homology.sh" \
	# 	-t 8 --id 0.85 --cov 0.7

	# # 4) Step 3 (cp phylogeny; coding supermatrix)
	# ./polap-bash-step3-cp-phylogeny -H /path/plastid_hmms.hmm -t 12 --mode coding
	#
	# # 5) Step 4 (gain/loss mapping)
	# ./polap-bash-step4-gain-loss --model all
	#
	# # 1) Model comparison (AICc, Dollo-like flag)
	# ./polap-bash-upgrade-model-compare.sh -b man/analysis --nsim 200 --qratio-thr 0.1
	#
	# # 2) Bootstrap robustness across ufboot trees
	# ./polap-bash-upgrade-bootstrap-asr.sh -b man/analysis --B 200 --nsim 10 --model auto
	#
	# # 3) Erosion spectrum (quantile slopes)
	# ./polap-bash-extras-erosion-spectrum.sh -b man/analysis -o man/analysis/results \
	# 	--taus 0.5,0.9 --boots 500 \
	# 	--clades man/analysis/metadata/species_to_clade.tsv # optional
	#
	# # 4) IR standardization for Mauve path
	# ./polap-bash-upgrade-ir-standardize.sh -b man/analysis -H /path/plastid_hmms.hmm --min-ir 5000

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_mtpt-tree {
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

	local MT="${_arg_mt_ref}"
	local PT="${_arg_pt_ref}"
	local OUT="${_arg_outdir}/mtpt"

	mkdir -p "$OUT/mtpt"

	# python3 "$_POLAPLIB_DIR/polap-py-mtpt-verify-prepare.py" \
	# 	--fasta "$MT" \
	# 	--reads "${_arg_long_reads}" \
	# 	--mtpt-tsv "$OUT/mtpt/mtpt.tsv" \
	# 	--flank 800 \
	# 	--threads 12 \
	# 	--preset map-ont \
	# 	--out "$OUT/mtpt/verify"
	#
	# python3 "$_POLAPLIB_DIR/polap-py-mtpt-verify-batch.py" \
	# 	--fasta "$MT" \
	# 	--reads "${_arg_long_reads}" \
	# 	--mtpt-tsv "$OUT/mtpt/mtpt.tsv" \
	# 	--out "$OUT/mtpt" \
	# 	--mode both

	# ----------------------------------------------------------------------------
	# Run per-species
	# ----------------------------------------------------------------------------
	# 1) init one species (repeat per species)
	# bash "$_POLAPLIB_DIR/polap-bash-init-analysis.sh" \
	# 	-s "${_arg_outdir}" \
	# 	--copy \
	# 	--mt "$MT" \
	# 	--pt "$PT"

	# 2) Step 1 per species (MTPT detector)
	# bash "$_POLAPLIB_DIR/polap-bash-step1-mtpt-detector.sh" \
	# 	-s "${_arg_outdir}" \
	# 	-H OatkDB/v20230921/embryophyta_pltd.fam \
	# 	-t 8 \
	# 	--min-len 150 \
	# --min-pid 0.85

	# 3) Step 2 (pool all species’ MTPTs)
	# bash "$_POLAPLIB_DIR/polap-bash-step2-mtpt-homology.sh" \
	# 	-t 8 --id 0.85 --cov 0.7

	# 4) Step 3 (cp phylogeny; coding supermatrix)
	# bash "$_POLAPLIB_DIR/polap-bash-step3-pt-phylogeny.sh" \
	# 	-H OatkDB/v20230921/embryophyta_pltd.fam \
	# 	-t 12 --mode coding

	# 5) Step 4 (gain/loss mapping)
	bash "$_POLAPLIB_DIR/polap-bash-step4-gain-loss.sh" \
		-t "man/analysis/pt_tree/concat/partitions.nex.treefile" \
		--outpdf "man/analysis/results/tree.pdf" \
		--model all

	# # 1) Model comparison (AICc, Dollo-like flag)
	# ./polap-bash-upgrade-model-compare.sh -b man/analysis --nsim 200 --qratio-thr 0.1
	#
	# # 2) Bootstrap robustness across ufboot trees
	# ./polap-bash-upgrade-bootstrap-asr.sh -b man/analysis --B 200 --nsim 10 --model auto
	#
	# # 3) Erosion spectrum (quantile slopes)
	# ./polap-bash-extras-erosion-spectrum.sh -b man/analysis -o man/analysis/results \
	# 	--taus 0.5,0.9 --boots 500 \
	# 	--clades man/analysis/metadata/species_to_clade.tsv # optional
	#
	# # 4) IR standardization for Mauve path
	# ./polap-bash-upgrade-ir-standardize.sh -b man/analysis -H /path/plastid_hmms.hmm --min-ir 5000

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_mtpt_old {
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

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_lib_conda-ensure_conda_env polap-evo || exit 1

	local MT="${_arg_mt_ref}"
	local CP="${_arg_pt_ref}"
	local OUT="${_arg_outdir}"

	mkdir -p "$OUT/mtpt"
	makeblastdb -in "$MT" -dbtype nucl -out "$OUT/mtpt/mt"
	blastn -task megablast \
		-db "$OUT/mtpt/mt" \
		-query "$CP" \
		-evalue 1e-5 -dust no -soft_masking false \
		-perc_identity 75 -word_size 11 \
		-outfmt "6 qseqid sseqid pident length qstart qend sstart send qcovs" \
		>"$OUT/mtpt/raw.blast6.tsv"

	python3 "$_POLAPLIB_DIR/scripts/polap_py_mtpt_scan.py" \
		--blast6 "$OUT/mtpt/raw.blast6.tsv" \
		--mt-fasta "$MT" \
		--min-len 150 \
		--recent 97 --intermediate 90 \
		--out-tsv "$OUT/mtpt/mtpt.tsv" \
		--out-bed "$OUT/mtpt/mtpt.bed"

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
