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

function _run_polap_polish-ont {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<EOF
Name:
  polap polish - polishes assembly using either long-read or short-read data.

Synopsis:
  polap polish [options]

Description:
  polap polish uses tools to polish sequences in fasta format using long-read or short-read data.

Options:
  -l FASTQ
    reads data file

Examples:
  Get organelle genome sequences:
    polap polish-ont -l l.fq -p pt.unpolished.fa -f pt.polished.fa

TODO:
  racon conda installation does not work. Install it from the source.

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

	bash "${_POLAPLIB_DIR}/polap-bash-polish-ptdna-ont.sh" \
		--threads "${_arg_threads}" \
		--outdir "${_arg_outdir}/polish-ont" \
		"${_arg_unpolished_fasta}" \
		"${_arg_long_reads}"

	cp "${_arg_outdir}/polish-ont/pt.best.fa" "${_arg_final_assembly}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# Rewrite of _run_polap_polish3 → calls polap-bash-polish.sh
# - Supports ONT-only or Hybrid polishing (auto-detect)
# - Accepts ONE short-read file via -s ( interleaved/merged R1+R2 or single-end )
#   • If two files are available, still accepts -1/-2 (takes precedence over -s)
# - Optionally preserves Flye GFA topology with --gfa-preserve
# - Keeps your help/man + logging style
#
# Requires: ${_POLAPLIB_DIR}/polap-bash-polish.sh  (from our previous step)
#
# References (for chosen polishers and graph preservation):
#  - Racon/minimap2 long-read polishing (no model dependence): Kolmogorov et al. 2019; Li 2018.
#  - FMLRC2 assembly polishing with short reads (alignment-free): Wang et al. 2021.
#  - Polypolish (repeat-aware with all alignments): Wick & Holt 2022.
#  - GFA extraction/reinjection pattern via gfatools to preserve paths/links: Li 2020.

# Version: v0.5.1
function _run_polap_polish2 { # polish organelle genome sequences (ONT-only or Hybrid) via backend
	# Enable debugging if _POLAP_DEBUG is set
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo "$FUNCNAME" | sed 's/_run_polap_//')"
	local polap_cmd="${FUNCNAME##*_}"

	# Verbosity routing
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose:-0}" -ge "${_polap_var_function_verbose:-2}" ] && _polap_output_dest="/dev/stderr"

	# Common vars
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# ---------- Help ----------
	local help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - polish a draft organelle assembly (FASTA or GFA) using ONT and/or short reads

Synopsis:
  polap ${polap_cmd} [options]

Description:
  Front-end driver that delegates to 'polap-bash-polish.sh' (backend v0.5.1).
  Automatically selects ONT-only or Hybrid polishing, optionally preserves Flye GFA topology (extract → polish → reinject).

Options (mapped to backend):
  -p FILE      Draft assembly (FASTA or GFA). [default: ${_arg_unpolished_fasta}]
  -f FILE      Output polished assembly (FASTA/GFA). [default: ${_arg_final_assembly}]
  -l FILE      ONT reads (FASTQ[.gz]) for Racon.
  -1 FILE      Short-read R1 (FASTQ[.gz]).
  -2 FILE      Short-read R2 (FASTQ[.gz]).
  -m MODE      Mode: auto|ont|hybrid  (default: ${_arg_mode:-auto})
  -r INT       Racon rounds (ONT-only). (default: ${_arg_racon_round:-3})
  -t INT       Threads. (default: ${_arg_threads:-16})
  -g           Preserve GFA topology (if input is .gfa).
  --polisher   Hybrid: auto|fmlrc2|polypolish (default: ${_arg_polisher:-auto})

Advanced speed/QC (pass-through):
  --speed normal|fast|turbo
  --recruit off|draft|bait|auto   (default backend: off → use ALL reads)
  --bait BAITS.fa
  --ont-cap X                     (0=disable; if enabling, selection = len>3000 & id>0.75, MAPQ ignored)
  --ont-cap-method rasusa|seqtk
  --minimap-preset map-ont|map-hifi
  --ab-test                       (A=ALL vs B=recruit+cap)
  --qc-report
  --eval                          (Merqury QV if SR available)
  --kmer K                        (default: 31)
EOF
	)
	if [[ ${_arg_menu[1]:-} == "help" || "${_arg_help:-off}" == "on" ]]; then
		local manfile
		manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]:-polish}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# ---------- Resolve inputs ----------
	local in_asm="${_arg_unpolished_fasta:-}"
	local out_final="${_arg_final_assembly:-}"
	local ont_fq="${_arg_long_reads:-}" # -l
	if [[ "${_arg_short_read1_is}" == "on" ]]; then
		local sr_one="${_arg_short_read1:-}" # maps to -s (single/interleaved)
		local sr_r1="${_arg_short_read1:-}"  # maps to -1
	else
		local sr_one=""
		local sr_r1=""
	fi
	if [[ "${_arg_short_read2_is}" == "on" ]]; then
		local sr_r2="${_arg_short_read2:-}" # maps to -1
	else
		local sr_r2=""
	fi
	_polap_log0 "s1: $sr_r1"

	local racon_rounds="${_arg_racon_round:-3}"
	local threads="${_arg_threads:-16}"
	local mode="${_arg_mode:-auto}"                # auto|ont|hybrid
	local gfa_preserve="${_arg_gfa_preserve:-off}" # on/off
	local do_eval="${_arg_eval:-off}"              # on/off
	local kmer="${_arg_kmer:-31}"
	local polisher="${_arg_polisher:-auto}" # auto|fmlrc2|polypolish

	# New speed/QC extras (optional in caller)
	local speed="${_arg_speed:-normal}"             # normal|fast|turbo
	local recruit="${_arg_recruit:-off}"            # off|draft|bait|auto
	local bait_fa="${_arg_bait_fasta:-}"            # path or empty
	local ont_cap="${_arg_ont_cap:-}"               # empty -> do not pass, backend default applies
	local ont_cap_method="${_arg_ont_cap_method:-}" # rasusa|seqtk
	local ab_test="${_arg_ab_test:-off}"            # on/off
	local qc_report="${_arg_qc_report:-off}"        # on/off
	local minimap_preset="${_arg_minimap_preset:-}" # map-ont|map-hifi

	# ---------- Conda env ----------
	_polap_lib_conda-ensure_conda_env polap-polish || return $EXIT_ERROR

	# ---------- Validate ----------
	if [[ -z "$in_asm" || ! -s "$in_asm" ]]; then
		_polap_log0 "ERROR: missing/empty draft assembly (-p): [$in_asm]"
		conda deactivate || true
		return $EXIT_ERROR
	fi
	if [[ -z "$out_final" ]]; then
		_polap_log0 "ERROR: missing final output path (-f)."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	# Mode resolution if auto
	_polap_log0 "mode1: $mode"
	if [[ "$mode" == "auto" ]]; then
		if [[ -n "$sr_r1" && -n "$sr_r2" ]]; then
			mode="hybrid"
		elif [[ -n "$ont_fq" ]]; then
			mode="ont"
		else
			_polap_log0 "ERROR: mode=auto needs -l (ONT) or -1/-2 (short reads)."
			conda deactivate || true
			return $EXIT_ERROR
		fi
	fi

	_polap_log0 "mode2: $mode"

	[[ "$mode" == "ont" && -z "$ont_fq" ]] && {
		_polap_log0 "ERROR: mode=ont requires -l."
		conda deactivate || true
		return $EXIT_ERROR
	}
	[[ "$mode" == "hybrid" && (-z "$sr_r1" || -z "$sr_r2") ]] && {
		_polap_log0 "ERROR: mode=hybrid requires -1 and -2."
		conda deactivate || true
		return $EXIT_ERROR
	}

	# ---------- Backend ----------
	local backend="${_POLAPLIB_DIR}/polap-bash-polish.sh"
	if [[ ! -s "$backend" ]]; then
		_polap_log0 "ERROR: backend not found: $backend"
		conda deactivate || true
		return $EXIT_ERROR
	fi

	local prefix
	prefix="$(basename "${out_final%.*}")"
	local args=("polish" "$in_asm" "--threads" "$threads" "--prefix" "$prefix" "--mode" "$mode")

	# ONT
	if [[ -n "$ont_fq" ]]; then
		args+=("--ont" "$ont_fq" "--racon-rounds" "$racon_rounds")
		[[ -n "$speed" ]] && args+=("--speed" "$speed")
		[[ -n "$recruit" ]] && args+=("--recruit" "$recruit")
		[[ -n "$bait_fa" ]] && args+=("--bait" "$bait_fa")
		[[ -n "$ont_cap" ]] && args+=("--ont-cap" "$ont_cap")
		[[ -n "$ont_cap_method" ]] && args+=("--ont-cap-method" "$ont_cap_method")
		[[ -n "$minimap_preset" ]] && args+=("--minimap-preset" "$minimap_preset")
	fi

	# Hybrid
	if [[ "$mode" == "hybrid" ]]; then
		args+=("--sr1" "$sr_r1" "--sr2" "$sr_r2" "--polisher" "$polisher")
	fi

	# GFA-preserve
	if [[ "${gfa_preserve}" =~ ^(on|true|1)$ ]]; then
		args+=("--gfa-preserve")
	fi

	# AB-test / QC report
	if [[ "${ab_test}" =~ ^(on|true|1)$ ]]; then
		args+=("--ab-test")
		# The backend will force QC when ab-test is on; we still allow explicit --qc-report for single-run too
	fi
	if [[ "${qc_report}" =~ ^(on|true|1)$ ]]; then
		args+=("--qc-report")
	fi

	# Merqury eval (will only be effective when SR are present; backend guards it)
	if [[ "${do_eval}" =~ ^(on|true|1)$ ]]; then
		args+=("--eval-merqury" "--kmer" "$kmer")
	fi

	_polap_log1 "INFO: backend cmd: $(printf "%q " "$backend" "${args[@]}")" >&2

	# ---------- Execute ----------
	if ! bash "$backend" "${args[@]}" >"${_polap_output_dest}" 2>&1; then
		_polap_log0 "ERROR: polishing failed (see backend log under out/${prefix}.polish_*/run.log)."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	# ---------- Locate and copy final artifact ----------
	local outdir produced
	outdir="$(ls -1dt out/"$prefix".polish_* 2>/dev/null | head -n1)"
	if [[ -z "$outdir" || ! -d "$outdir" ]]; then
		_polap_log0 "ERROR: cannot find backend output directory for prefix=$prefix"
		conda deactivate || true
		return $EXIT_ERROR
	fi

	if compgen -G "${outdir}/*.polished.gfa" >/dev/null; then
		produced="$(ls -1 ${outdir}/*.polished.gfa | head -n1)"
	elif compgen -G "${outdir}/*.polished.fa" >/dev/null; then
		produced="$(ls -1 ${outdir}/*.polished.fa | head -n1)"
	fi

	if [[ -z "$produced" || ! -s "$produced" ]]; then
		_polap_log0 "ERROR: backend did not produce a polished file."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	mkdir -p "$(dirname "$out_final")"
	cp -f -- "$produced" "$out_final"
	_polap_log1 "INFO: polished output → $out_final"

	conda deactivate || true
	_polap_log3 "Function end: $(echo "$FUNCNAME" | sed 's/_run_polap_//')"
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
	return 0
}

# Version: v0.3.0
function v3_run_polap_polish2 { # polish organelle genome sequences (ONT-only or Hybrid)
	# Enable debugging if _POLAP_DEBUG is set
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo "$FUNCNAME" | sed 's/_run_polap_//')"
	local polap_cmd="${FUNCNAME##*_}"

	# Verbosity routing
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose:-0}" -ge "${_polap_var_function_verbose:-2}" ] && _polap_output_dest="/dev/stderr"

	# Common vars
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# ---------- CLI help ----------
	local help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - polish a draft organelle assembly (FASTA or GFA) using ONT and/or short reads

Synopsis:
  polap ${polap_cmd} [options]

Description:
  Automatically selects an ONT-only or Hybrid (short-reads) polishing path and optionally preserves
  Flye GFA topology (circular paths) via extract → polish → reinject.

Options:
  -p FILE      Draft assembly (FASTA or GFA). [default: ${_arg_unpolished_fasta}]
  -f FILE      Output polished assembly (FASTA or GFA if --gfa-preserve). [default: ${_arg_final_assembly}]
  -l FILE      ONT reads (FASTQ[.gz]) for ONT-only polishing (Racon rounds).
  -s FILE      ONE short-read file (FASTQ[.gz]); supports interleaved/merged or single-end polishing.
  -1 FILE      Short-read R1 (FASTQ[.gz]); if provided with -2, takes precedence over -s.
  -2 FILE      Short-read R2 (FASTQ[.gz]).
  -g           Preserve GFA topology (if input -p is .gfa): extract S-seqs → polish → reinject.
  -m MODE      Mode: auto|ont|hybrid  (default: auto)
  -r INT       Racon rounds for ONT-only (default: ${_arg_racon_round:-3})
  -t INT       Threads (default: ${_arg_threads:-16})
  -k INT       k-mer for Merqury QV if short reads exist and --eval (default: 31)
  --eval       Run Merqury QV if short reads are available (writes QV report alongside outputs)
  --polisher   Hybrid polisher: auto|fmlrc2|polypolish (default: auto → fmlrc2)

Speed presets (passed to backend; ONT-only path):
  --speed {normal|fast|turbo}   normal=100×,3r,-w1000 ; fast=80×,2r,-w1200 ; turbo=60×,1r,-w1500
  --ont-cap INT                 override depth cap (×); 0 disables capping
  --ont-cap-method {rasusa|seqtk}  prefer rasusa; seqtk as fallback
  --minimap-preset {map-ont|map-hifi}  choose mapper profile for overlaps

Examples:
  ONT-only polishing (Racon×3):
    polap ${polap_cmd} -p pt.fa -l ont.fastq.gz -f pt.polished.fa

  Hybrid with a single short-read file (interleaved):
    polap ${polap_cmd} -p mt.fa -s illumina.fastq.gz -f mt.polished.fa --eval -k 31

  Hybrid with paired short reads (takes precedence over -s):
    polap ${polap_cmd} -p mt.fa -1 R1.fq.gz -2 R2.fq.gz -f mt.polished.fa --polisher polypolish

  Flye graph, preserve circular paths in GFA:
    polap ${polap_cmd} -p flye/assembly_graph.gfa -l ont.fastq.gz -g -f mito.polished.gfa

See Also:
  polap polish (backend): polap-bash-polish.sh
EOF
	)
	if [[ ${_arg_menu[1]:-} == "help" || "${_arg_help:-off}" == "on" ]]; then
		local manfile
		manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]:-polish}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# ---------- Resolve inputs (backward-compatible with existing _arg_* names) ----------
	local in_asm="${_arg_unpolished_fasta:-}"
	local out_final="${_arg_final_assembly:-}"
	local ont_fq="${_arg_long_reads:-}" # -l
	local sr_one=""                     # -s
	local sr_r1="" sr_r2=""             # -1/-2
	[[ "${_arg_short_read1_is:-off}" == "on" ]] && sr_one="${_arg_short_read1:-}" && sr_r1="${_arg_short_read1:-}"
	[[ "${_arg_short_read2_is:-off}" == "on" ]] && sr_r2="${_arg_short_read2:-}"

	# Aliases (compat)
	[[ -n "${_arg_l_fastq:-}" ]] && ont_fq="${_arg_l_fastq}"
	[[ -n "${_arg_s_fastq:-}" ]] && sr_one="${_arg_s_fastq}"
	[[ -n "${_arg_1_fastq:-}" ]] && sr_r1="${_arg_1_fastq}"
	[[ -n "${_arg_2_fastq:-}" ]] && sr_r2="${_arg_2_fastq}"
	[[ -n "${_arg_p_fasta:-}" ]] && in_asm="${_arg_p_fasta}"
	[[ -n "${_arg_f_fasta:-}" ]] && out_final="${_arg_f_fasta}"

	local racon_rounds="${_arg_racon_round:-3}"
	local threads="${_arg_threads:-16}"
	local mode="${_arg_mode:-auto}"                # auto|ont|hybrid
	local gfa_preserve="${_arg_gfa_preserve:-off}" # on|off
	local do_eval="${_arg_eval:-off}"              # on|off
	local kmer="${_arg_kmer:-31}"
	local polisher="${_arg_polisher:-auto}" # auto|fmlrc2|polypolish

	# New speed/preset knobs (optional in caller)
	local speed="${_arg_speed:-normal}"             # normal|fast|turbo
	local ont_cap="${_arg_ont_cap:-}"               # empty -> backend default
	local ont_cap_method="${_arg_ont_cap_method:-}" # rasusa|seqtk
	local minimap_preset="${_arg_minimap_preset:-}" # map-ont|map-hifi

	# ---------- Conda env ----------
	_polap_lib_conda-ensure_conda_env polap-polish || return $EXIT_ERROR

	# ---------- Validate ----------
	if [[ -z "$in_asm" || ! -s "$in_asm" ]]; then
		_polap_log0 "ERROR: missing/empty draft assembly (-p): [$in_asm]"
		conda deactivate || true
		return $EXIT_ERROR
	fi
	if [[ -z "$out_final" ]]; then
		_polap_log0 "ERROR: missing final output path (-f)."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	# Resolve mode if auto
	if [[ "$mode" == "auto" ]]; then
		if [[ -n "$sr_r1" && -n "$sr_r2" ]]; then
			mode="hybrid"
		elif [[ -n "$sr_one" ]]; then
			mode="hybrid"
		elif [[ -n "$ont_fq" ]]; then
			mode="ont"
		else
			_polap_log0 "ERROR: mode=auto needs -l (ONT) or -s / -1 -2 (short reads)."
			conda deactivate || true
			return $EXIT_ERROR
		fi
	fi
	[[ "$mode" == "ont" && -z "$ont_fq" ]] && {
		_polap_log0 "ERROR: mode=ont requires -l."
		conda deactivate || true
		return $EXIT_ERROR
	}
	[[ "$mode" == "hybrid" && -z "$sr_one" && (-z "$sr_r1" || -z "$sr_r2") ]] && {
		_polap_log0 "ERROR: mode=hybrid requires -s OR -1/-2."
		conda deactivate || true
		return $EXIT_ERROR
	}

	# ---------- Build backend command ----------
	local backend="${_POLAPLIB_DIR}/polap-bash-polish.sh"
	if [[ ! -s "$backend" ]]; then
		_polap_log0 "ERROR: backend not found: $backend"
		conda deactivate || true
		return $EXIT_ERROR
	fi

	local prefix
	prefix="$(basename "${out_final%.*}")"
	local args=("polish" "$in_asm" "--threads" "$threads" "--prefix" "$prefix" "--mode" "$mode")

	# ONT path
	if [[ -n "$ont_fq" ]]; then
		args+=("--ont" "$ont_fq" "--racon-rounds" "$racon_rounds")
		# speed/preset knobs → backend (only add if set to avoid shadowing backend defaults)
		[[ -n "$speed" ]] && args+=("--speed" "$speed")
		[[ -n "$ont_cap" ]] && args+=("--ont-cap" "$ont_cap")
		[[ -n "$ont_cap_method" ]] && args+=("--ont-cap-method" "$ont_cap_method")
		[[ -n "$minimap_preset" ]] && args+=("--minimap-preset" "$minimap_preset")
	fi

	# Hybrid path
	if [[ "$mode" == "hybrid" ]]; then
		if [[ -n "$sr_r1" && -n "$sr_r2" ]]; then
			args+=("--sr1" "$sr_r1" "--sr2" "$sr_r2")
		elif [[ -n "$sr_one" ]]; then
			args+=("--sr1" "$sr_one") # backend treats single-file as interleaved/single-end
		fi
		args+=("--polisher" "$polisher")
	fi

	# GFA-preserve
	if [[ "${gfa_preserve}" =~ ^(on|true|1)$ ]]; then
		args+=("--gfa-preserve")
	fi

	# Merqury eval (only meaningful if SR present; backend will sanity-check)
	if [[ "${do_eval}" =~ ^(on|true|1)$ ]]; then
		args+=("--eval-merqury" "--kmer" "$kmer")
	fi

	_polap_log1 "INFO: backend cmd: $(printf "%q " "$backend" "${args[@]}")" >&2

	# ---------- Execute ----------
	# honor verbosity routing but keep user-facing logs when verbose
	if ! bash "$backend" "${args[@]}" >"${_polap_output_dest}" 2>&1; then
		_polap_log0 "ERROR: polishing failed (see backend log under out/${prefix}.polish_*/run.log)."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	# ---------- Locate and copy final artifact ----------
	local outdir produced
	outdir="$(ls -1dt out/"$prefix".polish_* 2>/dev/null | head -n1)"
	if [[ -z "$outdir" || ! -d "$outdir" ]]; then
		_polap_log0 "ERROR: cannot find backend output directory for prefix=$prefix"
		conda deactivate || true
		return $EXIT_ERROR
	fi

	if compgen -G "${outdir}/*.polished.gfa" >/dev/null; then
		produced="$(ls -1 ${outdir}/*.polished.gfa | head -n1)"
	elif compgen -G "${outdir}/*.polished.fa" >/dev/null; then
		produced="$(ls -1 ${outdir}/*.polished.fa | head -n1)"
	fi

	if [[ -z "$produced" || ! -s "$produced" ]]; then
		_polap_log0 "ERROR: backend did not produce a polished file."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	mkdir -p "$(dirname "$out_final")"
	cp -f -- "$produced" "$out_final"
	_polap_log1 "INFO: polished output → $out_final"

	conda deactivate || true
	_polap_log3 "Function end: $(echo "$FUNCNAME" | sed 's/_run_polap_//')"
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
	return 0
}

function v1_run_polap_polish2 { # polish organelle genome sequences (ONT-only or Hybrid)
	# Enable debugging if _POLAP_DEBUG is set
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo "$FUNCNAME" | sed 's/_run_polap_//')"
	local polap_cmd="${FUNCNAME##*_}"

	# Verbosity routing
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose:-0}" -ge "${_polap_var_function_verbose:-2}" ] && _polap_output_dest="/dev/stderr"

	# Common vars
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# ---------- CLI help ----------
	local help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - polish a draft organelle assembly (FASTA or GFA) using ONT and/or short reads

Synopsis:
  polap ${polap_cmd} [options]

Description:
  Automatically selects an ONT-only or Hybrid (short-reads) polishing path and optionally preserves
  Flye GFA topology (circular paths) via extract, polish, and reinject.

Options:
  -p FILE      Draft assembly (FASTA or GFA). [default: ${_arg_unpolished_fasta}]
  -f FILE      Output polished assembly (FASTA or GFA if --gfa-preserve). [default: ${_arg_final_assembly}]
  -l FILE      ONT reads (FASTQ[.gz]) for ONT-only polishing (Racon rounds).
  -s FILE      ONE short-read file (FASTQ[.gz]); supports interleaved/merged or single-end polishing.
  -1 FILE      Short-read R1 (FASTQ[.gz]); if provided together with -2, takes precedence over -s.
  -2 FILE      Short-read R2 (FASTQ[.gz]).
  -g           Preserve GFA topology (if input -p is .gfa): extract S-seqs →  polish →  reinject.
  -m MODE      Mode: auto|ont|hybrid  (default: auto)
  -r INT       Racon rounds for ONT-only (default: ${_arg_racon_round:-3})
  -t INT       Threads (default: ${_arg_threads:-16})
  -k INT       k-mer for Merqury QV if short reads exist and --eval (default: 31)
  --eval       Run Merqury QV if short reads are available (writes QV report alongside outputs)
  --polisher   Hybrid polisher: auto|fmlrc2|polypolish (default: auto → fmlrc2)

Examples:
  ONT-only polishing (Racon×3):
    polap ${polap_cmd} -p pt.fa -l ont.fastq.gz -f pt.polished.fa

  Hybrid with a single short-read file (interleaved):
    polap ${polap_cmd} -p mt.fa -s illumina.fastq.gz -f mt.polished.fa --eval -k 31

  Hybrid with paired short reads (takes precedence over -s):
    polap ${polap_cmd} -p mt.fa -1 R1.fq.gz -2 R2.fq.gz -f mt.polished.fa --polisher polypolish

  Flye graph, preserve circular paths in GFA:
    polap ${polap_cmd} -p flye/assembly_graph.gfa -l ont.fastq.gz -g -f mito.polished.gfa

See Also:
  polap polish (backend): polap-bash-polish.sh

Author:
  Sang Chul Choi (c) 2025
EOF
	)

	if [[ ${_arg_menu[1]:-} == "help" || "${_arg_help:-off}" == "on" ]]; then
		local manfile
		manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]:-polish}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# ---------- Resolve inputs from your global args (backward-compatible) ----------
	# Defaults from your existing variable conventions (override-able via your CLI parser)
	local in_asm="${_arg_unpolished_fasta:-}"
	local out_final="${_arg_final_assembly:-}"
	local ont_fq="${_arg_long_reads:-}" # maps to -l
	if [[ "${_arg_short_read1_is}" == "on" ]]; then
		local sr_one="${_arg_short_read1:-}" # maps to -s (single/interleaved)
		local sr_r1="${_arg_short_read1:-}"  # maps to -1
	else
		local sr_one=""
		local sr_r1=""
	fi
	if [[ "${_arg_short_read2_is}" == "on" ]]; then
		local sr_r2="${_arg_short_read2:-}" # maps to -1
	else
		local sr_r2=""
	fi
	local racon_rounds="${_arg_racon_round:-3}"
	local threads="${_arg_threads:-16}"
	local mode="${_arg_mode:-auto}" # auto|ont|hybrid
	local gfa_preserve="${_arg_gfa_preserve:-off}"
	local do_eval="${_arg_eval:-off}"
	local kmer="${_arg_kmer:-31}"
	local polisher="${_arg_polisher:-auto}" # auto|fmlrc2|polypolish

	# Fallbacks to satisfy original template flags if they were passed:
	[[ -n "${_arg_l_fastq:-}" ]] && ont_fq="${_arg_l_fastq}"
	[[ -n "${_arg_s_fastq:-}" ]] && sr_one="${_arg_s_fastq}"
	[[ -n "${_arg_1_fastq:-}" ]] && sr_r1="${_arg_1_fastq}"
	[[ -n "${_arg_2_fastq:-}" ]] && sr_r2="${_arg_2_fastq}"

	# If caller provided -p / -f directly in this function’s context, prefer those
	[[ -n "${_arg_p_fasta:-}" ]] && in_asm="${_arg_p_fasta}"
	[[ -n "${_arg_f_fasta:-}" ]] && out_final="${_arg_f_fasta}"

	# ---------- Conda env ----------
	_polap_lib_conda-ensure_conda_env polap-polish || exit 1

	# ---------- Validate inputs ----------
	if [[ -z "$in_asm" || ! -s "$in_asm" ]]; then
		_polap_log0 "ERROR: missing/empty draft assembly (-p): [$in_asm]"
		return $EXIT_ERROR
	fi
	if [[ -z "$out_final" ]]; then
		_polap_log0 "ERROR: missing final output path (-f)."
		return $EXIT_ERROR
	fi

	# Resolve mode automatically if requested
	if [[ "$mode" == "auto" ]]; then
		if [[ -n "$sr_r1" && -n "$sr_r2" ]]; then
			mode="hybrid"
		elif [[ -n "$sr_one" ]]; then
			mode="hybrid"
		elif [[ -n "$ont_fq" ]]; then
			mode="ont"
		else
			_polap_log0 "ERROR: cannot resolve mode=auto; provide -l (ONT) or -s / -1 -2 (short reads)."
			return $EXIT_ERROR
		fi
	fi
	# _polap_log0 "mode=$mode"

	# If user forced ont, ensure ONT is there
	if [[ "$mode" == "ont" && -z "$ont_fq" ]]; then
		_polap_log0 "ERROR: mode=ont requires ONT reads (-l)."
		return $EXIT_ERROR
	fi

	# If user forced hybrid, ensure SR exist (either single -s or paired -1/-2)
	if [[ "$mode" == "hybrid" && -z "$sr_one" && (-z "$sr_r1" || -z "$sr_r2") ]]; then
		_polap_log0 "ERROR: mode=hybrid requires -s OR both -1 and -2."
		return $EXIT_ERROR
	fi

	# ---------- Build backend command for polap-bash-polish.sh ----------
	local backend="${_POLAPLIB_DIR}/polap-bash-polish.sh"
	if [[ ! -s "$backend" ]]; then
		_polap_log0 "ERROR: backend not found: $backend"
		return $EXIT_ERROR
	fi

	local args
	args=("polish" "$in_asm" "--threads" "$threads" "--prefix" "$(basename "${out_final%.*}")")
	# choose output container type by suffix; actual file written by backend:
	#   - If input is .gfa and --gfa-preserve, backend outputs .polished.gfa
	#   - Else outputs .polished.fa; we’ll copy/move to -f after completion

	# mode
	args+=("--mode" "$mode")

	# ONT (optional)
	if [[ -n "$ont_fq" ]]; then
		args+=("--ont" "$ont_fq" "--racon-rounds" "$racon_rounds")
	fi

	# Short reads: prefer paired; else single (-s)
	if [[ "$mode" == "hybrid" ]]; then
		if [[ -n "$sr_r1" && -n "$sr_r2" ]]; then
			args+=("--sr1" "$sr_r1" "--sr2" "$sr_r2")
		elif [[ -n "$sr_one" ]]; then
			# For the backend, we pass single file via --sr1 and leave --sr2 empty,
			# polap-bash-polish.sh interprets that as single/interleaved mode in the chosen polisher.
			args+=("--sr1" "$sr_one")
		fi
		args+=("--polisher" "$polisher")
	fi

	# GFA-preserve
	if [[ "${gfa_preserve}" == "on" || "${gfa_preserve}" == "true" || "${gfa_preserve}" == "1" ]]; then
		args+=("--gfa-preserve")
	fi

	# Merqury eval (only if SR exist)
	if [[ "${do_eval}" == "on" || "${do_eval}" == "true" || "${do_eval}" == "1" ]]; then
		args+=("--eval-merqury" "--kmer" "$kmer")
	fi

	_polap_log1 "INFO: running backend: $(printf "%q " "$backend" "${args[@]}")" >&2

	# ---------- Execute ----------
	if ! bash "$backend" "${args[@]}" >"${_polap_output_dest}" 2>&1; then
		_polap_log0 "ERROR: polishing failed (see log)."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	# ---------- Move result to -f ----------
	# Backend places outputs under out/<prefix>.polish_<timestamp>/
	# We detect whether the final artifact is .polished.gfa or .polished.fa and then copy to -f.
	local outdir
	outdir="$(ls -1dt out/"$(basename "${out_final%.*}")".polish_* 2>/dev/null | head -n1)"
	if [[ -z "$outdir" || ! -d "$outdir" ]]; then
		_polap_log0 "ERROR: cannot find backend output directory for prefix=$(basename "${out_final%.*}")"
		conda deactivate || true
		return $EXIT_ERROR
	fi

	local produced=""
	if compgen -G "${outdir}/*.polished.gfa" >/dev/null; then
		produced="$(ls -1 ${outdir}/*.polished.gfa | head -n1)"
	elif compgen -G "${outdir}/*.polished.fa" >/dev/null; then
		produced="$(ls -1 ${outdir}/*.polished.fa | head -n1)"
	fi

	if [[ -z "$produced" || ! -s "$produced" ]]; then
		_polap_log0 "ERROR: backend did not produce a polished file."
		conda deactivate || true
		return $EXIT_ERROR
	fi

	mkdir -p "$(dirname "$out_final")"
	cp -f "$produced" "$out_final"

	_polap_log1 "INFO: polished output → $out_final"

	conda deactivate || true

	_polap_log3 "Function end: $(echo "$FUNCNAME" | sed 's/_run_polap_//')"
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
	return 0
}

# NOTE: the following disassemble has its own polishing script, which is using
# Version 1.
# Look for these lines in polap-cmd-disassemble.sh
# command time -v bash "${_POLAPLIB_DIR}"/polap-bash-build-msbwt.sh \
# command time -v fmlrc -p "${_arg_threads}" \
# The two lines above are for subsampling-polishing.
# The simple-polishing of disassemble uses _run_polap_polish the Version 1.
# Have a separate msbwt folder like msbwt2: polap-variables-common.sh
# local _polap_var_outdir_msbwt2_dir="${_polap_var_outdir}/msbwt2"
# local _polap_var_outdir_msbwt2="${_polap_var_outdir}/msbwt2/comp_msbwt2.npy"
# local _polap_var_outdir_msbwt2_tar_gz="${_polap_var_outdir}/msbwt2.tar.gz"

################################################################################
# 2025-05-27
#
# function _run_polap_prepare-polishing { # prepare the polishing using FMLRC
#   -> Version 1 has been used since the beginning of the polap development.
#
# function _run_polap_prepare-polishing_v2 { # prepare the polishing using FMLRC
#   -> Version 2 needs testing.
#
# function _run_polap_polish { # polish organelle genome sequences using FMLRC
#   -> Version 1 has been used since the beginning of the polap development.
#
# function _run_polap_polish_v2 { # polish organelle genome sequences using FMLRC
#   -> Version 2 needs testing.
#
# function _run_polap_prepare-polishing_sort_efficient { # prepare the polishing using FMLRC
#   -> Version 1 with sort changed.
#   -> The code is more original than Version 1 itself except the sort line of code.
#
################################################################################

################################################################################
# Prepares the polishing using FMLRC.
#
# 2025-05-16: add the memory and CPU tracking
#
# Arguments:
#   -a s1.fq
#   -b s2.fq
# Inputs:
#   s1.fq
#   s2.fq
# Outputs:
#   $${_arg_outdir}/msbwt
################################################################################

summarize_polish_memlog_file() {
	local _memlog_file="$1"
	local _summary_base="$2" # e.g., outdir/summary-polish.txt
	local _start_ts="${3:-}"
	local _end_ts="${4:-}"
	local _pipeline_log="${5:-}"

	local _timestamp
	_timestamp=$(date +"%Y%m%d_%H%M")
	local _summary_file="${_summary_base}.${_timestamp}.txt"

	local _memlog_basename
	_memlog_basename="$(basename "$_memlog_file")"

	awk -F',' -v start="$_start_ts" -v end="$_end_ts" -v file="$_memlog_basename" '
    function format_duration(seconds) {
      h = int(seconds / 3600)
      m = int((seconds % 3600) / 60)
      s = seconds % 60
      return sprintf("%d:%02d:%02d (%.2f h)", h, m, s, seconds / 3600)
    }

    NR == 1 { next }
    $2 > max {
      max = $2
      ts = $1
      load = $4
      disk = $5
      cmd = $6
    }
    END {
      cmd_str = "date -d @" ts " +\"%Y-%m-%d %H:%M:%S\""
      cmd_str | getline ts_human
      close(cmd_str)

      print "📊 Peak total PSS during polishing:"
      print "File:             " file
      print "Timestamp:        " ts " (" ts_human ")"
      printf "Peak PSS:         %d KB (%.2f GB)\n", max, max / 1048576

      if (start && end) {
        elapsed = end - start
        elapsed_str = format_duration(elapsed)
        print "Elapsed time:     " elapsed_str
      } else {
        print "Elapsed time:     (not available)"
      }

      print "CPU load (1min):  " load
      print "Disk free (GB):   " disk
      print "Max PSS command:  " cmd
    }
  ' "$_memlog_file" | tee "$_summary_file"

	if [[ -n "$_pipeline_log" && -r "$_pipeline_log" ]]; then
		echo -e "\n📋 Full pipeline log:" >>"$_summary_file"
		cat "$_pipeline_log" >>"$_summary_file"
	fi

	# Create/update symlink to point to the latest summary
	ln -sf "$(basename "$_summary_file")" "$_summary_base"

	echo "[INFO] Summary written to: $_summary_file"
	echo "[INFO] Symlink created:    $_summary_base -> $(basename "$_summary_file")"
}

# The version without memory-usage tracking
################################################################################
# Prepares the polishing using FMLRC.
#
# 2025-05-16: add the memory and CPU tracking
#
# Arguments:
#   -a s1.fq
#   -b s2.fq
# Inputs:
#   s1.fq
#   s2.fq
# Outputs:
#   $${_arg_outdir}/msbwt
################################################################################
function _run_polap_prepare-polishing { # prepare the polishing using FMLRC
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local filesize
	_polap_set-variables-short-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<EOF
Name:
  polap prepare-polishing - prepares the polishing using FMLRC.

Synopsis:
  polap prepare-polishing [options]

Description:
  polap prepare-polishing uses FMLRC to prepare the polishing of sequences in FASTQ format using short-read data.

Options:
  -a FASTQ
    short-read data file 1 [default: ${_arg_short_read1}]

  -b FASTQ
    short-read data file 2 [default: ${_arg_short_read2}]

Examples:
  Prepare for short-read polishing:
    polap prepare-polishing -a s1.fq -b s2.fq -o outdir

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
		ls -l "${_arg_outdir}/msbwt" >&2
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
		# Get the file size in bytes
		filesize=$(stat --format=%s "${_polap_var_outdir_msbwt}")
	else
		filesize=0
	fi

	if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
		_polap_log1_file "${_polap_var_outdir_msbwt_tar_gz}"
		if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
			_polap_log1_file "${_polap_var_outdir_msbwt}"
			_polap_log1 "  skipping the short-read polishing preparation."
		else
			tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
		fi
	elif ((filesize > 1024)); then
		# Check if the file size is greater than 100 KB (100 * 1024 bytes)
		_polap_log1_file "${_polap_var_outdir_msbwt}"
		_polap_log1 "  skipping the short-read polishing preparation."
	else

		# Initialize Conda
		# _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
		# source $HOME/miniconda3/etc/profile.d/conda.sh
		# conda activate polap-fmlrc
		#
		# source $HOME/miniconda3/bin/activate polap-fmlrc

		_polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
		if ! run_check2; then
			echoerr "ERROR: change your conda environment to polap-fmlrc."
			echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
			echoerr "INFO: (base) $ conda activate polap-fmlrc"
			exit $EXIT_ERROR
		fi

		# check_file_existence "${_arg_short_read1}"

		_polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"

		# one or two files
		# fq or gzipped fq
		# exist or not
		#!/bin/bash

		_polap_log1 "  decompressing short-read data if necessary..."
		if [[ -s "${_arg_short_read1}" ]]; then
			_arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
		else
			_arg_short_read1=""
		fi
		if [[ -s "${_arg_short_read2}" ]]; then
			_arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
		else
			_arg_short_read2=""
		fi
		if [[ -n "${_arg_short_read1}" ]]; then
			_polap_log1 "    short-read1: ${_arg_short_read1}"
		else
			_polap_log1 "    short-read1: no such data file"
		fi
		if [[ -n "${_arg_short_read2}" ]]; then
			_polap_log1 "    short-read2: ${_arg_short_read2}"
		else
			_polap_log1 "    short-read2: no such data file"
		fi

		if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
			_polap_log1 "    short-read1 file: ${_arg_short_read1}"
		else
			_polap_log1 "    short-read1: no fastq or fq file: ${_arg_short_read1}"
		fi

		if [[ ${_arg_short_read2} = *.fastq || ${_arg_short_read2} = *.fq ]]; then
			_polap_log1 "    short-read2 file: ${_arg_short_read2}"
		else
			_polap_log1 "    short-read2: no fastq or fq file: ${_arg_short_read2}"
		fi

		rm -rf "${_arg_outdir}/msbwt"

		# 2025-05-27
		# Version 1
		# mkdir -p "${_polap_var_outdir_msbwt_dir}"
		if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
			_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
			cat "${_arg_short_read1}" "${_arg_short_read2:-/dev/null}" |
				awk 'NR % 4 == 2' | sort | tr NT TN |
				ropebwt2 -LR 2>"${_polap_output_dest}" |
				tr NT TN |
				msbwt convert "${_arg_outdir}"/msbwt \
					>/dev/null 2>&1
		elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
			_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
			zcat "${_arg_short_read1}" "${_arg_short_read2}" |
				awk 'NR % 4 == 2' | sort | tr NT TN |
				ropebwt2 -LR 2>"${_polap_output_dest}" |
				tr NT TN |
				msbwt convert "${_arg_outdir}"/msbwt \
					>/dev/null 2>&1
		else
			_polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
			_polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
		fi

		# Version 2
		# https://github.com/HudsonAlpha/rust-msbwt
		# mkdir -p "${_polap_var_outdir_msbwt_dir}"
		# if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
		# 	_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
		#
		# 	msbwt2-build \
		# 		-o "${_polap_var_outdir_msbwt}" \
		# 		"${_arg_short_read1}" "${_arg_short_read2:-}"
		#
		# elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
		# 	_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
		#
		# 	msbwt2-build \
		# 		-o "${_polap_var_outdir_msbwt}" \
		# 		"${_arg_short_read1}" "${_arg_short_read2:-}"
		#
		# else
		# 	_polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
		# 	_polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
		# fi

		conda deactivate
	fi

	_polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_prepare-polishing_v2 { # prepare the polishing using FMLRC
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local filesize
	_polap_set-variables-short-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Prepares the polishing using FMLRC.
# Arguments:
#   -a ${_arg_short_read1}
#   -b ${_arg_short_read2}
#   or
#   --bioproject ${_arg_bioproject}
# Inputs:
#   ${_arg_short_read1}
#   ${_arg_short_read2}
# Outputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
# Precondition:
#   get-bioproject --bioproject ${_arg_bioproject}
Example: $(basename "$0") ${_arg_menu[0]} -a ${_arg_short_read1} [-b ${_arg_short_read2}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		ls -l "${_arg_outdir}/msbwt" >&2
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
		# Get the file size in bytes
		filesize=$(stat --format=%s "${_polap_var_outdir_msbwt}")
	else
		filesize=0
	fi

	if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
		_polap_log1_file "${_polap_var_outdir_msbwt_tar_gz}"
		if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
			_polap_log1_file "${_polap_var_outdir_msbwt}"
			_polap_log1 "  skipping the short-read polishing preparation."
		else
			tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
		fi
	elif ((filesize > 1024)); then
		# Check if the file size is greater than 100 KB (100 * 1024 bytes)
		_polap_log1_file "${_polap_var_outdir_msbwt}"
		_polap_log1 "  skipping the short-read polishing preparation."
	else

		# File preparation
		# check_file_existence "${_arg_short_read1}"

		# one or two files
		# fq or gzipped fq
		# exist or not
		#!/bin/bash

		_polap_log1 "  decompressing short-read data if necessary..."
		if [[ -s "${_arg_short_read1}" ]]; then
			_arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
		else
			_arg_short_read1=""
		fi
		if [[ -s "${_arg_short_read2}" ]]; then
			_arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
		else
			_arg_short_read2=""
		fi
		if [[ -n "${_arg_short_read1}" ]]; then
			_polap_log1 "    short-read1: ${_arg_short_read1}"
		else
			_polap_log1 "    short-read1: no such data file"
		fi
		if [[ -n "${_arg_short_read2}" ]]; then
			_polap_log1 "    short-read2: ${_arg_short_read2}"
		else
			_polap_log1 "    short-read2: no such data file"
		fi

		if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
			_polap_log1 "    short-read1 file: ${_arg_short_read1}"
		else
			_polap_log1 "    short-read1: no fastq or fq file: ${_arg_short_read1}"
		fi

		if [[ ${_arg_short_read2} = *.fastq || ${_arg_short_read2} = *.fq ]]; then
			_polap_log1 "    short-read2 file: ${_arg_short_read2}"
		else
			_polap_log1 "    short-read2: no fastq or fq file: ${_arg_short_read2}"
		fi

		rm -rf "${_arg_outdir}/msbwt"

		# Initialize Conda
		_polap_lib_conda-ensure_conda_env polap-fmlrc2 || exit 1
		if ! run_check2_v2; then
			echoerr "ERROR: change your conda environment to polap-fmlrc."
			echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
			echoerr "INFO: (base) $ conda activate polap-fmlrc"
			exit $EXIT_ERROR
		fi

		_polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"

		# Version 2
		# https://github.com/HudsonAlpha/rust-msbwt
		mkdir -p "${_polap_var_outdir_msbwt_dir}"
		if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
			_polap_log1 "  creating ${_arg_outdir}/msbwt ..."

			msbwt2-build \
				-o "${_polap_var_outdir_msbwt}" \
				"${_arg_short_read1}" "${_arg_short_read2:-}"

		elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
			_polap_log1 "  creating ${_arg_outdir}/msbwt ..."

			msbwt2-build \
				-o "${_polap_var_outdir_msbwt}" \
				"${_arg_short_read1}" "${_arg_short_read2:-}"

		else
			_polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
			_polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
		fi

		conda deactivate
	fi

	_polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Polishes using FMLRC.
#
# Arguments:
#   -p mt.0.fasta
#   -f mt.1.fa
# Inputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
################################################################################
function _run_polap_polish { # polish organelle genome sequences using FMLRC
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"
	local polap_cmd="${FUNCNAME##*_}"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - polish a draft sequence

Synopsis:
  polap ${polap_cmd} [options]

Description:
  polap ${polap_cmd} uses either short-read data or long-read data to polish a draft genome assembly sequence in FASTA format.

Options:
  -p FASTA [default: ${_arg_unpolished_fasta}]
    draft genome assembly sequence file

  -f FASTA [default: ${_arg_final_assembly}]
    final genome assembly sequence file

  -l FASTQ
    reads data file

Examples:
  Get organelle genome sequences:
    polap ${polap_cmd} -a s1.fq -b s2.fq -p mt.0.fasta -f mt.1.fa

See Also:
  polap prepare-polishing - prepares the polishing using FMLRC.

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

	# Initialize Conda
	_polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
	if ! run_check2; then
		_polap_log0 "ERROR: change your conda environment to polap-fmlrc."
		_polap_log0 "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
		_polap_log0 "INFO: (base) $ conda activate polap-fmlrc"
		exit $EXIT_ERROR
	fi

	if [[ ! -s "${_polap_var_outdir_msbwt}" ]]; then
		_polap_log0 "ERROR: no msbwt at ${_polap_var_outdir_msbwt}"
		_polap_log0 "HINT: $(basename "$0") prepare-polishing [-a s1.fq] [-b s2.fq]"
		exit $EXIT_ERROR
	fi

	_polap_log1 "INFO: executing fmlrc on the draft sequence ${_arg_unpolished_fasta} ... be patient!"
	if [[ -s "${_arg_unpolished_fasta}" ]]; then
		# fmlrc -p "${_arg_threads}" "${_polap_var_outdir_msbwt}" "${_arg_unpolished_fasta}" "${_arg_final_assembly}" >/dev/null 2>&1
		# fmlrc -p "${_arg_threads}" "${_polap_var_outdir_msbwt}" "${_arg_unpolished_fasta}" "${_arg_final_assembly}" >/dev/null 2>&1

		# 2025-05-27
		# Version 1
		_polap_log3_pipe "fmlrc -p ${_arg_threads_fmlrc} ${_polap_var_outdir_msbwt} ${_arg_unpolished_fasta} ${_arg_final_assembly} >/dev/null 2>&1"

	else
		_polap_log0 "ERROR: no unpolished fasta file: [${_arg_unpolished_fasta}]"
		exit $EXIT_ERROR
	fi

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Polishes using FMLRC2.
#
# Arguments:
#   -p mt.0.fasta
#   -f mt.1.fa
# Inputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
################################################################################
function _run_polap_polish_v2 { # polish organelle genome sequences using FMLRC
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Polish a draft sequence using FMLRC2.
#
# Arguments:
#   -p ${_arg_unpolished_fasta}: a long-read draft genome assembly
#   -f ${_arg_final_assembly}: a final genome assembly sequence name
# Inputs:
#   ${_polap_var_outdir_msbwt}
#   ${_arg_unpolished_fasta}
# Outputs:
#   ${_arg_final_assembly}
Example: $(basename "$0") ${_arg_menu[0]} -p ${_arg_unpolished_fasta} -f ${_arg_final_assembly}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Initialize Conda
	_polap_lib_conda-ensure_conda_env polap-fmlrc2 || exit 1
	if ! run_check2_v2; then
		_polap_log0 "ERROR: change your conda environment to polap-fmlrc."
		_polap_log0 "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
		_polap_log0 "INFO: (base) $ conda activate polap-fmlrc"
		exit $EXIT_ERROR
	fi

	if [[ ! -s "${_polap_var_outdir_msbwt}" ]]; then
		_polap_log0 "ERROR: no msbwt at ${_polap_var_outdir_msbwt}"
		_polap_log0 "HINT: $(basename "$0") prepare-polishing [-a s1.fq] [-b s2.fq]"
		exit $EXIT_ERROR
	fi

	_polap_log1 "INFO: executing fmlrc2 on the draft sequence ${_arg_unpolished_fasta} ... be patient!"
	if [[ -s "${_arg_unpolished_fasta}" ]]; then

		# Version 2
		# https://github.com/HudsonAlpha/fmlrc2
		fmlrc2 \
			-t "${_arg_threads_fmlrc}" \
			${_polap_var_outdir_msbwt} \
			"${_arg_unpolished_fasta}" \
			"${_arg_final_assembly}"

	else
		_polap_log0 "ERROR: no unpolished fasta file: [${_arg_unpolished_fasta}]"
		exit $EXIT_ERROR
	fi

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# 2025-05-27
# This used to be original version 1 with sort changed in the main script.
# The other versions are modified so that we could keep two versions.
# sort --buffer-size=1G --temporary-directory=/tmp |
function _run_polap_prepare-polishing_sort_efficient { # prepare the polishing using FMLRC
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local filesize
	_polap_set-variables-short-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Prepares the polishing using FMLRC.
# Arguments:
#   -a ${_arg_short_read1}
#   -b ${_arg_short_read2}
#   or
#   --bioproject ${_arg_bioproject}
# Inputs:
#   ${_arg_short_read1}
#   ${_arg_short_read2}
# Outputs:
#   ${_arg_outdir}/msbwt/comp_msbwt.npy
# Precondition:
#   get-bioproject --bioproject ${_arg_bioproject}
Example: $(basename "$0") ${_arg_menu[0]} -a ${_arg_short_read1} [-b ${_arg_short_read2}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		ls -l "${_arg_outdir}/msbwt" >&2
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
		# Get the file size in bytes
		filesize=$(stat --format=%s "${_polap_var_outdir_msbwt}")
	else
		filesize=0
	fi

	if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
		_polap_log1_file "${_polap_var_outdir_msbwt_tar_gz}"
		if [[ -s "${_polap_var_outdir_msbwt}" ]]; then
			_polap_log1_file "${_polap_var_outdir_msbwt}"
			_polap_log1 "  skipping the short-read polishing preparation."
		else
			tar -zxf "${_polap_var_outdir_msbwt_tar_gz}" -C "${_arg_outdir}"
		fi
	elif ((filesize > 1024)); then
		# Check if the file size is greater than 100 KB (100 * 1024 bytes)
		_polap_log1_file "${_polap_var_outdir_msbwt}"
		_polap_log1 "  skipping the short-read polishing preparation."
	else

		# Initialize Conda
		source $HOME/miniconda3/etc/profile.d/conda.sh
		conda activate polap-fmlrc
		# source $HOME/miniconda3/bin/activate polap-fmlrc

		if ! run_check2; then
			echoerr "ERROR: change your conda environment to polap-fmlrc."
			echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
			echoerr "INFO: (base) $ conda activate polap-fmlrc"
			exit $EXIT_ERROR
		fi

		# check_file_existence "${_arg_short_read1}"

		_polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"

		# one or two files
		# fq or gzipped fq
		# exist or not
		#!/bin/bash

		_polap_log1 "  decompressing short-read data if necessary..."
		if [[ -s "${_arg_short_read1}" ]]; then
			_arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
		else
			_arg_short_read1=""
		fi
		if [[ -s "${_arg_short_read2}" ]]; then
			_arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
		else
			_arg_short_read2=""
		fi
		if [[ -n "${_arg_short_read1}" ]]; then
			_polap_log1 "    short-read1: ${_arg_short_read1}"
		else
			_polap_log1 "    short-read1: no such data file"
		fi
		if [[ -n "${_arg_short_read2}" ]]; then
			_polap_log1 "    short-read2: ${_arg_short_read2}"
		else
			_polap_log1 "    short-read2: no such data file"
		fi

		if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
			_polap_log1 "    short-read1 file: ${_arg_short_read1}"
		else
			_polap_log1 "    short-read1: no fastq or fq file: ${_arg_short_read1}"
		fi

		if [[ ${_arg_short_read2} = *.fastq || ${_arg_short_read2} = *.fq ]]; then
			_polap_log1 "    short-read2 file: ${_arg_short_read2}"
		else
			_polap_log1 "    short-read2: no fastq or fq file: ${_arg_short_read2}"
		fi

		rm -rf "${_arg_outdir}/msbwt"
		if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
			_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
			cat "${_arg_short_read1}" "${_arg_short_read2:-/dev/null}" |
				awk 'NR % 4 == 2' |
				sort --buffer-size=1G --temporary-directory=/tmp |
				tr NT TN |
				ropebwt2 -LR 2>"${_polap_output_dest}" |
				tr NT TN |
				msbwt convert "${_arg_outdir}"/msbwt \
					>/dev/null 2>&1
		elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
			_polap_log1 "  creating ${_arg_outdir}/msbwt ..."
			zcat "${_arg_short_read1}" "${_arg_short_read2}" |
				awk 'NR % 4 == 2' |
				sort --buffer-size=1G --temporary-directory=/tmp |
				tr NT TN |
				ropebwt2 -LR 2>"${_polap_output_dest}" |
				tr NT TN |
				msbwt convert "${_arg_outdir}"/msbwt \
					>/dev/null 2>&1
		else
			_polap_log0 "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
			_polap_log0 "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
		fi

		conda deactivate
	fi

	_polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
