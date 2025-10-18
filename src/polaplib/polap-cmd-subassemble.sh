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

# Version: v0.1.0
# File   : polaplib/polap-subcmd-ptdna-subassemble.sh
# Role   : Polap dispatcher shim -> recruit → subsample → (optional) Flye for cpDNA
_run_polap_ptdna-subassemble() {
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: ${FUNCNAME[0]}"

	local polap_cmd="${FUNCNAME##*_}"
	local seeds="" reads=""
	local outdir="${_arg_outdir}"
	local cov="150" gsize="" method="auto" minlen="2000"
	local threads="${_POLAP_THREADS:-16}" window="2000" asm_cov="" no_assemble=0 keep_intermediate=0
	local mapper_opts=""

	local help_message=$(
		cat <<'USAGE'
Name:
  polap ptdna-subassemble — recruit ONT reads to cpDNA seeds, coverage-aware subsample, optional Flye.

Synopsis:
  polap ptdna-subassemble -i seeds.fa -j ont.fq[.gz] -o outdir
                          [--cov 150] [--gsize N] [--method auto|rasusa|filtlong|longest]
                          [--minlen 2000] [--threads 16] [--window 2000]
                          [--asm-coverage C] [--no-assemble] [--keep-intermediate]
                          [--mapper-opts "…"]

Description:
  1) minimap2 recruit (keep primary+secondary+supplementary) → 2) subsample to target coverage
  (rasusa/filtlong/longest fallback) → 3) pre/post coverage QC (mosdepth if present) →
  4) optional Flye assembly. Writes a single CSV of metrics.

Options:
  -i, --seeds         cpDNA seed contigs (FASTA).
  -j, --reads         ONT reads (FASTQ[.gz]).
  -o, --outdir        Output directory.
  --cov           Target coverage for subsampling (default 150).
  --gsize         Plastome size; if unset, estimated as sum(seeds) × 1.9 (IR heuristic).
  --method        auto|rasusa|filtlong|longest (default auto).
  --minlen        Min read length for selection (default 2000).
  --threads       Threads for mapping/assembly (default 16).
  --window        Mosdepth window size (default 2000 bp).
  --asm-coverage  Flye --asm-coverage (default: equals --cov).
  --no-assemble   Only recruit+subsample+QC; skip Flye.
  --keep-intermediate  Keep tmp BAMs and intermediates.
  --mapper-opts   Extra args appended to minimap2.

USAGE
	)

	set -- "${_arg_unknown_opts[@]}"

	# helpers
	_die() {
		_polap_log0 "[ERROR] $*"
		echo "$help_message" >&2
		return 2
	}
	_needval() { # usage: _needval <opt> <next>
		local _opt="$1" _next="${2-}"
		[[ $# -ge 2 ]] || _die "$_opt requires a value"
		[[ -n "$_next" && "$_next" != -- && ! "$_next" =~ ^- ]] || _die "$_opt requires a non-option value"
	}

	# arrays for pass-through (positional/unknown), if you want to keep them
	positionals=()
	unknown_opts=()

	while [[ $# -gt 0 ]]; do
		_polap_log0 "[INFO] arg: $1 (argc=$#)"
		case "$1" in
		view | help | -h | --help)
			echo "$help_message"
			return 0
			;;

		redo) shift ;; # explicit no-op

		--)
			shift
			positionals+=("$@")
			break
			;;

		-i | --seeds)
			_needval "$1" "${2-}"
			seeds="$2"
			shift 2
			;;

		-j | --reads)
			_needval "$1" "${2-}"
			reads="$2"
			shift 2
			;;

		-o | --outdir)
			_needval "$1" "${2-}"
			outdir="$2"
			shift 2
			;;

		--cov)
			_needval "$1" "${2-}"
			cov="$2"
			shift 2
			;;

		--gsize)
			_needval "$1" "${2-}"
			gsize="$2"
			shift 2
			;;

		--method)
			_needval "$1" "${2-}"
			method="$2"
			shift 2
			;;

		--minlen)
			_needval "$1" "${2-}"
			minlen="$2"
			shift 2
			;;

		--threads)
			_needval "$1" "${2-}"
			threads="$2"
			shift 2
			;;

		--window)
			_needval "$1" "${2-}"
			window="$2"
			shift 2
			;;

		--asm-coverage)
			_needval "$1" "${2-}"
			asm_cov="$2"
			shift 2
			;;

		--no-assemble)
			no_assemble=1
			shift
			;;

		--keep-intermediate)
			keep_intermediate=1
			shift
			;;

		--mapper-opts)
			_needval "$1" "${2-}"
			mapper_opts="$2"
			shift 2
			;;

		# unknown option starting with '-'  → error (no silent shift)
		-*)
			_die "Unknown option: $1"
			;;

		# positional arg
		*)
			positionals+=("$1")
			shift
			;;
		esac
	done

	# hard validation (avoid “silent” failures later)
	[[ -n "${seeds-}" ]] || _die "Missing required option: --seeds"
	[[ -n "${reads-}" ]] || _die "Missing required option: --reads"
	[[ -n "${outdir-}" ]] || _die "Missing required option: --outdir"

	# optional: existence checks with clear errors
	[[ -e "$seeds" ]] || _die "--seeds not found: $seeds"
	[[ -e "$reads" ]] || _die "--reads not found: $reads"

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	[[ -z "$seeds" || -z "$reads" || -z "$outdir" ]] && {
		echo "$help_message"
		return 2
	}

	_polap_log0 "1"

	local wrapper="${_POLAPLIB_DIR}/polap-bash-ptdna-recruit-subsample-flye.sh"
	if [[ ! -x "$wrapper" && ! -f "$wrapper" ]]; then
		_polap_log0 "[ERROR] wrapper missing: $wrapper"
		return 3
	fi

	_polap_log0 "2"

	local args=(--seeds "$seeds" --reads "$reads" --outdir "$outdir" --cov "$cov" --minlen "$minlen"
		--threads "$threads" --window "$window")
	[[ -n "$gsize" ]] && args+=(--gsize "$gsize")
	[[ -n "$asm_cov" ]] && args+=(--asm-coverage "$asm_cov")
	[[ "$no_assemble" -eq 1 ]] && args+=(--no-assemble)
	[[ "$keep_intermediate" -eq 1 ]] && args+=(--keep-intermediate)
	[[ -n "$mapper_opts" ]] && args+=(--mapper-opts "$mapper_opts")
	args+=(--method "$method")

	_polap_log1 "[RUN] $wrapper ${args[*]}"

	_polap_lib_conda-ensure_conda_env polap || exit 1
	_polap_log0 bash "$wrapper" "${args[@]}"
	bash "$wrapper" "${args[@]}"
	conda deactivate
}
