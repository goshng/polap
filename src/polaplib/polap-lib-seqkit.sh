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
# Convert numbers between different units.
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

# _polap_lib_seqkit_stats  (note: Bash identifiers cannot contain '-', so we use underscores)
# Version: v0.1.6
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Minimal, screen-only progress wrapper for:
#   seqkit stats -Ta
#
# Features:
#   - Shows a pv progress bar to the terminal (stderr) with a real ETA using the
#     *target* file size (follows symlinks).
#   - Optional --threads N to pass through to `seqkit stats -j N`.
#   - Clean CLI options; no logs, no GNU time.
#
# Examples:
#   _polap_lib_seqkit_stats \
#     --in  "${_brg_tmpdir}/s_1.fq.gz" \
#     --out "${_brg_rundir}/s_1.fq.seqkit.stats.ta.txt" \
#     --label "s_1.fq.gz" \
#     --threads 8
#
#   # Disable the progress bar (no pv):
#   _polap_lib_seqkit_stats --in s_1.fq.gz --out s_1.stats.txt --no-pv
#
_polap_lib_seqkit-stats() {
	local FQ="" OUT="" LABEL=""
	local NO_PV=0
	local THREADS=4

	_pls_usage_seqkit_stats() {
		cat <<'EOF' 1>&2
Usage:
  _polap_lib_seqkit_stats --in FILE.fq[.gz] --out stats.ta.txt
                          [--label STR] [--threads N] [--no-pv] [-h|--help]

Options:
  --in FILE      Input FASTQ(.gz); may be a symlink (size is resolved on target)
  --out FILE     Output table path for seqkit (-o FILE)
  --label STR    pv progress label (default: basename of --in)
  --threads N    Pass -j N to 'seqkit stats' (0 = omit -j; default 0)
  --no-pv        Disable pv (run seqkit directly; no progress bar)
  -h, --help     Show this help
EOF
	}

	# ---- parse options ----
	while (($#)); do
		case "$1" in
		--in)
			FQ="${2:?}"
			shift 2
			;;
		--out)
			OUT="${2:?}"
			shift 2
			;;
		--label)
			LABEL="${2:?}"
			shift 2
			;;
		--threads)
			THREADS="${2:?}"
			shift 2
			;;
		--no-pv)
			NO_PV=1
			shift
			;;
		-h | --help)
			_pls_usage_seqkit_stats
			return 0
			;;
		*)
			echo "[_polap_lib_seqkit_stats] ERROR: unknown option: $1" >&2
			_pls_usage_seqkit_stats
			return 2
			;;
		esac
	done

	# ---- validate & defaults ----
	[[ -n "$FQ" ]] || {
		echo "[_polap_lib_seqkit_stats] --in required" >&2
		return 2
	}
	[[ -n "$OUT" ]] || {
		echo "[_polap_lib_seqkit_stats] --out required" >&2
		return 2
	}
	[[ -e "$FQ" ]] || {
		echo "[_polap_lib_seqkit_stats] input not found: $FQ" >&2
		return 2
	}
	[[ -r "$FQ" ]] || {
		echo "[_polap_lib_seqkit_stats] input not readable: $FQ" >&2
		return 2
	}
	[[ -n "$LABEL" ]] || LABEL="$(basename -- "$FQ")"
	[[ "$THREADS" =~ ^[0-9]+$ ]] || {
		echo "[_polap_lib_seqkit_stats] --threads must be an integer" >&2
		return 2
	}

	# ---- resolve symlink to target for accurate pv -s ----
	_realpath_portable() {
		local p="$1"
		if command -v realpath >/dev/null 2>&1; then
			realpath -- "$p"
		elif command -v python3 >/dev/null 2>&1; then
			python3 - "$p" <<'PY'
import os,sys
print(os.path.realpath(sys.argv[1]))
PY
		elif readlink -f / >/dev/null 2>&1; then
			readlink -f -- "$p"
		else
			# fallback: resolve a few hops by hand
			local target="$p" link dir
			for _ in 1 2 3 4 5; do
				if link=$(readlink "$target" 2>/dev/null); then
					dir=$(cd "$(dirname -- "$target")" && pwd -P)
					target="${dir%/}/${link}"
				else
					break
				fi
			done
			printf '%s' "$(cd "$(dirname -- "$target")" && pwd -P)/$(basename -- "$target")"
		fi
	}

	local FQ_REAL
	FQ_REAL="$(_realpath_portable "$FQ" 2>/dev/null || printf '%s' "$FQ")"

	# ---- compute size of the target (portable) ----
	local SIZE=""
	if stat --version >/dev/null 2>&1; then
		SIZE="$(stat -Lc %s -- "$FQ_REAL" 2>/dev/null || true)" # GNU stat
	else
		SIZE="$(stat -f%z -L -- "$FQ_REAL" 2>/dev/null || stat -f%z -- "$FQ_REAL" 2>/dev/null || true)" # BSD/macOS
	fi

	# ---- build seqkit args (add -j if requested) ----
	local -a SK_ARGS=(stats -Ta)
	if ((THREADS > 0)); then
		SK_ARGS+=(-j "$THREADS")
	fi

	# ---- run with pv (or without) ----
	if ((NO_PV == 0)) && command -v pv >/dev/null 2>&1; then
		if [[ -n "$SIZE" ]]; then
			pv -N "$LABEL" -petbr -s "$SIZE" -- "$FQ" | seqkit "${SK_ARGS[@]}" - -o "$OUT"
		else
			pv -N "$LABEL" -petbr -- "$FQ" | seqkit "${SK_ARGS[@]}" - -o "$OUT"
		fi
	else
		# no progress bar
		seqkit "${SK_ARGS[@]}" "$FQ" -o "$OUT"
	fi
}
