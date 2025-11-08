#!/usr/bin/env bash
# polaplib/polap-bash-test-crash-r.sh
# Version: v0.1.4 (stack-on-by-default)
#
# Parent launcher for R crash tests. Enables your Bash failsafe, then executes
# the R child with plain `Rscript --vanilla`. If the child exits non-zero, your
# v0.3.0 failsafe prints a concise call-site crash line (and, if you have the
# env/expansion variant, the line will include the conda env + expanded command).

# Re-exec under bash if called via /bin/sh
had_u=0
case $- in *u*) had_u=1 ;; esac
set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Enable your v0.3.0 failsafe in THIS process (prints call-site if child fails)
# shellcheck disable=SC1091
source "${POLAPLIB}/polap-lib-failsafe.sh"
polap_enable_failsafe

# ---------- Arguments we forward to the R child ----------
message="intentional test failure (R)"
code=7
infile=""
outfile=""
want_stack="on" # default ON; pass --no-stack to silence R stack

while [[ $# -gt 0 ]]; do
	case "$1" in
	--message)
		message="${2:-$message}"
		shift 2
		;;
	--code | --exit-code)
		code="${2:-$code}"
		shift 2
		;;
	--in)
		infile="${2:-$infile}"
		shift 2
		;;
	--out)
		outfile="${2:-$outfile}"
		shift 2
		;;
	--no-stack)
		want_stack="off"
		shift
		;;
	-h | --help)
		cat <<'EOF'
Usage: polaplib/polap-bash-test-crash-r.sh [--message STR] [--code N] [--in FILE] [--out FILE] [--no-stack]
  --no-stack   do NOT request a stack from the R child (stack is ON by default)
EOF
		exit 0
		;;
	--)
		shift
		break
		;;
	-*)
		echo "[ERROR] unknown option: $1" >&2
		exit 2
		;;
	*) break ;;
	esac
done

child_r="${POLAPLIB}/scripts/test_crash.R"

# Optional parent context line (useful for grep; your failsafe will print a nicer one)
printf '[%s:%s] parent(R): message=%s code=%s in=%s out=%s stack=%s\n' \
	"${BASH_SOURCE[0]##*/}" "${LINENO}" "$message" "$code" "${infile:-}" "${outfile:-}" "$want_stack" >&2

# Build child argv safely; use --flag=value to avoid quoting pitfalls for optparse
args=("--message=${message}" "--exit-code=${code}")
[[ -n "$infile" ]] && args+=("--in=${infile}")
[[ -n "$outfile" ]] && args+=("--out=${outfile}")
if [[ "$want_stack" != "off" ]]; then args+=("--stack"); fi

# Run child; if it exits non-zero, ERR trap prints our call-site header
Rscript --vanilla "$child_r" "${args[@]}"
