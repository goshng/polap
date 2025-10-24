#!/usr/bin/env bash
set -euo pipefail

# mm2_pv.sh — minimap2 with a pv progress bar for FASTQ inputs.
# Requirements: minimap2, pv; (optional) seqkit for fast/accurate read counting.

usage() {
	cat <<'USAGE'
Usage:
  mm2_pv.sh -r REF.fa|REF.mmi -i READS.fq[.gz] -o OUT.paf [options]

Options:
  -r REF     Reference FASTA or prebuilt index (.mmi)
  -i READS   FASTQ (gz or plain). If "-", read from stdin (progress w/o total)
  -o OUT     Output file (PAF by default; add -a via --sam)
  -x PRESET  minimap2 preset (default: map-ont). Examples: map-hifi, map-pb, ava-ont, asm20
  -t N       Threads (default: 8)
  --sam      Output SAM/BAM (adds -a; you can pipe to samtools yourself)
  --args "…" Extra minimap2 args to append (quoted string)
  --log LOG  Save minimap2 stderr to LOG (default: OUT.log)
  --est N    Estimated total reads (override counting; useful for stdin)
  -h|--help  This help

Notes:
* Progress uses pv on the decompressed FASTQ stream.
* If seqkit is available, read counting uses seqkit stats; otherwise counts '+' lines.
* When input is "-", no pre-count is possible, so pv shows rate/ETA but no percent.
USAGE
}

# Defaults
PRESET="map-ont"
THREADS=8
OUT=""
REF=""
IN=""
SAM=0
EXTRA_ARGS=""
EST_READS=""
LOG_OVERRIDE=""

# Parse args
while (($#)); do
	case "$1" in
	-r)
		REF="${2:?}"
		shift 2
		;;
	-i)
		IN="${2:?}"
		shift 2
		;;
	-o)
		OUT="${2:?}"
		shift 2
		;;
	-x)
		PRESET="${2:?}"
		shift 2
		;;
	-t)
		THREADS="${2:?}"
		shift 2
		;;
	--sam)
		SAM=1
		shift
		;;
	--args)
		EXTRA_ARGS="${2:-}"
		shift 2
		;;
	--log)
		LOG_OVERRIDE="${2:-}"
		shift 2
		;;
	--est)
		EST_READS="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown option: $1" >&2
		usage
		exit 2
		;;
	esac
done

# Validate
[[ -n "$REF" && -n "${IN:-}" && -n "$OUT" ]] || {
	usage
	exit 2
}
command -v minimap2 >/dev/null 2>&1 || {
	echo "[ERR] minimap2 not found" >&2
	exit 127
}
command -v pv >/dev/null 2>&1 || {
	echo "[ERR] pv not found" >&2
	exit 127
}

LOG="${LOG_OVERRIDE:-${OUT}.log}"

# Decompressor function
decompress() {
	local f="$1"
	if [[ "$f" == "-" ]]; then
		cat # already stdin
	elif [[ "$f" =~ \.gz$ ]]; then
		# Prefer gzip -dc; fall back to zcat/gunzip -c
		if command -v gzip >/dev/null 2>&1; then
			gzip -dc -- "$f"
		elif command -v zcat >/dev/null 2>&1; then
			zcat -- "$f"
		else gunzip -c -- "$f"; fi
	else
		cat -- "$f"
	fi
}

# Count total reads (FASTQ)
count_reads() {
	local f="$1"
	# If user provided an estimate, trust it.
	if [[ -n "${EST_READS:-}" ]]; then
		echo "$EST_READS"
		return 0
	fi

	# Can't pre-count stdin
	if [[ "$f" == "-" ]]; then
		echo ""
		return 0
	fi

	# Use seqkit if present (fast and robust)
	if command -v seqkit >/dev/null 2>&1; then
		# seqkit stats -T prints a TSV; column "num_seqs" is usually 4th
		seqkit stats -T --fq -- "$f" 2>/dev/null | awk 'NR==2{print $4}'
		return 0
	fi

	# Fallback: count '+' lines in FASTQ (works for plain or gz)
	if [[ "$f" =~ \.gz$ ]]; then
		zgrep -c '^+$' -- "$f" || true
	else
		grep -c '^+$' -- "$f" || true
	fi
}

# Build minimap2 command
MM2_ARGS=(-t "$THREADS" -x "$PRESET")
if ((SAM)); then MM2_ARGS+=(-a); fi
# If REF is .mmi or .fa, both acceptable to minimap2
MM2_ARGS+=("$REF" -)

# Add any user-specified extra args (split safely)
if [[ -n "$EXTRA_ARGS" ]]; then
	# shellcheck disable=SC2206
	EXTRA_ARR=($EXTRA_ARGS)
	MM2_ARGS=("${MM2_ARGS[@]}" "${EXTRA_ARR[@]}")
fi

# Pre-count reads for progress bar (if possible)
TOTAL_READS="$(count_reads "$IN" || true)"
# pv config: line mode (-l). For FASTQ, 4 lines per read -> set -s to 4*reads
PV_OPTS=(-l --eta --rate --progress)

if [[ -n "$TOTAL_READS" && "$TOTAL_READS" =~ ^[0-9]+$ ]]; then
	TOTAL_LINES=$((TOTAL_READS * 4))
	PV_OPTS+=(-s "$TOTAL_LINES")
	echo "[INFO] Total reads: $TOTAL_READS (FASTQ lines: $TOTAL_LINES)" >&2
else
	echo "[INFO] Total reads unknown (stdin or counting unavailable). Showing live rate only." >&2
fi

# Run
# - pv runs on the decompressed FASTQ stream
# - minimap2 stderr -> LOG; pv progress goes to stderr automatically
# - output -> OUT
{
	decompress "$IN" |
		pv "${PV_OPTS[@]}" |
		minimap2 "${MM2_ARGS[@]}"
} >"$OUT" 2>"$LOG"

echo "[OK] Done. Output: $OUT  |  minimap2 log: $LOG" >&2
