#!/usr/bin/env bash
# polap-bash-pt-remove-only.sh
# Version: v0.1.0
#
# One-shot wrapper to run **PT removal only** (no nuclear deplete).
# It forwards the relevant knobs to polap-bash-remove-ptdna-reads.sh and
# intentionally omits the --nuclear-deplete path.
#
# You can optionally override the auto identity cutoff with --identity-min.
# If you have weak labels (--pt-origin, --nuc-origin), the underlying
# helper will auto-tune the identity cutoff from your PAFs.
#
# Required:
#   -r, --reads      FASTQ(.gz) long reads
#   -p, --pt-ref     plastid reference (FASTA, single contig)
#
# Optional:
#   -o, --outdir     output dir                     [ptfilter.out]
#   -t, --threads    threads                        [16]
#   --pt-origin      PT-origin read IDs (weak labels)
#   --nuc-origin     Nuclear-origin read IDs (weak labels)
#   --alen-min       aligned length guard           [1000]
#   --identity-min   override auto identity cutoff  (e.g., 0.93)
#   --fpr            nuclear FPR target for ident   [0.01]
#   --tpr            PT TPR target for ident        [0.95]
#   --dry            log only; no execution
#   -v, --verbose    verbose notes
#   --quiet          quiet
#   --profiling      write per-step timings (profile.tsv)
#
# Outputs (in OUTDIR):
#   reads.nonpt.fq.gz
#   pt.ids, pt_thresh.vars, pt_thresh.diag.tsv
#
set -euo pipefail

# ─────────── logging ───────────
VERBOSE=1
QUIET=0
DRY=0
LOG_FILE=""
note() {
	[[ $QUIET -eq 1 ]] && return 0
	local ts src ln
	ts="$(date +'%F %T')"
	src="${BASH_SOURCE[1]##*/}"
	ln="${BASH_LINENO[0]}"
	printf "[%s][%s:%s] %s\n" "$ts" "$src" "$ln" "$*" |
		{ if [[ -n "$LOG_FILE" ]]; then tee -a "$LOG_FILE"; else cat; fi; } >&2
}

# ─────────── CLI ───────────
READS=""
PT_REF=""
OUTDIR="ptfilter.out"
THREADS=16
PT_ORIGIN=""
NUC_ORIGIN=""
ALEN_MIN=1000
IDENTITY_MIN=""
FPR=0.01
TPR=0.95
PROF=0

usage() {
	cat <<EOF
polap-bash-pt-remove-only.sh v0.1.0
Usage:
  $0 -r reads.fq.gz -p plastid.fa -o outdir [options]

Required:
  -r, --reads          FASTQ(.gz) long reads
  -p, --pt-ref         plastid reference (FASTA)

Optional:
  -o, --outdir         output directory            [${OUTDIR}]
  -t, --threads        threads                     [${THREADS}]
  --pt-origin          PT-origin IDs (weak labels)
  --nuc-origin         Nuclear-origin IDs (weak labels)
  --alen-min           aligned length guard        [${ALEN_MIN}]
  --identity-min       override identity cutoff    (e.g., 0.93)
  --fpr                nuclear FPR target          [${FPR}]
  --tpr                PT TPR target               [${TPR}]
  --dry                log only; no execution
  -v, --verbose        verbose notes
  --quiet              quiet
  --profiling          per-step timings to log/profile.tsv
EOF
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	-r | --reads)
		READS="$2"
		shift 2
		;;
	-p | --pt-ref)
		PT_REF="$2"
		shift 2
		;;
	-o | --outdir)
		OUTDIR="$2"
		shift 2
		;;
	-t | --threads)
		THREADS="$2"
		shift 2
		;;
	--pt-origin)
		PT_ORIGIN="$2"
		shift 2
		;;
	--nuc-origin)
		NUC_ORIGIN="$2"
		shift 2
		;;
	--alen-min)
		ALEN_MIN="$2"
		shift 2
		;;
	--identity-min)
		IDENTITY_MIN="$2"
		shift 2
		;;
	--fpr)
		FPR="$2"
		shift 2
		;;
	--tpr)
		TPR="$2"
		shift 2
		;;
	--dry)
		DRY=1
		shift
		;;
	-v | --verbose)
		VERBOSE=$((VERBOSE + 1))
		shift
		;;
	--quiet)
		QUIET=1
		shift
		;;
	--profiling)
		PROF=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		note "ERR unknown arg: $1"
		usage
		exit 2
		;;
	esac
done

[[ -s "$READS" ]] || {
	note "ERR missing --reads"
	exit 2
}
[[ -s "$PT_REF" ]] || {
	note "ERR missing --pt-ref"
	exit 2
}
mkdir -p "$OUTDIR" "$OUTDIR/log"
LOG_FILE="$OUTDIR/log/pt_remove_only.log"
((DRY)) && note "--dry on (no execution)"

# ─────────── locate the main script ───────────
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${_POLAPLIB_DIR}/polap-bash-remove-ptdna-reads.sh"
[[ -x "$MAIN" || -s "$MAIN" ]] || {
	note "ERR missing main: $MAIN"
	exit 127
}

# ─────────── build arg list for the main (PT removal only) ───────────
args=("$MAIN" -r "$READS" -p "$PT_REF" -o "$OUTDIR" -t "$THREADS")
((DRY)) && args+=(--dry)
((VERBOSE)) && args+=(-v)
((QUIET)) && args+=(--quiet)
((PROF)) && args+=(--profiling)

# PT removal knobs
[[ -n "$PT_ORIGIN" ]] && args+=(--pt-origin "$PT_ORIGIN")
[[ -n "$NUC_ORIGIN" ]] && args+=(--nuc-origin "$NUC_ORIGIN")
args+=(--alen-min "$ALEN_MIN")
[[ -n "$IDENTITY_MIN" ]] && args+=(--identity-min "$IDENTITY_MIN")
args+=(--fpr "$FPR" --tpr "$TPR")

# Intentionally DO NOT add --nuclear-deplete
# (No BUSCO/miniprot/overlapness invoked by this wrapper)

# ─────────── run ───────────
note "PT-only wrapper -> running: ${args[*]}"
"${args[@]}"

note "PT-only done. Output: ${OUTDIR}/reads.nonpt.fq.gz (and pt.ids, pt_thresh.*)"
