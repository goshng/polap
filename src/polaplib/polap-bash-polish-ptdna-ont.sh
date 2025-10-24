#!/usr/bin/env bash
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

# polap-bash-polish-ptdna-ont.sh
export LC_ALL=C
set -euo pipefail

# Defaults
OUTDIR=""
THREADS=32
ROUNDS=3
VERBOSE=0
QUIET=0
DRYRUN=0

usage() {
	cat <<EOF
Usage: $(basename "$0") [options] <ptdna.fa> <reads.fastq[.gz]>

Required positional arguments:
  ptdna.fa           Plastid assembly FASTA
  reads.fastq[.gz]   ONT reads (FASTQ; gz OK)

Options:
  -o, --outdir DIR     Output directory (default: basename of ptdna.fa without extension)
  -t, --threads INT    Threads [default: 32]
  -n, --rounds INT     Racon rounds [default: 3]
  -v, --verbose        Verbose logging to stdout
  -q, --quiet          Suppress logging except errors
      --dry-run        Print planned commands and exit
  -h, --help           Show this help

Examples:
  $(basename "$0") pt.fa ont.fq.gz
  $(basename "$0") -o polish_out -t 48 -n 2 --dry-run pt.fa ont.fq.gz
EOF
}

# Parse options
ARGS=()
while [[ $# -gt 0 ]]; do
	case "$1" in
	-o | --outdir)
		OUTDIR="$2"
		shift 2
		;;
	-t | --threads)
		THREADS="$2"
		shift 2
		;;
	-n | --rounds)
		ROUNDS="$2"
		shift 2
		;;
	-v | --verbose)
		VERBOSE=1
		shift
		;;
	-q | --quiet)
		QUIET=1
		shift
		;;
	--dry-run)
		DRYRUN=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	--)
		shift
		break
		;;
	-*)
		echo "Unknown option: $1" >&2
		usage
		exit 1
		;;
	*)
		ARGS+=("$1")
		shift
		;;
	esac
done
# Append any remaining as positional
if [[ $# -gt 0 ]]; then ARGS+=("$@"); fi

# Positional args
if [[ ${#ARGS[@]} -lt 2 ]]; then
	echo "Error: need <ptdna.fa> and <reads.fastq[.gz]>." >&2
	usage
	exit 1
fi
ASM_IN="${ARGS[0]}"
READS="${ARGS[1]}"

# Default outdir from FASTA basename (strip common fasta extensions)
if [[ -z "$OUTDIR" ]]; then
	base="$(basename -- "$ASM_IN")"
	OUTDIR="${base%.*}"
	# Also strip a second extension if present (e.g., .fa.gz -> .fa)
	case "$OUTDIR" in
	*.fa | *.fna | *.fasta) OUTDIR="${OUTDIR%.*}" ;;
	esac
fi

# Logging
log() {
	local msg="[$(date '+%F %T')] $*"
	if [[ $QUIET -eq 0 ]]; then
		if [[ $VERBOSE -eq 1 ]]; then
			echo "$msg"
		else
			echo "$msg" >>"$OUTDIR/logs/run.log"
		fi
	fi
}

# Tools check
require() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing: $1" >&2
	exit 1
}; }

# Medaka model panel (covers common public ONT sets)
MEDAKA_MODELS=(
	r941_min_hac_g507
	r941_min_sup_g507
	r1041_e82_260bps_hac_v4.2.0
	r1041_e82_260bps_sup_v4.2.0
	r1041_e82_400bps_hac_v4.2.0
	r1041_e82_400bps_sup_v4.2.0
)

# If dry-run: print plan and exit
if [[ $DRYRUN -eq 1 ]]; then
	cat <<PLAN
[D R Y  R U N]
Input FASTA      : $ASM_IN
Input reads      : $READS
Output directory : $OUTDIR
Threads          : $THREADS
Racon rounds     : $ROUNDS
Medaka models    : ${MEDAKA_MODELS[*]}

Plan:
  mkdir -p "$OUTDIR"/{tmp,logs}
  cp "$ASM_IN" "$OUTDIR/tmp/pt.r0.fa"
  for i in 1..$ROUNDS:
    minimap2 -x map-ont -t $THREADS \$ASM_CUR $READS > \$PAF
    racon/build/bin/racon -t $THREADS $READS \$PAF \$ASM_CUR > \$NEXT_ASM
  medaka_consensus for each model -> pick best by samtools stats error rate
  final remap + samtools stats + optional VCF stats

Outputs:
  $OUTDIR/pt.racon.fa
  $OUTDIR/pt.best.fa
  $OUTDIR/logs/medaka_model_scores.tsv
  $OUTDIR/logs/samtools_stats.txt
PLAN
	exit 0
fi

# Real run
mkdir -p "$OUTDIR"/{tmp,logs}

for tool in minimap2 samtools racon medaka_consensus medaka_variant bcftools; do require "$tool"; done

map_bam() {
	local asm=$1 bam=$2
	minimap2 -x map-ont -a -t "$THREADS" "$asm" "$READS" |
		samtools sort -@ "$((THREADS / 4))" -o "$bam"
	samtools index "$bam"
}

racon_round() {
	local asm=$1 paf=$2 out=$3
	racon/build/bin/racon -t "$THREADS" "$READS" "$paf" "$asm" >"$out"
}

error_rate() {
	local bam=$1
	samtools stats "$bam" 2>/dev/null |
		awk -F'\t' '/^error rate:/ {print $3; exit}' |
		sed 's/%//g' | awk '{printf "%.6f\n", $1/100}'
}

try_medaka() {
	local asm=$1 model=$2 outdir=$3
	rm -rf "$outdir"
	mkdir -p "$outdir"
	set +e
	medaka_consensus -i "$READS" -d "$asm" -o "$outdir" -t "$THREADS" -m "$model" \
		>"$OUTDIR/logs/medaka_${model}.log" 2>&1
	local rc=$?
	set -e
	[[ $rc -ne 0 || ! -s "$outdir/consensus.fasta" ]] && echo "" || echo "$outdir/consensus.fasta"
}

error_rate() {
	local bam=$1
	local stats er
	stats=$(samtools stats "$bam" 2>/dev/null) || {
		echo "1"
		return
	}

	# Try direct "error rate:" first; if missing, compute from SN lines.
	er=$(awk -F'\t' '
    BEGIN{found=0; bm=0; mm=ins=del=0}
    /^error rate:/ {
      x=$3; gsub(/%/,"",x);
      if (x ~ /^[0-9.]+$/) {printf("%.6f\n", x/100); found=1}
    }
    /^SN\tbases mapped \(cigar\):/ {gsub(/,/,"",$3); bm=$3}
    /^SN\tnumber of mismatches:/   {gsub(/,/,"",$3); mm=$3}
    /^SN\tnumber of insertions:/   {gsub(/,/,"",$3); ins=$3}
    /^SN\tnumber of deletions:/    {gsub(/,/,"",$3); del=$3}
    END{
      if (found) exit;
      if (bm+0 > 0) {
        e=(mm+ins+del)/bm;
        if (e<0 || e!=e) e=1;   # NaN/negative guard
        printf("%.6f\n", e);
      }
    }' <<<"$stats")

	[[ -z "$er" ]] && er="1"
	echo "$er"
}

ASM_CUR="$OUTDIR/tmp/pt.r0.fa"
cp "$ASM_IN" "$ASM_CUR"

log "Starting $ROUNDS Racon round(s)…"
for i in $(seq 1 "$ROUNDS"); do
	log "Racon round $i"
	PAF="$OUTDIR/tmp/r$i.paf"
	minimap2 -x map-ont -t "$THREADS" "$ASM_CUR" "$READS" >"$PAF"
	RACON_OUT="$OUTDIR/tmp/pt.r$i.fa"
	racon_round "$ASM_CUR" "$PAF" "$RACON_OUT"
	ASM_CUR="$RACON_OUT"
done
cp "$ASM_CUR" "$OUTDIR/pt.racon.fa"

log "Medaka model sweep…"
BEST_ERR=1
BEST_ASM=""
BEST_MODEL=""
: >"$OUTDIR/logs/medaka_model_scores.tsv" || true
for model in "${MEDAKA_MODELS[@]}"; do
	log "Trying Medaka model: $model"
	CONS=$(try_medaka "$ASM_CUR" "$model" "$OUTDIR/tmp/medaka_$model")
	if [[ -z "$CONS" ]]; then
		log "  -> $model unavailable/failed; skipping."
		continue
	fi

	BAM="$OUTDIR/tmp/${model}.bam"
	map_bam "$CONS" "$BAM"

	ERR="$(error_rate "$BAM")"
	# Sanitize: if empty or non-numeric, force to 1
	if ! awk -v x="$ERR" 'BEGIN{exit (x ~ /^[0-9.]+$/)?0:1}'; then ERR="1"; fi

	log "  -> error rate: $ERR"
	echo -e "$model\t$ERR" >>"$OUTDIR/logs/medaka_model_scores.tsv"

	# Safe numeric comparison (no bare shell math on floats)
	if awk -v err="$ERR" -v best="$BEST_ERR" 'BEGIN{exit !(err+0 < best+0)}'; then
		BEST_ERR="$ERR"
		BEST_ASM="$CONS"
		BEST_MODEL="$model"
	fi
done

if [[ -z "$BEST_ASM" ]]; then
	log "No Medaka success; using Racon output"
	cp "$OUTDIR/pt.racon.fa" "$OUTDIR/pt.best.fa"
	FINAL="$OUTDIR/pt.best.fa"
else
	log "Best Medaka model: $BEST_MODEL (error rate=$BEST_ERR)"
	cp "$BEST_ASM" "$OUTDIR/pt.best.fa"
	FINAL="$OUTDIR/pt.best.fa"
fi

map_bam "$FINAL" "$OUTDIR/final.bam"
samtools stats "$OUTDIR/final.bam" >"$OUTDIR/logs/samtools_stats.txt" || true

log "Done. Final polished assembly: $FINAL"
log "Model scores: $OUTDIR/logs/medaka_model_scores.tsv"
log "Samtools stats: $OUTDIR/logs/samtools_stats.txt"
