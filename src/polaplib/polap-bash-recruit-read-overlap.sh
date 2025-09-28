#!/usr/bin/env bash
# polap-bash-recruit-read-overlap.sh  v0.2.0
set -euo pipefail

VERBOSE=0
QUIET=0
READS=""
PANEL=""
OUT="ovlp-recruit.out"
THREADS=8
MIN_ID=0.70    # nmatch/alnlen
MIN_TSPAN=1000 # aligned span on target (tend-tstart)

usage() {
	cat <<EOF
read-overlap recruiter via minimap2 (iterative, 5% incremental stop)
v0.2.0

Usage:
  bash polap-bash-recruit-read-overlap.sh -l reads.fq.gz -p panel.fa -o out -t 16 -v

Options:
  -l FILE   reads FASTQ[.gz]
  -p FILE   panel FASTA (mito seeds, etc.)
  -o DIR    outdir [${OUT}]
  -t INT    threads [${THREADS}]
  -v        verbose
  --quiet   quiet
EOF
}
while [[ $# -gt 0 ]]; do
	case "$1" in
	-l)
		READS="$2"
		shift 2
		;;
	-p)
		PANEL="$2"
		shift 2
		;;
	-o)
		OUT="$2"
		shift 2
		;;
	-t)
		THREADS="$2"
		shift 2
		;;
	-v)
		VERBOSE=$((VERBOSE + 1))
		shift
		;;
	--quiet)
		QUIET=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERR: $1"
		usage
		exit 2
		;;
	esac
done
[[ -s "$READS" && -s "$PANEL" ]] || {
	echo "need -l and -p"
	exit 2
}

mkdir -p "$OUT/tmp"
log() { [[ ${QUIET:-0} -eq 0 ]] && echo "[$(date +%T)] $*" >&2; }
logv() { [[ ${VERBOSE:-0} -gt 0 && ${QUIET:-0} -eq 0 ]] && echo "[$(date +%T)] $*" >&2; }

CUR="$OUT/tmp/round0.fq"
cp -f "$READS" "$CUR"

round=0
while :; do
	round=$((round + 1))
	PRE_N=$(zcat -f -- "$CUR" | awk 'END{print NR/4}')
	[[ -z "$PRE_N" ]] && PRE_N=0

	PAF="$OUT/tmp/r${round}.paf"
	log "round ${round}: minimap2 overlap vs panel"
	minimap2 -t "$THREADS" -x map-ont --secondary=no "$PANEL" "$CUR" >"$PAF"

	# filter: id >= MIN_ID and target span >= MIN_TSPAN
	awk -v ID="$MIN_ID" -v TS="$MIN_TSPAN" '
    ($10>0 && $11>0) {
      id=$10/$11; tspan=$9-$8; if(tspan<0) tspan=-tspan;
      if(id>=ID && tspan>=TS) print $1
    }' "$PAF" | sort -u >"$OUT/tmp/r${round}.keep.ids"

	if [[ -s "$OUT/tmp/r${round}.keep.ids" ]]; then
		seqtk subseq "$READS" "$OUT/tmp/r${round}.keep.ids" >"$OUT/tmp/r${round}.keep.fq"
	else
		: >"$OUT/tmp/r${round}.keep.fq"
	fi

	# union with previous
	zcat -f -- "$OUT/tmp/r${round}.keep.fq" | awk 'NR%4==1{print substr($0,2)}' | sort -u >"$OUT/tmp/r${round}.new.ids"
	zcat -f -- "$CUR" | awk 'NR%4==1{print substr($0,2)}' | sort -u >"$OUT/tmp/r${round}.prev.ids"
	cat "$OUT/tmp/r${round}.new.ids" "$OUT/tmp/r${round}.prev.ids" | sort -u >"$OUT/tmp/r${round}.union.ids"
	seqtk subseq "$READS" "$OUT/tmp/r${round}.union.ids" >"$OUT/tmp/r${round}.union.fq"
	mv -f "$OUT/tmp/r${round}.union.fq" "$CUR"

	POST_N=$(zcat -f -- "$CUR" | awk 'END{print NR/4}')
	[[ -z "$POST_N" ]] && POST_N=0
	added=$((POST_N - PRE_N))
	ratio=$(awk -v a="$added" -v p="$POST_N" 'BEGIN{if(p>0) printf "%.3f", (a/p)*100; else print "0.000"}')
	logv "round${round}: prev=${PRE_N} post=${POST_N} added=${added} (+${ratio}%)"
	awk -v r="$ratio" 'BEGIN{exit (r<5.0)?0:1}' && {
		log "stop: <5% new reads"
		break
	}
done

cp -f "$CUR" "$OUT/reads.recruited.fq"
log "overlap recruiter DONE: $OUT/reads.recruited.fq"
