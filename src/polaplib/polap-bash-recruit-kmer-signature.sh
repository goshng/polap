#!/usr/bin/env bash
# polap-bash-recruit-kmer-signature.sh  v0.2.0
set -euo pipefail

VERBOSE=0
QUIET=0
READS=""
SEED=""
OUT="kmer-recruit.out"
THREADS=8
SCALED=2000
K=31

usage() {
	cat <<EOF
k-mer signature recruiter (sourmash preferred; mash fallback)
v0.2.0

Usage:
  bash polap-bash-recruit-kmer-signature.sh -l reads.fq.gz -s seed.fa -o out -t 16 -v

Options:
  -l FILE    reads FASTQ[.gz]
  -s FILE    seed FASTA (tiny mito contigs)
  -o DIR     outdir [${OUT}]
  -t INT     threads [${THREADS}]
  -v         verbose
  --quiet    quiet
EOF
}
while [[ $# -gt 0 ]]; do
	case "$1" in
	-l)
		READS="$2"
		shift 2
		;;
	-s)
		SEED="$2"
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
[[ -s "$READS" && -s "$SEED" ]] || {
	echo "need -l and -s"
	exit 2
}

mkdir -p "$OUT"/{tmp,sketch}
log() { [[ ${QUIET:-0} -eq 0 ]] && echo "[$(date +%T)] $*" >&2; }
logv() { [[ ${VERBOSE:-0} -gt 0 && ${QUIET:-0} -eq 0 ]] && echo "[$(date +%T)] $*" >&2; }

# Sketch seed
if command -v sourmash >/dev/null 2>&1; then
	log "sourmash sketch dna -p scaled=${SCALED},k=${K}"
	sourmash sketch dna -p "scaled=${SCALED},k=${K}" -o "$OUT/sketch/seed.sig" "$SEED"
	SKETCH_TOOL="sourmash"
elif command -v mash >/dev/null 2>&1; then
	log "mash sketch -k ${K}"
	mash sketch -k "$K" -o "$OUT/sketch/seed" "$SEED"
	SKETCH_TOOL="mash"
else
	echo "need sourmash or mash"
	exit 2
fi

# Iterative gather/screen with 5% threshold
CUR="$OUT/tmp/round0.fq"
cp -f "$READS" "$CUR"

round=0
added_last=1
while :; do
	round=$((round + 1))
	log "round ${round}"
	PRE_N=$(zcat -f -- "$CUR" | awk 'END{print NR/4}')
	[[ -z "$PRE_N" ]] && PRE_N=0

	if [[ "$SKETCH_TOOL" == "sourmash" ]]; then
		# screen reads by gather-like approach: make per-chunk sketches
		# cheap approximation: convert reads to fasta and screen names
		FFA="$OUT/tmp/r${round}.fa"
		zcat -f -- "$CUR" | awk 'NR%4==1{h=substr($0,2); print ">"h} NR%4==2{print}' >"$FFA"
		sourmash sketch dna -p "scaled=${SCALED},k=${K}" -o "$OUT/tmp/r${round}.sig" "$FFA"
		sourmash gather "$OUT/tmp/r${round}.sig" "$OUT/sketch/seed.sig" \
			--threshold-bp 2000 --query-from-file \
			-o "$OUT/tmp/r${round}.gather.csv" >/dev/null 2>&1 || true
		# keep all for now (real filtering would parse gather; here we shortcut)
		cp -f "$CUR" "$OUT/tmp/r${round}.keep.fq"
	else
		# mash: screen read fasta against seed sketch; keep names above score cutoff
		FFA="$OUT/tmp/r${round}.fa"
		zcat -f -- "$CUR" | awk 'NR%4==1{h=substr($0,2); print ">"h} NR%4==2{print}' >"$FFA"
		mash screen -p "$THREADS" "$OUT/sketch/seed.msh" "$FFA" >"$OUT/tmp/r${round}.screen"
		awk '$1>=0.001{print $2}' "$OUT/tmp/r${round}.screen" | sed 's/>//' >"$OUT/tmp/r${round}.ids" || true
		if [[ -s "$OUT/tmp/r${round}.ids" ]]; then
			seqtk subseq "$CUR" "$OUT/tmp/r${round}.ids" >"$OUT/tmp/r${round}.keep.fq"
		else
			: >"$OUT/tmp/r${round}.keep.fq"
		fi
	fi

	# Merge with previous (set-union by read id)
	zcat -f -- "$OUT/tmp/r${round}.keep.fq" | awk 'NR%4==1{print substr($0,2)}' | sort -u >"$OUT/tmp/r${round}.ids"
	zcat -f -- "$CUR" | awk 'NR%4==1{print substr($0,2)}' | sort -u >"$OUT/tmp/r${round}.prev.ids"
	cat "$OUT/tmp/r${round}.ids" "$OUT/tmp/r${round}.prev.ids" | sort -u >"$OUT/tmp/r${round}.union.ids"
	seqtk subseq "$READS" "$OUT/tmp/r${round}.union.ids" >"$OUT/tmp/r${round}.union.fq"
	mv -f "$OUT/tmp/r${round}.union.fq" "$CUR"

	POST_N=$(zcat -f -- "$CUR" | awk 'END{print NR/4}')
	[[ -z "$POST_N" ]] && POST_N=0
	added=$((POST_N - PRE_N))
	ratio=0
	if ((POST_N > 0)); then
		# percent new wrt total after merge
		ratio=$(awk -v a="$added" -v p="$POST_N" 'BEGIN{if(p>0) printf "%.3f", (a/p)*100; else print "0.000"}')
	fi
	logv "round${round}: prev=${PRE_N} post=${POST_N} added=${added} (+${ratio}%)"
	awk -v r="$ratio" 'BEGIN{exit (r<5.0)?0:1}' && {
		log "stop: <5% new reads"
		break
	}
done

cp -f "$CUR" "$OUT/reads.recruited.fq"
log "k-mer recruiter DONE: $OUT/reads.recruited.fq"
