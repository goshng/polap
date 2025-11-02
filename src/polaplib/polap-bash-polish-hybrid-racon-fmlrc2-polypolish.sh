#!/usr/bin/env bash
# Version: v0.3.0
# Name: polap-bash-polish-hybrid-racon-fmlrc2-polypolish.sh
# Purpose: Polishing pipeline that handles three cases with one script:
#   1) ONT-only:          Racon xN using ONT overlaps
#   2) SR-only:           fmlrc2 (FM-index) -> Polypolish (Bowtie2 pileup)
#   3) ONT + SR (hybrid): Racon xN -> fmlrc2 -> Polypolish
#
# Final output is always: <outdir>/polished.fa
#
# Usage:
#   polap-bash-polish-hybrid-racon-fmlrc2-polypolish.sh \
#     --fasta draft.fa \
#     [--ont ont.fq.gz] \
#     [--sr1 R1.fq.gz --sr2 R2.fq.gz] \
#     --outdir outdir \
#     [--threads 32] [--rounds 2] \
#     [--ont-preset map-ont|map-hifi] \
#     [--min-ident 0.91] [--min-alen 3000] \
#     [--bt2-index PREFIX] [-v ...]
#
# Dependencies (conditional by stage):
#   - Racon stage (ONT):       minimap2, racon, python3
#   - fmlrc2 stage (SR):       ropebwt2, fmlrc2-convert, fmlrc2
#   - Polypolish stage (SR):   bowtie2, samtools, polypolish
# Common: bash, coreutils, awk, gzip

set -Eeuo pipefail

# ---------------- defaults ----------------
FASTA=""
ONT="" # optional
SR1="" # optional (requires SR2 if given)
SR2="" # optional
OUT=""
THREADS=32
ROUNDS=2
ONT_PRESET="map-ont" # or map-hifi
MIN_IDENT=0.91
MIN_ALEN=3000
BT2_INDEX=""
VERB=1
ROTATE_BETWEEN_ROUNDS=1 # off by default
ROTATE_POS=""           # optional explicit POS

POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAF_FILTER="${POLAPLIB_DIR}/scripts/polap-py-paf-filter.py" # S points to your scripts/ folder

# ---------------- logging ----------------
log() {
	local lvl="$1"
	shift
	((VERB > lvl)) && printf '[polish] %s\n' "$*" >&2 || true
}
die() {
	printf '[ERROR] %s\n' "$*" >&2
	exit 2
}

have_seqkit() { command -v seqkit >/dev/null 2>&1; }

is_single_contig() (
	awk 'BEGIN{n=0} /^>/{n++} END{print (n==1)?"YES":"NO"}' "$1"
)

contig_length() (
	# total length for a single-contig FASTA; headers ignored
	awk '/^>/{next}{L+=length($0)}END{print L+0}' "$1"
)

# ---------------- cli ----------------
print_help() {
	cat <<EOF
Usage: $0 --fasta draft.fa --outdir out/ [--ont ont.fq.gz] [--sr1 R1.fq.gz --sr2 R2.fq.gz] [options]
  --threads N         Default 32
  --rounds N          Racon rounds when ONT present (default 2)
  --ont-preset PRE    map-ont (default) or map-hifi
  --min-ident F       PAF filter nmatch/alen cutoff (default 0.91)
  --min-alen N        PAF filter aligned length cutoff (default 3000)
  --bt2-index PREF    Reuse a bowtie2 index (optional)
  -v                  Increase verbosity (repeatable)
  -h|--help           Show help
EOF
}

while (("$#")); do
	case "$1" in
	--rotate-between-rounds)
		ROTATE_BETWEEN_ROUNDS=1
		shift
		;;
	--rotate-pos)
		ROTATE_POS="${2:?}"
		shift 2
		;;
	--fasta)
		FASTA="${2:?}"
		shift 2
		;;
	--ont)
		ONT="${2:?}"
		shift 2
		;;
	--sr1)
		SR1="${2:?}"
		shift 2
		;;
	--sr2)
		SR2="${2:?}"
		shift 2
		;;
	--outdir)
		OUT="${2:?}"
		shift 2
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;
	--rounds)
		ROUNDS="${2:?}"
		shift 2
		;;
	--ont-preset)
		ONT_PRESET="${2:?}"
		shift 2
		;;
	--min-ident)
		MIN_IDENT="${2:?}"
		shift 2
		;;
	--min-alen)
		MIN_ALEN="${2:?}"
		shift 2
		;;
	--bt2-index)
		BT2_INDEX="${2:?}"
		shift 2
		;;
	-v)
		VERB=$((VERB + 1))
		shift
		;;
	-h | --help)
		print_help
		exit 0
		;;
	*) die "Unknown option: $1" ;;
	esac
done

# --------- validate I/O and decide mode ---------
[[ -n "$FASTA" && -n "$OUT" ]] || {
	print_help
	die "Required: --fasta and --outdir"
}
[[ -s "$FASTA" ]] || die "FASTA not found: $FASTA"

HAVE_ONT=0
[[ -n "$ONT" ]] && {
	[[ -s "$ONT" ]] || die "ONT not found: $ONT"
	HAVE_ONT=1
}

HAVE_SR=0
if [[ -n "$SR1" || -n "$SR2" ]]; then
	[[ -n "$SR1" && -n "$SR2" ]] || die "Both --sr1 and --sr2 are required for short-read mode"
	[[ -s "$SR1" && -s "$SR2" ]] || die "Short-read files not found"
	HAVE_SR=1
fi

((HAVE_ONT == 1 || HAVE_SR == 1)) || die "Provide at least one input: --ont or --sr1/--sr2"

# --------- dependency checks by stage ---------
if ((HAVE_ONT == 1)); then
	for t in minimap2 racon python3; do command -v "$t" >/dev/null 2>&1 || die "Missing dependency: $t"; done
fi
if ((HAVE_SR == 1)); then
	for t in ropebwt2 fmlrc2-convert fmlrc2 bowtie2 samtools polypolish; do
		command -v "$t" >/dev/null 2>&1 || die "Missing dependency: $t"
	done
fi

# --------- setup ---------
mkdir -p "$OUT"/{stage1,stage2,stage3,final,logs}
FINAL_OUT="$OUT/polished.fa"

gz_or_cat() { [[ "$1" =~ \.gz$ ]] && gzip -cd "$1" || cat "$1"; }

# Start from the provided FASTA
CUR="$FASTA"

# =========================
# Stage 1: Racon xN with ONT (optional)
# =========================
if ((HAVE_ONT == 1)); then
	log 1 "Stage 1: Racon x${ROUNDS} with ONT preset=${ONT_PRESET}"
	for r in $(seq 1 "$ROUNDS"); do
		RAW="$OUT/stage1/r${r}.raw.paf"
		FLT="$OUT/stage1/r${r}.flt.paf"
		OUTFA="$OUT/stage1/polished.r${r}.fa"

		log 2 "  Round $r: minimap2 overlaps (primary-only) -> PAF"
		minimap2 -x "$ONT_PRESET" -t "$THREADS" -c --secondary=no "$CUR" "$ONT" \
			>"$RAW" 2>"$OUT/logs/r${r}.minimap2.err"

		log 2 "  Round $r: filter PAF (ident>=$MIN_IDENT, alen>=$MIN_ALEN)"
		python3 "$PAF_FILTER" "$MIN_IDENT" "$MIN_ALEN" <"$RAW" >"$FLT" \
			2>>"$OUT/logs/r${r}.filter.err"

		if [[ ! -s "$FLT" ]]; then
			log 0 "Minimap2 -> PAF filter: $FLT is empty (raw=$(wc -l <"$RAW" 2>/dev/null || echo 0))"
			# optional: fallback to RAW to avoid empty-overlap crash
			# cp -f "$RAW" "$FLT"
		fi

		log 2 "  Round $r: Racon consensus"
		racon -t "$THREADS" "$ONT" "$FLT" "$CUR" >"$OUTFA" 2>"$OUT/logs/r${r}.racon.err"

		# ----- Optional: rotate for next round -----
		if ((ROTATE_BETWEEN_ROUNDS)) && ((r < ROUNDS)); then
			if ! have_seqkit; then
				log 0 "seqkit not found; skipping rotation"
				CUR="$OUTFA"
			else
				if [[ "$(is_single_contig "$OUTFA")" == "YES" ]]; then
					if [[ -n "$ROTATE_POS" ]]; then
						POS="$ROTATE_POS"
					else
						L=$(contig_length "$OUTFA")
						POS=$(((L / 2) + 1)) # 1-based; rotate ~half the length
					fi
					ROTFA="$OUT/stage1/polished.r${r}.rot.fa"
					log 1 "  Round $r: rotating contig to POS=$POS with seqkit restart"
					if seqkit restart -i "$POS" "$OUTFA" >"$ROTFA" 2>>"$OUT/logs/r${r}.rotate.err"; then
						CUR="$ROTFA"
					else
						log 0 "seqkit restart failed; using unrotated $OUTFA"
						CUR="$OUTFA"
					fi
				else
					log 0 "Multiple contigs detected; skipping rotation for round $r"
					CUR="$OUTFA"
				fi
			fi
		else
			CUR="$OUTFA"
		fi
		# ----- end rotation block -----
	done
	log 2 "  Racon output: $CUR"
else
	log 1 "Stage 1: ONT not provided; skipping Racon"
fi

# no code
if false; then
	log 1 "Stage 1: Racon x${ROUNDS} with ONT preset=${ONT_PRESET}"
	for r in $(seq 1 "$ROUNDS"); do
		RAW="$OUT/stage1/r${r}.raw.paf"
		FLT="$OUT/stage1/r${r}.flt.paf"
		OUTFA="$OUT/stage1/polished.r${r}.fa"

		log 2 "  Round $r: minimap2 overlaps (primary-only) -> PAF"
		minimap2 -x "$ONT_PRESET" -t "$THREADS" -c --secondary=no "$CUR" "$ONT" \
			>"$RAW" 2>"$OUT/logs/r${r}.minimap2.err"

		log 2 "  Round $r: filter PAF (ident>=$MIN_IDENT, alen>=$MIN_ALEN)"
		log 2 "  Round $r: RAW: $RAW"

		python3 "$PAF_FILTER" "$MIN_IDENT" "$MIN_ALEN" <"$RAW" >"$FLT" 2>>"$OUT/logs/r${r}.filter.err"

		if [[ ! -s "$FLT" ]]; then
			log 0 "Minimap2 -> PAF filter: $FLT: empty"
		fi

		log 2 "  Round $r: racon consensus"
		racon -t "$THREADS" "$ONT" "$FLT" "$CUR" >"$OUTFA" 2>"$OUT/logs/r${r}.racon.err"
		CUR="$OUTFA"
	done
	log 2 "  Racon output: $CUR"
else
	log 1 "Stage 1: ONT not provided; skipping Racon"
fi

# =========================
# Stage 2: fmlrc2 with short reads (optional)
# =========================
if ((HAVE_SR == 1)); then
	log 1 "Stage 2: fmlrc2 short-read correction (threads=$THREADS)"
	log 2 "  phase msbwt"
	# Stage 2: Build MSBWT (FM-index) from short reads for fmlrc2
	# Notes:
	#   - Extract only sequence lines from FASTQ (NR%4==2)
	#   - Swap N↔T (upper & lower) before and after ropebwt2
	#   - Use -L since input is one sequence per line
	#   - No -t option (ropebwt2 r187 is single-threaded)
	{
		gz_or_cat "$SR1"
		gz_or_cat "$SR2"
	} |
		awk 'NR%4==2' |
		tr 'NTnt' 'TNtn' |
		ropebwt2 -LR 2>"$OUT/logs/ropebwt2.err" |
		tr 'NTnt' 'TNtn' |
		fmlrc2-convert "$OUT/stage2/comp_msbwt.npy" \
			1>"$OUT/logs/fmlrc2-convert.out" \
			2>"$OUT/logs/fmlrc2-convert.err"

	# ropebwt2 -LR -t "$THREADS" 1>"$OUT/logs/ropebwt2.out" 2>"$OUT/logs/ropebwt2.err" |

	# if [[ ! -s "$OUT/stage2/comp_msbwt.npy" || $(stat -c%s "$OUT/stage2/comp_msbwt.npy") -lt 1024 ]]; then
	# 	die "MSBWT build failed or incomplete: $OUT/stage2/comp_msbwt.npy"
	# fi

	log 2 "  phase fmlrc2"
	FMLRC2_OUT="$OUT/stage2/segments.fmlrc2.fa"
	fmlrc2 -t "$THREADS" "$OUT/stage2/comp_msbwt.npy" "$CUR" "$FMLRC2_OUT" \
		1>"$OUT/logs/fmlrc2.out" 2>"$OUT/logs/fmlrc2.err"
	[[ -s "$FMLRC2_OUT" ]] || die "fmlrc2 output missing: $FMLRC2_OUT"
	CUR="$FMLRC2_OUT"
	log 2 "  fmlrc2 output: $CUR"
else
	log 1 "Stage 2: short reads not provided; skipping fmlrc2"
fi

# =========================
# Stage 3: Polypolish (requires SR)
# =========================
if ((HAVE_SR == 1)); then
	log 1 "Stage 3: Polypolish alignment-based refinement"

	# Index must be built from the *same* FASTA you give polypolish
	IDX="$OUT/stage3/bt2idx/asm"
	mkdir -p "$(dirname "$IDX")"
	bowtie2-build "$OUT/stage2/segments.fmlrc2.fa" "$IDX" \
		>"$OUT/logs/bt2-build.out" 2>"$OUT/logs/bt2-build.err"

	# Map R1/R2 directly; report ALL alignments; name-sort for polypolish
	bowtie2 -x "$IDX" \
		-1 "$SR1" -2 "$SR2" \
		-a -p "$THREADS" \
		2>"$OUT/logs/polypolish.bt2.err" |
		samtools sort -n -O BAM -@ "$THREADS" \
			-o "$OUT/stage3/short.name.bam" - \
			2>"$OUT/logs/polypolish.sortn.err"

	POLY_OUT="$OUT/stage3/polished.polypolish.fa"

	samtools view "$OUT/stage3/short.name.bam" |
		polypolish polish "$CUR" /dev/stdin \
			>"$POLY_OUT" 2>"$OUT/logs/polypolish.err"
	[[ -s "$POLY_OUT" ]] || die "Polypolish output missing: $POLY_OUT"
	CUR="$POLY_OUT"
	log 2 "  Polypolish output: $CUR"
else
	log 1 "Stage 3: short reads not provided; skipping Polypolish"
fi

# =========================
# Final
# =========================
cp -f "$CUR" "$FINAL_OUT"
log 1 "DONE: $FINAL_OUT"
