#!/usr/bin/env bash
# polap-bash-oatk-ont-s2.sh
# Stage-2 driver:
#  1) Assemble backbone (k-ladder) using polap-bash-oatk-core-ont.sh (HPC optional)
#  2) Bait ONT subset from backbone (bbduk|meryl)
#  3) Recover:
#      - LIFT (RLE) + 1× polish   OR
#      - POLISH (racon, medaka)
#  (Annotation/Pathfinder are optional and OFF by default here)
set -euo pipefail

: "${POLAP_LOG_LEVEL:=1}"
log() {
	local l=$1
	shift
	[[ $POLAP_LOG_LEVEL -ge $l ]] && echo "[INFO]" "$@" >&2
}
die() {
	echo "[ERR]" "$@" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

# ───────── CLI defaults ─────────
READS_PRE=""
OUT=""
OATKDB=""
CLADE="magnoliopsida"
THREADS=32
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORE="${CORE:-"$SCRIPT_DIR/polap-bash-oatk-core-ont.sh"}"

# Stage-1 (backbone) knobs
HPC_ENABLE=1
HPC_OUT=""
K_LIST_A="251,151,121,91"
# K_LIST_A="91,71,51"
SMER=31
AARC=0.25
WEAKX=0.20
UNZIP=6
NO_READ_EC=0
MAX_BUB="100000"
MAX_TIP=""
COV=""
MT_SIZE_EST=0

# Bait
BAIT_METHOD="bbduk" # bbduk | meryl
BBDUK_BIN="bbduk.sh"
MERYL_K=21

# Recovery (Stage-2)
USE_LIFT=0
LIFT_BIN="${LIFT_BIN:-oatk-rle-lift.py}"
LIFT_MIN_MAPQ=10
LIFT_MIN_ALN=100
POLISH_AFTER_LIFT=1 # 1× racon after lift
RACON_ROUNDS=2
NO_MEDAKA=0
MEDAKA_MODEL="r104_e81_sup_g615"
MEDAKA_BIN=""

# Mapping speed
MAP_CHUNKS=4
MM2_EXTRA=""
MM2_BATCH="2g"

# Optional annotate/pathfinder (OFF by default)
DO_ANNOT=0
HMMBIN="hmmannot"
NHMMSCAN="nhmmscan"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --reads-pre OUT/pre/reads.pre.fq --out OUTDIR --oatkdb OATKDB [options]

Required:
  --reads-pre FILE        preprocessed reads from Stage-1
  --out DIR               output folder
  --oatkdb DIR            OatkDB root (must contain v20230921/)

General:
  --threads INT           (default ${THREADS})
  --clade NAME            fam prefix (default ${CLADE})
  -v|--verbose | -q|--quiet

Backbone (Stage-1 via core):
  --hpc / --no-hpc        use seqtk hpc for assembly input (default: on)
  --hpc-out PATH          precomputed HPC file
  --k-list "251,151,121,91"
  --smer INT              (default ${SMER})
  --a FLOAT               (default ${AARC})
  --weak-cross FLOAT      (default ${WEAKX})
  --unzip-round INT       (default ${UNZIP})
  --no-read-ec            pass --no-read-ec to core
  --mt-size-est INT       optional, tighten judge

Bait:
  --bait bbduk|meryl      (default bbduk)
  --meryl-k INT           (default ${MERYL_K})

Recovery (Stage-2):
  --lift                  do RLE lift + 1× polish (racon); else full polish path
  --lift-bin PATH         lifter (default ${LIFT_BIN})
  --lift-min-mapq INT     (default ${LIFT_MIN_MAPQ})
  --lift-min-aln  INT     (default ${LIFT_MIN_ALN})
  --racon-rounds INT      (default ${RACON_ROUNDS})
  --no-medaka             skip medaka
  --medaka-bin PATH
  --medaka-model STR      (default ${MEDAKA_MODEL})

Mapping speed:
  --map-chunks INT        split reads for mm2 (default ${MAP_CHUNKS})
  --mm2-extra "FLAGS"     extra minimap2 flags; batch -K uses ${MM2_BATCH}

Optional annotate (OFF):
  --annot                 run hmmannot (mt fam) after recovery (requires DBs)
EOF
}

# ───────── parse CLI ─────────
while [[ $# -gt 0 ]]; do
	case "$1" in
	--reads-pre)
		READS_PRE="$2"
		shift 2
		;;
	--out)
		OUT="$2"
		shift 2
		;;
	--oatkdb)
		OATKDB="$2"
		shift 2
		;;
	--threads)
		THREADS="$2"
		shift 2
		;;
	--clade)
		CLADE="$2"
		shift 2
		;;
	--hpc)
		HPC_ENABLE=1
		shift
		;;
	--no-hpc)
		HPC_ENABLE=0
		shift
		;;
	--hpc-out)
		HPC_OUT="$2"
		shift 2
		;;
	--k-list)
		K_LIST_A="$2"
		shift 2
		;;
	--c)
		COV="$2"
		shift 2
		;;
	--max-bubble)
		MAX_BUB="$2"
		shift 2
		;;
	--max-tip)
		MAX_TIP="$2"
		shift 2
		;;
	--smer)
		SMER="$2"
		shift 2
		;;
	--a)
		AARC="$2"
		shift 2
		;;
	--weak-cross)
		WEAKX="$2"
		shift 2
		;;
	--unzip-round)
		UNZIP="$2"
		shift 2
		;;
	--no-read-ec)
		NO_READ_EC=1
		shift
		;;
	--mt-size-est)
		MT_SIZE_EST="$2"
		shift 2
		;;
	--bait)
		BAIT_METHOD="$2"
		shift 2
		;;
	--meryl-k)
		MERYL_K="$2"
		shift 2
		;;
	--lift)
		USE_LIFT=1
		shift
		;;
	--lift-bin)
		LIFT_BIN="$2"
		shift 2
		;;
	--lift-min-mapq)
		LIFT_MIN_MAPQ="$2"
		shift 2
		;;
	--lift-min-aln)
		LIFT_MIN_ALN="$2"
		shift 2
		;;
	--racon-rounds)
		RACON_ROUNDS="$2"
		shift 2
		;;
	--no-medaka)
		NO_MEDAKA=1
		shift
		;;
	--medaka-bin)
		MEDAKA_BIN="$2"
		shift 2
		;;
	--medaka-model)
		MEDAKA_MODEL="$2"
		shift 2
		;;
	--map-chunks)
		MAP_CHUNKS="$2"
		shift 2
		;;
	--mm2-extra)
		MM2_EXTRA="$2"
		shift 2
		;;
	--annot)
		DO_ANNOT=1
		shift
		;;
	-v | --verbose)
		POLAP_LOG_LEVEL=2
		shift
		;;
	-q | --quiet)
		POLAP_LOG_LEVEL=0
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] unknown arg: $1" >&2
		usage
		exit 1
		;;
	esac
done
[[ -z "$READS_PRE" || -z "$OUT" || -z "$OATKDB" ]] && {
	usage
	exit 1
}

# ───────── prep & deps ─────────
READS_PRE="$(readlink -f "$READS_PRE")"
OUT="$(readlink -f "$OUT")"
OATKDB="$(readlink -f "$OATKDB")"
mkdir -p "$OUT/stage2"
cd "$OUT/stage2"

need seqkit
[[ -s "$CORE" ]] || die "core script not found: $CORE"
[[ $HPC_ENABLE -eq 1 ]] && need seqtk
[[ "$BAIT_METHOD" == "bbduk" ]] && need "$BBDUK_BIN" || { [[ "$BAIT_METHOD" == "meryl" ]] && {
	need meryl
	need meryl-lookup
	need seqtk
}; }

# ───────── helpers ─────────
mm2_map() {
	# mm2_map ref reads out.paf [mode paf|sam] [extra...]
	local ref="$1" reads="$2" out="$3" mode="${4:-paf}"
	shift 4 || true
	local extra=("$@")
	local xpreset=("-x" "map-ont" "--secondary=no" "-N" "1" "-K" "$MM2_BATCH")
	if ((MAP_CHUNKS > 1)); then
		mkdir -p map_chunks
		seqkit split2 -p "$MAP_CHUNKS" -O map_chunks "$reads" >/dev/null
		rm -f "$out"
		local nth=$((THREADS / (MAP_CHUNKS > 0 ? MAP_CHUNKS : 1)))
		((nth < 1)) && nth=1
		for f in map_chunks/*; do
			if [[ "$mode" == "sam" ]]; then
				minimap2 -a -t "$nth" "${xpreset[@]}" "${extra[@]}" "$ref" "$f" >>"$out"
			else minimap2 -t "$nth" "${xpreset[@]}" "${extra[@]}" "$ref" "$f" >>"$out"; fi
		done
	else
		if [[ "$mode" == "sam" ]]; then
			minimap2 -a -t "$THREADS" "${xpreset[@]}" "${extra[@]}" "$ref" "$reads" >"$out"
		else minimap2 -t "$THREADS" "${xpreset[@]}" "${extra[@]}" "$ref" "$reads" >"$out"; fi
	fi
}

# ───────── Stage-1: backbone (core ladder) ─────────
CUR="$READS_PRE"
CUR_ASM="$CUR"
if [[ $HPC_ENABLE -eq 1 ]]; then
	HPC_FILE="${HPC_OUT:-$OUT/stage2/reads.hpc.fa}"
	log 1 "[hpc] seqtk hpc -> $HPC_FILE"
	if [[ "$CUR" =~ \.gz$ ]]; then
		seqtk seq -A <(gzip -dc "$CUR") | seqtk hpc - >"$HPC_FILE"
	else
		seqtk seq -A "$CUR" | seqtk hpc - >"$HPC_FILE"
	fi
	CUR_ASM="$HPC_FILE"
else
	log 1 "[hpc] disabled for assembly input"
fi

log 1 "[stage2] backbone via core (k=${K_LIST_A}; s=${SMER}; a=${AARC})"
rm -rf "$OUT/ont_stage1"
bash "$CORE" \
	--reads "$CUR_ASM" \
	--out "$OUT/ont_stage1" \
	--threads "$THREADS" \
	--c "$COV" \
	--max-bubble "$MAX_BUB" \
	--max-tip "$MAX_TIP" \
	--k-list "$K_LIST_A" \
	--smer "$SMER" \
	--a "$AARC" \
	--weak-cross "$WEAKX" \
	--unzip-round "$UNZIP" \
	$([[ $NO_READ_EC -eq 1 ]] && echo --no-read-ec) \
	$([[ $MT_SIZE_EST -gt 0 ]] && echo --mt-size-est "$MT_SIZE_EST")

BACKBONE="$OUT/ont_stage1/k1/unitigs.fa"
[[ -s "$BACKBONE" ]] || die "backbone not found: $BACKBONE"

# ───────── Bait: backbone -> ONT subset ─────────
mkdir -p bait
SUBSET="$OUT/stage2/bait/ont.mt.fq"
mkdir -p "$(dirname "$SUBSET")"
if [[ "$BAIT_METHOD" == "bbduk" ]]; then
	log 1 "[bait] bbduk: ref=$(basename "$BACKBONE") -> $SUBSET"
	"$BBDUK_BIN" in="$READS_PRE" outm="$SUBSET" outu="$OUT/stage2/bait/nonmt.fq" ref="$BACKBONE" k=31 hdist=1 threads="$THREADS"
else
	log 1 "[bait] meryl k=$MERYL_K"
	meryl k="$MERYL_K" count "$BACKBONE" output "$OUT/stage2/bait/seeds.k$MERYL_K.meryl"
	# Convert FASTQ to FASTA for meryl-lookup; then grep IDs from original FASTQ
	if [[ "$READS_PRE" =~ \.gz$ ]]; then
		seqtk seq -A <(gzip -dc "$READS_PRE") >"$OUT/stage2/bait/reads.pre.fa"
	else
		seqtk seq -A "$READS_PRE" >"$OUT/stage2/bait/reads.pre.fa"
	fi
	meryl-lookup -existence -sequence "$OUT/stage2/bait/reads.pre.fa" \
		"$OUT/stage2/bait/seeds.k$MERYL_K.meryl" |
		awk '/^>/{print substr($0,2)}' >"$OUT/stage2/bait/mt.ids"
	seqkit grep -f "$OUT/stage2/bait/mt.ids" "$READS_PRE" -o "$SUBSET"
fi
log 1 "[bait] subset bytes: $( (wc -c <"$SUBSET") 2>/dev/null || echo 0)"

# ───────── Stage-2 Recovery: LIFT or POLISH ─────────
mkdir -p recov
if [[ $USE_LIFT -eq 1 ]]; then
	need minimap2
	command -v "$LIFT_BIN" >/dev/null 2>&1 || die "missing lifter: $LIFT_BIN"
	log 1 "[recover] LIFT (RLE) on backbone + 1× polish"
	# map RAW -> HPC contigs (cs tag helpful)
	mm2_map "$BACKBONE" "$SUBSET" recov/raw_vs_hpc.paf paf "--cs=long" $MM2_EXTRA
	"$LIFT_BIN" --hpc-fa "$BACKBONE" --paf recov/raw_vs_hpc.paf \
		--out-fa recov/unitigs.lifted.fa \
		--min-mapq "$LIFT_MIN_MAPQ" --min-aln "$LIFT_MIN_ALN"
	cp -f "$BACKBONE" recov/unitigs.hpc.fa
	cp -f recov/unitigs.lifted.fa recov/unitigs.fa

	# 1× racon
	need racon
	log 1 "[recover] racon x1"
	mm2_map recov/unitigs.fa "$SUBSET" recov/map.r1.paf paf $MM2_EXTRA
	racon -t "$THREADS" "$SUBSET" recov/map.r1.paf recov/unitigs.fa >recov/unitigs.polished.fa
	ln -sf unitigs.polished.fa recov/unitigs.final.fa 2>/dev/null || cp -f recov/unitigs.polished.fa recov/unitigs.final.fa

else
	# Full polish path
	need minimap2
	DRAFT="$BACKBONE"
	if [[ "$RACON_ROUNDS" -gt 0 ]]; then
		need racon
		for ((i = 1; i <= RACON_ROUNDS; i++)); do
			log 1 "[recover] racon round $i"
			mm2_map "$DRAFT" "$SUBSET" "recov/map.r${i}.paf" paf $MM2_EXTRA
			racon -t "$THREADS" "$SUBSET" "recov/map.r${i}.paf" "$DRAFT" >"recov/racon.r${i}.fa"
			DRAFT="recov/racon.r${i}.fa"
		done
	fi
	if [[ "$NO_MEDAKA" -eq 0 ]]; then
		local MED="$MEDAKA_BIN"
		if [[ -z "$MED" ]]; then
			if command -v medaka_consensus >/dev/null 2>&1; then
				MED="medaka_consensus"
			elif command -v medaka >/dev/null 2>&1; then MED="medaka"; fi
		fi
		if [[ -n "$MED" ]]; then
			log 1 "[recover] medaka model=$MEDAKA_MODEL"
			mm2_map "$DRAFT" "$SUBSET" recov/map.medaka.sam sam $MM2_EXTRA
			samtools sort -@ "$THREADS" -o recov/map.medaka.bam recov/map.medaka.sam
			samtools index recov/map.medaka.bam
			"$MED" -i "$SUBSET" -d "$DRAFT" -o recov/medaka -t "$THREADS" -m "$MEDAKA_MODEL" >/dev/null 2>&1 || true
			if [[ -s recov/medaka/consensus.fasta ]]; then DRAFT="recov/medaka/consensus.fasta"; fi
		else
			log 1 "[recover] medaka not found; skipping"
		fi
	fi
	ln -sf "$(basename "$DRAFT")" recov/unitigs.final.fa 2>/dev/null || cp -f "$DRAFT" recov/unitigs.final.fa
fi

# (optional) annotation
if [[ $DO_ANNOT -eq 1 ]]; then
	need "$HMMBIN"
	need "$NHMMSCAN"
	FAM_M="$OATKDB/v20230921/${CLADE}_mito.fam"
	[[ -s "$FAM_M" ]] || die "mt fam not found: $FAM_M"
	log 1 "[annot] $HMMBIN vs mito fam"
	"$HMMBIN" -t "$THREADS" --nhmmscan "$NHMMSCAN" -o recov/unitigs.hmm.txt "$FAM_M" recov/unitigs.final.fa || true
fi

log 1 "[Stage-2 done]"
log 1 "  backbone : $BACKBONE"
log 1 "  subset   : $SUBSET"
log 1 "  final    : $OUT/stage2/recov/unitigs.final.fa"
