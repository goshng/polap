#!/usr/bin/env bash
# polap-bash-oatk-ont.sh
# ONT mtDNA assembly with Oatk components, step-by-step, with:
#  - k-ladder assemble
#  - optional RLE lift (--lift)
#  - polish (racon/medaka)
#  - annotation + pathfinder
#  - summary (PF stats + identity)
#  - NEW: --reads-map to override mapping reads for lift/polish (e.g., mt-baited subset)
#  - NEW: chunked minimap2 helper for large read sets (--map-chunks)

set -euo pipefail

: "${POLAP_LOG_LEVEL:=1}" # 0=quiet, 1=info, 2=verbose
log() {
	local lvl="$1"
	shift
	[[ "$POLAP_LOG_LEVEL" -ge "$lvl" ]] && echo "$@" >&2
}
die() {
	echo "[ERR]" "$@" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

require_seqkit_210() {
	need seqkit
	local req="2.10.0" ver
	ver="$(seqkit version 2>/dev/null | awk '{print $2}' | sed 's/^v//')"
	[[ -n "$ver" ]] || die "could not read 'seqkit version' output"
	if ! printf '%s\n%s\n' "$req" "$ver" | sort -V | head -n1 | grep -qx "$req"; then
		die "seqkit >= $req required (found v$ver). Please upgrade seqkit."
	fi
}

# ---- CLI & defaults ----
READS=""     # raw ONT reads (used for assemble by default)
READS_MAP="" # NEW: mapping reads override (used in lift/polish); defaults to READS
OUT=""
OATKDB=""
CLADE="magnoliopsida"
THREADS=32

K_LIST="401,251,151,121"
SMER=31
NO_READ_EC=0
AARC=""
MAX_BUB="200000"
MAX_TIP=""
WEAK_X="0.20"
UNZIP_ROUND="4"
COV=""

DO_POLISH=1
RACON_ROUNDS=2
MEDAKA_ENABLE=1
MEDAKA_BIN=""
MEDAKA_MODEL="r104_e81_sup_g615"

LIFT_ENABLE=0
LIFT_BIN="${LIFT_BIN:-oatk-rle-lift.py}"
LIFT_MIN_MAPQ=10
LIFT_MIN_ALN=100
POLISH_AFTER_LIFT=0

# NEW: chunked mapping options
MAP_CHUNKS=0   # 0 = no chunking; N = split reads-map into N parts
MM2_EXTRA=""   # extra flags to mm2 mapping (e.g., --cs=short)
MM2_BATCH="2g" # -K batch size

HMMBIN="hmmannot"
NHMMSCAN="nhmmscan"
STEP_REQ="all"

usage() {
	cat <<EOF
Usage: $(basename "$0") --reads READS.fq.gz --out OUTDIR --oatkdb /path/OatkDB [options]

Required:
  --reads FILE                 ONT reads for assembly (fa/fq[.gz])
  --out   DIR                  output directory
  --oatkdb DIR                 OatkDB root (must contain v20230921/)

Optional:
  --reads-map FILE             NEW: mapping reads override (used in lift/polish). If omitted, uses --reads.
  --map-chunks INT             NEW: split mapping reads into INT chunks for parallel minimap2 (default 0=no chunk)
  --mm2-extra "FLAGS"          NEW: extra minimap2 flags for mapping (e.g., --cs=short); default ""

Databases:
  --clade NAME                 fam prefix under v20230921 (default ${CLADE})
  --hmm-bin PATH               hmmannot binary (default ${HMMBIN})
  --nhmmscan PATH              nhmmscan (default ${NHMMSCAN})

Steps:
  --step, -s LIST              qc,assemble,lift,polish,annotate,pathfinder,summary or "all" (default all)

syncasm (assemble):
  --k-list "401,251,151,121"   k ladder (first success wins)
  --smer INT                   syncasm -s (<=31; default ${SMER})
  --no-read-ec                 pass --no-read-ec
  --a FLOAT                    -a min arc coverage
  --max-bubble INT             --max-bubble (default ${MAX_BUB})
  --max-tip INT                --max-tip
  --weak-cross FLOAT           --weak-cross (default ${WEAK_X})
  --unzip-round INT            --unzip-round (default ${UNZIP_ROUND})
  --c INT                      -c min k-mer coverage (use carefully for ONT)

RLE lift (Option B):
  --lift                       enable lifter (expand HPC contigs to DNA)
  --lift-bin PATH              lifter script (default ${LIFT_BIN})
  --lift-min-mapq INT          min MAPQ for lift PAF (default ${LIFT_MIN_MAPQ})
  --lift-min-aln  INT          min aln len for lift PAF (default ${LIFT_MIN_ALN})
  --polish-after-lift          1× racon after lifting (uses --reads-map if set, else --reads)

Polish (Option A):
  --polish-off                 disable polish in 'all'
  --racon-rounds INT           racon rounds (default ${RACON_ROUNDS})
  --no-medaka                  skip medaka
  --medaka-bin PATH            medaka binary
  --medaka-model STR           medaka model (default ${MEDAKA_MODEL})

General:
  --threads INT                threads (default ${THREADS})
  -v, --verbose                verbose logs
  -q, --quiet                  quiet logs
  -h, --help                   help
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
	case "$1" in
	--reads)
		READS="$2"
		shift 2
		;;
	--reads-map)
		READS_MAP="$2"
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
	--clade)
		CLADE="$2"
		shift 2
		;;
	--threads)
		THREADS="$2"
		shift 2
		;;
	--step | -s)
		STEP_REQ="$2"
		shift 2
		;;

	--k-list)
		K_LIST="$2"
		shift 2
		;;
	--smer)
		SMER="$2"
		shift 2
		;;
	--no-read-ec)
		NO_READ_EC=1
		shift
		;;
	--a)
		AARC="$2"
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
	--weak-cross)
		WEAK_X="$2"
		shift 2
		;;
	--unzip-round)
		UNZIP_ROUND="$2"
		shift 2
		;;
	--c)
		COV="$2"
		shift 2
		;;

	--lift)
		LIFT_ENABLE=1
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
	--polish-after-lift)
		POLISH_AFTER_LIFT=1
		shift
		;;

	--polish-off)
		DO_POLISH=0
		shift
		;;
	--racon-rounds)
		RACON_ROUNDS="$2"
		shift 2
		;;
	--no-medaka)
		MEDAKA_ENABLE=0
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

	--hmm-bin)
		HMMBIN="$2"
		shift 2
		;;
	--nhmmscan)
		NHMMSCAN="$2"
		shift 2
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
[[ -z "$READS" || -z "$OUT" || -z "$OATKDB" ]] && {
	usage
	exit 1
}

# ---- deps & paths ----
need syncasm
if ! command -v "$HMMBIN" >/dev/null 2>&1; then
	if command -v hmm_annotation >/dev/null 2>&1; then HMMBIN="hmm_annotation"; else die "missing hmmannot/hmm_annotation"; fi
fi
need pathfinder
need seqkit
need awk
need grep
require_seqkit_210
[[ "$LIFT_ENABLE" -eq 1 ]] && {
	need minimap2
	command -v "$LIFT_BIN" >/dev/null 2>&1 || die "missing lifter: $LIFT_BIN"
}

READS="$(readlink -f "$READS")"
[[ -n "$READS_MAP" ]] && READS_MAP="$(readlink -f "$READS_MAP")" || READS_MAP="$READS"
OUT="$(readlink -f "$OUT")"
OATKDB="$(readlink -f "$OATKDB")"

mkdir -p "$OUT"
cd "$OUT"
FAM_M="$OATKDB/v20230921/${CLADE}_mito.fam"
[[ -s "$FAM_M" ]] || die "mito fam not found: $FAM_M"

# ---- helper: chunked minimap2 mapping ----
mm2_map() {
	# mm2_map ref reads out.paf [mode paf|sam] [extra flags...]
	local ref="$1" reads="$2" out="$3" mode="${4:-paf}"
	shift 4 || true
	local extra=("$@")
	local xpreset=("-x" "map-ont" "--secondary=no" "-N" "1" "-K" "$MM2_BATCH")
	if ((MAP_CHUNKS > 1)); then
		need seqkit
		mkdir -p map_chunks
		# split by records (~equal sized parts)
		seqkit split2 -p "$MAP_CHUNKS" -O map_chunks "$reads" >/dev/null
		rm -f "$out"
		local nthreads_per=$((THREADS / (MAP_CHUNKS > 0 ? MAP_CHUNKS : 1)))
		((nthreads_per < 1)) && nthreads_per=1
		for f in map_chunks/*; do
			if [[ "$mode" == "sam" ]]; then
				minimap2 -a -t "$nthreads_per" "${xpreset[@]}" "${extra[@]}" "$ref" "$f" >>"$out"
			else
				minimap2 -t "$nthreads_per" "${xpreset[@]}" "${extra[@]}" "$ref" "$f" >>"$out"
			fi
		done
	else
		if [[ "$mode" == "sam" ]]; then
			minimap2 -a -t "$THREADS" "${xpreset[@]}" "${extra[@]}" "$ref" "$reads" >"$out"
		else
			minimap2 -t "$THREADS" "${xpreset[@]}" "${extra[@]}" "$ref" "$reads" >"$out"
		fi
	fi
}

# ---- steps ----
want_step() {
	local n="$1"
	[[ "$STEP_REQ" == "all" ]] && return 0
	IFS=',' read -r -a arr <<<"$STEP_REQ"
	for s in "${arr[@]}"; do [[ "${s// /}" == "$n" ]] && return 0; done
	return 1
}

step_qc() {
	mkdir -p qc
	log 1 "[qc] seqkit stats -Ta ${READS##*/}"
	local line
	line="$(seqkit stats -Ta "$READS" | awk 'NR==2{print}')"
	log 2 "$line"
	printf "%s\n" "$line" >qc/seqkit.stats.line.txt
	local AVGQ
	AVGQ="$(awk '{print $(NF-2)}' <<<"$line")"
	log 1 "[qc] AvgQual=$AVGQ (ONT often <25; proceed)"
}

run_syncasm_try() {
	local k="$1" outd="k1_k${k}"
	mkdir -p "$outd"
	local cmd=(syncasm -k "$k" -s "$SMER" -t "$THREADS" -o "$outd/syncasm.asm")
	[[ $NO_READ_EC -eq 1 ]] && cmd+=(--no-read-ec)
	[[ -n "$AARC" ]] && cmd+=(-a "$AARC")
	[[ -n "$MAX_BUB" ]] && cmd+=(--max-bubble "$MAX_BUB")
	[[ -n "$MAX_TIP" ]] && cmd+=(--max-tip "$MAX_TIP")
	[[ -n "$WEAK_X" ]] && cmd+=(--weak-cross "$WEAK_X")
	[[ -n "$UNZIP_ROUND" ]] && cmd+=(--unzip-round "$UNZIP_ROUND")
	[[ -n "$COV" ]] && cmd+=(-c "$COV")
	log 1 "[assemble] ${cmd[*]} ${READS##*/}"
	/usr/bin/time -f "TIME(s) %E  RSS(KB) %M" "${cmd[@]}" "$READS" 2>"$outd/run.log" || true

	local fa1="$outd/syncasm.asm.fa" fa2="$outd/syncasm.asm.fasta"
	local gfaF="$outd/syncasm.asm.utg.final.gfa" gfaI="$outd/syncasm.asm.utg.gfa" gfa="$outd/syncasm.asm.gfa"
	if [[ -s "$fa1" ]]; then
		echo "$fa1"
		return 0
	elif [[ -s "$fa2" ]]; then
		echo "$fa2"
		return 0
	elif [[ -s "$gfaF" ]] && awk 'BEGIN{s=0} /^S\t/{s++} END{exit (s==0)}' "$gfaF"; then
		awk '$1=="S"{print ">"$2"\n"$3}' "$gfaF" >"$outd/unitigs.fa"
		echo "$outd/unitigs.fa"
		return 0
	elif [[ -s "$gfaI" ]] && awk 'BEGIN{s=0} /^S\t/{s++} END{exit (s==0)}' "$gfaI"; then
		awk '$1=="S"{print ">"$2"\n"$3}' "$gfaI" >"$outd/unitigs.fa"
		echo "$outd/unitigs.fa"
		return 0
	elif [[ -s "$gfa" ]] && command -v gfatools >/dev/null 2>&1; then
		gfatools asm -u "$gfa" >"$outd/unitigs.gfa"
		if awk 'BEGIN{s=0} /^S\t/{s++} END{exit (s==0)}' "$outd/unitigs.gfa"; then
			awk '$1=="S"{print ">"$2"\n"$3}' "$outd/unitigs.gfa" >"$outd/unitigs.fa"
			echo "$outd/unitigs.fa"
			return 0
		fi
	fi
	return 1
}

step_assemble() {
	mkdir -p qc k1
	log 1 "[assemble] k ladder=${K_LIST} ; smer=${SMER} ; threads=${THREADS} ; EC=$((NO_READ_EC ? 0 : 1))"
	: >qc/assemble.out
	local K_USED="" UFA=""
	IFS=',' read -r -a KS <<<"$K_LIST"
	for k in "${KS[@]}"; do
		k="${k// /}"
		[[ -z "$k" ]] && continue
		if UFA="$(run_syncasm_try "$k")"; then
			K_USED="$k"
			echo -e "[ok]\tk=${k}\t${UFA}" | tee -a qc/assemble.out
			ln -sf "../${UFA}" k1/unitigs.fa 2>/dev/null || cp -f "${UFA}" k1/unitigs.fa
			printf "%s\n" "$K_USED" >k1/USED_K.txt
			local N TBP
			N=$(grep -c '^>' k1/unitigs.fa)
			TBP=$(awk '/^>/{if(len){t+=len}len=0;next}{len+=length}END{t+=len;print t}' k1/unitigs.fa)
			printf "k_used\t%s\nunitigs\t%d\ntotal_bp\t%d\n" "$K_USED" "$N" "$TBP" | tee qc/assemble.quick.tsv
			return 0
		else
			echo -e "[fail]\tk=${k}" | tee -a qc/assemble.out
		fi
	done
	log 1 "[assemble] no unitigs for any k in [${K_LIST}]"
	return 2
}

step_lift() {
	[[ -s k1/unitigs.fa ]] || die "k1/unitigs.fa not found; run assemble"
	need minimap2
	command -v "$LIFT_BIN" >/dev/null 2>&1 || die "missing lifter $LIFT_BIN"
	mkdir -p lift
	log 1 "[lift] RAW→HPC PAF with cs (reads-map=${READS_MAP##*/})"
	mm2_map k1/unitigs.fa "$READS_MAP" lift/raw_vs_hpc.paf paf "--cs=long" $MM2_EXTRA
	log 1 "[lift] lifter on PAF"
	"$LIFT_BIN" --hpc-fa k1/unitigs.fa --paf lift/raw_vs_hpc.paf \
		--out-fa k1/unitigs.lifted.fa \
		--min-mapq "$LIFT_MIN_MAPQ" --min-aln "$LIFT_MIN_ALN"
	cp -f k1/unitigs.fa k1/unitigs.hpc.fa
	cp -f k1/unitigs.lifted.fa k1/unitigs.fa
	log 1 "[lift] k1/unitigs.fa is now DNA-expanded (lifted)"
	return 0
}

step_polish() {
	[[ -s k1/unitigs.fa ]] || die "k1/unitigs.fa not found; run assemble/lift"
	mkdir -p polish
	need minimap2
	local draft="k1/unitigs.fa"
	local reads="$READS_MAP"
	local i
	if [[ "$RACON_ROUNDS" -gt 0 ]]; then
		need racon
		for ((i = 1; i <= RACON_ROUNDS; i++)); do
			log 1 "[polish] racon round $i (reads-map=${reads##*/})"
			mm2_map "$draft" "$reads" "polish/map.r${i}.paf" paf $MM2_EXTRA
			racon -t "$THREADS" "$reads" "polish/map.r${i}.paf" "$draft" >"polish/racon.r${i}.fa"
			draft="polish/racon.r${i}.fa"
		done
	else
		log 1 "[polish] skipping racon (rounds=0)"
	fi
	if [[ "$MEDAKA_ENABLE" -eq 1 ]]; then
		local medaka_cmd=""
		if [[ -n "$MEDAKA_BIN" ]]; then
			need "$MEDAKA_BIN"
			medaka_cmd="$MEDAKA_BIN"
		elif command -v medaka_consensus >/dev/null 2>&1; then
			medaka_cmd="medaka_consensus"
		elif command -v medaka >/dev/null 2>&1; then medaka_cmd="medaka"; fi
		if [[ -n "$medaka_cmd" ]]; then
			[[ -n "$MEDAKA_MODEL" ]] || die "set --medaka-model"
			log 1 "[polish] medaka ($medaka_cmd) model=$MEDAKA_MODEL"
			# build BAM efficiently on the subset
			mm2_map "$draft" "$reads" polish/map.medaka.sam sam $MM2_EXTRA
			samtools sort -@ "$THREADS" -o polish/map.medaka.bam polish/map.medaka.sam
			samtools index polish/map.medaka.bam
			"$medaka_cmd" -i "$reads" -d "$draft" -o polish/medaka -t "$THREADS" -m "$MEDAKA_MODEL" >/dev/null 2>&1 || true
			if [[ -s polish/medaka/consensus.fasta ]]; then draft="polish/medaka/consensus.fasta"; fi
		else
			log 1 "[polish] medaka not found; skipping"
		fi
	fi
	cp -f k1/unitigs.fa k1/unitigs.fa.bak
	cp -f "$draft" k1/unitigs.polished.fa
	ln -sf "unitigs.polished.fa" k1/unitigs.fa 2>/dev/null || cp -f k1/unitigs.polished.fa k1/unitigs.fa

	mkdir -p qc
	mm2_map k1/unitigs.fa "$reads" qc/raw_vs_polished.paf paf $MM2_EXTRA
	awk '{if($11>0){id=($10/$11)*100; print id}}' qc/raw_vs_polished.paf | sort -n >qc/ident.pct.txt || true
	MEDID=$(awk 'NF{a[NR]=$1} END{if(NR==0){print "NA"} else if (NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}' qc/ident.pct.txt)
	echo "$MEDID" >qc/ident.pct.median.txt
	log 1 "[polish] median mapping identity: ${MEDID}%"
	return 0
}

step_annotate() {
	[[ -s k1/unitigs.fa ]] || die "k1/unitigs.fa not found"
	log 1 "[annotate] $HMMBIN -t $THREADS --nhmmscan $NHMMSCAN -o k1.hmmannot.txt $FAM_M k1/unitigs.fa"
	"$HMMBIN" -t "$THREADS" --nhmmscan "$NHMMSCAN" -o k1.hmmannot.txt "$FAM_M" k1/unitigs.fa || return 3
	grep -v '^#' k1.hmmannot.txt | awk '{print $1}' | sort -u >qc/genes.hit.names
	wc -l <qc/genes.hit.names | tee k1/GENE_HIT_COUNT.txt >/dev/null
	log 1 "[annotate] distinct mt genes: $(cat k1/GENE_HIT_COUNT.txt)"
	return 0
}

step_pathfinder() {
	[[ -s k1/unitigs.fa ]] || die "k1/unitigs.fa not found"
	[[ -s k1.hmmannot.txt ]] || die "k1.hmmannot.txt not found"
	log 1 "[pathfinder] pathfinder -t $THREADS -o pathfinder k1/unitigs.fa k1.hmmannot.txt"
	pathfinder -t "$THREADS" -o pathfinder k1/unitigs.fa k1.hmmannot.txt || return 3
	return 0
}

eval_pathfinder() {
	mkdir -p qc
	local pf="pathfinder" out="qc/pathfinder.summary.tsv"
	echo -e "file\tn_seqs\ttotal_bp\tn_circular\tbp_circular\tnote" >"$out"
	fa_stats() { awk 'BEGIN{n=0;bp=0;nc=0;bc=0;hdr=""}
    /^>/{n++; if(hdr~/circular|circ|CIRC/){nc++;bc+=len} hdr=$0; len=0; next}
    {len+=length}
    END{ if(hdr~/circular|circ|CIRC/){nc++;bc+=len}; print n, bp+len, nc, bc }' "$1" 2>/dev/null || echo "0 0 0 0"; }
	local total_n=0 total_bp=0 total_nc=0 total_bc=0 noted="."
	shopt -s nullglob
	local candidates=()
	for f in "$pf"/final*.fa* "$pf"/*consensus*.fa* "$pf"/*.fa "$pf"/*.fasta; do candidates+=("$f"); done
	if [[ ${#candidates[@]} -eq 0 ]]; then
		echo -e "-\t0\t0\t0\t0\tno_fasta_found" >>"$out"
	else
		for f in "${candidates[@]}"; do
			read -r n bp nc bc < <(fa_stats "$f")
			echo -e "$(basename "$f")\t$n\t$bp\t$nc\t$bc\t-" >>"$out"
			total_n=$((total_n + n))
			total_bp=$((total_bp + bp))
			total_nc=$((total_nc + nc))
			total_bc=$((total_bc + bc))
		done
	fi
	echo -e "TOTAL\t${total_n}\t${total_bp}\t${total_nc}\t${total_bc}\t${noted}" >>"$out"
	log 1 "[summary] PF summary → $out"
}

step_summary() {
	log 1 "[summary]"
	[[ -s qc/seqkit.stats.line.txt ]] && log 1 "[EMIT] QC line         : $OUT/qc/seqkit.stats.line.txt"
	[[ -s qc/assemble.out ]] && log 1 "[EMIT] assemble ladder : $OUT/qc/assemble.out"
	[[ -s qc/assemble.quick.tsv ]] && log 1 "[EMIT] assemble quick  : $OUT/qc/assemble.quick.tsv"
	[[ -s k1/USED_K.txt ]] && log 1 "[EMIT] used k          : $(cat k1/USED_K.txt)"
	[[ -s k1/unitigs.hpc.fa ]] && log 1 "[EMIT] unitigs (HPC)   : $OUT/k1/unitigs.hpc.fa"
	[[ -s k1/unitigs.lifted.fa ]] && log 1 "[EMIT] unitigs (lifted): $OUT/k1/unitigs.lifted.fa"
	[[ -s k1/unitigs.fa.bak ]] && log 1 "[EMIT] unitigs.bak     : $OUT/k1/unitigs.fa.bak"
	[[ -s k1/unitigs.fa ]] && log 1 "[EMIT] unitigs (final) : $OUT/k1/unitigs.fa"
	[[ -s qc/ident.pct.median.txt ]] && log 1 "[EMIT] median identity : $(cat qc/ident.pct.median.txt)%"
	[[ -s k1.hmmannot.txt ]] && log 1 "[EMIT] annotate table  : $OUT/k1.hmmannot.txt (hits=$(cat k1/GENE_HIT_COUNT.txt 2>/dev/null || echo 0))"
	[[ -d pathfinder ]] && {
		eval_pathfinder
		log 1 "[EMIT] PF summary     : $OUT/qc/pathfinder.summary.tsv"
	}
}

# ---- run orchestration ----
mkdir -p qc
require_seqkit_210

# Default ‘all’: swap polish for lift if requested
if [[ "$STEP_REQ" == "all" ]]; then
	if [[ "$LIFT_ENABLE" -eq 1 ]]; then
		STEP_REQ=$([[ "$POLISH_AFTER_LIFT" -eq 1 ]] && echo "qc,assemble,lift,polish,annotate,pathfinder,summary" || echo "qc,assemble,lift,annotate,pathfinder,summary")
	else
		STEP_REQ=$([[ "$DO_POLISH" -eq 1 ]] && echo "qc,assemble,polish,annotate,pathfinder,summary" || echo "qc,assemble,annotate,pathfinder,summary")
	fi
fi

want_step() {
	local n="$1"
	[[ "$STEP_REQ" == "all" ]] && return 0
	IFS=',' read -r -a arr <<<"$STEP_REQ"
	for s in "${arr[@]}"; do [[ "${s// /}" == "$n" ]] && return 0; done
	return 1
}

rc=0
if want_step qc; then step_qc || rc=$?; fi
if want_step assemble; then step_assemble || rc=$?; fi
if want_step lift; then step_lift || rc=$?; fi
if want_step polish; then step_polish || rc=$?; fi
if want_step annotate; then step_annotate || rc=$?; fi
if want_step pathfinder; then step_pathfinder || rc=$?; fi
if want_step summary; then step_summary || rc=$?; fi
exit "$rc"
