#!/usr/bin/env bash
################################################################################
# polap-bash-polish.sh
# Version: v0.7.2
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Purpose:
#   Hybrid/ONT polishing with optional GFA reinjection.
#   Short-read inputs are auto-streamed (no temp FASTQs):
#     • polypolish / pilon  → interleaved PE stream (Bowtie2 --interleaved)
#     • racon / nextpolish  → single-end concatenated stream (R1 then R2)
#
# Notes:
#   - Default streaming uses process substitution <(...)>.
#   - If a tool refuses /dev/fd/*, use: --stream fifo  (named pipe; zero disk growth)
#   - This script avoids eval-strings; pipelines are executed directly to stay
#     neovim/shfmt/shellcheck friendly and reduce grammar hazards.
################################################################################
set -euo pipefail
IFS=$'\n\t'

_show_usage() {
	cat <<'USAGE'
Usage:
  polap-bash-polish.sh polish <assembly.fa|assembly.gfa>
                         [--ont ONT.fq[.gz]] [--sr1 R1.fq[.gz] --sr2 R2.fq[.gz]]
                         [--mode auto|ont|hybrid]
                         [--gfa-preserve]
                         [--polisher auto|fmlrc2|polypolish|pilon|racon|nextpolish]
                         [--aligner bowtie2|bwa] [--ref-index PREFIX]
                         [--racon-rounds N] [--threads N]
                         [--minimap-preset map-ont|map-hifi]
                         [--speed normal|fast|turbo]
                         [--recruit auto|draft|bait|off] [--bait BAITS.fa]
                         [--ont-cap X] [--ont-cap-method rasusa|seqtk]
                         [--eval-merqury] [--kmer K] [--prefix NAME]
                         [--ab-test] [--qc-report]
                         [--workdir DIR]
                         [--stream ps|fifo]               # select process substitution or FIFO
                         [--bt2-cache-dir DIR] [--bt2-no-cache]
                         [--keep-bam] [--keep-paf] [--dry-run]
                         [-h|--help]

Behavior (short reads):
  • Polypolish / Pilon → interleaved stream to Bowtie2 (--interleaved).
    Polypolish needs all-per-read alignments (-a).
  • Racon / NextPolish → single-end concatenated stream to minimap2 (-x sr).

Index cache:
  • Default cache dir: <workdir>/bt2cache (created automatically)
  • Cache key: sha1(path|size|mtime) of the FASTA used.

Examples:
  # Polypolish with interleaved stream + cache into species dir
  polap-bash-polish.sh polish asm.fa --sr1 s1.fq.gz --sr2 s2.fq.gz \
    --mode hybrid --polisher polypolish --bt2-cache-dir species/.bt2

  # If a tool refuses /dev/fd/*, use FIFO:
  polap-bash-polish.sh polish asm.fa --sr1 s1.fq.gz --sr2 s2.fq.gz \
    --mode hybrid --polisher polypolish --stream fifo
USAGE
}

# ---------- subcommand ----------
if [[ "${1:-}" != "polish" || $# -lt 2 ]]; then
	_show_usage
	exit 2
fi
shift

# ---------- defaults ----------
asm="${1:?}"
shift
ont=""
sr1=""
sr2=""
mode="auto" # auto|ont|hybrid
gfa_preserve=0
polisher="auto" # auto|fmlrc2|polypolish|pilon|racon|nextpolish
racon_rounds=3
threads=16
minimap_preset="map-ont" # map-ont|map-hifi
speed="normal"           # normal|fast|turbo
minimap_K="2g"

recruit="off" # off|draft|bait|auto
bait_fa=""
ont_cap=0
ont_cap_method="rasusa" # rasusa|seqtk

eval_merqury=0
kmer=31
prefix="polap"
ab_test=0
qc_report=0

aligner="bowtie2"
ref_index=""

workdir="polish.work"
stream_mode="ps" # ps|fifo
keep_bam=0
keep_paf=0
dry_run=0
bt2_no_cache=0
bt2_cache_dir=""

# ---------- parse args ----------
while (($#)); do
	case "$1" in
	--ont)
		ont="${2:?}"
		shift 2
		;;
	--sr1)
		sr1="${2:?}"
		shift 2
		;;
	--sr2)
		sr2="${2:-}"
		shift 2
		;;
	--mode)
		mode="${2:?}"
		shift 2
		;;
	--gfa-preserve)
		gfa_preserve=1
		shift
		;;
	--polisher)
		polisher="${2:?}"
		shift 2
		;;
	--racon-rounds)
		racon_rounds="${2:?}"
		shift 2
		;;
	--threads)
		threads="${2:?}"
		shift 2
		;;
	--minimap-preset)
		minimap_preset="${2:?}"
		shift 2
		;;
	--speed)
		speed="${2:?}"
		shift 2
		;;
	--recruit)
		recruit="${2:?}"
		shift 2
		;;
	--bait)
		bait_fa="${2:?}"
		shift 2
		;;
	--ont-cap)
		ont_cap="${2:?}"
		shift 2
		;;
	--ont-cap-method)
		ont_cap_method="${2:?}"
		shift 2
		;;
	--eval-merqury)
		eval_merqury=1
		shift
		;;
	--kmer)
		kmer="${2:?}"
		shift 2
		;;
	--prefix)
		prefix="${2:?}"
		shift 2
		;;
	--ab-test)
		ab_test=1
		shift
		;;
	--qc-report)
		qc_report=1
		shift
		;;
	--aligner)
		aligner="${2:?}"
		shift 2
		;;
	--ref-index)
		ref_index="${2:?}"
		shift 2
		;;
	--workdir)
		workdir="${2:?}"
		shift 2
		;;
	--stream)
		stream_mode="${2:?}"
		shift 2
		;; # ps|fifo
	--keep-bam)
		keep_bam=1
		shift
		;;
	--keep-paf)
		keep_paf=1
		shift
		;;
	--dry-run)
		dry_run=1
		shift
		;;
	--bt2-cache-dir)
		bt2_cache_dir="${2:?}"
		shift 2
		;;
	--bt2-no-cache)
		bt2_no_cache=1
		shift
		;;
	-h | --help)
		_show_usage
		exit 0
		;;
	*)
		echo "[ERROR] Unknown argument: $1" >&2
		_show_usage
		exit 2
		;;
	esac
done

# ---------- helpers ----------
log() { printf "[polish] %s\n" "$*" >&2; }
die() {
	printf "[ERROR] %s\n" "$*" >&2
	exit 1
}
have() { command -v "$1" >/dev/null 2>&1; }

# ---------- apply speed preset (ONT only) ----------
case "$speed" in
normal) : ;;
fast)
	[[ "$ont_cap" -eq 0 ]] && ont_cap=80
	[[ "$racon_rounds" == "3" ]] && racon_rounds=2
	minimap_K="2g"
	;;
turbo)
	[[ "$ont_cap" -eq 0 ]] && ont_cap=60
	[[ "$racon_rounds" == "3" ]] && racon_rounds=1
	minimap_K="4g"
	;;
*) echo "[WARN] Unknown --speed '$speed' (using 'normal')" ;;
esac

# ---------- setup ----------
run="$(date +%Y%m%d_%H%M%S)"
odir="out/${prefix}.polish_${run}"
mkdir -p "$odir" "$workdir"
if ((dry_run == 0)); then exec > >(tee -i "$odir/run.log") 2>&1; fi

log "assembly: $asm"
log "ONT:      ${ont:-'(none)'}"
log "SR1:      ${sr1:-'(none)'}"
log "SR2:      ${sr2:-'(none)'}"
log "mode:     $mode    polisher: $polisher   aligner: $aligner"
log "gfa-preserve: $gfa_preserve"
log "speed:    $speed   racon-rounds: $racon_rounds"
log "recruit:  $recruit   bait=${bait_fa:-'(none)'}"
log "ont-cap:  $ont_cap ($ont_cap_method)   minimap: -x $minimap_preset -K $minimap_K"
log "workdir:  $workdir  stream: $stream_mode"

[[ -f "$asm" ]] || die "assembly not found: $asm"
ext="${asm##*.}"
is_gfa=0
[[ "$ext" =~ ^([Gg][Ff][Aa])$ ]] && is_gfa=1

# ---------- resolve mode ----------
if [[ "$mode" == "auto" ]]; then
	if [[ -n "$sr1" && -n "$sr2" ]]; then
		mode="hybrid"
	elif [[ -n "$ont" ]]; then
		mode="ont"
	else
		die "auto mode found no usable reads; provide --ont or --sr1/--sr2"
	fi
fi
[[ "$mode" == "ont" && -z "$ont" ]] && die "--mode ont requires --ont"
[[ "$mode" == "hybrid" && (-z "$sr1" || -z "$sr2") ]] && die "--mode hybrid requires --sr1/--sr2"

# ---------- local refs ----------
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scripts_dir="$_here/scripts"
mkdir -p "$scripts_dir"

# ---------- fast helpers ----------
quiet_zcat() {
	set +e
	zcat -f "$1" 2>/dev/null || cat "$1"
	set -e
}
fa_len_total() { awk '/^>/{next}{L+=length($0)}END{print (L?L:1)}' "$1"; }

# ---------- GFA export (optional) ----------
target_fa="$asm"
if ((is_gfa == 1)); then
	if ((gfa_preserve == 1)); then
		log "GFA preserve ON → extracting S-sequences"
		if ((dry_run)); then
			echo "[DRYRUN] gfatools gfa2fa '$asm' > '$odir/segments.fa'" >&2
		else
			gfatools gfa2fa "$asm" >"$odir/segments.fa"
		fi
		target_fa="$odir/segments.fa"
	else
		log "Input is GFA without --gfa-preserve; polishing exported FASTA only."
		if ((dry_run)); then
			echo "[DRYRUN] gfatools gfa2fa '$asm' > '$odir/export.fa'" >&2
		else
			gfatools gfa2fa "$asm" >"$odir/export.fa"
		fi
		target_fa="$odir/export.fa"
	fi
fi

# ---------- ONT selection/capping ----------
recruit_for_capping() {
	local target="$1" reads="$2" outfq="$3" preset="$4" th="$5"
	local paf="$odir/recruit.paf" ids="$odir/recruit.ids"
	log "[RECRUIT] preset=$preset (len>3000 & id>0.75; MAPQ ignored)"
	if ((dry_run)); then
		echo "[DRYRUN] minimap2 -x $preset --secondary=no -K $minimap_K -t $th '$target' '$reads' > '$paf'" >&2
		echo "[DRYRUN] awk filter + seqkit grep → '$outfq'" >&2
		return
	fi
	minimap2 -x "$preset" --secondary=no -K "$minimap_K" -t "$th" "$target" "$reads" >"$paf"
	awk '$11+0>3000 {pid=$10/$11; if(pid>0.75) print $1}' "$paf" | sort -u >"$ids"
	if [[ ! -s "$ids" ]]; then
		log "[RECRUIT] no candidates; using full read set"
		cp -f "$reads" "$outfq"
	else
		seqkit grep -f "$ids" "$reads" | gzip -dc >"$outfq"
	fi
}

cap_ont_to() {
	local draft="$1" reads="$2" targetx="$3" outfq="$4"
	if [[ "$targetx" -le 0 ]]; then
		if ((dry_run)); then echo "[DRYRUN] cp '$reads' '$outfq'"; else cp -f "$reads" "$outfq"; fi
		return
	fi
	local asm_len
	asm_len=$(fa_len_total "$draft")
	if have rasusa; then
		local bases=$((targetx * asm_len))
		log "[CAP] rasusa ~${targetx}x (bases=$bases)"
		if ((dry_run)); then echo "[DRYRUN] rasusa -s 42 -i '$reads' -o '$outfq' -b $bases"; else rasusa -s 42 -i "$reads" -o "$outfq" -b "$bases"; fi
	else
		local total
		total=$(quiet_zcat "$reads" | awk 'NR%4==2{n+=length($0)}END{print (n?n:1)}')
		local p
		p=$(
			python3 - <<PY
asm_len=$asm_len; total=$total; targetx=$targetx
p = min(1.0, (targetx*asm_len)/total) if total>0 else 1.0
print(f"{p:.6f}")
PY
		)
		log "[CAP] seqtk sample p=$p to ~${targetx}x"
		if ((dry_run)); then echo "[DRYRUN] seqtk sample -s42 '$reads' $p > '$outfq'"; else seqtk sample -s42 "$reads" "$p" >"$outfq"; fi
	fi
}

# ---------- streaming views ----------
: "${POLAP_SHORT_VIEWS:=${_here}/polap-bash-short-views.sh}"
declare -a READS_OPT=() READS_SRC=() _POLAP_BG_PIDS=()
POLAP_FIFO_TO_CLEAN=""

cleanup_stream() {
	if ((${#_POLAP_BG_PIDS[@]})); then kill "${_POLAP_BG_PIDS[@]}" 2>/dev/null || true; fi
	if [[ -n "$POLAP_FIFO_TO_CLEAN" && -p "$POLAP_FIFO_TO_CLEAN" ]]; then rm -f "$POLAP_FIFO_TO_CLEAN"; fi
}
trap cleanup_stream EXIT

stream_init() {
	local mode="$1" label="$2" # se|ilv
	READS_OPT=()
	READS_SRC=()
	_POLAP_BG_PIDS=()
	POLAP_FIFO_TO_CLEAN=""
	[[ -x "$POLAP_SHORT_VIEWS" ]] || die "Missing helper: $POLAP_SHORT_VIEWS"
	local -a view=("$POLAP_SHORT_VIEWS" -1 "$sr1")
	[[ -n "$sr2" ]] && view+=(-2 "$sr2")
	view+=(--mode "$mode" --threads "$threads")

	case "$stream_mode" in
	fifo)
		local fifo="$workdir/s.${label}.fq"
		[[ "$mode" == "se" ]] && fifo+=".gz"
		if ((dry_run)); then
			echo "[DRYRUN] mkfifo '$fifo'" >&2
			echo "[DRYRUN] ${view[*]} --fifo '$fifo' &" >&2
		else
			rm -f "$fifo"
			mkfifo -m 600 "$fifo"
			POLAP_FIFO_TO_CLEAN="$fifo"
			"${view[@]}" --fifo "$fifo" &
			_POLAP_BG_PIDS+=("$!")
		fi
		[[ "$mode" == "ilv" ]] && READS_OPT=(--interleaved)
		READS_SRC=("$fifo")
		;;
	ps | *)
		[[ "$mode" == "ilv" ]] && READS_OPT=(--interleaved)
		if ((dry_run)); then
			echo "[DRYRUN] <( ${view[*]} --stdout )" >&2
			READS_SRC=("/dev/fd/63")
		else
			# shellcheck disable=SC2207
			READS_SRC=(<("${view[@]}" --stdout))
		fi
		;;
	esac
}

# ---------- Bowtie2 index with cache ----------
bt2_index_prefix() {
	if [[ -n "$ref_index" ]]; then
		printf "%s" "$ref_index"
		return
	fi
	if ((bt2_no_cache)); then
		local pref="$workdir/bt2idx/asm"
		if [[ ! -e "${pref}.1.bt2" && ! -e "${pref}.1.bt2l" ]]; then
			mkdir -p "$(dirname "$pref")"
			if ((dry_run)); then echo "[DRYRUN] bowtie2-build '$target_fa' '$pref'"; else bowtie2-build "$target_fa" "$pref"; fi
		fi
		printf "%s" "$pref"
		return
	fi
	[[ -n "$bt2_cache_dir" ]] || bt2_cache_dir="$workdir/bt2cache"
	mkdir -p "$bt2_cache_dir"
	local ap sig pref
	ap="$(readlink -f "$target_fa")"
	if ((dry_run)); then
		sig="<sha1>"
	else
		sig="$(stat -c '%n|%s|%Y' "$ap" | sha1sum | awk '{print $1}')"
	fi
	pref="${bt2_cache_dir}/${sig}"
	if ((dry_run)); then
		echo "[DRYRUN] (bt2 cache prefix) $pref" >&2
	else
		if [[ ! -e "${pref}.1.bt2" && ! -e "${pref}.1.bt2l" ]]; then
			log "[BT2] cache miss → build index: $pref"
			bowtie2-build "$target_fa" "$pref"
		else
			log "[BT2] cache hit: $pref"
		fi
	fi
	printf "%s" "$pref"
}

# ---------- polishing ----------
polished="$target_fa"
final_out=""

if [[ "$mode" == "ont" ]]; then
	log "ONT-only polishing with Racon (rounds=$racon_rounds)"
	cur="$target_fa"
	use_reads="$ont"
	label="ALL"

	if ((ab_test)); then
		log "AB test: A(ALL) vs B(recruit+cap)."
		if ((dry_run)); then
			echo "[DRYRUN] bash '$0' polish '$asm' --ont '$ont' --mode ont --threads $threads --speed $speed --recruit off --ont-cap 0 --prefix ${prefix}_A" >&2
			echo "[DRYRUN] bash '$0' polish '$asm' --ont '$ont' --mode ont --threads $threads --speed $speed --recruit draft --ont-cap 200 --prefix ${prefix}_B" >&2
		else
			bash "$0" polish "$asm" --ont "$ont" --mode ont --threads "$threads" --speed "$speed" --recruit off --ont-cap 0 --prefix "${prefix}_A" >/dev/null
			bash "$0" polish "$asm" --ont "$ont" --mode ont --threads "$threads" --speed "$speed" --recruit draft --ont-cap 200 --prefix "${prefix}_B" >/dev/null
		fi
		qc_report=1
		if ((dry_run == 0)); then
			A_DIR=$(ls -1dt out/"${prefix}_A".polish_* 2>/dev/null | head -n1 || true)
			[[ -n "${A_DIR:-}" ]] && polished="$(ls -1 "$A_DIR"/*polished.* 2>/dev/null | head -n1 || true)"
			[[ -z "${polished:-}" ]] && polished="$target_fa"
		fi
	else
		if [[ "$recruit" != "off" ]]; then
			recruited="$odir/ont.recruit.fq"
			case "$recruit" in
			draft) recruit_for_capping "$cur" "$ont" "$recruited" "$minimap_preset" "$threads" ;;
			bait)
				if [[ -n "$bait_fa" && -s "$bait_fa" ]]; then
					recruit_for_capping "$bait_fa" "$ont" "$recruited" "$minimap_preset" "$threads"
				else
					log "[RECRUIT] --bait missing; using all reads"
					if ((dry_run)); then echo "[DRYRUN] cp '$ont' '$recruited'"; else cp -f "$ont" "$recruited"; fi
				fi
				;;
			auto) recruit_for_capping "$cur" "$ont" "$recruited" "$minimap_preset" "$threads" ;;
			off | *) if ((dry_run)); then echo "[DRYRUN] cp '$ont' '$recruited'"; else cp -f "$ont" "$recruited"; fi ;;
			esac
			capped="$odir/ont.capped.fq"
			cap_ont_to "$cur" "$recruited" "$ont_cap" "$capped"
			use_reads="$capped"
			label="RECRUIT+CAP"
		fi

		# Racon rounds
		racon_w=1000
		[[ "$speed" == "fast" ]] && racon_w=1200
		[[ "$speed" == "turbo" ]] && racon_w=1500
		for i in $(seq 1 "$racon_rounds"); do
			log "[minimap2][$label] round $i"
			if ((dry_run)); then
				echo "[DRYRUN] minimap2 -x $minimap_preset --secondary=no -K $minimap_K -t $threads '$cur' '$use_reads' > '$odir/r${i}.paf'" >&2
				echo "[DRYRUN] racon -t $threads -w $racon_w '$use_reads' '$odir/r${i}.paf' '$cur' > '$odir/segments.racon${i}.fa'" >&2
			else
				minimap2 -x "$minimap_preset" --secondary=no -K "$minimap_K" -t "$threads" "$cur" "$use_reads" >"$odir/r${i}.paf"
				racon -t "$threads" -w "$racon_w" "$use_reads" "$odir/r${i}.paf" "$cur" >"$odir/segments.racon${i}.fa"
			fi
			cur="$odir/segments.racon${i}.fa"
		done
		polished="$cur"
		if ((keep_paf == 0)) && ((dry_run == 0)); then rm -f "$odir"/r*.paf || true; fi
	fi

elif [[ "$mode" == "hybrid" ]]; then
	[[ "$polisher" == "auto" ]] && polisher="fmlrc2"
	log "Hybrid polishing using $polisher"

	case "$polisher" in
	fmlrc2)
		log "[FMLRC2] build BWT & polish"
		if ((dry_run)); then
			echo "[DRYRUN] (zcat R1/R2 | ropebwt2 -LR | fmlrc2-convert comp.npy)" >&2
			echo "[DRYRUN] fmlrc2 -t $threads comp.npy '$target_fa' '$odir/segments.fmlrc2.fa'" >&2
		else
			(quiet_zcat "$sr1" "$sr2" || cat "$sr1" "$sr2") | awk 'NR%4==2' | ropebwt2 -LR | fmlrc2-convert "$odir/comp_msbwt.npy"
			fmlrc2 -t "$threads" "$odir/comp_msbwt.npy" "$target_fa" "$odir/segments.fmlrc2.fa"
		fi
		polished="$odir/segments.fmlrc2.fa"
		;;

	polypolish)
		[[ "$aligner" == "bowtie2" ]] || die "Polypolish path requires bowtie2"
		idx="$(bt2_index_prefix)"
		stream_init ilv ilv
		log "[BT2] interleaved -a → name-sorted BAM"
		if ((dry_run)); then
			echo "[DRYRUN] bowtie2 -x '$idx' ${READS_OPT[*]} ${READS_SRC[*]} -a -p $threads | samtools sort -n -o '$workdir/short.name.bam'" >&2
			echo "[DRYRUN] polypolish polish '$target_fa' '$workdir/short.name.bam' > '$odir/segments.polypolish.fa'" >&2
		else
			bowtie2 -x "$idx" "${READS_OPT[@]}" "${READS_SRC[@]}" -a -p "$threads" |
				samtools sort -n -o "$workdir/short.name.bam"
			polypolish polish "$target_fa" "$workdir/short.name.bam" >"$odir/segments.polypolish.fa"
		fi
		polished="$odir/segments.polypolish.fa"
		if ((keep_bam == 0)) && ((dry_run == 0)); then rm -f "$workdir/short.name.bam" || true; fi
		;;

	pilon)
		[[ "$aligner" == "bowtie2" ]] || die "Pilon path wired for bowtie2"
		idx="$(bt2_index_prefix)"
		stream_init ilv ilv
		log "[BT2] paired → coord-sorted BAM"
		if ((dry_run)); then
			echo "[DRYRUN] bowtie2 -x '$idx' ${READS_OPT[*]} ${READS_SRC[*]} -p $threads | samtools sort -o '$workdir/short.coord.bam'" >&2
			echo "[DRYRUN] samtools index '$workdir/short.coord.bam'" >&2
			echo "[DRYRUN] pilon --genome '$target_fa' --frags '$workdir/short.coord.bam' --threads $threads --output '$odir/pilon'" >&2
		else
			bowtie2 -x "$idx" "${READS_OPT[@]}" "${READS_SRC[@]}" -p "$threads" |
				samtools sort -o "$workdir/short.coord.bam"
			samtools index "$workdir/short.coord.bam"
			pilon --genome "$target_fa" --frags "$workdir/short.coord.bam" --threads "$threads" --output "$odir/pilon"
		fi
		polished="$odir/pilon.fasta"
		if ((keep_bam == 0)) && ((dry_run == 0)); then rm -f "$workdir/short.coord.bam" "$workdir/short.coord.bam.bai" || true; fi
		;;

	racon)
		stream_init se se
		log "[racon] SE stream via minimap2 -x sr → PAF"
		if ((dry_run)); then
			echo "[DRYRUN] minimap2 -x sr -t $threads '$target_fa' ${READS_SRC[*]} > '$workdir/short.paf'" >&2
			echo "[DRYRUN] racon ${READS_SRC[*]} '$workdir/short.paf' '$target_fa' > '$odir/segments.racon.fa'" >&2
		else
			minimap2 -x sr -t "$threads" "$target_fa" "${READS_SRC[@]}" >"$workdir/short.paf"
			racon "${READS_SRC[@]}" "$workdir/short.paf" "$target_fa" >"$odir/segments.racon.fa"
		fi
		polished="$odir/segments.racon.fa"
		if ((keep_paf == 0)) && ((dry_run == 0)); then rm -f "$workdir/short.paf" || true; fi
		;;

	nextpolish)
		stream_init se se
		log "[NextPolish] prepare SE BAM"
		if ((dry_run)); then
			echo "[DRYRUN] minimap2 -x sr -t $threads '$target_fa' ${READS_SRC[*]} | samtools sort -o '$workdir/sr.bam'" >&2
			echo "[DRYRUN] samtools index '$workdir/sr.bam'" >&2
		else
			minimap2 -x sr -t "$threads" "$target_fa" "${READS_SRC[@]}" | samtools sort -o "$workdir/sr.bam"
			samtools index "$workdir/sr.bam"
		fi
		if ((dry_run == 0)); then cp -f "$target_fa" "$odir/segments.nextpolish.fa"; fi
		polished="$odir/segments.nextpolish.fa"
		if ((keep_bam == 0)) && ((dry_run == 0)); then rm -f "$workdir/sr.bam" "$workdir/sr.bam.bai" || true; fi
		;;

	*) die "unsupported --polisher '$polisher' (use auto|fmlrc2|polypolish|pilon|racon|nextpolish)" ;;
	esac
fi

# ---------- reinject into GFA or finalize FASTA ----------
if ((is_gfa == 1)) && ((gfa_preserve == 1)); then
	log "Reinserting polished S-sequences into GFA"
	if ((dry_run)); then
		echo "[DRYRUN] python3 '$scripts_dir/gfa_inject_polished.py' '$asm' '$polished' > '$odir/${prefix}.polished.gfa'" >&2
	else
		python3 "$scripts_dir/gfa_inject_polished.py" "$asm" "$polished" >"$odir/${prefix}.polished.gfa"
	fi
	final_out="$odir/${prefix}.polished.gfa"
else
	if ((dry_run)); then
		echo "[DRYRUN] cp '$polished' '$odir/${prefix}.polished.fa'" >&2
	else
		cp -f "$polished" "$odir/${prefix}.polished.fa"
	fi
	final_out="$odir/${prefix}.polished.fa"
fi

# ---------- QC ----------
if ((qc_report == 1)); then
	mkdir -p "$odir/qc"
	_ont_qc() {
		local ref="$1" tag="$2"
		local paf="$odir/qc/${tag}.paf" tsv="$odir/qc/${tag}.pafstats.tsv"
		if ((dry_run)); then
			echo "[DRYRUN] minimap2 -x $minimap_preset --secondary=no -K $minimap_K -t $threads '$ref' '$ont' > '$paf'" >&2
			echo "[DRYRUN] awk -f '$scripts_dir/paf_stats.awk' '$paf' > '$tsv'" >&2
		else
			minimap2 -x "$minimap_preset" --secondary=no -K "$minimap_K" -t "$threads" "$ref" "$ont" >"$paf"
			awk -f "$scripts_dir/paf_stats.awk" "$paf" >"$tsv" || true
			if ((keep_paf == 0)); then rm -f "$paf" || true; fi
		fi
	}
	if ((ab_test == 1)); then
		if ((dry_run)); then
			echo "[DRYRUN] QC for A/B outputs" >&2
		else
			A_DIR=$(ls -1dt out/"${prefix}_A".polish_* 2>/dev/null | head -n1 || true)
			B_DIR=$(ls -1dt out/"${prefix}_B".polish_* 2>/dev/null | head -n1 || true)
			A_OUT=$(ls -1 "$A_DIR"/*polished.* 2>/dev/null | head -n1 || true)
			B_OUT=$(ls -1 "$B_DIR"/*polished.* 2>/dev/null | head -n1 || true)
			[[ -n "${A_OUT:-}" ]] && _ont_qc "$A_OUT" "A_all_reads"
			[[ -n "${B_OUT:-}" ]] && _ont_qc "$B_OUT" "B_recruit_cap"
		fi
	else
		_ont_qc "$final_out" "single_run"
	fi
	if ((eval_merqury == 1)) && [[ -n "${sr1:-}" && -n "${sr2:-}" ]]; then
		if ((dry_run)); then
			echo "[DRYRUN] meryl k=$kmer count '$sr1' '$sr2' output '$odir/qc/reads.meryl'" >&2
			echo "[DRYRUN] merqury.sh '$odir/qc/reads.meryl' '$final_out' '$odir/qc/SINGLE'" >&2
		else
			meryl k="$kmer" count "$sr1" "$sr2" output "$odir/qc/reads.meryl"
			merqury.sh "$odir/qc/reads.meryl" "$final_out" "$odir/qc/SINGLE" || true
		fi
	fi
	if ((dry_run == 0)); then
		echo "[QC] Files:"
		ls -1 "$odir/qc" | sed "s|^|[QC]   $odir/qc/|"
	fi
fi

echo "[DONE] Output: $final_out"
