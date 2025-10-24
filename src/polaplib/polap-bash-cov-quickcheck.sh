#!/usr/bin/env bash
# ==============================================================================
# polap-bash-cov-quickcheck.sh — v0.1.0
# Quick coverage sanity check: map reads -> BAM -> contig stats + windowed depth
# - Supports ONT/HiFi presets (minimap2)
# - Uses mosdepth if available (fast); otherwise bedtools genomecov fallback
# - Emits:
#     <out>/map.bam, map.bam.bai
#     <out>/contig_coverage.tsv            (samtools coverage per-contig)
#     <out>/window_coverage.bed.gz         (mosdepth or bedtools window coverage)
#     <out>/summary.txt                    (tiny summary)
# ==============================================================================

set -euo pipefail

VER="v0.1.0"

# ---------------------------- defaults & CLI ----------------------------------
outdir="covqc.out"
reads=""
assembly=""
threads="${POLAP_THREADS:-8}"
preset="map-ont" # map-ont | map-hifi | map-pb
window=2000      # bp
verbose=0
quiet=0

usage() {
	cat <<EOF
polap-bash-cov-quickcheck.sh ${VER}

USAGE:
  bash \$0 -r reads.fastq[.gz] -a assembly.fasta [options]

REQUIRED:
  -r, --reads    FASTQ(.gz) input
  -a, --asm      Assembly FASTA

OPTIONS:
  -o, --outdir   Output directory [${outdir}]
  -x, --preset   minimap2 preset: map-ont | map-hifi | map-pb [${preset}]
  -w, --window   Window size for depth track (bp) [${window}]
  -t, --threads  Threads [${threads}]
  -v, --verbose  Verbose (repeatable)
      --quiet    Quiet
      --version  Print version and exit
  -h, --help     This help

OUTPUTS (in --outdir):
  map.bam, map.bam.bai
  contig_coverage.tsv
  window_coverage.bed.gz     (mosdepth: *.regions.bed.gz; fallback: bedtools .bed.gz)
  summary.txt
EOF
}

# parse
while [[ $# -gt 0 ]]; do
	case "$1" in
	-r | --reads)
		reads="$2"
		shift 2
		;;
	-a | --asm | --assembly)
		assembly="$2"
		shift 2
		;;
	-o | --outdir)
		outdir="$2"
		shift 2
		;;
	-x | --preset)
		preset="$2"
		shift 2
		;;
	-w | --window)
		window="$2"
		shift 2
		;;
	-t | --threads)
		threads="$2"
		shift 2
		;;
	-v | --verbose)
		verbose=$((verbose + 1))
		shift
		;;
	--quiet)
		quiet=1
		verbose=0
		shift
		;;
	--version)
		echo "$VER"
		exit 0
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERROR] unknown arg: $1" >&2
		usage
		exit 2
		;;
	esac
done

# ------------------------------- logging --------------------------------------
log() { [[ $quiet -eq 0 ]] && echo -e "$@" >&2; }
vlog() { [[ $verbose -gt 0 && $quiet -eq 0 ]] && echo -e "$@" >&2; }

# ----------------------------- sanity checks ----------------------------------
[[ -s "$reads" ]] || {
	echo "[ERROR] missing --reads" >&2
	exit 2
}
[[ -s "$assembly" ]] || {
	echo "[ERROR] missing --asm" >&2
	exit 2
}

need_exec() {
	local miss=0 x
	for x in "$@"; do
		command -v "$x" >/dev/null 2>&1 || {
			echo "[ERROR] missing executable: $x" >&2
			miss=1
		}
	done
	return $miss
}
need_exec minimap2 samtools bgzip || exit 2
have_mosdepth=0
command -v mosdepth >/dev/null 2>&1 && have_mosdepth=1
have_bedtools=0
command -v bedtools >/dev/null 2>&1 && have_bedtools=1

mkdir -p "$outdir"
# own log file (since this runner can be used independently of polap)
logfile="${outdir}/covqc.log"
exec > >(tee -a "$logfile") 2>&1

log "# polap-bash-cov-quickcheck.sh ${VER}"
vlog "# args: reads=$reads asm=$assembly preset=$preset window=$window threads=$threads"

# ----------------------------- 1) index FASTA ---------------------------------
if [[ ! -s "${assembly}.fai" ]]; then
	log "# faidx: $assembly"
	samtools faidx "$assembly"
fi

# ----------------------------- 2) map -> BAM -----------------------------------
bam="${outdir}/map.bam"
if [[ ! -s "$bam" ]]; then
	log "# map: minimap2 -x $preset"
	minimap2 -t "$threads" -x "$preset" "$assembly" "$reads" |
		samtools view -@ "$threads" -b |
		samtools sort -@ "$threads" -o "$bam" -
	samtools index -@ "$threads" "$bam"
else
	vlog "# reuse BAM: $bam"
fi

# -------------------------- 3) per-contig coverage ----------------------------
contig_tsv="${outdir}/contig_coverage.tsv"
log "# per-contig coverage: samtools coverage"
samtools coverage -A -o "$contig_tsv" "$bam"

# Tiny overall summary (mean depth across contigs, mapped %, etc.)
summary="${outdir}/summary.txt"
{
	echo "# summary (quick):"
	# mean depth weighted by contig length
	awk 'BEGIN{FS="\t"} NR==1{next} {cov+=($7*$3); len+=$3} END{if(len>0) printf("weighted_mean_depth\t%.2f\n", cov/len); else print "weighted_mean_depth\t0"}' "$contig_tsv"
	# mapped reads
	samtools flagstat "$bam" | sed -n '1,4p'
} >"$summary"
vlog "# wrote: $contig_tsv"
vlog "# wrote: $summary"

# ------------------------ 4) windowed depth track -----------------------------
# Prefer mosdepth (fast; produces .regions.bed.gz). Otherwise bedtools fallback.
if [[ $have_mosdepth -eq 1 ]]; then
	log "# windowed coverage: mosdepth (window=${window})"
	# mosdepth outputs: <prefix>.regions.bed.gz with mean depth per region
	mosdepth --threads "$threads" --by "$window" --fast-mode "${outdir}/mos" "$bam"
	# normalize filename to a stable name for downstream:
	if [[ -s "${outdir}/mos.regions.bed.gz" ]]; then
		mv "${outdir}/mos.regions.bed.gz" "${outdir}/window_coverage.bed.gz"
	else
		# new mosdepth (>=0.3) names: <prefix>.regions.bed.gz already
		gzip -t "${outdir}/mos.regions.bed.gz" 2>/dev/null || true
	fi
	# keep the per-base and per-chunk outputs if present (optional clean)
	vlog "# mosdepth outputs at: ${outdir}"
else
	if [[ $have_bedtools -eq 0 ]]; then
		echo "[WARN] neither mosdepth nor bedtools found: skip windowed track." >&2
	else
		log "# windowed coverage: bedtools (window=${window})"
		# make windows from assembly.fai
		winbed="${outdir}/windows.w${window}.bed"
		awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "${assembly}.fai" |
			bedtools makewindows -w "$window" -b - >"$winbed"
		# mean depth per window
		bedtools coverage -mean -a "$winbed" -b "$bam" |
			bgzip -c >"${outdir}/window_coverage.bed.gz"
		rm -f "$winbed"
	fi
fi

# ------------------------------ 5) final notes --------------------------------
log "# done."
log "# outputs:"
log "#   - $bam (and .bai)"
log "#   - $contig_tsv"
[[ -s "${outdir}/window_coverage.bed.gz" ]] && log "#   - ${outdir}/window_coverage.bed.gz"
log "#   - $summary"
