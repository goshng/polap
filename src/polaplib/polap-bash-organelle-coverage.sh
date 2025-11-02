#!/usr/bin/env bash
set -euo pipefail

usage() {
	cat <<'USAGE'
polap-bash-organelle-coverage.sh \
  -a assembly.fasta [-r reads.fq[.gz] | -b aln.bam] \
  [--preset map-ont|map-hifi|lr:hq] [--mq 20] [--bq 0] \
  [--bin 500] [--count-dels] [-t 8] [-o outdir] [-s SAMPLE]

Outputs: BAM, per-base depth (zeros kept), binned summaries, metrics, outliers, and plots.

Options:
  -a, --assembly    Polished organelle assembly FASTA
  -r, --reads       Long-read FASTQ(.gz) (ONT/PacBio)
  -b, --bam         Use existing sorted BAM (skips mapping)
  --preset          minimap2 preset (default: map-ont)
  --mq              min mapping quality for depth (default: 20)
  --bq              min base quality for depth (default: 0)
  --bin             bin size for summaries (bp; default: 500)
  --count-dels      include deletions in depth (samtools depth -J)
  -t, --threads     threads (default: 8)
  -o, --outdir      output directory (default: cov_eval)
  -s, --sample      sample name prefix (default: sample)
USAGE
}

POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"

# defaults
preset="map-ont"
threads=8
mq=20
bq=0
binsize=500
count_dels=0
outdir="cov_eval"
sample="sample"
reads=""
ref=""
bam_in=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-a | --assembly)
		ref="$2"
		shift 2
		;;
	-r | --reads)
		reads="$2"
		shift 2
		;;
	-b | --bam)
		bam_in="$2"
		shift 2
		;;
	--preset)
		preset="$2"
		shift 2
		;;
	--mq)
		mq="$2"
		shift 2
		;;
	--bq)
		bq="$2"
		shift 2
		;;
	--bin)
		binsize="$2"
		shift 2
		;;
	--count-dels)
		count_dels=1
		shift 1
		;;
	-t | --threads)
		threads="$2"
		shift 2
		;;
	-o | --outdir)
		outdir="$2"
		shift 2
		;;
	-s | --sample)
		sample="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown arg: $1"
		usage
		exit 1
		;;
	esac
done

if [[ -z "$ref" ]]; then
	echo "[ERR] --assembly required"
	exit 1
fi
mkdir -p "$outdir"/{bam,depth,bins,metrics,plots,logs}

# (1) Map if no BAM provided
bam="$outdir/bam/${sample}.sorted.bam"
if [[ -n "$bam_in" ]]; then
	echo "[INFO] Using existing BAM: $bam_in"
	bam="$bam_in"
else
	if [[ -z "$reads" ]]; then
		echo "[ERR] Provide --reads or --bam"
		exit 1
	fi
	echo "[INFO] Mapping with minimap2 ($preset)..."
	minimap2 -t "$threads" -ax "$preset" "$ref" "$reads" |
		samtools sort -@ "$threads" -o "$bam" -
	samtools index -@ "$threads" "$bam"
fi

# (2) Per-base depth incl. zeros
depth_tsv="$outdir/depth/${sample}.depth.tsv"
opt="-aa -Q ${mq} -q ${bq}"
[[ $count_dels -eq 1 ]] && opt="$opt -J"
echo "[INFO] samtools depth $opt ..."
samtools depth $opt "$bam" >"$depth_tsv"

# (3) Contig lengths (FAI)
if [[ ! -s "${ref}.fai" ]]; then samtools faidx "$ref"; fi
cp "${ref}.fai" "$outdir/metrics/${sample}.contigs.fai"

# (4) Bin coverage with AWK
binned_tsv="$outdir/bins/${sample}.bins.${binsize}bp.tsv"
awk -v B="$binsize" -f "${POLAPLIB_DIR}/scripts/bin_depth_by_contig.awk" "$depth_tsv" >"$binned_tsv"

# (5) Uniformity metrics
metrics_tsv="$outdir/metrics/${sample}.uniformity.tsv"
python3 "${POLAPLIB_DIR}/scripts/compute_coverage_uniformity_metrics.py" "$depth_tsv" "$metrics_tsv"

# (6) Outlier segments (spikes/troughs) on binned means
outliers_bed="$outdir/metrics/${sample}.depth_outliers.bed"
python3 "${POLAPLIB_DIR}/scripts/detect_depth_outliers.py" "$binned_tsv" "$outliers_bed"

# (7) Plots: flatness lines, uniformity ECDF, Lorenz+Gini
Rscript "${POLAPLIB_DIR}/scripts/make_uniformity_plots.R" \
	"$depth_tsv" "$binned_tsv" "$metrics_tsv" "$outdir/plots" "$binsize"

echo "[DONE] Outputs in: $outdir"
