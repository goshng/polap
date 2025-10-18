#!/usr/bin/env bash
# polap-bash-ptdna-recruit-subsample-flye.sh
# Version: v0.3.0
# License: GPL-3.0+
#
# Pipeline (cpDNA):
#   minimap2 (recruit) → mosdepth (pre) → subsample (rasusa|filtlong|longest) →
#   mosdepth (post) → optional Flye → aggregate metrics CSV
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'EOF'
Usage:
  polap-bash-ptdna-recruit-subsample-flye.sh \
      --seeds seeds.fa --reads ont.fq[.gz] --outdir out \
      [--cov 150] [--gsize N] [--method auto|rasusa|filtlong|longest]
      [--minlen 2000] [--threads 16] [--window 2000]
      [--asm-coverage C] [--no-assemble] [--keep-intermediate]
      [--mapper-opts "extra minimap2 options"]

Outputs:
  out/tmp/recruited.fq.gz, ds.fq.gz, *stats.tsv, pre/post mosdepth files, PNGs (if R),
  and out/ptdna.metrics.csv (single-row summary).
EOF
}

# ---------- parse args
seeds="" reads="" outdir=""
cov=150 gsize="" method="auto" minlen=2000 threads=16 window=2000 asm_cov=""
no_assemble=0 keep_intermediate=0 mapper_opts=""

while (($#)); do
	case "$1" in
	--seeds)
		seeds="${2:?}"
		shift 2
		;;
	--reads)
		reads="${2:?}"
		shift 2
		;;
	--outdir)
		outdir="${2:?}"
		shift 2
		;;
	--cov)
		cov="${2:?}"
		shift 2
		;;
	--gsize)
		gsize="${2:?}"
		shift 2
		;;
	--method)
		method="${2:?}"
		shift 2
		;;
	--minlen)
		minlen="${2:?}"
		shift 2
		;;
	--threads)
		threads="${2:?}"
		shift 2
		;;
	--window)
		window="${2:?}"
		shift 2
		;;
	--asm-coverage)
		asm_cov="${2:?}"
		shift 2
		;;
	--no-assemble)
		no_assemble=1
		shift
		;;
	--keep-intermediate)
		keep_intermediate=1
		shift
		;;
	--mapper-opts)
		mapper_opts="${2:-}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown arg: $1" >&2
		exit 2
		;;
	esac
done

[[ -z "$seeds" || -z "$reads" || -z "$outdir" ]] && {
	usage
	exit 2
}
mkdir -p "$outdir/tmp"
tmp="$outdir/tmp"

have() { command -v "$1" >/dev/null 2>&1; }
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py_extract="${script_dir}/scripts/extract_fastq_by_readnames.py"
py_stats="${script_dir}/scripts/fastq_stats.py"
py_longest="${script_dir}/scripts/sample_longest_to_target.py"
py_agg="${script_dir}/scripts/aggregate_ptdna_metrics.py"
r_cov="${script_dir}/scripts/mosdepth_regions_to_csv_and_plot.R"

for f in "$py_extract" "$py_stats" "$py_longest" "$py_agg"; do
	[[ -f "$f" ]] || {
		echo "[ERR] Missing helper: $f" >&2
		exit 3
	}
done

# ---------- estimate genome size if not provided (IR heuristic ×1.9)
if [[ -z "$gsize" ]]; then
	base_len=$(awk '/^>/ {next} {s+=length($0)} END{print s+0}' "$seeds")
	gsize=$(awk -v b="$base_len" 'BEGIN{printf "%d", (b*1.9)}')
	echo "[INFO] Estimated plastome size from seeds: ${base_len} → gsize≈${gsize}"
fi

# ---------- compute target bases
target_bases=$(awk -v c="$cov" -v g="$gsize" 'BEGIN{printf "%d", (c*g)}')
echo "[INFO] Target bases: cov(${cov}) × gsize(${gsize}) = ${target_bases}"

# ---------- 1) Recruit with minimap2 → PAF (keep all mapped: primary+secondary)
echo "[STEP] minimap2 recruit (PAF)"
minimap2 -x map-ont -t "$threads" $mapper_opts "$seeds" "$reads" >"$tmp/seeds.map.paf"

# (Optional) mapQ filter: keep all mappings (including secondaries) with MAPQ >= 0 (default all)
# If you want to raise the bar a bit, e.g. MAPQ ≥ 5:
# awk '$12 >= 5 {print $1}' "$tmp/seeds.map.paf" | sort -u > "$tmp/recruited.names.txt"
awk '{print $1}' "$tmp/seeds.map.paf" | sort -u >"$tmp/recruited.names.txt"

# Extract recruited FASTQ using seqkit (preferred) or seqtk fallback
echo "[STEP] extract recruited FASTQ (seqkit/seqtk)"
if command -v seqkit >/dev/null 2>&1; then
	# -n match by read name; -f file of names; -o output gz file
	seqkit grep -f "$tmp/recruited.names.txt" -o "$tmp/recruited.fq.gz" "$reads"
elif command -v seqtk >/dev/null 2>&1; then
	# seqtk needs uncompressed stdout; gzip afterward
	seqtk subseq "$reads" "$tmp/recruited.names.txt" | gzip >"$tmp/recruited.fq.gz"
else
	# fallback to Python extractor if neither tool is available
	python3 "$py_extract" "$reads" "$tmp/recruited.names.txt" "$tmp/recruited.fq.gz"
fi

# Quick stats (you can keep py_stats for uniform columns)
python3 "$py_stats" "$tmp/recruited.fq.gz" >"$tmp/recruited.stats.tsv"
# Or, if you prefer seqkit’s built-in stats (different columns):
# seqkit stats -a "$tmp/recruited.fq.gz" > "$tmp/recruited.seqkit.stats.txt"

# ---------- 2) QC pre-downsample (mosdepth if available)
pre_bed="-"
if have mosdepth; then
	echo "[STEP] mosdepth pre"
	minimap2 -ax map-ont -t "$threads" "$seeds" "$tmp/recruited.fq.gz" |
		samtools view -b -F 4 -@ "$threads" |
		samtools sort -@ "$threads" -o "$tmp/pre.cov.bam"
	samtools index -@ "$threads" "$tmp/pre.cov.bam"
	mosdepth -n --by "$window" "$tmp/pre.cov" "$tmp/pre.cov.bam"
	pre_bed="$tmp/pre.cov.regions.bed.gz"
	if have Rscript; then
		Rscript "$r_cov" "$pre_bed" "$tmp/pre.cov.csv" "$tmp/pre.cov.png" "Pre-downsample coverage (window=${window})"
	fi
fi

# ---------- 3) Subsample
choose_method() {
	case "$method" in
	rasusa | filtlong | longest) echo "$method" ;;
	auto) if have rasusa; then echo rasusa; elif have filtlong; then echo filtlong; else echo longest; fi ;;
	*) echo longest ;;
	esac
}
m=$(choose_method)
echo "[STEP] subsample method: $m"

case "$m" in
rasusa)
	tmpfq="$tmp/ds.tmp.fq"
	rasusa reads \
		--coverage "$cov" \
		--genome-size "$gsize" \
		-o "$tmpfq" \
		"$tmp/recruited.fq.gz"
	gzip -c "$tmpfq" >"$tmp/ds.fq.gz"
	rm -f "$tmpfq"
	;;
filtlong)
	filtlong --min_length "$minlen" --target_bases "$target_bases" "$tmp/recruited.fq.gz" | gzip >"$tmp/ds.fq.gz"
	;;
longest)
	python3 "$py_longest" --minlen "$minlen" "$tmp/recruited.fq.gz" "$target_bases" "$tmp/ds.fq.gz"
	;;
esac

python3 "$py_stats" "$tmp/ds.fq.gz" >"$tmp/ds.stats.tsv"

# ---------- 4) QC post (mosdepth if available)
post_bed="-"
if have mosdepth; then
	echo "[STEP] mosdepth post"
	minimap2 -ax map-ont -t "$threads" "$seeds" "$tmp/ds.fq.gz" |
		samtools view -b -F 4 -@ "$threads" |
		samtools sort -@ "$threads" -o "$tmp/post.cov.bam"
	samtools index -@ "$threads" "$tmp/post.cov.bam"
	mosdepth -n --by "$window" "$tmp/post.cov" "$tmp/post.cov.bam"
	post_bed="$tmp/post.cov.regions.bed.gz"
	if have Rscript; then
		Rscript "$r_cov" "$post_bed" "$tmp/post.cov.csv" "$tmp/post.cov.png" "Post-downsample coverage (window=${window})"
	fi
fi

# ---------- 5) Optional Flye
if [[ "$no_assemble" -eq 0 ]]; then
	echo "[STEP] Flye"
	local_cov="${asm_cov:-$cov}"
	mkdir -p "$outdir/flye"
	flye --nano-raw "$tmp/ds.fq.gz" --asm-coverage "$local_cov" --genome-size 200000 --threads "$threads" --out-dir "$outdir"
	# flye --nano-raw "$tmp/ds.fq.gz" --threads "$threads" --out-dir "$outdir/flye"
fi

# ---------- 6) Aggregate a single CSV
python3 "$py_agg" \
	"$tmp/recruited.stats.tsv" "$tmp/ds.stats.tsv" "$pre_bed" "$post_bed" \
	"$outdir/ptdna.metrics.csv" \
	"outdir=$(basename "$outdir")" "cov=$cov" "gsize=$gsize" "method=$m" "minlen=$minlen" "window=$window"

[[ "$keep_intermediate" -eq 0 ]] && {
	rm -f "$tmp/recruited.names.txt" "$tmp/pre.cov.bam" "$tmp/pre.cov.bam.bai" "$tmp/post.cov.bam" "$tmp/post.cov.bam.bai" 2>/dev/null || true
}

echo "[OK] Done. Metrics: $outdir/ptdna.metrics.csv"
