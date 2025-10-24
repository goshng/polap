#!/usr/bin/env bash
# polap-bash-mt-stoich-fractions.sh
# Version : v1.1.0
# Purpose : Per-repeat recombination usage (isomer stoichiometry) from ONT reads, step-by-step.
#
# Steps:
#   1) Repeat discovery (MUMmer/BLAST) -> repeats.tsv/.bed
#   2) Junction synthesis (parental & recombinant) -> junc_parent.fa/.tsv, junc_recomb.fa/.tsv
#   3) Junction-bridge counting (map ONT -> junctions) -> junc_parent_support.tsv, junc_recomb_support.tsv
#   4) Fraction quantification (and optional depth normalization) -> repeat_fractions.tsv
#   5) Summary PDF -> report/report.pdf
#
# Example:
#   bash /path/to/polaplib/polap-bash-mt-stoich-fractions.sh \
#     -f Anthoceros_agrestis/v6/0/polap-readassemble-1-miniasm/mt.1.fa \
#     -l Anthoceros_agrestis/tmp/l.fq \
#     -o Anthoceros_agrestis-0/stoich-S1/ \
#     --label S1 --threads 12 \
#     --min-repeat-len 1000 --min-repeat-pid 90 \
#     --flank 800 --min-mapq 20 --min-span 250 --prefer mummer
#
# Run a subset of steps:
#   --step 1,2          # only repeat discovery + junction synthesis
#   --step 3-5          # from mapping to PDF
#   --step all          # (default) run everything
#
# Inputs (required):
#   -f, --fasta   FASTA   Mitochondrial assembly (reference)
#   -l, --long    FASTQ   ONT reads (needed for steps 3–5)
#   -o, --out     DIR     Output directory
#
# Other options (defaults are sensible for ONT):
#   --label STR [sample]          Tag for outputs (TSV/PDF)
#   --threads INT [8]             Threads for minimap2/samtools
#   --min-repeat-len INT [1000]   Repeat length threshold (bp)
#   --min-repeat-pid INT [95]     Repeat identity threshold (%)
#   --flank INT [500]             Junction template flank per side (bp)
#   --min-mapq INT [20]           MAPQ cutoff for bridge counting
#   --min-span INT [200]          Required span (bp) on each side of boundary
#   --prefer mummer|blast [mummer] Backend for repeat discovery
#   --step SPEC [all]             Step selector: e.g., 1,2 or 3-5 or all
#
set -euo pipefail
set -o errtrace
IFS=$'\n\t'

# ---------- locate polaplib ----------
_POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
export _POLAPLIB_DIR POLAPLIB_DIR="${_POLAPLIB_DIR}"

# ---------- tiny logger ----------
_now() { date "+%Y-%m-%d %H:%M:%S"; }
logi() { printf "[%s %s:%s] %s\n" "$(_now)" "${BASH_SOURCE[1]##*/}" "${BASH_LINENO[0]}" "$*" >&1; }
loge() { printf "[%s %s:%s] ERROR: %s\n" "$(_now)" "${BASH_SOURCE[1]##*/}" "${BASH_LINENO[0]}" "$*" >&2; }

# ---------- defaults ----------
reads="" out="" asm=""
label="sample"
threads=8
min_rep_len=1000
min_rep_pid=95
flank=500
min_mapq=20
min_span=200
prefer="mummer"
step_spec="all"

# ---------- parse args ----------
print_help() {
	sed -n '1,120p' "$0" | sed -n '1,120p' >/dev/null 2>&1 || true
	grep -E '^#( |$)' "$0" | sed 's/^# \{0,1\}//'
}

# long opts with values
while [[ $# -gt 0 ]]; do
	case "$1" in
	-f | --fasta)
		asm="$2"
		shift 2
		;;
	-l | --long)
		reads="$2"
		shift 2
		;;
	-o | --out)
		out="$2"
		shift 2
		;;
	--label)
		label="$2"
		shift 2
		;;
	--threads)
		threads="$2"
		shift 2
		;;
	--min-repeat-len)
		min_rep_len="$2"
		shift 2
		;;
	--min-repeat-pid)
		min_rep_pid="$2"
		shift 2
		;;
	--flank)
		flank="$2"
		shift 2
		;;
	--min-mapq)
		min_mapq="$2"
		shift 2
		;;
	--min-span)
		min_span="$2"
		shift 2
		;;
	--prefer)
		prefer="$2"
		shift 2
		;;
	--step)
		step_spec="$2"
		shift 2
		;;
	-h | --help)
		print_help
		exit 0
		;;
	*)
		loge "Unknown option: $1"
		print_help
		exit 1
		;;
	esac
done

[[ -z "$asm" || -z "$out" ]] && {
	loge "Missing required -f/--fasta and -o/--out"
	print_help
	exit 1
}
# Reads are required for steps >=3
if [[ "$step_spec" != "1" && "$step_spec" != "1,2" && "$step_spec" != "1-2" && "$step_spec" != "2" && "$step_spec" != "all" ]]; then
	[[ -z "${reads}" ]] && {
		loge "Steps 3–5 need -l/--long ONT reads"
		exit 1
	}
fi

# ---------- step selector ----------
# supports: all | "1" | "1,2" | "3-5" | mixed like "1,3-5"
_expand_steps() {
	local spec="$1"
	local -a outA=()
	IFS=',' read -r -a parts <<<"$spec"
	for p in "${parts[@]}"; do
		if [[ "$p" =~ ^([0-9]+)-([0-9]+)$ ]]; then
			local a="${BASH_REMATCH[1]}" b="${BASH_REMATCH[2]}"
			local i
			for ((i = a; i <= b; i++)); do outA+=("$i"); done
		elif [[ "$p" =~ ^[0-9]+$ ]]; then
			outA+=("$p")
		fi
	done
	printf "%s\n" "${outA[@]}" | sort -n | uniq
}
declare -A _want=()
if [[ "$step_spec" == "all" ]]; then
	for s in 1 2 3 4 5; do _want[$s]=1; done
else
	while read -r s; do _want[$s]=1; done < <(_expand_steps "$step_spec")
fi
step_enabled() { [[ "${_want[$1]:-0}" -eq 1 ]]; }

# ---------- deps ----------
require_cmd() { command -v "$1" >/dev/null 2>&1 || {
	loge "Missing '$1'"
	exit 2
}; }
require_cmd python3
require_cmd Rscript
require_cmd samtools
require_cmd minimap2

mode=""
if [[ "$prefer" == "mummer" ]] && command -v nucmer >/dev/null 2>&1; then
	mode="mummer"
elif [[ "$prefer" == "blast" ]] && command -v blastn >/dev/null 2>&1; then
	mode="blast"
elif command -v nucmer >/dev/null 2>&1; then
	mode="mummer"
elif command -v blastn >/dev/null 2>&1; then
	mode="blast"
else
	loge "Need MUMmer4 (nucmer/show-coords/delta-filter) or BLAST+ (blastn/makeblastdb)"
	exit 2
fi
[[ "$mode" == "mummer" ]] && {
	require_cmd nucmer
	require_cmd delta-filter
	require_cmd show-coords
}
[[ "$mode" == "blast" ]] && {
	require_cmd blastn
	require_cmd makeblastdb
}

# ---------- layout ----------
mkdir -p "$out"/{align,reads,report}
rep_tsv="$out/repeats.tsv"
rep_bed="$out/repeats.bed"
# parents & recombinants
par_fa="$out/junc_parent.fa"
par_meta="$out/junc_parent.tsv"
par_supp="$out/junc_parent_support.tsv"
rec_fa="$out/junc_recomb.fa"
rec_meta="$out/junc_recomb.tsv"
rec_supp="$out/junc_recomb_support.tsv"
# depth files
mtref_bam="$out/reads/mtref.bam"
mtref_depth="$out/reads/mtref.depth.tsv"
fra_tsv="$out/repeat_fractions.tsv"
pdf="$out/report/report.pdf"

# ========== STEP 1: Repeat discovery ==========
if step_enabled 1; then
	logi "Step 1/5: Repeat discovery via ${mode} (L>=$min_rep_len, PID>=$min_rep_pid)"
	if [[ "$mode" == "mummer" ]]; then
		nucmer --maxmatch -p "$out/align/self" "$asm" "$asm" >/dev/null
		delta-filter -i "$min_rep_pid" -l "$min_rep_len" "$out/align/self.delta" >"$out/align/self.filt.delta"
		# keep header; no -H
		show-coords -T -r -c -l "$out/align/self.filt.delta" >"$out/align/self.coords.tsv"
		python3 "$_POLAPLIB_DIR/scripts/polap_py_mt_selfrepeats.py" \
			--mode mummer --coords "$out/align/self.coords.tsv" \
			--assembly "$asm" \
			--min-len "$min_rep_len" --min-pid "$min_rep_pid" \
			--out-tsv "$rep_tsv" --out-bed "$rep_bed" --log-level INFO
	else
		makeblastdb -in "$asm" -dbtype nucl -out "$out/align/mtdb" >/dev/null 2>&1 || true
		blastn -task megablast -db "$out/align/mtdb" -query "$asm" -evalue 1e-20 -dust no -soft_masking false \
			-perc_identity "$min_rep_pid" -word_size 28 \
			-outfmt "6 qseqid sseqid pident length qstart qend sstart send qcovs qlen slen" \
			>"$out/align/self.blast6.tsv"
		python3 "$_POLAPLIB_DIR/scripts/polap_py_mt_selfrepeats.py" \
			--mode blast --blast6 "$out/align/self.blast6.tsv" \
			--assembly "$asm" \
			--min-len "$min_rep_len" --min-pid "$min_rep_pid" \
			--out-tsv "$rep_tsv" --out-bed "$rep_bed" --log-level INFO
	fi
	logi "Step 1 done. Expect: $rep_tsv (tab) and $rep_bed (bed)."
fi

# ========== STEP 2: Junction synthesis (parental/recombinant) ==========
if step_enabled 2; then
	logi "Step 2/5: Junction synthesis (flank=${flank}bp) -> parental & recombinant"
	python3 "$_POLAPLIB_DIR/scripts/polap_py_synthesize_junctions_full.py" \
		--assembly "$asm" --repeats "$rep_tsv" --flank "$flank" \
		--out-parent-fasta "$par_fa" --out-parent-meta "$par_meta" \
		--out-recomb-fasta "$rec_fa" --out-recomb-meta "$rec_meta"
	logi "Step 2 done. Expect: $par_fa/$par_meta and $rec_fa/$rec_meta."
fi

# ========== STEP 3: Junction-bridge counting ==========
if step_enabled 3; then
	logi "Step 3/5: Mapping ONT reads to junctions (MAPQ>=$min_mapq, span>=$min_span)"
	mkdir -p "$out/reads"
	# Parentals
	minimap2 -x map-ont -a -t "$threads" "$par_fa" "$reads" |
		samtools view -bh -@ "$threads" - |
		samtools sort -@ "$threads" -o "$out/reads/parent.bam"
	samtools index -@ "$threads" "$out/reads/parent.bam"
	python3 "$_POLAPLIB_DIR/scripts/polap_py_count_junction_support.py" \
		--bam "$out/reads/parent.bam" --fasta "$par_fa" --meta "$par_meta" \
		--min-mapq "$min_mapq" --min-span "$min_span" --out-tsv "$par_supp"
	# Recombinants
	minimap2 -x map-ont -a -t "$threads" "$rec_fa" "$reads" |
		samtools view -bh -@ "$threads" - |
		samtools sort -@ "$threads" -o "$out/reads/recomb.bam"
	samtools index -@ "$threads" "$out/reads/recomb.bam"
	python3 "$_POLAPLIB_DIR/scripts/polap_py_count_junction_support.py" \
		--bam "$out/reads/recomb.bam" --fasta "$rec_fa" --meta "$rec_meta" \
		--min-mapq "$min_mapq" --min-span "$min_span" --out-tsv "$rec_supp"
	logi "Step 3 done. Expect: $par_supp and $rec_supp with per-junction bridging counts."
fi

# ========== STEP 4: Fraction quantification ==========
if step_enabled 4; then
	logi "Step 4/5: Quantifying per-repeat recombinant fraction (f_hat = R/(R+P))"
	# optional normalization by flank depth (map reads to assembly)
	minimap2 -x map-ont -a -t "$threads" "$asm" "$reads" |
		samtools view -bh -@ "$threads" - |
		samtools sort -@ "$threads" -o "$mtref_bam"
	samtools index -@ "$threads" "$mtref_bam"
	samtools depth -aa "$mtref_bam" >"$mtref_depth"

	python3 "$_POLAPLIB_DIR/scripts/polap_py_quantify_repeat_fraction.py" \
		--repeats "$rep_tsv" \
		--parent "$par_supp" --parent-meta "$par_meta" \
		--recomb "$rec_supp" --recomb-meta "$rec_meta" \
		--depth "$mtref_depth" --flank "$flank" --label "$label" \
		--out-tsv "$fra_tsv"
	logi "Step 4 done. Expect: $fra_tsv with P/R raw & normalized counts, f_hat and 95% CI."
fi

# ========== STEP 5: Report PDF ==========
if step_enabled 5; then
	logi "Step 5/5: Rendering summary PDF"
	Rscript --vanilla "$_POLAPLIB_DIR/scripts/polap_r_stoich_plots.R" \
		"$fra_tsv" "$rep_tsv" "$par_meta" "$rec_meta" "$pdf" "$label"
	logi "Step 5 done. Expect: $pdf"
fi

logi "All requested steps completed."
