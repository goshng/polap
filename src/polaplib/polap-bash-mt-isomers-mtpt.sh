#!/usr/bin/env bash
# polap-bash-mt-isomers-mtpt.sh
# Version : v1.1.0
# Purpose : (A) Mitochondrial isomer evidence (repeats → junctions → ONT bridges)
#           (B) MTPT turnover (cp → mt tracts, recency bins)
#           with step-by-step execution control.
#
# STEPS (select with --step):
#   1) Repeat discovery (MUMmer/BLAST)                      → repeats.tsv / repeats.bed
#   2) Junction synthesis (DIR/INV; flank)                  → junctions.fasta / junctions.tsv / network.gfa
#   3) Junction-bridge counting (ONT → junction templates)  → junction_support.tsv
#   4) MTPT scan (cp → mt BLAST; collapse; bin by PID)        → mtpt/mtpt.tsv / mtpt/mtpt.bed
#   5) MTPT–repeat enrichment (permutation; ±window)        → enrichment/mtpt_repeat_enrichment.tsv
#   6) Report PDF (compact overview)                        → report/report.pdf
#
# REQUIRED INPUTS (flags):
#   -f, --fasta   <FASTA>   mitochondrial assembly (reference)
#   -o, --out     <DIR>     output directory
#   -l, --long    <FASTQ>   ONT reads (needed for step ≥3 and for depth in report)
#   -c, --cp      <FASTA>   chloroplast FASTA (needed for step ≥4/5)
#
# OTHER OPTIONS:
#   --label STR            [sample]     tag recorded in TSV/PDF
#   --threads INT          [8]
#   --min-repeat-len INT   [1000]       repeat threshold (bp)
#   --min-repeat-pid INT   [95]         repeat threshold (% identity)
#   --flank INT            [500]        junction flank per side (bp)
#   --min-mapq INT         [20]         bridge MAPQ filter
#   --min-span INT         [200]        required span per side (bp)
#   --prefer mummer|blast  [mummer]     backend for Step 1
#   --mtpt-min-len INT     [100]        minimum cp → mt tract length (bp)
#   --mtpt-min-pid INT     [85]         minimum cp → mt PID for hits
#   --mtpt-recent INT      [97]         PID bin: recent ≥ this
#   --mtpt-intermediate INT[90]         PID bin: intermediate ≥ this
#   --enrich-window INT    [1000]       ±bp window for MTPT–repeat enrichment
#   --enrich-perm INT      [1000]       permutations for enrichment p-value
#   --step SPEC            [all]        e.g., 1,2    or 3-6   or 1-2,4,6   or all
#
# USAGE EXAMPLE:
#   bash "$_POLAPLIB_DIR/polap-bash-mt-isomers-mtpt.sh" \
#     -f Anthoceros_agrestis/v6/0/polap-readassemble-1-miniasm/mt.1.fa \
#     -l Anthoceros_agrestis/tmp/l.fq \
#     -c Anthoceros_agrestis/pt.fa \
#     -o Anthoceros_agrestis-0/isomers-mtpt-S1/ \
#     --label S1 --threads 12 \
#     --min-repeat-len 1000 --min-repeat-pid 90 \
#     --flank 800 --min-mapq 20 --min-span 250 --prefer mummer \
#     --mtpt-min-len 150 --mtpt-min-pid 85 --mtpt-recent 97 --mtpt-intermediate 90 \
#     --enrich-window 1000 --enrich-perm 2000 --step all
#
set -euo pipefail
set -o errtrace
IFS=$'\n\t'

# -------- locate polaplib ----------
_POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
export _POLAPLIB_DIR POLAPLIB_DIR="${_POLAPLIB_DIR}"

# -------- tiny logger --------------
_now() { date "+%Y-%m-%d %H:%M:%S"; }
logi() { printf "[%s %s:%s] %s\n" "$(_now)" "${BASH_SOURCE[1]##*/}" "${BASH_LINENO[0]}" "$*" >&1; }
loge() { printf "[%s %s:%s] ERROR: %s\n" "$(_now)" "${BASH_SOURCE[1]##*/}" "${BASH_LINENO[0]}" "$*" >&2; }

# -------- defaults -----------------
asm="" out="" reads="" cpfa=""
label="sample"
threads=8
min_rep_len=1000
min_rep_pid=95
flank=500
min_mapq=20
min_span=200
prefer="mummer"
mtpt_min_len=100
mtpt_min_pid=85
mtpt_recent=97
mtpt_inter=90
enrich_win=1000
enrich_perm=1000
step_spec="all"

print_help() {
	grep -E '^#( |$)' "$0" | sed 's/^# \{0,1\}//'
}

# -------- parse args ---------------
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
	-c | --cp)
		cpfa="$2"
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
	--mtpt-min-len)
		mtpt_min_len="$2"
		shift 2
		;;
	--mtpt-min-pid)
		mtpt_min_pid="$2"
		shift 2
		;;
	--mtpt-recent)
		mtpt_recent="$2"
		shift 2
		;;
	--mtpt-intermediate)
		mtpt_inter="$2"
		shift 2
		;;
	--enrich-window)
		enrich_win="$2"
		shift 2
		;;
	--enrich-perm)
		enrich_perm="$2"
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
	loge "Need -f/--fasta and -o/--out"
	print_help
	exit 1
}

# reads required if running step ≥3; cpfa required if running step ≥4/5
needs_reads() { [[ "$1" -ge 3 ]]; }
needs_cp() { [[ "$1" -ge 4 ]]; }

# step expansion: "all" or "1,2" or "3-6" etc.
_expand_steps() {
	local spec="$1"
	local -a a=()
	IFS=',' read -r -a parts <<<"$spec"
	for p in "${parts[@]}"; do
		if [[ "$p" =~ ^([0-9]+)-([0-9]+)$ ]]; then
			local i
			for ((i = ${BASH_REMATCH[1]}; i <= ${BASH_REMATCH[2]}; i++)); do a+=("$i"); done
		elif [[ "$p" =~ ^[0-9]+$ ]]; then
			a+=("$p")
		fi
	done
	printf "%s\n" "${a[@]}" | sort -n | uniq
}
declare -A WANT=()
if [[ "$step_spec" == "all" ]]; then
	for s in 1 2 3 4 5 6; do WANT[$s]=1; done
else while read -r s; do WANT[$s]=1; done < <(_expand_steps "$step_spec"); fi
step_enabled() { [[ "${WANT[$1]:-0}" -eq 1 ]]; }

# -------- deps ---------------------
require_cmd() { command -v "$1" >/dev/null 2>&1 || {
	loge "Missing '$1'"
	exit 2
}; }
require_cmd python3
require_cmd Rscript
require_cmd minimap2
require_cmd samtools

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
	loge "Need MUMmer4 or BLAST+"
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

# -------- layout -------------------
mkdir -p "$out"/{align,reads,mtpt,enrichment,report}
rep_tsv="$out/repeats.tsv"
rep_bed="$out/repeats.bed"
junc_fa="$out/junctions.fasta"
junc_meta="$out/junctions.tsv"
net_gfa="$out/network.gfa"
junc_bam="$out/reads/junc.bam"
junc_supp="$out/junction_support.tsv"
mtpt_tsv="$out/mtpt/mtpt.tsv"
mtpt_bed="$out/mtpt/mtpt.bed"
enr_tsv="$out/enrichment/mtpt_repeat_enrichment.tsv"
pdf="$out/report/report.pdf"

# =========================================================
# STEP 1: Repeat discovery
# =========================================================
if step_enabled 1; then
	logi "Step 1/6: Repeat discovery via ${mode} (L>=$min_rep_len, PID>=$min_rep_pid)"
	if [[ "$mode" == "mummer" ]]; then
		nucmer --maxmatch -p "$out/align/self" "$asm" "$asm" >/dev/null
		delta-filter -i "$min_rep_pid" -l "$min_rep_len" "$out/align/self.delta" >"$out/align/self.filt.delta"
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
	logi "Step 1 done → expect: $rep_tsv , $rep_bed"
fi

# =========================================================
# STEP 2: Junction synthesis
# =========================================================
if step_enabled 2; then
	logi "Step 2/6: Junction synthesis (flank=${flank}bp)"
	python3 "$_POLAPLIB_DIR/scripts/polap_py_predict_mt_isomers.py" \
		--assembly "$asm" --repeats "$rep_tsv" --flank "$flank" \
		--out-fasta "$junc_fa" --out-meta "$junc_meta" --out-gfa "$net_gfa"
	logi "Step 2 done → expect: $junc_fa , $junc_meta , $net_gfa"
fi

# =========================================================
# STEP 3: Junction-bridge counting (ONT → junctions)
# =========================================================
if step_enabled 3; then
	[[ -z "$reads" ]] && {
		loge "Step 3 needs ONT reads (-l/--long)"
		exit 1
	}
	logi "Step 3/6: Mapping ONT to junction templates (MAPQ>=$min_mapq, span>=$min_span)"
	minimap2 -x map-ont -a -t "$threads" "$junc_fa" "$reads" |
		samtools view -bh -@ "$threads" - |
		samtools sort -@ "$threads" -o "$junc_bam"
	samtools index -@ "$threads" "$junc_bam"
	python3 "$_POLAPLIB_DIR/scripts/polap_py_count_junction_support.py" \
		--bam "$junc_bam" --fasta "$junc_fa" --meta "$junc_meta" \
		--min-mapq "$min_mapq" --min-span "$min_span" --out-tsv "$junc_supp"
	logi "Step 3 done → expect: $junc_supp (per-junction bridging counts)"
fi

# =========================================================
# STEP 4: MTPT (cp → mt) scan & bins
# =========================================================
if step_enabled 4; then
	[[ -z "$cpfa" ]] && {
		loge "Step 4 needs cp FASTA (-c/--cp)"
		exit 1
	}
	logi "Step 4/6: MTPT scan (cp → mt) with min_len=$mtpt_min_len, min_pid=$mtpt_min_pid"
	makeblastdb -in "$asm" -dbtype nucl -out "$out/mtpt/mt" >/dev/null 2>&1 || true
	blastn -task megablast -db "$out/mtpt/mt" -query "$cpfa" -evalue 1e-20 -dust no -soft_masking false \
		-perc_identity "$mtpt_min_pid" -word_size 28 \
		-outfmt "6 qseqid sseqid pident length qstart qend sstart send qcovs" \
		>"$out/mtpt/raw.blast6.tsv"
	python3 "$_POLAPLIB_DIR/scripts/polap_py_mtpt_scan.py" \
		--blast6 "$out/mtpt/raw.blast6.tsv" \
		--mt-fasta "$asm" \
		--min-len "$mtpt_min_len" \
		--recent "$mtpt_recent" --intermediate "$mtpt_inter" \
		--out-tsv "$mtpt_tsv" --out-bed "$mtpt_bed"
	logi "Step 4 done → expect: $mtpt_tsv , $mtpt_bed"
fi

# =========================================================
# STEP 5: MTPT–repeat enrichment (optional)
# =========================================================
if step_enabled 5; then
	[[ -s "$mtpt_bed" && -s "$rep_bed" ]] || { logi "Step 5 skipped (need $mtpt_bed and $rep_bed)"; }
	if [[ -s "$mtpt_bed" && -s "$rep_bed" ]]; then
		logi "Step 5/6: Enrichment test (MTPT within ±${enrich_win} bp of repeats; ${enrich_perm} perms)"
		python3 "$_POLAPLIB_DIR/scripts/polap_py_mtpt_repeat_enrichment.py" \
			--sample "$label" \
			--mtpt-bed "$mtpt_bed" \
			--repeat-bed "$rep_bed" \
			--mt-fasta "$asm" \
			--window "$enrich_win" \
			--permutations "$enrich_perm" \
			--out-tsv "$enr_tsv"
		logi "Step 5 done → expect: $enr_tsv"
	fi
fi

# =========================================================
# STEP 6: Report PDF
# =========================================================
if step_enabled 6; then
	logi "Step 6/6: Rendering report PDF"
	Rscript --vanilla "$_POLAPLIB_DIR/scripts/polap_r_mt_isomer_mtpt_report.R" \
		"$rep_tsv" "$junc_meta" "$junc_supp" "$mtpt_tsv" "$pdf"
	logi "Step 6 done → expect: $pdf"
fi

logi "Requested steps completed."
