#!/usr/bin/env bash
# polap-bash-evo-quick.sh
# Purpose : One-button evolutionary package for ONT mt/pt:
#           (i) per-repeat recombination usage (stoichiometry),
#           (ii) MTPT load & recency,
#           (iii) MTPT–repeat enrichment,
#           (iv) recombinant-junction vs CDS proximity,
#           (v) cross-sample summary + PDF report.
#
# Usage   : polap-bash-evo-quick.sh <manifest.tsv> <outdir> [options]
#
# Options :
#   --threads INT            [8]
#   --min-repeat-len INT     [1000]
#   --min-repeat-pid INT     [95]
#   --flank INT              [500]
#   --min-mapq INT           [20]
#   --min-span INT           [200]
#   --mtpt-min-len INT       [100]
#   --mtpt-min-pid INT       [85]
#   --mtpt-recent INT        [98]
#   --mtpt-intermediate INT  [90]
#   --enrich-window INT      [1000]   # MTPT–repeat proximity window (bp)
#   --enrich-perm INT        [1000]   # permutations for enrichment p-value
#   --cds-threshold INT      [500]    # "near gene" threshold (bp)
#   --cds-perm INT           [5000]   # permutations for CDS proximity p-value
#   --prefer mummer|blast    [mummer]
#
# Deps    : minimap2, samtools, MUMmer4 or BLAST+, python3 (biopython,pysam), R (ggplot2,data.table,gridExtra)
# Notes   : Expects the two launchers in PATH:
#             - polap-bash-mt-stoich-fractions.sh
#             - polap-bash-mt-isomers-mtpt.sh
set -euo pipefail
shopt -s nullglob

if [[ $# -lt 2 ]]; then
	echo "Usage: $0 <manifest.tsv> <outdir> [options]" >&2
	exit 1
fi

man="$1"
out="$2"
shift 2 || true

threads=8
minL=1000
minPID=95
flank=500
minq=20
mspan=200
mtpt_min_len=100
mtpt_min_pid=85
mtpt_recent=98
mtpt_inter=90
enrich_win=1000
enrich_perm=1000
cds_thr=500
cds_perm=5000
prefer="mummer"

while [[ $# -gt 0 ]]; do
	case "$1" in
	--threads)
		threads="$2"
		shift 2
		;;
	--min-repeat-len)
		minL="$2"
		shift 2
		;;
	--min-repeat-pid)
		minPID="$2"
		shift 2
		;;
	--flank)
		flank="$2"
		shift 2
		;;
	--min-mapq)
		minq="$2"
		shift 2
		;;
	--min-span)
		mspan="$2"
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
	--cds-threshold)
		cds_thr="$2"
		shift 2
		;;
	--cds-perm)
		cds_perm="$2"
		shift 2
		;;
	--prefer)
		prefer="$2"
		shift 2
		;;
	*)
		echo "Unknown option: $1" >&2
		exit 1
		;;
	esac
done

mkdir -p "$out"/{samples,summary,logs}

# loop over manifest
tail -n +2 "$man" | awk -F'\t' 'NF>=4 && $1!~/^#/ {print}' >"$out/summary/_rows.tsv" || true
if [[ ! -s "$out/summary/_rows.tsv" ]]; then
	echo "Manifest has no usable rows (need: sample, mt_fasta, cp_fasta, reads_fastq)." >&2
	exit 2
fi

# aggregate holders
>"$out/summary/stoich_all.tsv" # sample, pair_id, orient, repeat_len, pid, f_hat
echo -e "sample\tobs_overlap\tmtpt_total\texp_mean\tz\tp_emp\toverlap_frac\twindow_bp" >"$out/summary/enrichment_summary.tsv"
echo -e "sample\tN_junc\tmean_dist\tfrac_<${cds_thr}\texp_mean_dist\texp_frac\temp_p_mean\temp_p_frac" >"$out/summary/cds_proximity_summary.tsv"

while IFS=$'\t' read -r sample mtfa cpfa reads cds_gff group; do
	sdir="$out/samples/$sample"
	mkdir -p "$sdir"
	echo "[$(date +'%F %T')] >>> Sample: $sample" >&2

	# 1) Per-repeat stoichiometry
	sto="$sdir/stoich"
	polap-bash-mt-stoich-fractions.sh "$mtfa" "$reads" "$sto" \
		--label "$sample" --threads "$threads" \
		--min-repeat-len "$minL" --min-repeat-pid "$minPID" --flank "$flank" \
		--min-mapq "$minq" --min-span "$mspan" --prefer "$prefer" \
		>"$sdir/stoich.log" 2>&1 || true

	# 2) MTPT (and repeat catalog for enrichment)
	tmpman="$sdir/_one.tsv"
	{
		echo -e "sample\tmt_fasta\tcp_fasta\treads_fastq\tgroup"
		echo -e "$sample\t$mtfa\t$cpfa\t$reads\t${group:-novel}"
	} >"$tmpman"

	mtpt="$sdir/mtpt"
	polap-bash-mt-isomers-mtpt.sh "$tmpman" "$mtpt" \
		--threads "$threads" \
		--min-repeat-len "$minL" --min-repeat-pid "$minPID" --flank "$flank" \
		--min-mapq "$minq" --min-span "$mspan" \
		--mtpt-min-len "$mtpt_min_len" --mtpt-min-pid "$mtpt_min_pid" \
		--mtpt-recent "$mtpt_recent" --mtpt-intermediate "$mtpt_inter" \
		--prefer "$prefer" \
		>"$sdir/mtpt.log" 2>&1 || true

	# 3) MTPT–repeat enrichment
	# repeats from stoich stage: $sto/repeats.tsv -> bed at $sto/repeats.bed ; MTPT bed at $mtpt/$sample/mtpt/mtpt.bed
	mtpt_bed="$mtpt/$sample/mtpt/mtpt.bed"
	rep_bed="$sto/repeats.bed"
	enr="$sdir/enrichment"
	mkdir -p "$enr"
	python3 scripts/polap_py_mtpt_repeat_enrichment.py \
		--sample "$sample" \
		--mtpt-bed "$mtpt_bed" \
		--repeat-bed "$rep_bed" \
		--mt-fasta "$mtfa" \
		--window "$enrich_win" \
		--permutations "$enrich_perm" \
		--out-tsv "$enr/mtpt_repeat_enrichment.tsv" \
		>"$sdir/enrich.log" 2>&1 || true
	if [[ -s "$enr/mtpt_repeat_enrichment.tsv" ]]; then
		tail -n +2 "$enr/mtpt_repeat_enrichment.tsv" >>"$out/summary/enrichment_summary.tsv"
	fi

	# 4) Recombinant-junction vs CDS proximity (skip if CDS GFF missing)
	if [[ -n "${cds_gff:-}" && -s "$cds_gff" ]]; then
		cdsout="$sdir/cds"
		mkdir -p "$cdsout"
		# recombinant metadata made by stoich launcher: junc_recomb.tsv
		python3 scripts/polap_py_junction_vs_cds_distance.py \
			--sample "$sample" \
			--junctions "$sto/junc_recomb.tsv" \
			--gff "$cds_gff" \
			--threshold "$cds_thr" \
			--permutations "$cds_perm" \
			--out-tsv "$cdsout/junc_cds_distance.tsv" \
			>"$sdir/cds.log" 2>&1 || true
		if [[ -s "$cdsout/junc_cds_distance.tsv" ]]; then
			tail -n +2 "$cdsout/junc_cds_distance.tsv" >>"$out/summary/cds_proximity_summary.tsv"
		fi
	fi

	# 5) Collect per-repeat f_hat for cross-sample plot
	if [[ -s "$sto/repeat_fractions.tsv" && -s "$sto/repeats.tsv" ]]; then
		# sample  pair_id orient repeat_len pid f_hat
		awk -v S="$sample" 'BEGIN{FS=OFS="\t"}
      FNR==NR{if(NR>1) {rep[$1]=$4; pid[$1]=$5; next}} 
      FNR>1 {print S,$2,$3,rep[$2]+0,pid[$2]+0,$8+0.0}' \
			"$sto/repeats.tsv" "$sto/repeat_fractions.tsv" >>"$out/summary/stoich_all.tsv"
	fi

done <"$out/summary/_rows.tsv"

# 6) Build project-level summary & PDF
python3 scripts/polap_py_collect_evo_quick_summary.py \
	--manifest "$man" \
	--workdir "$out/samples" \
	--stoich-all "$out/summary/stoich_all.tsv" \
	--enrichment "$out/summary/enrichment_summary.tsv" \
	--cds-prox "$out/summary/cds_proximity_summary.tsv" \
	--out-tsv "$out/summary/evo_quick_summary.tsv"

Rscript scripts/polap_r_evo_quick_report.R \
	"$out/summary/evo_quick_summary.tsv" \
	"$out/summary/stoich_all.tsv" \
	"$out/summary/enrichment_summary.tsv" \
	"$out/summary/cds_proximity_summary.tsv" \
	"$out/summary/evo_quick_report.pdf"

echo "DONE. Key outputs:"
echo "  $out/summary/evo_quick_summary.tsv"
echo "  $out/summary/evo_quick_report.pdf"
