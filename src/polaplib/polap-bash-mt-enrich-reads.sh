#!/usr/bin/env bash
################################################################################
# polap-bash-mt-enrich-reads.sh
# Version : v0.1.3
# Author  : POLAP
# License : GPL-3.0+
#
# Purpose :
#   Aggressively remove plastid (ptDNA) reads and keep confident mitochondrial
#   (mtDNA) reads from ONT WGS using minimap2 -> pysam classifier -> seqtk, with
#   optional coverage downsampling (rasusa) and QC plots (R).
################################################################################

# ---- debug guard: refuse being sourced; require execution
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
	echo "[ERROR] This script must be executed, not sourced: use 'bash ${BASH_SOURCE[0]}'" >&2
	return 1 2>/dev/null || exit 1
fi

set -euo pipefail
IFS=$'\n\t'

# ---- version & tracing wiring (POLAP_DEBUG=1 -> set -x)
_VERSION="v0.1.3"
: "${POLAP_DEBUG:=0}"
if [[ "$POLAP_DEBUG" -eq 1 ]]; then
	set -x
fi

# ---- resolve script dir to find companion scripts robustly
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${_here}/scripts"

usage() {
	cat <<'EOF'
Usage:
  polap-bash-mt-enrich-reads.sh \
    --reads ont.fq.gz \
    --mt-ref mt.fasta \
    --pt-ref pt.fasta \
    --out outdir [options]

Options:
  -t, --threads N            [8] threads for minimap2/samtools
      --id-mt F             [0.75] min identity for mt keep
      --frac-mt F           [0.40] min aligned fraction (aln_len/read_len) for mt keep
      --id-pt F             [0.85] min identity for pt removal
      --frac-pt F           [0.30] min aligned fraction for pt removal
      --max-clip N          [4000] max soft+hard clipping (bp)
      --min-read N          [2000] min read length to consider (bp)
      --ambig drop|keep|none [drop] policy for ambiguous reads
      --emit-metrics        write per-read metrics TSV for QC plotting
      --downsample COV GENOMESIZE  e.g. "--downsample 80x 700000"

  # NEW (kept off by default):
      --primary-only        only consider primary alignments
      --min-mapq N          [0] ignore alignments with MAPQ < N
      --margin F            [0.10] tie-break margin (10% -> 1.10x)

  -h, --help
  -V, --version

Outputs (in OUT):
  mt.keep.fq.gz, pt.drop.fq.gz, ambig.fq.gz, none.fq.gz, mtpt.summary.tsv
  qc/ (plots if R available), tmp/{pt.bam,mt.bam,*.mmi (if built)}
EOF
}

# ---- quick flag handling (before heavy parsing)
if [[ "${1:-}" == "--version" || "${1:-}" == "-V" ]]; then
	echo "polap-bash-mt-enrich-reads.sh ${_VERSION}"
	exit 0
fi
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
	usage
	exit 0
fi

# ---------- defaults
reads="" mt_ref="" pt_ref="" out=""
threads=8
id_mt=0.75 frac_mt=0.40
id_pt=0.85 frac_pt=0.30
max_clip=4000 min_read=2000
ambig_policy="drop"
emit_metrics=0
down_cov="" down_gsize=""

# NEW defaults
primary_only=0
min_mapq=0
margin=0.10

# ---------- parse args
while (($#)); do
	case "$1" in
	--reads)
		reads="${2:?}"
		shift 2
		;;
	--mt-ref)
		mt_ref="${2:?}"
		shift 2
		;;
	--pt-ref)
		pt_ref="${2:?}"
		shift 2
		;;
	--out)
		out="${2:?}"
		shift 2
		;;
	-t | --threads)
		threads="${2:?}"
		shift 2
		;;
	--id-mt)
		id_mt="${2:?}"
		shift 2
		;;
	--frac-mt)
		frac_mt="${2:?}"
		shift 2
		;;
	--id-pt)
		id_pt="${2:?}"
		shift 2
		;;
	--frac-pt)
		frac_pt="${2:?}"
		shift 2
		;;
	--max-clip)
		max_clip="${2:?}"
		shift 2
		;;
	--min-read)
		min_read="${2:?}"
		shift 2
		;;
	--ambig)
		ambig_policy="${2:?}"
		shift 2
		;;
	--emit-metrics)
		emit_metrics=1
		shift
		;;
	--downsample)
		down_cov="${2:?}"
		down_gsize="${3:?}"
		shift 3
		;;
	# NEW
	--primary-only)
		primary_only=1
		shift
		;;
	--min-mapq)
		min_mapq="${2:?}"
		shift 2
		;;
	--margin)
		margin="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	-V | --version)
		echo "polap-bash-mt-enrich-reads.sh ${_VERSION}"
		exit 0
		;;
	*)
		echo "[ERROR] Unknown arg: $1" >&2
		usage
		exit 2
		;;
	esac
done

# ---------- checks
[[ -n "$reads" ]] || {
	echo "[ERROR] --reads is required" >&2
	usage
	exit 2
}
[[ -n "$mt_ref" ]] || {
	echo "[ERROR] --mt-ref is required" >&2
	usage
	exit 2
}
[[ -n "$pt_ref" ]] || {
	echo "[ERROR] --pt-ref is required" >&2
	usage
	exit 2
}
[[ -n "$out" ]] || {
	echo "[ERROR] --out is required" >&2
	usage
	exit 2
}
[[ -f "$reads" ]] || {
	echo "[ERROR] --reads not found: $reads" >&2
	exit 2
}
[[ -f "$mt_ref" ]] || {
	echo "[ERROR] --mt-ref not found: $mt_ref" >&2
	exit 2
}
[[ -f "$pt_ref" ]] || {
	echo "[ERROR] --pt-ref not found: $pt_ref" >&2
	exit 2
}

for tool in minimap2 samtools seqtk python3; do
	command -v "$tool" >/dev/null 2>&1 || {
		echo "[ERROR] missing tool: $tool" >&2
		exit 3
	}
done
python3 - <<'PY' >/dev/null 2>&1 || {
import pysam
PY
	echo "[ERROR] python 'pysam' missing (pip/conda install pysam)" >&2
	exit 3
}
R_OK=0
command -v Rscript >/dev/null 2>&1 && R_OK=1

mkdir -p "$out"/{tmp,qc}
[[ -f "${SCRIPTS_DIR}/filter_mtpt_from_bam.py" ]] || {
	echo "[ERROR] Missing ${SCRIPTS_DIR}/filter_mtpt_from_bam.py" >&2
	exit 3
}

# ---------- index (create .mmi on demand)
[[ -f "${mt_ref}.mmi" ]] || minimap2 -d "${mt_ref}.mmi" "$mt_ref"
[[ -f "${pt_ref}.mmi" ]] || minimap2 -d "${pt_ref}.mmi" "$pt_ref"

# ---------- map to pt and mt (SAM -> sorted BAM)
pt_bam="$out/tmp/pt.bam"
mt_bam="$out/tmp/mt.bam"

minimap2 -x map-ont -a -t "$threads" "${pt_ref}.mmi" "$reads" |
	samtools sort -@ "$threads" -o "$pt_bam"
samtools index -@ "$threads" "$pt_bam"

minimap2 -x map-ont -a -t "$threads" "${mt_ref}.mmi" "$reads" |
	samtools sort -@ "$threads" -o "$mt_bam"
samtools index -@ "$threads" "$mt_bam"

# ---------- classify
metrics_flag=()
[[ $emit_metrics -eq 1 ]] && metrics_flag=(--emit-metrics "$out/mtpt.metrics.tsv")

py_args=(
	--pt-bam "$pt_bam"
	--mt-bam "$mt_bam"
	--id-mt "$id_mt" --frac-mt "$frac_mt"
	--id-pt "$id_pt" --frac-pt "$frac_pt"
	--max-clip "$max_clip" --min-read "$min_read"
	--ambig "$ambig_policy"
	--margin "$margin"
	--min-mapq "$min_mapq"
	--out-prefix "$out/mtpt"
	"${metrics_flag[@]}"
)
[[ $primary_only -eq 1 ]] && py_args=(--primary-only "${py_args[@]}")

python3 "${SCRIPTS_DIR}/filter_mtpt_from_bam.py" "${py_args[@]}"

# ---------- materialize FASTQs
declare -A TAG2OUT=(["keep"]="mt.keep.fq.gz" ["drop"]="pt.drop.fq.gz"
	["ambig"]="ambig.fq.gz" ["none"]="none.fq.gz")

for tag in keep drop ambig none; do
	ids="$out/mtpt.${tag}_ids.txt"
	if [[ -s "$ids" ]]; then
		seqtk subseq "$reads" "$ids" | gzip -c >"$out/${TAG2OUT[$tag]}"
	else
		: >"$out/${TAG2OUT[$tag]}"
	fi
done

# ---------- optional downsample (requires rasusa)
if [[ -n "$down_cov" && -n "$down_gsize" ]]; then
	if command -v rasusa >/dev/null 2>&1; then
		rasusa reads -c "$down_cov" -g "$down_gsize" -o "$out/mt.keep.down.fq.gz" --input "$out/mt.keep.fq.gz"
	else
		echo "[WARN] rasusa not found; skipping --downsample" >&2
	fi
fi

# ---------- QC plots
if [[ $R_OK -eq 1 && -f "$out/mtpt.metrics.tsv" ]]; then
	if [[ -f "${SCRIPTS_DIR}/plot_mtpt_qc.R" ]]; then
		Rscript "${SCRIPTS_DIR}/plot_mtpt_qc.R" "$out/mtpt.metrics.tsv" "$out/qc"
	fi
fi

# ---------- summary (counts + thresholds)
count_or_zero() { [[ -s "$1" ]] && wc -l <"$1" || echo 0; }
cat >"$out/mtpt.summary.tsv" <<EOF
#key	value
n_keep	$(count_or_zero "$out/mtpt.keep_ids.txt")
n_drop	$(count_or_zero "$out/mtpt.drop_ids.txt")
n_ambig	$(count_or_zero "$out/mtpt.ambig_ids.txt")
n_none	$(count_or_zero "$out/mtpt.none_ids.txt")
id_mt	${id_mt}
frac_mt	${frac_mt}
id_pt	${id_pt}
frac_pt	${frac_pt}
max_clip	${max_clip}
min_read	${min_read}
ambig_policy	${ambig_policy}
primary_only	${primary_only}
min_mapq	${min_mapq}
margin	${margin}
EOF

# ---------- next step
cat <<EOF
[OK] mt-first enrichment complete in: $out
  - mt.keep.fq.gz
  - pt.drop.fq.gz
  - ambig.fq.gz
  - none.fq.gz
  - mtpt.summary.tsv
  - qc/ (if metrics & R)

Suggested assembly:
  flye --nano-raw "$out/mt.keep.fq.gz" -o "$out/flye_mt" -g 700k -t $threads
EOF
