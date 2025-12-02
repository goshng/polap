#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step6-a2a-mask.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -A asm.fa -O other.fa -t threads -i minIdent -l minLen -p pad -d mask_dir -Q run_mask_qc(0|1)
while getopts ":A:O:t:i:l:p:d:Q:" opt; do
	case "$opt" in
	A) ASM="$OPTARG" ;;
	O) OTH="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	i) MINID="$OPTARG" ;;
	l) MINLEN="$OPTARG" ;;
	p) PAD="$OPTARG" ;;
	d) MD="$OPTARG" ;;
	Q) RUNQC="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$ASM" && -s "$OTH" && -n "$MD" ]] || {
	echo "args?"
	exit 2
}
install -d -- "$MD"
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"
PAF_FILTER="${SCRIPTS_DIR}/polap_py_paf_filter.py"
MASK_FROM_PAF="${SCRIPTS_DIR}/polap_py_mask_from_paf.py"
BED_COMPLEMENT="${SCRIPTS_DIR}/polap_py_bed_complement.py"
PAF_TO_HITS="${SCRIPTS_DIR}/polap_py_paf_to_hits_tsv.py"
BED_ADD_LEN="${SCRIPTS_DIR}/polap_py_bed_add_len.py"
MASK_QC_R="${SCRIPTS_DIR}/polap_r_mask_qc.R"

A2A_RAW="$MD/a2a.raw.paf"
A2A_FILT="$MD/a2a.filt.paf"
MASK_BED="$MD/transfer.mask.bed"
ALLOW_BED="$MD/allow.bed"

minimap2 -x asm10 -t "$THREADS" -c --secondary=no "$ASM" "$OTH" >"$A2A_RAW" 2>"$MD/a2a.minimap2.err"
python3 "$PAF_FILTER" --min-ident "$MINID" --min-alen "$MINLEN" --in "$A2A_RAW" --out "$A2A_FILT" 2>>"$MD/a2a.filter.err"
python3 "$MASK_FROM_PAF" --paf "$A2A_FILT" --fasta "$ASM" --pad "$PAD" --out-bed "$MASK_BED" 2>>"$MD/mask_from_paf.err"
python3 "$BED_COMPLEMENT" --fasta "$ASM" --bed "$MASK_BED" --out-bed "$ALLOW_BED" 2>>"$MD/bed_complement.err"
python3 "$PAF_TO_HITS" --paf "$A2A_FILT" --out "$MD/transfer.mask.hits.tsv" 2>>"$MD/paf_to_hits.err"
python3 "$BED_ADD_LEN" --fasta "$ASM" --bed "$MASK_BED" --out "$MD/transfer.mask.merged.tsv" 2>>"$MD/bed_add_len.err"
if [[ "$RUNQC" == "1" ]] && command -v Rscript >/dev/null 2>&1 && [[ -s "$MASK_QC_R" ]]; then
	Rscript --vanilla "$MASK_QC_R" --fasta "$ASM" --bed "$MASK_BED" \
		--outpdf "$MD/mask_qc.pdf" --outtsv "$MD/mask_qc.tsv" >/dev/null 2>&1 || true
fi
