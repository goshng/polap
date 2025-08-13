#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./ont100k.sh -r reads.fq.gz -o outdir [-g mt_genes.fa] [-t 100000] [-N 5] [-T 16]
#
# Requires: minimap2, awk, python3; (optional) seqkit

reads=""
genes=""
out="ont100k_out"
target=100000
n_contigs=5
threads=${THREADS:-8}
min_id=0.90  # ONT identity threshold (fraction)
min_ovl=4000 # min overlap length (bp)

while getopts "r:g:o:t:N:T:I:L:" opt; do
	case $opt in
	r) reads="$OPTARG" ;;
	g) genes="$OPTARG" ;;
	o) out="$OPTARG" ;;
	t) target="$OPTARG" ;;
	N) n_contigs="$OPTARG" ;;
	T) threads="$OPTARG" ;;
	I) min_id="$OPTARG" ;;
	L) min_ovl="$OPTARG" ;;
	*)
		echo "Bad option"
		exit 1
		;;
	esac
done

[[ -z "$reads" ]] && {
	echo "ERROR: -r reads required"
	exit 1
}
mkdir -p "$out"

# 1) Optional recruitment of mt-like reads
inreads="$reads"
if [[ -n "${genes}" ]]; then
	echo "[INFO] Recruiting mt-like reads..."
	minimap2 -t "$threads" -x map-ont "$genes" "$reads" -c --cs=long |
		awk 'BEGIN{FS=OFS="\t"} $11>=400 && ($10/$11)>=0.75 {print $1}' |
		sort -u >"$out/mt_like.reads.txt"

	if command -v seqkit >/dev/null 2>&1; then
		seqkit grep -f "$out/mt_like.reads.txt" "$reads" >"$out/recr.fq"
	else
		python3 - "$reads" "$out/mt_like.reads.txt" "$out/recr.fq" <<'PY'
import sys, gzip
reads, keepf, outf = sys.argv[1:]
keep=set(x.strip() for x in open(keepf))
out = gzip.open(outf,'wt') if outf.endswith('.gz') else open(outf,'w')
rec=[]; 
def name(h): 
    h=h.split()[0]; 
    return h[1:] if h[:1] in '@>' else h
def flush():
    if not rec: return
    if name(rec[0]) in keep: out.write(''.join(rec))
with (gzip.open(reads,'rt') if reads.endswith('.gz') else open(reads)) as f:
    for ln in f:
        if ln[:1] in '@>':
            flush(); rec=[ln]
        else:
            rec.append(ln)
    flush()
out.close()
PY
	fi
	inreads="$out/recr.fq"
fi

# 2) All-vs-all overlaps
echo "[INFO] Computing overlaps (ava-ont)..."
paf="$out/overlaps.paf"
minimap2 -t "$threads" -x ava-ont "$inreads" "$inreads" >"$paf"

# 3) Greedy ~100 kb constructor
echo "[INFO] Assembling ~${target} bp unitigs (N=${n_contigs})..."
# python3 "$(dirname "$0")/paf_greedy_100k_ont.py" \
bashdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$bashdir/polap-py-paf-greedy-ont100k.py" \
	--reads "$inreads" --paf "$paf" \
	--min-id "$min_id" --min-ovl "$min_ovl" \
	--target "$target" --max-contigs "$n_contigs" \
	--topk 3 --branch_delta 0.10 \
	--out "$out/greedy_100k.fasta" \
	--paths "$out/greedy_paths.tsv" \
	>"$out/greedy_100k.log"

# 4) Optional annotation
if [[ -n "${genes}" ]]; then
	echo "[INFO] Annotating constructed contigs with mt genes..."
	minimap2 -t "$threads" -x map-ont "$genes" "$out/greedy_100k.fasta" -c --cs=long >"$out/mt_annot.paf"
	awk 'BEGIN{FS=OFS="\t"}{cov[$6]+=$11} END{for(k in cov) print k,cov[k]}' "$out/mt_annot.paf" |
		sort -k2,2nr >"$out/mt_annot.coverage.tsv"
fi

echo "[DONE] Output:"
echo "  ${out}/greedy_100k.fasta"
echo "  ${out}/greedy_paths.tsv"
[[ -f "$out/mt_annot.coverage.tsv" ]] && echo "  ${out}/mt_annot.coverage.tsv"
