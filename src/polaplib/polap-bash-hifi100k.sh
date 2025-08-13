#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./hifi100k.sh -r reads.fq.gz -o outdir [-g mt_genes.fa] [-t 100000] [-N 5]
#
# Requires: minimap2, gzip (or zcat), awk, python3

reads=""
genes=""
out="hifi100k_out"
target=100000 # ~100 kb
n_contigs=5   # how many ~100 kb constructs to attempt
threads=${THREADS:-8}
min_id=0.985 # HiFi-friendly identity filter
min_ovl=2000 # minimum bp overlap for extension

while getopts "r:g:o:t:N:T:" opt; do
	case $opt in
	r) reads="$OPTARG" ;;
	g) genes="$OPTARG" ;;
	o) out="$OPTARG" ;;
	t) target="$OPTARG" ;;
	N) n_contigs="$OPTARG" ;;
	T) threads="$OPTARG" ;;
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

# 1) Optional recruitment of likely-mt reads
recr="$out/recr.fq"
if [[ -n "${genes}" ]]; then
	echo "[INFO] Recruiting mt-like reads..."
	# HiFi long-read to gene sequences (nucleotide) mapping
	# Keep reads with any hit ≥500 bp and ≥85% identity (adjust as needed)
	minimap2 -t "$threads" -x map-hifi "$genes" "$reads" -c --cs=long |
		awk 'BEGIN{FS=OFS="\t"} $11>=500 && ($10/$11)>=0.85 {print $1}' |
		sort -u >"$out/mt_like.reads.txt"

	# Extract recruited reads (seqkit if available; otherwise awk fallback)
	if command -v seqkit >/dev/null 2>&1; then
		seqkit grep -f "$out/mt_like.reads.txt" "$reads" >"$recr"
	else
		# very simple extractor (FASTQ/FASTA)
		python3 - "$reads" "$out/mt_like.reads.txt" "$recr" <<'PY'
import sys, gzip
reads, keep, out = sys.argv[1:]
keep = set(x.strip() for x in open(keep))
op = gzip.open(out, 'wt') if out.endswith('.gz') else open(out,'w')
def put(rec): op.write(''.join(rec))
rec=[]; hdr=None
def flush():
    if not rec: return
    h = rec[0].split()[0]
    h = h[1:] if h.startswith('>') or h.startswith('@') else h
    if h in keep: put(rec)
with (gzip.open(reads,'rt') if reads.endswith('.gz') else open(reads)) as f:
    for line in f:
        if line[:1] in '>@' and rec:
            flush(); rec=[line]
        else:
            rec.append(line)
    flush()
op.close()
PY
	fi
	inreads="$recr"
else
	inreads="$reads"
fi

# 2) All-vs-all overlaps (PAF)
echo "[INFO] Computing overlaps..."
paf="$out/overlaps.paf"
minimap2 -t "$threads" -x ava-pb "$inreads" "$inreads" >"$paf"

# 3) Greedy grow to ~100 kb
echo "[INFO] Assembling ~${target} bp unitigs (N=${n_contigs})..."
# python3 "$(dirname "$0")/paf_greedy_100k.py" \
bashdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$bashdir/polap-py-paf-greedy-hifi100k.py" \
	--reads "$inreads" --paf "$paf" \
	--min-id "$min_id" --min-ovl "$min_ovl" \
	--target "$target" --max-contigs "$n_contigs" \
	--out "$out/greedy_100k.fasta" \
	>"$out/greedy_100k.log"

# 4) Optional annotation of the constructed contigs against mt genes
if [[ -n "${genes}" ]]; then
	echo "[INFO] Annotating constructed contigs with mt genes..."
	minimap2 -t "$threads" -x map-hifi "$genes" "$out/greedy_100k.fasta" -c --cs=long >"$out/mt_annot.paf"
	# A tiny summary: per contig total aligned bp
	awk 'BEGIN{FS=OFS="\t"}{len=$11; cov[$6]+=len} END{for(k in cov) print k,cov[k]}' "$out/mt_annot.paf" |
		sort -k2,2nr >"$out/mt_annot.coverage.tsv"
fi

echo "[DONE] Output:"
echo "  ${out}/greedy_100k.fasta"
[[ -f "$out/mt_annot.coverage.tsv" ]] && echo "  ${out}/mt_annot.coverage.tsv"
