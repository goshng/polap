#!/usr/bin/env bash
set -euo pipefail

# polap-bash-recruit-competitive.sh  v0.2.1
# Competitive recruit of reads against: mito seeds ∪ pt (2 isomers) ∪ optional nuclear decoy.
# Picks the single best target per read, then filters by identity, query coverage, and target span.

show_help() {
	cat <<'EOF'
Usage:
  polap-bash-recruit-competitive.sh \
    --threads INT \
    --reads READS.fq[.gz] \
    --mito mito.fa \
    [--pt1 plastid_form1.fa] [--pt2 plastid_form2.fa] \
    [--nuc nuclear_decoy.fa] \
    --id-min FLOAT --qcov-min FLOAT --tspan-min INT \
    -o OUTDIR

Notes
- Produces: OUTDIR/competitive.paf (raw), OUTDIR/competitive.best.tsv (best per read),
            OUTDIR/mapped.ids, OUTDIR/mapped.fastq.gz, OUTDIR/unmapped.fastq.gz
- Requires: minimap2, seqtk, python3
EOF
}

# ---------- defaults
THREADS=8
READS=""
MITO=""
PT1=""
PT2=""
NUC=""
ID_MIN=0.80
QCOV_MIN=0.50
TSPAN_MIN=1000
OUTDIR=""

# ---------- parse CLI
while [[ $# -gt 0 ]]; do
	case "$1" in
	--threads)
		THREADS="$2"
		shift 2
		;;
	--reads)
		READS="$2"
		shift 2
		;;
	--mito)
		MITO="$2"
		shift 2
		;;
	--pt1)
		PT1="$2"
		shift 2
		;;
	--pt2)
		PT2="$2"
		shift 2
		;;
	--nuc)
		NUC="$2"
		shift 2
		;;
	--id-min)
		ID_MIN="$2"
		shift 2
		;;
	--qcov-min)
		QCOV_MIN="$2"
		shift 2
		;;
	--tspan-min)
		TSPAN_MIN="$2"
		shift 2
		;;
	-o | --out)
		OUTDIR="$2"
		shift 2
		;;
	-h | --help)
		show_help
		exit 0
		;;
	*)
		echo "[ERROR] unknown option: $1" >&2
		show_help
		exit 2
		;;
	esac
done

# ---------- sanity
[[ -n "$READS" && -s "$READS" ]] || {
	echo "[ERROR] --reads missing" >&2
	exit 2
}
[[ -n "$MITO" && -s "$MITO" ]] || {
	echo "[ERROR] --mito missing" >&2
	exit 2
}
[[ -n "$OUTDIR" ]] || {
	echo "[ERROR] -o/--out missing" >&2
	exit 2
}
command -v minimap2 >/dev/null || {
	echo "[ERROR] minimap2 not found" >&2
	exit 2
}
command -v seqtk >/dev/null || {
	echo "[ERROR] seqtk not found" >&2
	exit 2
}
command -v python3 >/dev/null || {
	echo "[ERROR] python3 not found" >&2
	exit 2
}

mkdir -p "$OUTDIR"

# ---------- build panel FASTA
PANEL="$OUTDIR/panel.fasta"
: >"$PANEL"
cat "$MITO" >>"$PANEL"
[[ -n "$PT1" && -s "$PT1" ]] && cat "$PT1" >>"$PANEL"
[[ -n "$PT2" && -s "$PT2" ]] && cat "$PT2" >>"$PANEL"
[[ -n "$NUC" && -s "$NUC" ]] && cat "$NUC" >>"$PANEL"

# ---------- map all reads competitively (keep sec hits so best can be chosen)
PAF_RAW="$OUTDIR/competitive.paf"
minimap2 -t "$THREADS" -x map-ont -k14 --secondary=yes -N 50 -p 0.5 -c \
	"$PANEL" "$READS" >"$PAF_RAW"

# ---------- helper python (prefix rule)
PY="$OUTDIR/polap-bash-recruit-competitive.filter_best.py"
cat >"$PY" <<'PY'
#!/usr/bin/env python3
import sys, math, gzip

# Inputs
paf = sys.argv[1]
id_min = float(sys.argv[2])
qcov_min = float(sys.argv[3])
tspan_min = int(sys.argv[4])
out_tsv = sys.argv[5]

# PAF cols (1-based doc): 1 qname 2 qlen 3 qstart 4 qend 6 tname 8 tstart 9 tend 10 nmatch 11 alnlen 12 mapq
best = {}  # qname -> (score_tuple, fields)
with open(paf, 'r') as f:
    for line in f:
        if not line.strip() or line.startswith('#'): continue
        a = line.rstrip('\n').split('\t')
        if len(a) < 12: continue
        qname = a[0]; qlen = int(a[1]); qstart = int(a[2]); qend = int(a[3])
        tname = a[5]; tstart = int(a[7]); tend = int(a[8])
        nmatch = int(a[9]); alnlen = int(a[10]); mapq = int(a[11])
        if qlen <= 0 or alnlen <= 0: continue

        ident = nmatch / float(alnlen)
        qcov  = (qend - qstart) / float(qlen)
        tspan = abs(tend - tstart)

        # keep everything for best-pick; threshold later
        # score: prefer higher ident, then qcov, then tspan, then mapq
        score = (ident, qcov, tspan, mapq)
        prev = best.get(qname)
        if (prev is None) or (score > prev[0]):
            best[qname] = (score, (qname, qlen, qstart, qend, tname, tstart, tend, nmatch, alnlen, mapq, ident, qcov, tspan))

# write best per read, then apply thresholds
with open(out_tsv, 'w') as w:
    w.write("qname\tqlen\tqstart\tqend\ttname\ttstart\ttend\tnmatch\talnlen\tmapq\tident\tqcov\ttspan\tpass\n")
    for qname,(score,rec) in best.items():
        ident = rec[10]; qcov = rec[11]; tspan = rec[12]
        ok = (ident >= id_min) and (qcov >= qcov_min) and (tspan >= tspan_min)
        w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6f}\t{:.6f}\t{}\t{}\n".format(
            rec[0], rec[1], rec[2], rec[3], rec[4], rec[5], rec[6],
            rec[7], rec[8], rec[9], ident, qcov, tspan, "PASS" if ok else "FAIL"
        ))
PY

chmod +x "$PY"

BEST_TSV="$OUTDIR/competitive.best.tsv"
python3 "$PY" "$PAF_RAW" "$ID_MIN" "$QCOV_MIN" "$TSPAN_MIN" "$BEST_TSV"

# ---------- extract mapped/unmapped IDs and FASTQ
MAPPED_IDS="$OUTDIR/mapped.ids"
UNMAPPED_IDS="$OUTDIR/unmapped.ids"
awk -F'\t' 'NR>1 && $14=="PASS"{print $1}' "$BEST_TSV" | sort -u >"$MAPPED_IDS"

# all read IDs (ID only, no desc)
ALL_IDS="$OUTDIR/all.ids"
if [[ "$READS" =~ \.gz$ ]]; then
	zcat -- "$READS" | awk 'NR%4==1{split($0,a,/ /); print substr(a[1],2)}' | sort -u >"$ALL_IDS"
else
	awk 'NR%4==1{split($0,a,/ /); print substr(a[1],2)}' "$READS" | sort -u >"$ALL_IDS"
fi
comm -23 "$ALL_IDS" "$MAPPED_IDS" >"$UNMAPPED_IDS"

seqtk subseq "$READS" "$MAPPED_IDS" | gzip -c >"$OUTDIR/mapped.fastq.gz"
seqtk subseq "$READS" "$UNMAPPED_IDS" | gzip -c >"$OUTDIR/unmapped.fastq.gz"

echo "[OK] competitive recruit:"
echo "  raw PAF:      $PAF_RAW"
echo "  best table:   $BEST_TSV"
echo "  mapped ids:   $MAPPED_IDS"
echo "  mapped fastq: $OUTDIR/mapped.fastq.gz"
echo "  unmapped:     $OUTDIR/unmapped.fastq.gz"
