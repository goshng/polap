#!/usr/bin/env bash
# polap-bash-orthodb-viridi-busco-cds.sh
# Build Viridiplantae BUSCO CDS FASTA from OrthoDB v12 flat files.
# Inputs:
#   --busco-links FILE   BUSCO links_to_ODB12.txt (1st col = OG id like 10041at33090)
#   --cds FILE           Path to OrthoDB v12 CDS dump: odb12v1_cds_fasta.gz (~54 GB)
#   -o DIR               Output dir [./orthodb_viridi_busco]
#   --threads INT        Threads for gzip [4]
#   --clade TAXID        NCBI taxon id (Viridiplantae=33090) [33090]
#   --download           Download needed TSVs (levels/species/level2species/OG2genes/genes)
# Usage example:
#   bash polap-bash-orthodb-viridi-busco-cds.sh \
#     --busco-links links_to_ODB12.txt \
#     --cds /mnt/db/orthodb/odb12v1_cds_fasta.gz \
#     -o viridi_busco
set -euo pipefail

OUT=./orthodb_viridi_busco
T=4
CLADE=33090
DL=0
BUSCO_LINKS=""
CDS=""
while [[ $# -gt 0 ]]; do
	case "$1" in
	--busco-links)
		BUSCO_LINKS="$2"
		shift 2
		;;
	--cds)
		CDS="$2"
		shift 2
		;;
	-o | --outdir)
		OUT="$2"
		shift 2
		;;
	--threads)
		T="$2"
		shift 2
		;;
	--clade)
		CLADE="$2"
		shift 2
		;;
	--download)
		DL=1
		shift 1
		;;
	-h | --help)
		grep -E '^# ' "$0" | sed 's/^# //'
		exit 0
		;;
	*)
		echo "[ERR] Unknown arg: $1" >&2
		exit 2
		;;
	esac
done

[[ -s "$BUSCO_LINKS" ]] || {
	echo "[ERR] --busco-links required" >&2
	exit 2
}

[[ -s "$CDS" ]] || {
	echo "[ERR] --cds odb12v1_cds_fasta.gz required" >&2
	exit 2
}

mkdir -p "$OUT"/{raw,tmp,logs,scripts}

BASE="https://data.orthodb.org/current/download"
declare -A F=(
	[levels]="odb12v1_levels.tab.gz"
	[species]="odb12v1_species.tab.gz"
	[l2s]="odb12v1_level2species.tab.gz"
	[og2genes]="odb12v1_OG2genes.tab.gz"
	[genes]="odb12v1_genes.tab.gz"
)

if [[ $DL -eq 1 ]]; then
	for k in levels species l2s og2genes genes; do
		f="${F[$k]}"
		if [[ -s "$OUT/raw/$f" ]]; then
			continue
		fi
		echo "[info] downloading $f ..."
		echo wget -q -O "$OUT/raw/$f" "$BASE/$f" >>1.sh
		# wget -q -O "$OUT/raw/$f" "$BASE/$f"
	done
else
	for k in levels species l2s og2genes genes; do
		f="$OUT/raw/${F[$k]}"
		[[ -s "$f" ]] || {
			echo "[ERR] Missing $f (re-run with --download)" >&2
			exit 3
		}
	done
fi
exit

# 1) Parse BUSCO OG IDs (1st column of links_to_ODB12.txt)
awk 'NF{print $1}' "$BUSCO_LINKS" | sort -u >"$OUT/tmp/busco_og_ids.txt"
echo "[info] BUSCO OG count: $(wc -l <"$OUT/tmp/busco_og_ids.txt")"

# 2) Build Viridiplantae organism id list (from level2species)
gunzip -c "$OUT/raw/${F[l2s]}" | awk -v c="$CLADE" -F'\t' '$1==c{print $2}' |
	sort -u >"$OUT/tmp/viridi.organism_ids.txt"
echo "[info] Viridi organism ids: $(wc -l <"$OUT/tmp/viridi.organism_ids.txt")"

# 3) Restrict OG->gene map to BUSCO OGs
#    OG2genes columns: <OGid>\t<gene_id>
#    Keep rows with OG in our list, then restrict to Viridi by organism prefix of gene_id.
#    OrthoDB gene_id format begins with organism id like 3702_0:XXXX
gunzip -c "$OUT/raw/${F[og2genes]}" |
	awk 'NF==2' |
	LC_ALL=C join -t $'\t' -1 1 -2 1 \
		<(sort "$OUT/tmp/busco_og_ids.txt") - |
	cut -f2 >"$OUT/tmp/busco_og_genes.all.txt"

# 4) Filter to Viridi gene_ids by organism prefix match
awk 'NR==FNR{v[$1]; next} {split($1,a,":"); if(a[1] in v) print $1}' \
	"$OUT/tmp/viridi.organism_ids.txt" \
	"$OUT/tmp/busco_og_genes.all.txt" |
	sort -u >"$OUT/viridiplantae_busco_gene_ids.txt"

echo "[info] Viridi BUSCO gene ids: $(wc -l <"$OUT/viridiplantae_busco_gene_ids.txt")"

# 5) Stream-filter the CDS FASTA by gene IDs
#    FASTA headers start with OrthoDB gene id (per OrthoDB docs)
if [[ ! -s "$OUT/scripts/fasta_select_by_geneids.py" ]]; then
	cat >"$OUT/scripts/fasta_select_by_geneids.py" <<'PY'
#!/usr/bin/env python3
import gzip, argparse
def open_any(p, m='rt'):
    return gzip.open(p, m) if p.endswith('.gz') else open(p, m)
ap=argparse.ArgumentParser()
ap.add_argument("--ids", required=True)
ap.add_argument("--fasta", required=True)
ap.add_argument("--out", required=True)
args=ap.parse_args()
ids=set(x.strip() for x in open(args.ids) if x.strip())
with open_any(args.fasta,'rt') as fi, gzip.open(args.out,'wt') as fo:
    keep=False
    for line in fi:
        if line.startswith('>'):
            hdr=line[1:].strip().split()[0]
            keep = hdr in ids
        if keep: fo.write(line)
PY
	chmod +x "$OUT/scripts/fasta_select_by_geneids.py"
fi

"$OUT/scripts/fasta_select_by_geneids.py" \
	--ids "$OUT/viridiplantae_busco_gene_ids.txt" \
	--fasta "$CDS" \
	--out "$OUT/viridiplantae_busco_cds.fa.gz"

echo "[done] Wrote $OUT/viridiplantae_busco_cds.fa.gz"
