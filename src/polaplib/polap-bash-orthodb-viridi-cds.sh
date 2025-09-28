#!/usr/bin/env bash
# polap-bash-orthodb-viridi-cds.sh
# Pull CDS sequences for Viridiplantae (taxon 33090) from OrthoDB v12 flat files.
# Inputs:
#   -o, --outdir DIR          Output directory [default: ./orthodb_viridi]
#   --clade  TAXID            NCBI taxon for clade (default 33090 = Viridiplantae)
#   --cds   FILE              Path to odb12v1_cds_fasta.gz (53.6 GB). If absent, we can wget it.
#   --wget-cds                Also download the big CDS file (requires ~54 GB + network)
#   --threads INT             pigz threads for (re)compression [4]
# Usage:
#   bash polap-bash-orthodb-viridi-cds.sh --cds /path/odb12v1_cds_fasta.gz
set -euo pipefail

OUTDIR="./orthodb_viridi"
CLADE="33090"
CDS=""
WGET_CDS=0
T=4
while [[ $# -gt 0 ]]; do
	case "$1" in
	-o | --outdir)
		OUTDIR="$2"
		shift 2
		;;
	--clade)
		CLADE="$2"
		shift 2
		;;
	--cds)
		CDS="$2"
		shift 2
		;;
	--wget-cds)
		WGET_CDS=1
		shift 1
		;;
	--threads)
		T="$2"
		shift 2
		;;
	-h | --help)
		grep -E '^# ' "$0" | sed 's/^# //'
		exit 0
		;;
	*)
		echo "Unknown arg: $1" >&2
		exit 1
		;;
	esac
done

mkdir -p "$OUTDIR"/{raw,tmp,logs}

# 1) Small TSVs from OrthoDB v12 "current" flat-file mirror
# They are small enough to (re)download each run safely.
base="https://data.orthodb.org/current/download"
declare -A files=(
	[levels]="odb12v1_levels.tab.gz"
	[species]="odb12v1_species.tab.gz"
	[l2s]="odb12v1_level2species.tab.gz"
	[genes]="odb12v1_genes.tab.gz"
	[og2genes]="odb12v1_OG2genes.tab.gz"
)

for k in levels species l2s genes; do
	f="${files[$k]}"
	if [[ ! -s "$OUTDIR/raw/$f" ]]; then
		echo "[info] downloading $f ..."
		wget -q -O "$OUTDIR/raw/$f" "$base/$f"
	fi
done

# Optionally download the big CDS FASTA (53.6 GB)
if [[ $WGET_CDS -eq 1 ]]; then
	f="odb12v1_cds_fasta.gz"
	if [[ ! -s "$OUTDIR/raw/$f" ]]; then
		echo "[info] downloading $f (this is ~54 GB) ..."
		wget -O "$OUTDIR/raw/$f" "$base/$f"
	fi
	CDS="$OUTDIR/raw/$f"
fi

if [[ -z "${CDS:-}" ]]; then
	echo "[info] Using user-supplied --cds: $CDS"
	[[ -s "$CDS" ]] || {
		echo "[err] --cds file not found: $CDS" >&2
		exit 2
	}
fi

# 2) Build Viridiplantae organism list (OrthoDB organism ids, e.g. 3702_0)
#    level2species: columns = <level_taxid>\t<organism_id>
gunzip -c "$OUTDIR/raw/${files[l2s]}" |
	awk -v clade="$CLADE" -F'\t' '$1==clade{print $2}' |
	sort -u >"$OUTDIR/tmp/viridi.organism_ids.txt"

# 3) Map organism -> gene_ids using genes.tab
#    genes.tab columns (v10 README; v12 similar): includes organism id and gene id.
gunzip -c "$OUTDIR/raw/${files[genes]}" |
	awk -F'\t' '{
      # Heuristic column mapping:
      # Expect columns like: ncbi_taxid, organism_id, gene_id (public), lengths, etc.
      # OrthoDB changes column count across releases; we detect organism_id as field with "_"
      oid=""; gid="";
      for(i=1;i<=NF;i++){
        if(oid=="" && $i ~ /_[0-9]+$/) oid=$i;
        if(gid=="" && $i ~ /^[0-9]+_[0-9]+:[0-9A-Za-z._-]+$/) gid=$i; # e.g. 3702_0:00ABCd
      }
      if(oid!="" && gid!="") print oid"\t"gid;
    }' |
	LC_ALL=C sort -S 1G -T "$OUTDIR/tmp" >"$OUTDIR/tmp/genes_organism_gene.tsv"

# 3b) Join to Viridiplantae organism ids → Viridi gene ids
LC_ALL=C join -t $'\t' -o 2.2 \
	<(LC_ALL=C sort "$OUTDIR/tmp/viridi.organism_ids.txt") \
	<(cut -f1,2 "$OUTDIR/tmp/genes_organism_gene.tsv" | LC_ALL=C sort) |
	sort -u >"$OUTDIR/viridiplantae_gene_ids.txt"

echo "[info] Viridiplantae gene IDs: $(wc -l <"$OUTDIR/viridiplantae_gene_ids.txt")"

# 4) Stream-filter the big CDS FASTA by OrthoDB gene IDs in FASTA headers
#    FASTA headers begin with the OrthoDB gene id (per OrthoDB docs).
python3 "$(dirname "$0")/scripts/fasta_select_by_geneids.py" \
	--ids "$OUTDIR/viridiplantae_gene_ids.txt" \
	--fasta "$CDS" \
	--out "$OUTDIR/viridiplantae_cds.fa.gz" \
	--threads "$T"

echo "[done] CDS written: $OUTDIR/viridiplantae_cds.fa.gz"
