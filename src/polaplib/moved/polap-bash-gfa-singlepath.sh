#!/usr/bin/env bash
# polap-bash-gfa-singlepath.sh
# Version: v0.1.0
# Goal: From a (normalized) organelle subgraph GFA, emit ONE connected FASTA:
#   1) Try gfatk linear
#   2) If still multiple sequences, auto-pick a single best path (by ec:i) and run gfatk path → fasta
#
# Usage:
#   polap-bash-gfa-singlepath.sh \
#     --gfa mtpt/sample.mito.gfa \
#     --outdir mtpt \
#     --prefix sample.mito \
#     [--default-ovl 0] [--greedy] [--min-len 1000]
#
# Inputs:
#   --gfa : organelle subgraph (e.g., from gfatk extract-mito/chloro), with numeric IDs, ll:f, ec:i, simple CIGARs
# Outputs:
#   outdir/prefix.linear.fa     (if gfatk linear succeeds with 1 record)
#   outdir/prefix.bestpath.txt  (comma path "1+,2-,3+")
#   outdir/prefix.single.fa     (one connected FASTA from best path)
#
set -euo pipefail
IFS=$'\n\t'

gfa="" outdir="out" prefix="organelle"
default_ovl=0 greedy=false minlen=1000

_show() {
	cat <<'USAGE'
Usage:
  polap-bash-gfa-singlepath.sh --gfa in.gfa --outdir out --prefix NAME [--default-ovl 0] [--greedy] [--min-len 1000]
USAGE
}

while (($#)); do
	case "$1" in
	--gfa)
		gfa="$2"
		shift 2
		;;
	--outdir)
		outdir="$2"
		shift 2
		;;
	--prefix)
		prefix="$2"
		shift 2
		;;
	--default-ovl)
		default_ovl="$2"
		shift 2
		;;
	--greedy)
		greedy=true
		shift
		;;
	--min-len)
		minlen="$2"
		shift 2
		;;
	-h | --help)
		_show
		exit 0
		;;
	*)
		echo "ERR: unknown arg $1" >&2
		exit 2
		;;
	esac
done

[[ -s "$gfa" ]] || {
	echo "ERR: --gfa missing/empty: $gfa" >&2
	exit 2
}
mkdir -p "$outdir"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "ERR: missing tool: $1" >&2
	exit 127
}; }
need gfatk
need python3

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pick_py="${script_dir}/scripts/polap_py_gfa_pick_best_path.py"

# 1) Try gfatk linear first
linfa="${outdir}/${prefix}.linear.fa"
if gfatk linear "$gfa" >"$linfa" 2>/dev/null; then
	nrec=$(grep -c '^>' "$linfa" || true)
	if ((nrec == 1)); then
		echo "[OK] gfatk linear produced a single sequence: $linfa"
		exit 0
	else
		echo "[INFO] gfatk linear produced $nrec sequences; will pick a single best path."
	fi
else
	echo "[INFO] gfatk linear failed or not applicable; will pick a single best path."
fi

# 2) Auto-pick one best path by ec:i and build a single path FASTA via gfatk path
pathfile="${outdir}/${prefix}.bestpath.txt"
python3 "$pick_py" --gfa "$gfa" --out-path "$pathfile" \
	--default-ovl "$default_ovl" $($greedy && echo --greedy) --min-len "$minlen"

single="${outdir}/${prefix}.single.fa"
gfatk path "$gfa" "$pathfile" >"${outdir}/${prefix}.single.gfa"
gfatk fasta "${outdir}/${prefix}.single.gfa" >"$single"

# Verify single record
nrec=$(grep -c '^>' "$single" || true)
if ((nrec != 1)); then
	echo "[WARN] expected 1 record in $single but saw $nrec; keep the best we could generate."
else
	echo "[OK] Wrote single connected FASTA: $single"
fi
