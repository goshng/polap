#!/usr/bin/env bash
# polap-bash-gfa-extract-with-gfatk.sh
# Version: v0.6.0
# Normalize Flye GFA (IDs→1..N, ensure H header, 0M/clean CIGARs, ll:f on S, ec:i on L,
# and add missing reciprocal L links), then gfatk extract-{chloro,mito} (+fasta,+linear).

set -euo pipefail
IFS=$'\n\t'

gfa=""
outdir="gfatk_out"
target="both"
linearize=false
prefix="sample"
ec_mode="mean"
ec_const="1"
ec_scale="1.0"
ec_round="round"
oatk_path=""
oatk_path_file=""

usage() {
	grep -E '^# ' "$0" | sed 's/^# //'
	exit 2
}
[[ $# -eq 0 ]] && usage

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
	--target)
		target="$2"
		shift 2
		;;
	--linearize)
		linearize=true
		shift
		;;
	--prefix)
		prefix="$2"
		shift 2
		;;
	--ec-mode)
		ec_mode="$2"
		shift 2
		;;
	--ec-const)
		ec_const="$2"
		shift 2
		;;
	--ec-scale)
		ec_scale="$2"
		shift 2
		;;
	--ec-round)
		ec_round="$2"
		shift 2
		;;
	--oatk-path)
		oatk_path="$2"
		shift 2
		;;
	--oatk-path-file)
		oatk_path_file="$2"
		shift 2
		;;
	-h | --help) usage ;;
	*)
		echo "ERR: unknown arg $1" >&2
		exit 2
		;;
	esac
done

[[ -s "$gfa" ]] || {
	echo "ERR: --gfa missing: $gfa" >&2
	exit 2
}
mkdir -p "$outdir"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "ERR: missing tool: $1" >&2
	exit 127
}; }
need python3
need gfatk

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
norm_py="${here}/scripts/polap_py_gfa_normalize_ec_bidirectional.py"
check_py="${here}/scripts/polap_py_gfa_check_ec.py"
map_py="${here}/scripts/polap_py_map_path_ids.py"

norm="${outdir}/${prefix}.norm.gfa"
idmap="${outdir}/${prefix}.idmap.tsv"

python3 "$norm_py" \
	--gfa "$gfa" \
	--out "$norm" \
	--idmap "$idmap" \
	--ec-mode "$ec_mode" \
	--ec-const "$ec_const" \
	--ec-scale "$ec_scale" \
	--ec-round "$ec_round"

python3 "$check_py" --gfa "$norm" # verifies ec:i on every L

run_one() {
	local which="$1" sub outg segf linf
	case "$which" in
	pltd) sub="extract-chloro" ;;
	mito) sub="extract-mito" ;;
	*)
		echo "ERR: which must be pltd|mito" >&2
		exit 2
		;;
	esac
	outg="${outdir}/${prefix}.${which}.gfa"
	segf="${outdir}/${prefix}.${which}.segs.fa"
	linf="${outdir}/${prefix}.${which}.linear.fa"

	gfatk "$sub" "$norm" >"$outg"
	gfatk fasta "$outg" >"$segf"
	$linearize && gfatk linear "$outg" >"$linf" || true
	echo "[OK] $which -> $outg ; $segf ${linearize:+; $linf}"
}

case "$target" in
mito) run_one mito ;;
pltd) run_one pltd ;;
both)
	run_one pltd
	run_one mito
	;;
*)
	echo "ERR: --target mito|pltd|both" >&2
	exit 2
	;;
esac

if [[ -n "$oatk_path_file" ]]; then oatk_path="$(tr -d '\r\n\t ' <"$oatk_path_file")"; fi
if [[ -n "$oatk_path" ]]; then
	mapped="${outdir}/${prefix}.mapped_path.txt"
	python3 "$map_py" --map "$idmap" --path "$oatk_path" >"$mapped"
	echo "[OK] mapped path -> $mapped"
	echo "Tip: path_to_fasta $norm \"\$(cat $mapped)\" > ${outdir}/${prefix}.path.fa" # OATK
fi

echo "[DONE] norm=$norm ; idmap=$idmap"
# Docs: required tags & CIGAR rule: gfatk README; GFA1 format.  (see Sources)
