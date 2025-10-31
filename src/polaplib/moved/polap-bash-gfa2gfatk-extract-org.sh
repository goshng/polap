#!/usr/bin/env bash
# polap-bash-gfa-extract-with-gfatk.sh
# Version: v0.4.0
# Normalize Flye GFA → numeric IDs, simple overlaps, add ll:f, AND add ec:i, then run gfatk.

set -euo pipefail
IFS=$'\n\t'

_show() {
	cat <<'USAGE'
Usage:
  polap-bash-gfa-extract-with-gfatk.sh \
    --gfa assembly_graph.gfa \
    --outdir out \
    [--target mito|pltd|both] \
    [--linearize] \
    [--prefix NAME] \
    [--ec-mode const|min|max|mean] \
    [--ec-const 1] \
    [--ec-scale 1.0] \
    [--ec-round round|ceil|floor] \
    [--oatk-path "n1+,n2-,n3+"] [--oatk-path-file file] [--oatk-map-out path.txt]

Notes:
  - Renumbers S IDs to 1..N and rewrites L endpoints.
  - L overlap "*" → "0M".
  - S adds ll:f from dp:i/f if missing (defaults ll:f:1.0).
  - L adds ec:i from node coverages via --ec-* (default: const 1).
  - Then runs: gfatk extract-{chloro,mito} + fasta (+ linear if asked).
USAGE
}

gfa="" outdir="gfatk_out" target="both" linearize=false prefix="sample"
oatk_path="" oatk_path_file="" oatk_map_out=""
ec_mode="const" ec_const="1" ec_scale="1.0" ec_round="round"

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
	--oatk-map-out)
		oatk_map_out="$2"
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
	echo "ERR: --gfa missing or empty: $gfa" >&2
	exit 2
}
mkdir -p "$outdir"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "ERR: missing: $1" >&2
	exit 127
}; }
need awk
need gfatk

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_map_awk="${script_dir}/scripts/gfa_build_idmap.awk"
apply_map_awk="${script_dir}/scripts/gfa_apply_idmap_fix_and_ec.awk"
map_path_py="${script_dir}/scripts/map_path_ids.py"

# 1) Build ID map old->new
map_tsv="${outdir}/${prefix}.idmap.tsv"
awk -f "$build_map_awk" "$gfa" >"$map_tsv"

# 2) Apply map, fix overlaps, add ll:f, add ec:i
norm_gfa="${outdir}/${prefix}.norm.gfa"
awk -v EC_MODE="$ec_mode" -v EC_CONST="$ec_const" -v EC_SCALE="$ec_scale" -v EC_ROUND="$ec_round" \
	-f "$apply_map_awk" "$map_tsv" "$gfa" >"$norm_gfa"

# 3) gfatk extract
run_one() {
	local which="$1"
	local sub outg segf linf
	case "$which" in
	pltd) sub="extract-chloro" ;;
	mito) sub="extract-mito" ;;
	*)
		echo "ERR: which must be pltd|mito"
		exit 2
		;;
	esac
	outg="${outdir}/${prefix}.${which}.gfa"
	segf="${outdir}/${prefix}.${which}.segs.fa"
	linf="${outdir}/${prefix}.${which}.linear.fa"

	echo "[*] gfatk ${sub} …"
	gfatk "$sub" "$norm_gfa" >"$outg"

	echo "[*] gfatk fasta …"
	gfatk fasta "$outg" >"$segf"

	if $linearize; then
		echo "[*] gfatk linear …"
		if ! gfatk linear "$outg" >"$linf"; then
			echo "[!] linear failed for $which (segments FASTA kept)" >&2
		fi
	fi
	echo "[OK] $which → $outg ; $segf ${linearize:+; $linf}"
}

case "$target" in
mito) run_one mito ;;
pltd) run_one pltd ;;
both)
	run_one pltd
	run_one mito
	;;
*)
	echo "ERR: --target mito|pltd|both"
	exit 2
	;;
esac

# 4) Optional: map an original-ID path for OATK path_to_fasta
if [[ -n "$oatk_path_file" ]]; then
	oatk_path="$(tr -d '\r\n\t ' <"$oatk_path_file")"
fi
if [[ -n "$oatk_path" ]]; then
	[[ -n "$oatk_map_out" ]] || oatk_map_out="${outdir}/${prefix}.mapped_path.txt"
	python3 "$map_path_py" --map "$map_tsv" --path "$oatk_path" >"$oatk_map_out"
	echo "[OK] Mapped OATK path → $oatk_map_out"
	echo "Use: path_to_fasta ${norm_gfa} \"\$(cat $oatk_map_out)\" > ${outdir}/${prefix}.path.fa"
fi

echo "[DONE] ID map: $map_tsv"
echo "[DONE] Normalized GFA: $norm_gfa"
