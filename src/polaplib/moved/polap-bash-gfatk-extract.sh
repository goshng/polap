#!/usr/bin/env bash
# polap-bash-gfatk-extract.sh
# Version: v0.5.0
# Normalize → (optional ec-first) → gfatk extract-{chloro,mito} → fasta → (try linear) → best single path
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'USAGE'
Usage:
  polap-bash-gfatk-extract.sh --gfa assembly_graph.gfa --outdir out --prefix SAMPLE \
    [--target mito|pltd|both] [--ec-mode mean|min|max|const] [--ec-const 1] \
    [--ec-scale 1.0] [--ec-round round|ceil|floor] [--ec-first] [--linearize] [--greedy]
USAGE
}

gfa="" outdir="gfatk_out" prefix="sample" target="both"
ec_mode="mean" ec_const="1" ec_scale="1.0" ec_round="round"
ec_first=false linearize=false greedy=false

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
	--target)
		target="$2"
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
	--ec-first)
		ec_first=true
		shift
		;;
	--linearize)
		linearize=true
		shift
		;;
	--greedy)
		greedy=true
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERR: unknown arg $1" >&2
		exit 2
		;;
	esac
done

[[ -s "$gfa" ]] || {
	echo "ERR: --gfa missing/empty" >&2
	exit 2
}
mkdir -p "$outdir"

command -v python3 >/dev/null || {
	echo "ERR: need python3" >&2
	exit 127
}
command -v gfatk >/dev/null || {
	echo "ERR: need gfatk in PATH" >&2
	exit 127
}

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
norm_sh="${here}/polap-bash-gfa-normalize.sh"
ecfirst_sh="${here}/polap-bash-gfa-ecfirst.sh"
pick_py="${here}/scripts/polap_py_gfa_pick_best_path.py"

norm="${outdir}/${prefix}.norm.gfa"
idmap="${outdir}/${prefix}.idmap.tsv"
bash "$norm_sh" --gfa "$gfa" --out "$norm" --idmap "$idmap" \
	--ec-mode "$ec_mode" --ec-const "$ec_const" --ec-scale "$ec_scale" --ec-round "$ec_round"

if $ec_first; then
	norm2="${outdir}/${prefix}.norm.ecfirst.gfa"
	bash "$ecfirst_sh" "$norm" "$norm2" --drop-other-tags
	norm="$norm2"
fi

extract_one() {
	local which="$1"
	local sub outg segf linf single path
	case "$which" in
	pltd) sub="extract-chloro" ;;
	mito) sub="extract-mito" ;;
	*)
		echo "ERR: target must be pltd|mito" >&2
		exit 2
		;;
	esac
	outg="${outdir}/${prefix}.${which}.gfa"
	segf="${outdir}/${prefix}.${which}.segs.fa"
	linf="${outdir}/${prefix}.${which}.linear.fa"
	single="${outdir}/${prefix}.${which}.single.fa"
	path="${outdir}/${prefix}.${which}.bestpath.txt"

	gfatk "$sub" "$norm" >"$outg"
	gfatk fasta "$outg" >"$segf"

	if $linearize; then
		if gfatk linear "$outg" >"$linf" 2>/dev/null && [[ $(grep -c '^>' "$linf" || true) -eq 1 ]]; then
			cp "$linf" "$single"
			echo "[OK] $which single via linear: $single"
			return
		fi
	fi

	python3 "$pick_py" --gfa "$outg" --out-path "$path" --min-len 1 $($greedy && echo --greedy)
	gfatk path "$outg" "$path" >"${outg%.gfa}.single.gfa"
	gfatk fasta "${outg%.gfa}.single.gfa" >"$single"
	echo "[OK] $which single via best path: $single"
}

case "$target" in
pltd) extract_one pltd ;;
mito) extract_one mito ;;
both)
	extract_one pltd
	extract_one mito
	;;
*)
	echo "ERR: --target mito|pltd|both" >&2
	exit 2
	;;
esac
