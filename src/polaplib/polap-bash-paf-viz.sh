#!/usr/bin/env bash
# polap-bash-paf-viz.sh
# Version: v0.1.0
# Minimal PAF visualization pipeline: extracts core stats and makes plots.
# Usage:
#   bash polap-bash-paf-viz.sh -i in.paf[.gz] -o outdir [--min-alen 50] [--min-ident 0.0] [--prefix NAME]
# Deps: python3, Rscript, ggplot2, data.table (R); gzip support for .gz inputs.

set -euo pipefail
IFS=$'\n\t'

_show_usage() {
	cat <<'USAGE'
polap-bash-paf-viz.sh -i in.paf[.gz] -o outdir [--min-alen 50] [--min-ident 0.0] [--prefix NAME]

Outputs:
  outdir/NAME.core.csv                 # per-alignment rows
  outdir/plots/NAME_hist_alen.png
  outdir/plots/NAME_hist_identity.png
  outdir/plots/NAME_hex_alen_identity.png
  outdir/plots/NAME_hist_mapq.png
  outdir/plots/NAME_ecdf_q_aligned_fraction.png
USAGE
}

in="" out="" min_alen=50 min_ident=0.0 prefix="pafviz"

# parse args
while (($#)); do
	case "$1" in
	-i | --in | --paf)
		in="$2"
		shift 2
		;;
	-o | --outdir)
		out="$2"
		shift 2
		;;
	--min-alen)
		min_alen="$2"
		shift 2
		;;
	--min-ident)
		min_ident="$2"
		shift 2
		;;
	--prefix)
		prefix="$2"
		shift 2
		;;
	-h | --help)
		_show_usage
		exit 0
		;;
	*)
		echo "ERR: Unknown arg: $1" >&2
		_show_usage
		exit 2
		;;
	esac
done

[[ -n "$in" ]] || {
	echo "ERR: -i/--in required" >&2
	exit 2
}
[[ -n "$out" ]] || {
	echo "ERR: -o/--outdir required" >&2
	exit 2
}

mkdir -p "$out/plots"

corecsv="${out}/${prefix}.core.csv"

# 1) Extract core fields → CSV
python3 "$(dirname "$0")/scripts/paf_extract_core_to_csv.py" \
	--paf "$in" --out "$corecsv" \
	--min-alen "$min_alen" --min-ident "$min_ident"

# 2) Plot
Rscript "$(dirname "$0")/scripts/paf_plot_core.R" \
	"$corecsv" "$out/plots" "$prefix"

echo "Done. Plots in: $out/plots" >&2
