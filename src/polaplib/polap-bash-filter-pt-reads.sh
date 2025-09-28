#!/usr/bin/env bash
set -euo pipefail

################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# polap-bash-filter-pt-reads.sh
# Wrapper for filtering plastid-like contigs using gene density heuristics
# Calls: scripts/polap-r-filter-organelle-reads.R
################################################################################

# Default values
TABLE=""
OUTDIR="."
LENGTH=3e+7
DISP_PT=5000
DISP_MT=100000
DIFF=3
SEED=42
PLOT=""
CROP_PDF=false
CONVERT_PNG=false

usage() {
	cat <<EOF
Usage: $(basename "$0") -t <table.tsv> -o <output_dir> [options]

Required:
  -t, --table         Input TSV table (annotation + depth)
  -o, --outdir        Output directory

Options:
  -l, --length        Max cumulative length [default: ${LENGTH}]
  --disp-pt           Max dispersion_PT to retain contig [default: ${DISP_PT}]
  --disp-mt           Max dispersion_MT to retain contig [default: ${DISP_MT}]
  --seed              RNG seed; if <= 0, random from time [default: ${SEED}]
  --crop-pdf          Crop the output PDF with pdfcrop [default: off]
  --convert-png       Convert PDF plot to PNG using 'convert' [default: off]
  -h, --help          Show this help message
EOF
}

# Argument parsing
while [[ $# -gt 0 ]]; do
	case $1 in
	-t | --table)
		TABLE="$2"
		shift
		;;
	-o | --outdir)
		OUTDIR="$2"
		shift
		;;
	-l | --length)
		LENGTH="$2"
		shift
		;;
	--diff)
		DIFF="$2"
		shift
		;;
	--disp-pt)
		DISP_PT="$2"
		shift
		;;
	--disp-mt)
		DISP_MT="$2"
		shift
		;;
	--seed)
		SEED="$2"
		shift
		;;
	--crop-pdf) CROP_PDF=true ;;
	--convert-png) CONVERT_PNG=true ;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown option: $1" >&2
		usage
		exit 1
		;;
	esac
	shift
done

[[ -z "$TABLE" || -z "$OUTDIR" ]] && {
	echo "Error: --table and --outdir are required." >&2
	usage
	exit 1
}

mkdir -p "$OUTDIR"

# File naming
BASENAME=$(basename "$TABLE" .tsv)
OUT_IDS="$OUTDIR/${BASENAME}.pt.id.txt"
OUT_ALL="$OUTDIR/${BASENAME}.pt.filtered.all.txt"
OUT_PDF="$OUTDIR/${BASENAME}.scatter.pdf"
OUT_PNG="${OUT_PDF%.pdf}.png"

# Run the R script
Rscript --vanilla "$(dirname "$0")/polap-r-filter-organelle-reads.R" \
	--table "$TABLE" \
	--length "$LENGTH" \
	--diff "$DIFF" \
	--output "$OUT_IDS" \
	--output-all "$OUT_ALL" \
	--min-dispersion-pt "$DISP_PT" \
	--min-dispersion-mt "$DISP_MT" \
	--rng-seed "$SEED" \
	--plot "$OUT_PDF"

# Optional cropping and conversion
if $CROP_PDF; then
	echo "[INFO] Cropping PDF plot..."
	pdfcrop "$OUT_PDF" "$OUT_PDF" >/dev/null
fi

if $CONVERT_PNG; then
	echo "[INFO] Converting PDF to PNG..."
	convert -density 300 "$OUT_PDF" -quality 100 "$OUT_PNG"
fi

echo "[DONE] Selected contig IDs: $OUT_IDS"
