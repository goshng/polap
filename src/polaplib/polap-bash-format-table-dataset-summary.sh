#!/usr/bin/env bash
# polap-bash-format-table-dataset-summary.sh
# Version: v0.1.5
#
# Steps (no pipes; every step writes a file):
#   0) Generate raw TSV via Python (always with --markdown)
#   1) [qsv|xan] compute scaled columns -> step2_calc.tsv
#   2) [qsv|xan] select final columns -> step3_select.tsv
#   3) [qsv|xan] normalize delimiter to TSV -> step4_fmt.tsv
#   4) Move step4_fmt.tsv -> OUT_TSV
#   5) Optional: pandoc TSV -> markdown_mmd (OUT_TSV.md)

set -euo pipefail

# -------- Init --------
POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export POLAPLIB_DIR

FORMATTER="qsv"
MANIFEST=""
OUT_TSV=""
DO_MARKDOWN=0

die() {
	echo "[ERROR] $*" >&2
	exit 1
}
have() { command -v "$1" >/dev/null 2>&1; }

usage() {
	cat <<'USAGE'
Usage:
  polap-bash-format-table-dataset-summary.sh \
    --formatter xan|qsv \
    --manifest MANIFEST.json \
    --out TABLE.tsv \
    [--markdown]

Assumes input TSV header EXACTLY:
  Species<TAB>Bases<TAB>Reads<TAB>AvgLen<TAB>N50<TAB>AvgQ
USAGE
}

# -------- Args --------
while (($#)); do
	case "$1" in
	--formatter | --formater)
		shift
		FORMATTER="${1:-}"
		[[ -n "$FORMATTER" ]] || die "Missing value for --formatter"
		;;
	--manifest)
		shift
		MANIFEST="${1:-}"
		[[ -n "$MANIFEST" ]] || die "Missing value for --manifest"
		;;
	--out)
		shift
		OUT_TSV="${1:-}"
		[[ -n "$OUT_TSV" ]] || die "Missing value for --out"
		;;
	--markdown) DO_MARKDOWN=1 ;;
	-h | --help)
		usage
		exit 0
		;;
	*) die "Unknown argument: $1" ;;
	esac
	shift
done

[[ -f "$MANIFEST" ]] || die "Manifest not found: $MANIFEST"
[[ -n "$OUT_TSV" ]] || die "--out is required."

out_dir="$(dirname "$OUT_TSV")"
out_base="$(basename "$OUT_TSV" .tsv)"
mkdir -p "$out_dir"

# Intermediate files live next to OUT_TSV
STEP0_RAW="${OUT_TSV%.tsv}.step0_raw.tsv"
STEP1_HDRCHK="${OUT_TSV%.tsv}.step1_header.txt"
STEP2_CALC="${OUT_TSV%.tsv}.step2_calc.tsv"
STEP3_SELECT="${OUT_TSV%.tsv}.step3_select.tsv"
STEP4_FMT="${OUT_TSV%.tsv}.step4_fmt.tsv"

echo "[INFO] POLAPLIB_DIR = $POLAPLIB_DIR"
echo "[INFO] Formatter     = $FORMATTER"
echo "[INFO] Manifest      = $MANIFEST"
echo "[INFO] OUT_TSV       = $OUT_TSV"

# -------- Step 0: Run Python to generate raw TSV --------
echo "[STEP 0] Generating raw TSV via Python -> $STEP0_RAW"
python3 "${POLAPLIB_DIR}/scripts/make_table_dataset_summary.py" \
	--manifest "$MANIFEST" \
	--out "$STEP0_RAW" \
	--markdown

[[ -s "$STEP0_RAW" ]] || die "Python step produced no TSV at: $STEP0_RAW"

# -------- Step 1: Check header exactly --------
echo "[STEP 1] Checking header -> $STEP1_HDRCHK"
head -n1 "$STEP0_RAW" | tr -d '\r' >"$STEP1_HDRCHK"
expected=$'Species\tBases\tReads\tAvgLen\tN50\tAvgQ'
actual="$(cat "$STEP1_HDRCHK")"
if [[ "$actual" != "$expected" ]]; then
	echo "[ERROR] Unexpected header."
	echo "  Expected: [$expected]"
	echo "  Actual:   [$actual]"
	echo "  File: $STEP0_RAW"
	exit 1
fi

# ===== Formatter: qsv path =====
if [[ "$FORMATTER" == "qsv" ]]; then
	have qsv || die "qsv not found in PATH."

	echo "[STEP 2] qsv compute scaled columns -> $STEP2_CALC"
	qsv luau map -d $'\t' Bases_Gb,AvgLen_kb,N50_kb,AvgQ_int \
		"local function num(v)
       if v == nil then return 0 end
       if type(v) == 'number' then return v end
       if type(v) == 'string' then
         -- strip thousands separators if any
         local w = v:gsub(',', '')
         local n = tonumber(w)
         if n == nil then return 0 else return n end
       end
       return 0
     end
     local b  = num(col['Bases'])
     local al = num(col['AvgLen'])
     local n  = num(col['N50'])
     local q  = num(col['AvgQ'])
     -- 1 decimal format as string
     local gb  = string.format('%.1f', (b/1e9))
     local kb1 = string.format('%.1f', (al/1000))
     local kb2 = string.format('%.1f', (n/1000))
     local qi  = string.format('%d', math.floor(q + 0.5))
     return gb, kb1, kb2, qi" \
		"$STEP0_RAW" >"$STEP2_CALC"

	[[ -s "$STEP2_CALC" ]] || die "Empty after qsv luau map: $STEP2_CALC"

	echo "[STEP 3] qsv select final columns -> $STEP3_SELECT"
	qsv select -d $'\t' Species,Bases_Gb,AvgLen_kb,N50_kb,AvgQ_int \
		"$STEP2_CALC" >"$STEP3_SELECT"
	[[ -s "$STEP3_SELECT" ]] || die "Empty after qsv select: $STEP3_SELECT"

	echo "[STEP 4] qsv fmt (TSV) -> $STEP4_FMT"
	qsv fmt -d $'\t' "$STEP3_SELECT" >"$STEP4_FMT"
	[[ -s "$STEP4_FMT" ]] || die "Empty after qsv fmt: $STEP4_FMT"

# ===== Formatter: xan path =====
elif [[ "$FORMATTER" == "xan" ]]; then
	have xan || die "xan not found in PATH."

	# Detect a numeric cast function available in your xan build
	detect_cast() {
		for f in number to_number float; do
			if xan eval "${f}(\"1\")" >/dev/null 2>&1; then
				echo "$f"
				return 0
			fi
		done
		echo ""
		return 1
	}
	CAST_FN="$(detect_cast || true)"
	[[ -n "$CAST_FN" ]] || die "xan numeric cast function not found (tried: number,to_number,float)."

	# Step 2: compute scaled columns with explicit casts, one file per operation
	# 2.1 cast numeric columns
	echo "[STEP 2.1] xan cast numerics -> ${STEP2_CALC%.tsv}.cast.tsv"
	XAN_CAST="${STEP2_CALC%.tsv}.cast.tsv"
	xan -d $'\t' map "${CAST_FN}(Bases) as _Bases" "$STEP0_RAW" >"$XAN_CAST"
	xan -d $'\t' map "${CAST_FN}(AvgLen) as _AvgLen" "$XAN_CAST" >"${XAN_CAST%.tsv}.2.tsv"
	mv "${XAN_CAST%.tsv}.2.tsv" "$XAN_CAST"
	xan -d $'\t' map "${CAST_FN}(N50) as _N50" "$XAN_CAST" >"${XAN_CAST%.tsv}.3.tsv"
	mv "${XAN_CAST%.tsv}.3.tsv" "$XAN_CAST"
	xan -d $'\t' map "${CAST_FN}(AvgQ) as _AvgQ" "$XAN_CAST" >"${XAN_CAST%.tsv}.4.tsv"
	mv "${XAN_CAST%.tsv}.4.tsv" "$XAN_CAST"

	# 2.2 compute scaled columns
	echo "[STEP 2.2] xan math -> ${STEP2_CALC}"
	xan -d $'\t' map '(_Bases/1e9) as Bases_Gb' "$XAN_CAST" >"$STEP2_CALC"
	xan -d $'\t' map '(_AvgLen/1000) as AvgLen_kb' "$STEP2_CALC" >"${STEP2_CALC%.tsv}.2.tsv"
	mv "${STEP2_CALC%.tsv}.2.tsv" "$STEP2_CALC"
	xan -d $'\t' map '(_N50/1000) as N50_kb' "$STEP2_CALC" >"${STEP2_CALC%.tsv}.3.tsv"
	mv "${STEP2_CALC%.tsv}.3.tsv" "$STEP2_CALC"
	xan -d $'\t' map 'round(_AvgQ) as AvgQ_int' "$STEP2_CALC" >"${STEP2_CALC%.tsv}.4.tsv"
	mv "${STEP2_CALC%.tsv}.4.tsv" "$STEP2_CALC"

	# 2.3 format to 1 decimal where needed (becomes strings)
	echo "[STEP 2.3] xan fmt decimals -> ${STEP2_CALC%.tsv}.fmt.tsv"
	XAN_FMT="${STEP2_CALC%.tsv}.fmt.tsv"
	xan -d $'\t' map 'fmt(Bases_Gb, "%.1f") as Bases_Gb' "$STEP2_CALC" >"$XAN_FMT"
	xan -d $'\t' map 'fmt(AvgLen_kb, "%.1f") as AvgLen_kb' "$XAN_FMT" >"${XAN_FMT%.tsv}.2.tsv"
	mv "${XAN_FMT%.tsv}.2.tsv" "$XAN_FMT"
	xan -d $'\t' map 'fmt(N50_kb, "%.1f") as N50_kb' "$XAN_FMT" >"${XAN_FMT%.tsv}.3.tsv"
	mv "${XAN_FMT%.tsv}.3.tsv" "$XAN_FMT"

	# -------- Step 3: select final columns --------
	echo "[STEP 3] xan select -> $STEP3_SELECT"
	xan -d $'\t' select Species Bases_Gb AvgLen_kb N50_kb AvgQ_int "$XAN_FMT" >"$STEP3_SELECT"

	[[ -s "$STEP3_SELECT" ]] || die "Empty after xan select: $STEP3_SELECT"

	# -------- Step 4: force TSV delimiter (xan fmt) --------
	echo "[STEP 4] xan fmt (TSV) -> $STEP4_FMT"
	xan -d $'\t' fmt -d $'\t' "$STEP3_SELECT" >"$STEP4_FMT"

	[[ -s "$STEP4_FMT" ]] || die "Empty after xan fmt: $STEP4_FMT"

else
	die "--formatter must be 'xan' or 'qsv' (got: $FORMATTER)"
fi

# -------- Finalize --------
echo "[STEP 5] Move formatted TSV -> $OUT_TSV"
mv -f "$STEP4_FMT" "$OUT_TSV"
echo "[OK] Wrote: $OUT_TSV"

# -------- Optional: pandoc to markdown_mmd --------
if ((DO_MARKDOWN)); then
	[[ -n "${PANDOC_PATH:-}" ]] && PANDOC_BIN="$PANDOC_PATH" || PANDOC_BIN="pandoc"
	have "$PANDOC_BIN" || die "pandoc not found in PATH."
	MD_OUT="${OUT_TSV%.tsv}.md"
	echo "[STEP 6] pandoc TSV -> markdown_mmd -> $MD_OUT"
	"$PANDOC_BIN" -f csv --separator=$'\t' -t markdown_mmd "$OUT_TSV" -o "$MD_OUT"
	echo "[OK] Wrote: $MD_OUT"
fi

echo "[DONE]"
