#!/usr/bin/env bash
# polap-bash-paf-diagnose.sh  v0.0.1
set -euo pipefail

VER="v0.0.1"
_verbose=0
_quiet=0

usage() {
	cat <<EOF
polap-bash-paf-diagnose.sh ${VER}
Run PAF diagnostics (identity/qcov sweeps), write TSVs, PNGs, and an HTML report.

Usage:
  polap-bash-paf-diagnose.sh -i mt.paf -o out/prefix [options]

Required:
  -i, --paf FILE          Input PAF
  -o, --out-prefix STR    Output prefix (dir/prefix)

Options:
  --fixed-qcov FLOAT      qcov fixed when sweeping identity [0.50]
  --fixed-id   FLOAT      id   fixed when sweeping qcov     [0.80]
  --min-alnlen INT        minimal alignment length (PAF col11) prefilter [0]
  --primary-only          keep tp:A:P (primary) only [off]

  -v, --verbose           verbose logs
      --quiet             suppress info logs
      --version           print version
      --help              this help

Produces:
  <prefix>.sweep_identity.tsv/.png
  <prefix>.sweep_qcov.tsv/.png
  <prefix>.sweep_combined.tsv
  <prefix>.report.html
EOF
}

logi() {
	if [[ ${_quiet:-0} -eq 0 && ${_verbose:-0} -gt 0 ]]; then
		echo "[INFO]" "$@" >&2
	fi
}

loge() { echo "[ERROR]" "$@" >&2; }

# defaults
PAF=""
OUTP=""
FIXED_QCOV="0.50"
FIXED_ID="0.80"
MIN_ALNLEN="0"
PRIMARY_ONLY="false"

# parse
while [[ $# -gt 0 ]]; do
	case "$1" in
	-i | --paf)
		PAF="$2"
		shift 2
		;;
	-o | --out-prefix)
		OUTP="$2"
		shift 2
		;;
	--fixed-qcov)
		FIXED_QCOV="$2"
		shift 2
		;;
	--fixed-id)
		FIXED_ID="$2"
		shift 2
		;;
	--min-alnlen)
		MIN_ALNLEN="$2"
		shift 2
		;;
	--primary-only)
		PRIMARY_ONLY="true"
		shift
		;;
	-v | --verbose)
		_verbose=$((_verbose + 1))
		shift
		;;
	--quiet)
		_quiet=1
		_verbose=0
		shift
		;;
	--version)
		echo "polap-bash-paf-diagnose.sh ${VER}"
		exit 0
		;;
	--help)
		usage
		exit 0
		;;
	*)
		loge "unknown option: $1"
		usage
		exit 2
		;;
	esac
done

[[ -n "$PAF" && -s "$PAF" ]] || {
	loge "missing or empty --paf"
	usage
	exit 2
}
[[ -n "$OUTP" ]] || {
	loge "missing --out-prefix"
	usage
	exit 2
}

OUTDIR="$(dirname -- "$OUTP")"
mkdir -p "$OUTDIR"

# Resolve directory of this script (absolute path, even if symlinked)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# find companion R script (assume in PATH or alongside this bash)
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
R_DIAG="${R_DIAG:-"${SCRIPT_DIR}/polap-r-paf-diagnose.R"}"

command -v "$RSCRIPT_BIN" >/dev/null 2>&1 || {
	loge "Rscript not found"
	exit 2
}
# command -v "$R_DIAG" >/dev/null 2>&1 || {
# 	# try relative to this script
# 	HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 	if [[ -x "${HERE}/polap-r-paf-diagnose.R" ]]; then
# 		R_DIAG="${HERE}/polap-r-paf-diagnose.R"
# 	else
# 		loge "polap-r-paf-diagnose.R not found in PATH nor in ${HERE}"
# 		exit 2
# 	fi
# }
#

# set -x
# trap 'echo "ERR at line $LINENO (status=$?)"' ERR

logi "running R diagnose: $R_DIAG"
# printf '%q ' "$RSCRIPT_BIN" "--vanilla" "$R_DIAG" "$PAF" "$OUTP" "$FIXED_QCOV" "$FIXED_ID" "$MIN_ALNLEN" "$PRIMARY_ONLY" 2>&1

"$RSCRIPT_BIN" --vanilla "$R_DIAG" "$PAF" "$OUTP" "$FIXED_QCOV" "$FIXED_ID" "$MIN_ALNLEN" "$PRIMARY_ONLY"

ID_TSV="${OUTP}.sweep_identity.tsv"
QC_TSV="${OUTP}.sweep_qcov.tsv"
COMBO_TSV="${OUTP}.sweep_combined.tsv"
ID_PNG="${OUTP}.sweep_identity.png"
QC_PNG="${OUTP}.sweep_qcov.png"
HTML="${OUTP}.report.html"

# TSV->HTML table (minimal; replaces tabs with <td>)
tsv_to_table() {
	local tsv="$1"
	[[ -s "$tsv" ]] || {
		echo "<p><i>Missing: $(basename "$tsv")</i></p>"
		return
	}
	awk -v FN="$(basename "$tsv")" '
    BEGIN{
      print "<details open><summary><b>" FN "</b></summary>";
      print "<div style=\"overflow:auto\"><table border=1 cellpadding=4 cellspacing=0>";
    }
    NR==1{
      printf "<tr style=\"background:#f0f0f0\">"
      for(i=1;i<=NF;i++){printf "<th>%s</th>", $i}
      print "</tr>";
      next
    }
    {
      printf "<tr>"
      for(i=1;i<=NF;i++){gsub(/&/,"&amp;",$i); gsub(/</,"&lt;",$i); gsub(/>/,"&gt;",$i); printf "<td>%s</td>", $i}
      print "</tr>"
    }
    END{
      print "</table></div></details>"
    }
  ' OFS='\t' "$tsv"
}

# Build HTML report
logi "writing HTML report: $HTML"
{
	echo "<!doctype html><html><head><meta charset='utf-8'>"
	echo "<title>PAF Diagnose Report ($(basename "$PAF"))</title>"
	echo "<style>body{font-family:sans-serif;max-width:1000px;margin:2em auto;line-height:1.4} code,pre{background:#f7f7f7} img{max-width:100%}</style>"
	echo "</head><body>"
	echo "<h1>PAF Diagnose Report</h1>"
	echo "<p><b>Script</b>: polap-bash-paf-diagnose.sh ${VER}</p>"
	echo "<h2>Inputs</h2>"
	echo "<ul>"
	echo "<li>PAF: <code>$(basename "$PAF")</code></li>"
	echo "<li>Out prefix: <code>$OUTP</code></li>"
	echo "<li>fixed_qcov (for ID sweep): <code>$FIXED_QCOV</code></li>"
	echo "<li>fixed_id (for QCOV sweep): <code>$FIXED_ID</code></li>"
	echo "<li>min_alnlen (prefilter): <code>$MIN_ALNLEN</code></li>"
	echo "<li>primary_only: <code>$PRIMARY_ONLY</code></li>"
	echo "</ul>"

	echo "<h2>Plots</h2>"
	if [[ -s "$ID_PNG" ]]; then
		echo "<h3>Percent passing vs identity</h3><img src='$(basename "$ID_PNG")' alt='identity plot'/>"
	else
		echo "<p><i>Missing: $(basename "$ID_PNG")</i></p>"
	fi
	if [[ -s "$QC_PNG" ]]; then
		echo "<h3>Percent passing vs qcov</h3><img src='$(basename "$QC_PNG")' alt='qcov plot'/>"
	else
		echo "<p><i>Missing: $(basename "$QC_PNG")</i></p>"
	fi

	echo "<h2>Tables</h2>"
	tsv_to_table "$ID_TSV"
	tsv_to_table "$QC_TSV"
	tsv_to_table "$COMBO_TSV"

	echo "<hr><p><small>Generated by polap-bash-paf-diagnose.sh ${VER}</small></p>"
	echo "</body></html>"
} >"$HTML"

logi "done. Open: $HTML"
echo "$HTML"
