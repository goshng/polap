#!/usr/bin/env bash
# polap-bash-run-sim-report.sh
# Drive: simulate ONT reads with polap-bash-sim-ont.sh, then emit an HTML report.
#
# Requirements:
#   - polap-bash-sim-ont.sh (in PATH or same directory)
#   - seqkit
#   - coreutils, awk, sed
#
# Example:
#   ./polap-bash-run-sim-report.sh \
#       --mt mt.fa --pt pt.fa --out sim_clean --threads 8 \
#       --depth-nuc 10x --depth-mt 50x --depth-pt 500x
#
# Example (with NUMT/NUPT spiking):
#   ./polap-bash-run-sim-report.sh \
#       --mt mt.fa --pt pt.fa --nuc nuclear.fa --out sim_spiked \
#       --spike-numt 80 --spike-nupt 300 --threads 8
#
set -euo pipefail

SELF="$(readlink -f "${BASH_SOURCE[0]:-$0}")"
SELF_DIR="$(dirname "$SELF")"
SIM_BIN="${SIM_BIN:-polap-bash-sim-ont.sh}"

if ! command -v "$SIM_BIN" >/dev/null 2>&1; then
  # try alongside this script
  if [[ -x "${SELF_DIR}/polap-bash-sim-ont.sh" ]]; then
    SIM_BIN="${SELF_DIR}/polap-bash-sim-ont.sh"
  else
    echo "[ERR] polap-bash-sim-ont.sh not found (set SIM_BIN or place it in PATH)" >&2
    exit 1
  fi
fi

SIM_OUT="sim_ont"
ARGS=()

# pass-through arguments to the simulator, also capture --out if present
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) SIM_OUT="$2"; ARGS+=("$1" "$2"); shift 2;;
    *) ARGS+=("$1"); shift;;
  esac
done

mkdir -p "$SIM_OUT"
pushd "$SIM_OUT" >/dev/null

# Run simulation
echo "[sim] running: $SIM_BIN ${ARGS[*]}"
"$SIM_BIN" "${ARGS[@]}"

# Files we expect from simulator
ALL_FQ="all.sim.fq"
NUC_FQ="nuc.fq"
MT_FQ="mt.fq"
PT_FQ="pt.fq"
NUC_FA="nuclear.fa"

# Basic checks
for f in "$ALL_FQ" "$NUC_FQ" "$MT_FQ" "$PT_FQ" "$NUC_FA"; do
  [[ -s "$f" ]] || { echo "[ERR] missing output $f" >&2; exit 2; }
done

# Collect statistics with seqkit
seqkit stats -Ta "$ALL_FQ" "$NUC_FQ" "$MT_FQ" "$PT_FQ" > stats.tsv

# Detect spiking truth
TRUTH_BED=""
TRUTH_TSV=""
if [[ -s nuclear.spiked.inserts.bed ]]; then TRUTH_BED="nuclear.spiked.inserts.bed"; fi
if [[ -s nuclear.spiked.inserts.tsv ]]; then TRUTH_TSV="nuclear.spiked.inserts.tsv"; fi

# Summarize truth if present
TRUTH_SUMMARY=""
if [[ -n "$TRUTH_TSV" ]]; then
  # columns: chrom start end name strand src src_name src_start src_end req_len mut_len obs_sub obs_indel
  mapfile -t T_SUM < <(awk 'BEGIN{OFS="\t"} NR>1{n[$5]++; L+=$11; if($5=="NUMT") l_numt+=$11; if($5=="NUPT") l_nupt+=$11} END{
      print n["NUMT"]+0, n["NUPT"]+0, L+0, l_numt+0, l_nupt+0
  }' "$TRUTH_TSV")
  NUMT_N=${T_SUM[0]%$'\t'*}
  # safer split
  read -r NUMT_N NUPT_N TOT_L NUMT_L NUPT_L <<<"$(printf "%s\n" "${T_SUM[0]}")"
  TRUTH_SUMMARY=$(cat <<EOF
<p><b>Spiking truth:</b> NUMT inserts: ${NUMT_N}, NUPT inserts: ${NUPT_N}. Total inserted bp: ${TOT_L} (NUMT ${NUMT_L}, NUPT ${NUPT_L}).</p>
EOF
)
fi

# Helper: escape HTML
html_escape() {
  sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g'
}

# Convert TSV to HTML table
tsv_to_html_table() {
  local tsv="$1"
  local title="$2"
  echo "<h3>${title}</h3>"
  echo "<table>"
  local NR=0
  while IFS=$'\t' read -r -a F; do
    if (( NR==0 )); then
      echo "<thead><tr>"
      for c in "${F[@]}"; do
        printf "<th>%s</th>" "$(printf "%s" "$c" | html_escape)"
      done
      echo "</tr></thead><tbody>"
    else
      echo "<tr>"
      for c in "${F[@]}"; do
        printf "<td>%s</td>" "$(printf "%s" "$c" | html_escape)"
      done
      echo "</tr>"
    fi
    ((NR++))
  done < "$tsv"
  echo "</tbody></table>"
}

# Capture the exact simulator command/args into a block
SIM_CMD_HTML="<pre>$(printf "%q " "$SIM_BIN" "${ARGS[@]}" | html_escape)</pre>"

# Write HTML
REPORT="report.html"
cat > "$REPORT" <<'HTML'
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>POLAP ONT Simulation Report</title>
<style>
  body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 2rem; color: #111; }
  header { margin-bottom: 1rem; }
  h1 { font-size: 1.8rem; margin: 0 0 .5rem 0; }
  h2 { margin-top: 2rem; }
  h3 { margin-top: 1.2rem; }
  code, pre { background: #f6f8fa; padding: .4rem .6rem; border-radius: 6px; }
  table { border-collapse: collapse; margin: 1rem 0; width: 100%; }
  th, td { border: 1px solid #ddd; padding: .4rem .6rem; font-size: .95rem; }
  th { background: #fafafa; text-align: left; }
  .files a { text-decoration: none; color: #0366d6; }
  .badge { display: inline-block; padding: .2rem .5rem; border-radius: 999px; background: #eef; margin-left: .5rem; font-size: .8rem; }
</style>
</head>
<body>
<header>
  <h1>POLAP: ONT Simulation Report</h1>
  <p>This report summarizes outputs from <code>polap-bash-sim-ont.sh</code>.</p>
</header>

<section id="overview">
  <h2>Overview</h2>
  __SIM_CMD__
  __TRUTH_SUM__
</section>

<section id="stats">
  <h2>Read statistics (seqkit)</h2>
  __STATS_TABLE__
</section>

<section id="files">
  <h2>Output files</h2>
  <ul class="files">
    <li><a href="all.sim.fq">all.sim.fq</a> <span class="badge">combined</span></li>
    <li><a href="nuc.fq">nuc.fq</a> <span class="badge">nuclear</span></li>
    <li><a href="mt.fq">mt.fq</a> <span class="badge">mitochondrial</span></li>
    <li><a href="pt.fq">pt.fq</a> <span class="badge">plastid</span></li>
    <li><a href="nuclear.fa">nuclear.fa</a> <span class="badge">sim reference</span></li>
    <li>Spiking truth (if present):
      <ul>
        <li><a href="nuclear.spiked.inserts.bed">nuclear.spiked.inserts.bed</a></li>
        <li><a href="nuclear.spiked.inserts.tsv">nuclear.spiked.inserts.tsv</a></li>
      </ul>
    </li>
  </ul>
</section>

<footer>
  <p>Generated by <code>polap-bash-run-sim-report.sh</code>.</p>
</footer>
</body>
</html>
HTML

# Inject dynamic parts
# Insert the sim command block
sed -i "s|__SIM_CMD__|${SIM_CMD_HTML//|/\\|}|" "$REPORT"

# Insert truth summary (or blank)
if [[ -n "$TRUTH_SUMMARY" ]]; then
  esc=$(printf "%s" "$TRUTH_SUMMARY" | sed -e 's/[\/&]/\\&/g')
  sed -i "s|__TRUTH_SUM__|$esc|" "$REPORT"
else
  sed -i "s|__TRUTH_SUM__||" "$REPORT"
fi

# Insert stats table
TABLE_HTML="$(tsv_to_html_table stats.tsv "Seqkit stats for combined and per-origin reads")"
esc_table=$(printf "%s" "$TABLE_HTML" | sed -e 's/[\/&]/\\&/g')
sed -i "s|__STATS_TABLE__|$esc_table|" "$REPORT"

echo "[ok] Wrote $REPORT"
popd >/dev/null

# Print a local path to the report
echo "$(readlink -f "$SIM_OUT/$REPORT")"
