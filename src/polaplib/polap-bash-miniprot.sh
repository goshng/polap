#!/usr/bin/env bash
# polap-bash-miniprot.sh v0.1.0
# miniprot-delegate: miniprot wrapper with --progress
# Usage:
#   miniprot-delegate [miniprot options] <ref.fa> <query.faa> [...] [--progress]
#
# Notes:
# * Identical to miniprot CLI, plus:
#     --progress        Show a live progress meter on stderr
# * Progress = unique query IDs observed in output / total queries in input FASTA(s)
# * Supports default (PAF) and --gff/--gtf outputs.
# * If any query is "-", we cannot pre-count; progress falls back to "unknown total".
# * Requires: awk, grep, gzip (for .gz queries), and miniprot in PATH.

set -euo pipefail

# --- Parse argv: we only care about detecting:
#     - presence of --progress
#     - presence of --gff / --gtf (to choose parser)
#     - locating the first non-option as <ref> and the rest non-options as queries
progress=0
fmt="paf" # "paf" (default) | "gff" | "gtf"
declare -a passthru=()

# We’ll store positional non-options as they appear to detect ref + queries
declare -a positionals=()

for arg in "$@"; do
	case "$arg" in
	--progress) progress=1 ;; # swallow this; not passed to miniprot
	--gff)
		fmt="gff"
		passthru+=("$arg")
		;;
	--gtf)
		fmt="gtf"
		passthru+=("$arg")
		;;
	--*) passthru+=("$arg") ;;
	-*) passthru+=("$arg") ;;
	*) positionals+=("$arg") ;;
	esac
done

if ! command -v miniprot >/dev/null 2>&1; then
	echo "[ERR] miniprot not found in PATH." >&2
	exit 127
fi

# Need at least ref + one query
if ((${#positionals[@]} < 2)); then
	echo "[ERR] Usage: miniprot [options] <ref.fa> <query.faa> [...] [--progress]" >&2
	exit 2
fi

ref="${positionals[0]}"
queries=("${positionals[@]:1}")

# If no progress requested, just call through exactly
if ((progress == 0)); then
	exec miniprot "${passthru[@]}" "$ref" "${queries[@]}"
fi

# --- Compute total number of query sequences across all query files, if possible
#     If any query is "-", total becomes unknown (-1).
total=0
unknown_total=0

for q in "${queries[@]}"; do
	if [[ "$q" == "-" ]]; then
		unknown_total=1
		break
	elif [[ -f "$q" ]]; then
		if [[ "$q" == *.gz ]]; then
			# count header lines in gz FASTA
			c=$(gzip -cd -- "$q" | grep -c '^>' || true)
		else
			c=$(grep -c '^>' -- "$q" || true)
		fi
		total=$((total + c))
	else
		echo "[WARN] Query not found: $q (skipping from total count)" >&2
	fi
done

# --- Build AWK program to count unique queries seen in the stream.
# We pass through all output to stdout, but update progress to stderr.

# For PAF: query name is column 1.
read -r -d '' AWK_PAF <<'AWK'
BEGIN{FS="\t"; done=0; last=0}
{
  q=$1
  if(!(q in seen)) { seen[q]=1; done++ }
  now=systime()
  if (now-last>=1 || done%100==0) {
    if (TOTAL >= 0) {
      pct = (TOTAL>0)? (100.0*done/TOTAL) : 0
      printf("\r[miniprot] %5.1f%%  (%d/%d)", pct, done, TOTAL) > "/dev/stderr"
    } else {
      printf("\r[miniprot] %d done (total unknown)", done) > "/dev/stderr"
    }
    last=now
  }
  print $0
}
END{
  if (TOTAL >= 0) {
    printf("\r[miniprot] %5.1f%%  (%d/%d)\n", (TOTAL>0?100.0:0), done, TOTAL) > "/dev/stderr"
  } else {
    printf("\r[miniprot] %d done (total unknown)\n", done) > "/dev/stderr"
  }
}
AWK

# For GFF3: try Target=QNAME first; else Name=QNAME.
read -r -d '' AWK_GFF <<'AWK'
BEGIN{FS="\t"; OFS="\t"; done=0; last=0}
{
  if ($0 ~ /^#/) { print $0; next }
  attr=$9
  q=""
  if (match(attr, /Target=([^ ;\t]+)/, m)) q=m[1]
  else if (match(attr, /Name=([^;]+)/, n)) q=n[1]
  if (q != "" && !(q in seen)) { seen[q]=1; done++ }
  now=systime()
  if (now-last>=1 || done%50==0) {
    if (TOTAL >= 0) {
      pct = (TOTAL>0)? (100.0*done/TOTAL) : 0
      printf("\r[miniprot] %5.1f%%  (%d/%d)", pct, done, TOTAL) > "/dev/stderr"
    } else {
      printf("\r[miniprot] %d done (total unknown)", done) > "/dev/stderr"
    }
    last=now
  }
  print $0
}
END{
  if (TOTAL >= 0) {
    printf("\r[miniprot] %5.1f%%  (%d/%d)\n", (TOTAL>0?100.0:0), done, TOTAL) > "/dev/stderr"
  } else {
    printf("\r[miniprot] %d done (total unknown)\n", done) > "/dev/stderr"
  }
}
AWK

# For GTF: try transcript_id "QNAME"; otherwise gene_id "QNAME".
read -r -d '' AWK_GTF <<'AWK'
BEGIN{FS="\t"; OFS="\t"; done=0; last=0}
{
  if ($0 ~ /^#/) { print $0; next }
  attr=$9
  q=""
  if (match(attr, /transcript_id[[:space:]]+"([^"]+)"/, m)) q=m[1]
  else if (match(attr, /gene_id[[:space:]]+"([^"]+)"/, n)) q=n[1]
  if (q != "" && !(q in seen)) { seen[q]=1; done++ }
  now=systime()
  if (now-last>=1 || done%50==0) {
    if (TOTAL >= 0) {
      pct = (TOTAL>0)? (100.0*done/TOTAL) : 0
      printf("\r[miniprot] %5.1f%%  (%d/%d)", pct, done, TOTAL) > "/dev/stderr"
    } else {
      printf("\r[miniprot] %d done (total unknown)", done) > "/dev/stderr"
    }
    last=now
  }
  print $0
}
END{
  if (TOTAL >= 0) {
    printf("\r[miniprot] %5.1f%%  (%d/%d)\n", (TOTAL>0?100.0:0), done, TOTAL) > "/dev/stderr"
  } else {
    printf("\r[miniprot] %d done (total unknown)\n", done) > "/dev/stderr"
  }
}
AWK

# Decide awk program and TOTAL
if ((unknown_total == 1)); then
	TOTAL="-1"
else
	TOTAL="$total"
fi

# Run miniprot and pipe through the progress counter. Keep miniprot's own stderr intact.
case "$fmt" in
gff)
	miniprot "${passthru[@]}" "$ref" "${queries[@]}" 2> >(cat >&2) |
		awk -v TOTAL="$TOTAL" "$AWK_GFF"
	;;
gtf)
	miniprot "${passthru[@]}" "$ref" "${queries[@]}" 2> >(cat >&2) |
		awk -v TOTAL="$TOTAL" "$AWK_GTF"
	;;
*)
	miniprot "${passthru[@]}" "$ref" "${queries[@]}" 2> >(cat >&2) |
		awk -v TOTAL="$TOTAL" "$AWK_PAF"
	;;
esac
