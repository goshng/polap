#!/usr/bin/env bash
# polap-bash-convert-gfa2fasta.sh
# Convert GFA (v1/rGFA) to FASTA using gfatools, with optional ID filtering via seqkit/seqtk.
# Dependencies: gfatools (required), plus at least one of: seqkit or seqtk (required when --ids or --line-width used)
#
# Usage:
#   polap-bash-convert-gfa2fasta.sh [--ids LIST_OR_FILE] [--line-width N] [--stable] INPUT.gfa[.gz|-] OUTPUT.fa[.gz|-]
#
# Examples:
#   polap-bash-convert-gfa2fasta.sh graph.gfa out.fa
#   polap-bash-convert-gfa2fasta.sh --stable graph.rgfa out.fa
#   polap-bash-convert-gfa2fasta.sh --ids s1,s3,s7 graph.gfa out.fa.gz
#   zcat graph.gfa.gz | polap-bash-convert-gfa2fasta.sh --ids ids.txt - out.fa
#
# Notes:
# - --ids accepts a comma/space-separated list or a file (one ID per line). Trailing +/− are stripped.
# - --stable passes '-s' to gfatools gfa2fa (recommended for rGFA W-line/stable names).
# - When OUTPUT ends with .gz, the result is gzip-compressed.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  polap-bash-convert-gfa2fasta.sh [--ids LIST_OR_FILE] [--line-width N] [--stable] INPUT.gfa[.gz|-] OUTPUT.fa[.gz|-]

Options:
  --ids LIST_OR_FILE   Comma/space-separated IDs or a file with one ID per line (orientation suffix +/− is ignored).
  --line-width N       Wrap FASTA lines at N columns. If omitted, leave as produced by gfatools; 0 = single line.
  --stable             Use gfatools gfa2fa -s (rGFA stable FASTA).
  -h, --help           Show this help.

Examples:
  polap-bash-convert-gfa2fasta.sh graph.gfa out.fa
  polap-bash-convert-gfa2fasta.sh --ids ids.txt graph.gfa out.fa.gz
  zcat graph.gfa.gz | polap-bash-convert-gfa2fasta.sh --stable - out.fa
EOF
}

# --- parse args ---
ids_arg=""
width_opt=""   # empty means "don't touch wrapping"; integer means rewrap; 0 means single-line
stable=0
pos=()

while (($#)); do
  case "$1" in
    --ids)        ids_arg="${2:-}"; shift 2;;
    --line-width) width_opt="${2:-}"; shift 2;;
    --stable)     stable=1; shift;;
    -h|--help)    usage; exit 0;;
    --) shift; while (($#)); do pos+=("$1"); shift; done;;
    -*) echo "ERROR: unknown option: $1" >&2; usage; exit 2;;
    *)  pos+=("$1"); shift;;
  esac
done

[[ ${#pos[@]} -eq 2 ]] || { echo "ERROR: need INPUT and OUTPUT" >&2; usage; exit 2; }
in_gfa="${pos[0]}"
out_fa="${pos[1]}"

# --- deps ---
command -v gfatools >/dev/null 2>&1 || {
  echo "ERROR: gfatools not found. Try: mamba install -c bioconda gfatools" >&2
  exit 3
}
# seqkit/seqtk are only needed if filtering or rewrapping
need_filter=0
[[ -n "$ids_arg" ]] && need_filter=1
need_rewrap=0
[[ -n "$width_opt" ]] && need_rewrap=1

have_seqkit=0; command -v seqkit >/dev/null 2>&1 && have_seqkit=1
have_seqtk=0;  command -v seqtk  >/dev/null 2>&1 && have_seqtk=1

if (( need_filter || (need_rewrap && -n "$width_opt") )); then
  if (( ! have_seqkit && ! have_seqtk )); then
    echo "ERROR: --ids/--line-width require seqkit or seqtk. Install one (e.g., mamba install -c bioconda seqkit seqtk)." >&2
    exit 3
  fi
fi

# --- temps and cleanup ---
tmpfiles=()
cleanup() { ((${#tmpfiles[@]})) && rm -f "${tmpfiles[@]}" 2>/dev/null || true; }
trap cleanup EXIT

tmp_gfa=""
if [[ "$in_gfa" == "-" ]]; then
  tmp_gfa="$(mktemp -t gfa2fa.in.XXXXXX)"; tmpfiles+=("$tmp_gfa")
  cat > "$tmp_gfa"
  in_gfa="$tmp_gfa"
elif [[ "$in_gfa" =~ \.gz$ ]]; then
  tmp_gfa="$(mktemp -t gfa2fa.in.XXXXXX)"; tmpfiles+=("$tmp_gfa")
  gzip -dc -- "$in_gfa" > "$tmp_gfa"
  in_gfa="$tmp_gfa"
fi

tmp_fa="$(mktemp -t gfa2fa.fa.XXXXXX)"; tmpfiles+=("$tmp_fa")
tmp_sel="$(mktemp -t gfa2fa.sel.XXXXXX)"; tmpfiles+=("$tmp_sel")
tmp_out="$(mktemp -t gfa2fa.out.XXXXXX)"; tmpfiles+=("$tmp_out")

# --- 1) GFA -> FASTA via gfatools ---
#     -s for stable FASTA when rGFA is provided
if (( stable )); then
  gfatools gfa2fa -s "$in_gfa" > "$tmp_fa"
else
  gfatools gfa2fa "$in_gfa" > "$tmp_fa"
fi

# --- 2) Optional: subset by IDs ---
if [[ -n "$ids_arg" ]]; then
  # normalize ids to file (strip whitespace and trailing + or -)
  ids_file="$ids_arg"
  if [[ ! -f "$ids_arg" ]]; then
    ids_file="$tmp_sel"
    printf "%s\n" "$ids_arg" \
      | tr ',;\t ' '\n' \
      | sed -E 's/^[[:space:]]+|[[:space:]]+$//g; s/[+-]$//; /^$/d' \
      > "$ids_file"
  else
    # clean in-place copy to $tmp_sel
    sed -E 's/[[:space:]]+$//; s/[+-]$//; /^[[:space:]]*$/d' "$ids_arg" > "$tmp_sel"
    ids_file="$tmp_sel"
  fi

  if (( have_seqkit )); then
    # Match by name; preserves original FASTA header formatting. (seqkit is robust & fast)
    seqkit grep -n -f "$ids_file" "$tmp_fa" -o "$tmp_out"
  else
    # Fallback: seqtk subseq (IDs must match the name token)
    seqtk subseq "$tmp_fa" "$ids_file" > "$tmp_out"
  fi
  mv -f "$tmp_out" "$tmp_fa"
fi

# --- 3) Optional: rewrap line width ---
if [[ -n "${width_opt}" ]]; then
  width="${width_opt}"
  if (( have_seqkit )); then
    # seqkit: -w N (0 = no wrap)
    seqkit seq -w "$width" "$tmp_fa" -o "$tmp_out"
  else
    # seqtk: -l N (0 = single-line)
    seqtk seq -l "$width" "$tmp_fa" > "$tmp_out"
  fi
  mv -f "$tmp_out" "$tmp_fa"
fi

# --- 4) Write OUTPUT (support stdout and .gz) ---
if [[ "$out_fa" == "-" ]]; then
  cat "$tmp_fa"
elif [[ "$out_fa" =~ \.gz$ ]]; then
  gzip -c "$tmp_fa" > "$out_fa"
else
  mv -f "$tmp_fa" "$out_fa"
  # remove from tmpfiles list so trap doesn't try to rm a moved file
  tmpfiles=("${tmpfiles[@]/$tmp_fa}")
fi
