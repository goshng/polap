#!/usr/bin/env bash
# polap-bash-table-data.sh
# Version: v0.5.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a minimal data table from a manifest:
#   code2 | species | long_sra | long_total_bases | short_sra | short1_total_bases
# TSV -> raw numbers; MD -> Gb-formatted totals
#
# Prefers the Rust CLI **qsv** for CSV work (available via conda-forge).
# Falls back to Python if qsv is not found.
#
# Usage:
#   polap-bash-table-data.sh \
#     --manifest md/manifest.json \
#     --out-tsv   md/table-data.tsv \
#     --out-md    md/table-data.md \
#     [--out-csv  md/table-data.csv] \
#     [--sort code2|species]
#
set -euo pipefail
IFS=$'\n\t'

MANIFEST=""
OUT_TSV=""
OUT_MD=""
OUT_CSV=""
SORT="code2" # code2|species

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest MANIFEST.json --out-tsv OUT.tsv [--out-md OUT.md] [--out-csv OUT.csv] [--sort code2|species]

Description:
  Create a compact data table with columns:
    code2 | species | long_sra | long_total_bases | short_sra | short1_total_bases

  TSV has raw totals (bases). Markdown shows totals in Gb (e.g., "7.7 Gb").

Options:
  --manifest   Path to manifest JSON
  --out-tsv    Path to output TSV (required)
  --out-md     Path to output Markdown (optional)
  --out-csv    Path to output CSV (optional; will be created if --out-md is given)
  --sort       Sort key (code2|species) [default: code2]
  -h, --help   Show this help
EOF
}

while (($#)); do
	case "$1" in
	--manifest)
		MANIFEST="${2:?}"
		shift 2
		;;
	--out-tsv)
		OUT_TSV="${2:?}"
		shift 2
		;;
	--out-md)
		OUT_MD="${2:?}"
		shift 2
		;;
	--out-csv)
		OUT_CSV="${2:?}"
		shift 2
		;;
	--sort)
		SORT="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERROR] Unknown option: $1" >&2
		usage
		exit 2
		;;
	esac
done

[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || {
	echo "[ERROR] --manifest not found: $MANIFEST" >&2
	exit 2
}
[[ -n "$OUT_TSV" ]] || {
	echo "[ERROR] --out-tsv is required" >&2
	exit 2
}

mkdir -p "$(dirname "$OUT_TSV")"
[[ -n "${OUT_MD:-}" ]] && mkdir -p "$(dirname "$OUT_MD")"
[[ -n "${OUT_CSV:-}" ]] && mkdir -p "$(dirname "$OUT_CSV")"

have() { command -v "$1" >/dev/null 2>&1; }

# ---------- 1) Build TSV (raw) and CSV (Gb) via embedded Python ----------
TMP_CSV=""
if [[ -z "${OUT_CSV:-}" && -n "${OUT_MD:-}" ]]; then
	TMP_CSV="$(mktemp -t polap-table-data.XXXXXX.csv)"
	OUT_CSV="$TMP_CSV"
fi

python3 - "$MANIFEST" "$OUT_TSV" "${OUT_CSV:-}" "$SORT" <<'PY'
import sys, json, csv, math

def as_int_str(x):
    if x in ("", None): return ""
    try:
        f=float(x)
        return str(int(f)) if f.is_integer() else str(f)
    except: return str(x)

def gb_fmt(x):
    if x in ("", None): return ""
    try:
        v=float(x)/1e9
        return f"{v:.1f} Gb"
    except: return ""

def get(d,*ks):
    cur=d
    for k in ks:
        if not isinstance(cur,dict) or k not in cur: return None
        cur=cur[k]
    return cur

def main():
    if len(sys.argv)<3: sys.exit(2)
    manifest, out_tsv = sys.argv[1], sys.argv[2]
    out_csv = sys.argv[3] if len(sys.argv)>=4 and sys.argv[3] else None
    sortkey = sys.argv[4] if len(sys.argv)>=5 else "code2"

    with open(manifest,encoding="utf-8") as fh:
        man=json.load(fh)

    rows=[]
    for it in man.get("items",[]):
        code2   = it.get("code2","")
        species = it.get("species","")
        long_sra   = get(it,"data","sra_id") or ""
        long_total = get(it,"data","total_bases")
        short_sra  = get(it,"short1data","sra_id") or ""
        short1_tot = get(it,"short1data","total_bases")
        rows.append({
            "code2": code2, "species": species,
            "long_sra": long_sra, "long_total_bases": as_int_str(long_total),
            "short_sra": short_sra, "short1_total_bases": as_int_str(short1_tot),
            "long_total_bases_gb": gb_fmt(long_total),
            "short1_total_bases_gb": gb_fmt(short1_tot),
        })

    if sortkey=="species":
        rows.sort(key=lambda r:(r["species"], r["code2"]))
    else:
        rows.sort(key=lambda r:(r["code2"], r["species"]))

    # TSV
    with open(out_tsv,"w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["code2","species","long_sra","long_total_bases","short_sra","short1_total_bases"])
        for r in rows:
            w.writerow([r["code2"],r["species"],r["long_sra"],r["long_total_bases"],
                        r["short_sra"],r["short1_total_bases"]])

    # CSV (Gb) if requested
    if out_csv:
        with open(out_csv,"w",encoding="utf-8",newline="") as f:
            w=csv.writer(f)
            w.writerow(["code2","species","long_sra","long_total_bases_gb","short_sra","short1_total_bases_gb"])
            for r in rows:
                w.writerow([r["code2"],r["species"],r["long_sra"],r["long_total_bases_gb"],
                            r["short_sra"],r["short1_total_bases_gb"]])

if __name__=="__main__":
    main()
PY

echo "[OK] Wrote TSV: $OUT_TSV"
[[ -n "${OUT_CSV:-}" ]] && echo "[OK] Wrote CSV: $OUT_CSV"

# ---------- 2) Markdown: prefer qsv (+ inline awk) ; fallback to Python ----------
if [[ -n "${OUT_MD:-}" ]]; then
	if have qsv && [[ -n "${OUT_CSV:-}" && -s "$OUT_CSV" ]]; then
		# Use qsv to ensure clean CSV, then inline awk to render Markdown.
		# (qsv prints to stdout by default)
		echo "[INFO] Markdown via qsv + awk from CSV: $OUT_CSV"
		qsv input "$OUT_CSV" |
			awk -F',' '
      BEGIN {
        OFS="|"
      }
      NR==1 {
        # header
        printf("| %s | %s | %s | %s | %s | %s |\n", $1,$2,$3,$4,$5,$6)
        printf("| --- | --- | --- | --- | --- | --- |\n")
        next
      }
      {
        printf("| %s | %s | %s | %s | %s | %s |\n", $1,$2,$3,$4,$5,$6)
      }
    ' >"$OUT_MD"
		echo "[OK] Wrote MD : $OUT_MD"
	else
		echo "[WARN] qsv not found (or no CSV). Using Python fallback MD." >&2
		python3 - "$OUT_TSV" "$OUT_MD" <<'PY'
import sys,csv
if len(sys.argv)!=3: sys.exit(2)
tsv,out=sys.argv[1],sys.argv[2]
rows=[]
with open(tsv,encoding="utf-8") as f:
    r=csv.reader(f,delimiter="\t")
    hdr=next(r)
    idx={k:i for i,k in enumerate(hdr)}
    # header
    md=[]
    md.append("| code2 | species | long_sra | long_total_bases | short_sra | short1_total_bases |")
    md.append("| --- | --- | --- | --- | --- | --- |")
    for row in r:
        # Convert bases to Gb for MD presentation
        def gb(x):
            try: 
                v=float(x)/1e9
                return f"{v:.1f} Gb"
            except: return ""
        lt=gb(row[idx["long_total_bases"]])
        s1=gb(row[idx["short1_total_bases"]])
        md.append(f"| {row[idx['code2']]} | {row[idx['species']]} | {row[idx['long_sra']]} | {lt} | {row[idx['short_sra']]} | {s1} |")
with open(out,"w",encoding="utf-8") as fo:
    fo.write("\n".join(md)+"\n")
PY
		echo "[OK] Wrote MD : $OUT_MD"
	fi
fi

# Cleanup temp CSV if created
if [[ -n "${TMP_CSV:-}" && -f "$TMP_CSV" ]]; then
	rm -f "$TMP_CSV"
fi
