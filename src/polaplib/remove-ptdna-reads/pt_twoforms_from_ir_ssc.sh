#!/usr/bin/env bash
# pt_twoforms_from_ir_ssc.sh  v0.1.0
# INPUT : -r pt_ref.fa -o outdir -t threads
# OUTPUT: pt_isomerA.fa, pt_isomerB.fa, ir.pick.tsv
set -euo pipefail
pt_ref=""
outdir="pt_forms"
threads=8
usage() { echo "Usage: $0 -r ref.fa -o outdir [-t threads]"; }
while [[ $# -gt 0 ]]; do
	case "$1" in
	-r)
		pt_ref="$2"
		shift 2
		;;
	-o)
		outdir="$2"
		shift 2
		;;
	-t)
		threads="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] unknown: $1" >&2
		usage
		exit 2
		;;
	esac
done
[[ -s "$pt_ref" ]] || {
	echo "[ERR] -r pt_ref.fa missing" >&2
	exit 2
}
mkdir -p "$outdir"
refbase="${outdir}/pt.ref"
cp -f "$pt_ref" "${refbase}.fa"
nucmer --maxmatch -t "$threads" -p "${outdir}/pt_ir" \
	"${refbase}.fa" "${refbase}.fa"
delta-filter -r "${outdir}/pt_ir.delta" >"${outdir}/pt_ir.r.delta"
show-coords -THcl "${outdir}/pt_ir.r.delta" >"${outdir}/pt_ir.coords"

python - "$outdir/pt_ir.coords" "$outdir/ir.pick.tsv" <<'PY'
import sys
co, out = sys.argv[1], sys.argv[2]
best=(0,None)
with open(co) as f:
    for line in f:
        if not line.strip(): continue
        x=line.strip().split()
        s1,e1,s2,e2=int(x[0]),int(x[1]),int(x[2]),int(x[3])
        # x[12] is query label (1-based, same seq); skip near-diagonal
        if abs(s1 - s2) < 200: continue
        length=abs(e1-s1)+1
        if length>best[0]: best=(length,(s1,e1,s2,e2))
with open(out,'w') as w:
    if best[1]:
        s1,e1,s2,e2=best[1]
        w.write("s1\te1\ts2\te2\tlen\n")
        w.write(f"{s1}\t{e1}\t{s2}\t{e2}\t{best[0]}\n")
PY

[[ -s "${outdir}/ir.pick.tsv" ]] || {
	echo "[ERR] IR pick failed"
	exit 3
}
read s1 e1 s2 e2 len < <(awk 'NR==2{print $1,$2,$3,$4,$5}' \
	"${outdir}/ir.pick.tsv")
cat "${refbase}.fa" "${refbase}.fa" >"${outdir}/pt.double.fa"
extract() {
	local fa="$1" beg="$2" end="$3" lab="$4" outfa="$5"
	if ((beg <= end)); then
		seqkit subseq -r "${beg}:${end}" -w 0 -u "$fa" |
			sed "1s/>.*/>${lab}/" >"$outfa"
	else
		seqkit subseq -r "${beg}:-1" -w 0 -u "$fa" |
			sed "1s/>.*/>${lab}/" >"$outfa"
		seqkit subseq -r "1:${end}" -w 0 -u "$fa" |
			sed "1s/>.*/>${lab}_cont/" >>"$outfa"
	fi
}
IRa_beg=$s1
IRa_end=$e1
IRb_beg=$s2
IRb_end=$e2
extract "${outdir}/pt.double.fa" "$((IRb_end + 1))" "$((IRa_beg - 1))" "LSC" \
	"${outdir}/LSC.fa"
extract "${outdir}/pt.double.fa" "$((IRa_end + 1))" "$((IRb_beg - 1))" "SSC" \
	"${outdir}/SSC.fa"
extract "${outdir}/pt.double.fa" "$IRa_beg" "$IRa_end" "IRa" \
	"${outdir}/IRa.fa"
extract "${outdir}/pt.double.fa" "$IRb_beg" "$IRb_end" "IRb" \
	"${outdir}/IRb.fa"
# isomer A: LSC-IRb-SSC-IRa
cat "${outdir}/LSC.fa" "${outdir}/IRb.fa" "${outdir}/SSC.fa" \
	"${outdir}/IRa.fa" | seqkit seq -w 0 >"${outdir}/pt_isomerA.fa"
# isomer B: LSC-IRb-rev(SSC)-IRa
seqkit seq -r -p -w 0 "${outdir}/SSC.fa" >"${outdir}/SSC_rc.fa"
cat "${outdir}/LSC.fa" "${outdir}/IRb.fa" "${outdir}/SSC_rc.fa" \
	"${outdir}/IRa.fa" | seqkit seq -w 0 >"${outdir}/pt_isomerB.fa"
echo "${outdir}/pt_isomerA.fa"
echo "${outdir}/pt_isomerB.fa"
