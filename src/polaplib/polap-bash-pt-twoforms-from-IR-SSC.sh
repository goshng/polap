#!/usr/bin/env bash
# polap-bash-pt-twoforms-from-IR-SSC.sh  v0.2.0
# INPUT : -r <pt_ref.fasta>   (complete plastid circle)
# OUTPUT: <outdir>/pt_isomerA.fa, <outdir>/pt_isomerB.fa (+ IR/SSC coords tsv)

set -euo pipefail
rscript() { Rscript --vanilla "$@"; }

usage() {
	cat <<EOF
Usage: $0 -i pt_ref.fa -o outdir [-t threads]
EOF
}

pt_ref=""
outdir="pt_forms"
threads=8
while [[ $# -gt 0 ]]; do
	case "$1" in
	-i)
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
		echo "[ERROR] unknown arg: $1" >&2
		usage
		exit 2
		;;
	esac
done
[[ -s "$pt_ref" ]] || {
	echo "[ERROR] missing -r pt_ref.fa" >&2
	exit 2
}
mkdir -p "$outdir"

# 1) find IRs with nucmer
refbase="${outdir}/pt.ref"
cp -f "$pt_ref" "${refbase}.fa"
nucmer --maxmatch -p "${outdir}/pt_ir" "${refbase}.fa" "${refbase}.fa"
delta-filter -r "${outdir}/pt_ir.delta" >"${outdir}/pt_ir.r.delta"
show-coords -THcl "${outdir}/pt_ir.r.delta" >"${outdir}/pt_ir.coords"

# 2) pick the longest self-match off diagonal (IRa/IRb)
python - "$outdir/pt_ir.coords" "$outdir/ir.pick.tsv" <<'PY'
import sys
co, out = sys.argv[1], sys.argv[2]
best=(0,None)
with open(co) as f:
    for line in f:
        if not line.strip(): continue
        toks=line.strip().split()
        s1,e1,s2,e2,iden,qlen,qlab=[int(toks[0]),int(toks[1]),int(toks[2]),int(toks[3]),float(toks[6]),int(toks[7]),toks[12]]
        if qlab!='1' or abs(s1-s2)<200: continue       # skip diagonal, same copy
        length=abs(e1-s1)+1
        if length>best[0]: best=(length,(s1,e1,s2,e2))
with open(out,'w') as w:
    if best[1]: w.write("s1\te1\ts2\te2\tlen\n"+'\t'.join(map(str,best[1]))+f"\t{best[0]}\n")
PY
[[ -s "${outdir}/ir.pick.tsv" ]] || {
	echo "[ERROR] IR detection failed" >&2
	exit 3
}

# 3) compute LSC/SSC ranges; build two isomers with seqkit
read s1 e1 s2 e2 len < <(awk 'NR==2{print $1,$2,$3,$4,$5}' "${outdir}/ir.pick.tsv")
# make the reference doubled so we can rotate across origin
cat "${refbase}.fa" "${refbase}.fa" >"${outdir}/pt.double.fa"

# helper to extract (closed) intervals
extract() {
	local fa="$1" beg="$2" end="$3" lab="$4" outfa="$5"
	if ((beg <= end)); then
		seqkit subseq -r "${beg}:${end}" -w 0 -u "$fa" | sed "1s/>.*/>${lab}/" >"$outfa"
	else
		# wrap
		seqkit subseq -r "${beg}:-1" -w 0 -u "$fa" | sed "1s/>.*/>${lab}/" >"$outfa"
		seqkit subseq -r "1:${end}" -w 0 -u "$fa" | sed "1s/>.*/>${lab}_cont/" >>"$outfa"
	fi
}

# IRa as [s1..e1], IRb as [s2..e2] on the doubled circle
IRa_beg=$s1
IRa_end=$e1
IRb_beg=$s2
IRb_end=$e2
# LSC = between IRb_end -> IRa_beg ; SSC = between IRa_end -> IRb_beg
extract "${outdir}/pt.double.fa" "$((IRb_end + 1))" "$((IRa_beg - 1))" "LSC" "${outdir}/LSC.fa"
extract "${outdir}/pt.double.fa" "$((IRa_end + 1))" "$((IRb_beg - 1))" "SSC" "${outdir}/SSC.fa"
extract "${outdir}/pt.double.fa" "$IRa_beg" "$IRa_end" "IRa" "${outdir}/IRa.fa"
extract "${outdir}/pt.double.fa" "$IRb_beg" "$IRb_end" "IRb" "${outdir}/IRb.fa"

# isomer A: LSC-IRb-SSC-IRa
cat "${outdir}/LSC.fa" "${outdir}/IRb.fa" "${outdir}/SSC.fa" "${outdir}/IRa.fa" |
	seqkit seq -w 0 >"${outdir}/pt_isomerA.fa"

# isomer B: flip SSC
seqkit seq -r -p -w 0 "${outdir}/SSC.fa" >"${outdir}/SSC_rc.fa"
cat "${outdir}/LSC.fa" "${outdir}/IRb.fa" "${outdir}/SSC_rc.fa" "${outdir}/IRa.fa" |
	seqkit seq -w 0 >"${outdir}/pt_isomerB.fa"

echo -e "# wrote:\n${outdir}/pt_isomerA.fa\n${outdir}/pt_isomerB.fa"
