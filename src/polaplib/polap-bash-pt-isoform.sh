#!/usr/bin/env bash
# polap-bash-pt-isoform.sh
# Version: v0.7.0 (nucmer-based, origin-rotation normalization)
#
# Build plastid isoforms from a circular ptDNA reference by locating IRa/IRb
# with MUMmer4 self-alignment (nucmer → delta-filter → show-coords),
# then **rotate the circle** so the origin lies in LSC, and slice segments
# from a **single-record doubled** FASTA by (start,len).
# If IR detection fails, emit a SINGLE form (A only).
#
# Outputs (in --outdir):
#   pt_isomerA.fa
#   [pt_isomerB.fa]  (only when IR detection succeeds)
#   pt.rot.fa        (reference rotated to LSC start)
#   pt.double.fa     (single-record doubled from rotated)
#   ir.pick.tsv      (header + one row if success; header only if fallback)
#
# Requirements: nucmer, delta-filter, show-coords, seqkit, python3
#
# Example:
#   polap-bash-pt-isoform.sh -r plastid.fa -o pt_forms -t 16
#
set -euo pipefail

# ────────────────────────────────────────────────────────────────────
# Logging
# ────────────────────────────────────────────────────────────────────
VERBOSE=1
QUIET=0
LOG_FILE=""
note() {
	[[ $QUIET -eq 1 ]] && return 0
	local ts src ln
	ts="$(date +'%F %T')"
	src="${BASH_SOURCE[1]##*/}"
	ln="${BASH_LINENO[0]}"
	printf "[%s][%s:%s] %s\n" "$ts" "$src" "$ln" "$*" |
		{ if [[ -n "$LOG_FILE" ]]; then tee -a "$LOG_FILE"; else cat; fi; } >&2
}

# ────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────
pt_ref=""
outdir="pt_forms"
threads=16
min_ir_len=5000
max_ir_len=80000
len_rel_diff=0.20
near_diag_pad=200
usage() {
	cat <<EOF
polap-bash-pt-isoform.sh v0.7.0
Usage:
  $0 -r pt_ref.fa -o outdir [-t threads]
     [--min-ir-len N] [--max-ir-len N] [--len-rel-diff F]
     [--near-diag-pad N] [-q]

Options:
  -r FILE           plastid reference (circular, FASTA)   [required]
  -o DIR            output directory                      [${outdir}]
  -t INT            threads for nucmer                    [${threads}]
  --min-ir-len N    minimal IR length                     [${min_ir_len}]
  --max-ir-len N    maximal IR length                     [${max_ir_len}]
  --len-rel-diff F  max |lenA-lenB|/mean(lenA,lenB)       [${len_rel_diff}]
  --near-diag-pad N diagonal padding to ignore            [${near_diag_pad}]
  -q                 quiet
  -h                 help
EOF
}
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
	--min-ir-len)
		min_ir_len="$2"
		shift 2
		;;
	--max-ir-len)
		max_ir_len="$2"
		shift 2
		;;
	--len-rel-diff)
		len_rel_diff="$2"
		shift 2
		;;
	--near-diag-pad)
		near_diag_pad="$2"
		shift 2
		;;
	-q)
		QUIET=1
		shift
		;;
	-h)
		usage
		exit 0
		;;
	*)
		echo "[ERR] unknown arg: $1" >&2
		usage
		exit 2
		;;
	esac
done
[[ -s "$pt_ref" ]] || {
	echo "[ERR] -r pt_ref.fa missing" >&2
	exit 2
}
mkdir -p "$outdir" "$outdir/log"
LOG_FILE="$outdir/log/pt_isoform.log"
[[ $QUIET -eq 1 ]] && VERBOSE=0

# ────────────────────────────────────────────────────────────────────
# Paths & tools
# ────────────────────────────────────────────────────────────────────
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPDIR="${_POLAPLIB_DIR}/pt-isoform"
mkdir -p "$HELPDIR"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "[ERR] tool not in PATH: $1" >&2
	exit 127
}; }
need nucmer
need delta-filter
need show-coords
need seqkit
need python3

# ────────────────────────────────────────────────────────────────────
# Helpers (Python): pick IR from show-coords; emit/rotate utils
# ────────────────────────────────────────────────────────────────────
PY_PICK="${HELPDIR}/py_pick_ir_from_mummer.py"
if [[ ! -s "$PY_PICK" ]]; then
	cat >"$PY_PICK" <<'PY'
#!/usr/bin/env python3
# py_pick_ir_from_mummer.py  v0.2.1 (wrap-aware, origin-safe pairing)
import sys, argparse
def parse_coords(path):
    out=[]
    with open(path) as f:
        for ln in f:
            ln=ln.strip()
            if not ln: continue
            xs=ln.split()
            if len(xs) < 12: continue
            try:
                s1=int(xs[0]); e1=int(xs[1])
                s2=int(xs[2]); e2=int(xs[3])
                l1=int(xs[4]); l2=int(xs[5]); idy=float(xs[6])
            except ValueError:
                continue
            qn=xs[-2]; tn=xs[-1]
            out.append((s1,e1,s2,e2,l1,l2,idy,qn,tn))
    return out
def is_inv(s1,e1,s2,e2): return (s1<=e1 and s2>=e2) or (s1>=e1 and s2<=e2)
def near_diag(s1,e1,s2,e2,p): return abs(s1-s2)<p and abs(e1-e2)<p
def within(x,a,b): return a <= x <= b
def norm1(x,L): t=(x-1)%L; return t+1
def make_variant(a1,a2,L):
    if a1<=a2: s_d,e_d=a1,a2
    else:      s_d,e_d=a2,a1
    len_d=e_d-s_d+1
    len_w=L-len_d
    s_w=norm1(e_d-(len_w-1),L)
    e_w=norm1(s_w+(len_w-1),L)
    return (s_d,e_d,len_d),(s_w,e_w,len_w)
def circ_overlap(sA,lA,sB,lB,L):
    a0=(sA-1)%L; b0=(sB-1)%L
    def segs(s,l):
        e=(s+l-1)%L
        return [(s,e)] if s<=e else [(s,L-1),(0,e)]
    A=segs(a0,lA); B=segs(b0,lB)
    ov=0
    for x0,x1 in A:
        for y0,y1 in B:
            ov += max(0, min(x1,y1)-max(x0,y0)+1)
    return ov
def pick_pair(rows,L,minL,maxL,reld,diag):
    C=[]
    for s1,e1,s2,e2,l1,l2,idy,qn,tn in rows:
        if qn!=tn: continue
        if not is_inv(s1,e1,s2,e2): continue
        if near_diag(s1,e1,s2,e2,diag): continue
        A_d,A_w=make_variant(s1,e1,L)
        B_d,B_w=make_variant(s2,e2,L)
        C.append((A_d,A_w,B_d,B_w))
    best=None; best_score=(-1,-1)
    for A_d,A_w,B_d,B_w in C:
        for Ad in (A_d,A_w):
            sA,eA,lA=Ad
            if not within(lA,minL,maxL): continue
            for Bd in (B_d,B_w):
                sB,eB,lB=Bd
                if not within(lB,minL,maxL): continue
                mu=(lA+lB)/2.0
                if mu==0: continue
                if abs(lA-lB)/mu>reld: continue
                ov=circ_overlap(sA,lA,sB,lB,L)
                if ov>1000: continue
                score=((lA+lB)/2.0,-ov)
                if score>best_score:
                    best_score=score; best=(sA,eA,lA,sB,eB,lB)
    return best
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("coords"); ap.add_argument("--L",type=int,required=True)
    ap.add_argument("--min_ir_len",type=int,default=5000)
    ap.add_argument("--max_ir_len",type=int,default=80000)
    ap.add_argument("--len_rel_diff",type=float,default=0.20)
    ap.add_argument("--near_diag_pad",type=int,default=200)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    rows=parse_coords(a.coords)
    best=pick_pair(rows,a.L,a.min_ir_len,a.max_ir_len,a.len_rel_diff,a.near_diag_pad)
    with open(a.out,'w') as w:
        w.write("s1\te1\ts2\te2\tlen\n")
        if not best: return
        sA,eA,lA,sB,eB,lB=best
        # emit with A first by start
        if sA<=sB:
            w.write(f"{sA}\t{eA}\t{sB}\t{eB}\t{int((lA+lB)/2.0)}\n")
        else:
            w.write(f"{sB}\t{eB}\t{sA}\t{eA}\t{int((lA+lB)/2.0)}\n")
if __name__=="__main__":
    main()
PY
	chmod +x "$PY_PICK"
fi

PY_EMIT="${HELPDIR}/py_emit_ir_fields.py"
if [[ ! -s "$PY_EMIT" ]]; then
	cat >"$PY_EMIT" <<'PY'
#!/usr/bin/env python3
# py_emit_ir_fields.py  v0.1.0
import sys
p=sys.argv[1] if len(sys.argv)>1 else "-"
with open(p) as f:
    _=next(f,""); row=next(f,"").strip()
    if not row: exit(3)
    s=row.split(); print(s[0],s[1],s[2],s[3],s[4])
PY
	chmod +x "$PY_EMIT"
fi

# Rotate first FASTA record so new origin is at given 1-based position
PY_ROT="${HELPDIR}/py_rotate_fa_to_start.py"
if [[ ! -s "$PY_ROT" ]]; then
	cat >"$PY_ROT" <<'PY'
#!/usr/bin/env python3
# py_rotate_fa_to_start.py  v0.1.0
# Usage: py_rotate_fa_to_start.py in.fa start out.fa
import sys
if len(sys.argv)<4: sys.exit(2)
inp, start, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
name=None; seq=[]
with open(inp) as f:
    for ln in f:
        if ln.startswith('>'):
            if name is None: name=ln[1:].strip().split()[0]
            else: break
        else: seq.append(ln.strip())
S=''.join(seq)
L=len(S); start=((start-1)%L)+1
rot=S[start-1:]+S[:start-1]
with open(outp,'w') as w:
    w.write('>pt.rot\n'); w.write(rot+'\n')
PY
	chmod +x "$PY_ROT"
fi

# ────────────────────────────────────────────────────────────────────
# Prepare original reference
# ────────────────────────────────────────────────────────────────────
refbase="${outdir}/pt.ref"
note "Copy reference → ${refbase}.fa"
cp -f "$pt_ref" "${refbase}.fa"
L=$(seqkit fx2tab -l -n "${refbase}.fa" | cut -f2)
[[ "$L" =~ ^[0-9]+$ ]] || {
	echo "[ERR] failed to get genome length" >&2
	exit 3
}

# ────────────────────────────────────────────────────────────────────
# NUCmer self-alignment to find inverted repeats
# ────────────────────────────────────────────────────────────────────
note "Run nucmer self-alignment (repeats/inversions preserved)"
nucmer --maxmatch -c 100 -l 20 -p "${outdir}/pt_ir" \
	"${refbase}.fa" "${refbase}.fa" >>"$LOG_FILE" 2>&1
note "delta-filter by len/id (keep repeats)"
delta-filter -i 90 -l 500 "${outdir}/pt_ir.delta" >"${outdir}/pt_ir.keep.delta"
note "show-coords -THcl → coords"
show-coords -THcl "${outdir}/pt_ir.keep.delta" >"${outdir}/pt_ir.coords"

note "Pick IR (wrap-aware) → ir.pick.tsv"
python3 "$PY_PICK" "${outdir}/pt_ir.coords" \
	--L "$L" --min_ir_len "$min_ir_len" --max_ir_len "$max_ir_len" \
	--len_rel_diff "$len_rel_diff" --near_diag_pad "$near_diag_pad" \
	--out "${outdir}/ir.pick.tsv" || true

have_twoforms=1
if ! python3 "$PY_EMIT" "${outdir}/ir.pick.tsv" >/dev/null 2>&1; then
	have_twoforms=0
	note "WARN IR detection failed; fallback to single form"
fi

# ────────────────────────────────────────────────────────────────────
# Slicing helpers
# ────────────────────────────────────────────────────────────────────
norm1() {
	local x="$1" L="$2"
	local t=$(((x - 1) % L))
	echo $((t + 1))
}
# slice_by_start_len fa L start len label outfa  (uses single-record doubled)
slice_by_start_len() {
	local fa="$1" L="$2" start="$3" len="$4" lab="$5" outfa="$6"
	start=$(norm1 "$start" "$L")
	local stop=$((start + len - 1))
	if ((stop <= L)); then
		seqkit subseq -r "${start}:${stop}" -w 0 "$fa" |
			awk 'BEGIN{lab="'"$lab"'"} /^>/{print ">"lab; next} {print}' >"$outfa"
	else
		{
			echo ">$lab"
			seqkit subseq -r "${start}:${L}" -w 0 "$fa" | awk 'NR>1{print}'
			seqkit subseq -r "1:$((stop - L))" -w 0 "$fa" | awk 'NR>1{print}'
		} >"$outfa"
	fi
}

# ────────────────────────────────────────────────────────────────────
# Rotate origin to LSC start, rebuild doubled, and slice by (start,len)
# ────────────────────────────────────────────────────────────────────
if ((have_twoforms)); then
	read s1 e1 s2 e2 lenIR < <(python3 "$PY_EMIT" "${outdir}/ir.pick.tsv")
	# Ensure s<=e and compute lengths
	if ((e1 < s1)); then
		t=$s1
		s1=$e1
		e1=$t
	fi
	if ((e2 < s2)); then
		t=$s2
		s2=$e2
		e2=$t
	fi
	L_IRA=$((e1 - s1 + 1))
	L_IRB=$((e2 - s2 + 1))

	# LSC starts immediately after IRb: new_origin = IRb_end + 1
	LSC_START=$(((e2 % L) + 1))
	note "Rotate origin to LSC start: ${LSC_START}"

	# Rotate original reference so 1 == LSC_START
	ROT_FA="${outdir}/pt.rot.fa"
	python3 "$PY_ROT" "${refbase}.fa" "$LSC_START" "$ROT_FA"

	# Function to shift a coordinate after rotation by shift=(LSC_START-1)
	shift_coord() { # new = ((old - shift - 1) mod L) + 1
		local old="$1" shift="$2" L="$3"
		local val=$(((old - shift - 1) % L))
		((val < 0)) && val=$((val + L))
		echo $((val + 1))
	}
	SHIFT=$((LSC_START - 1))
	s1r=$(shift_coord "$s1" "$SHIFT" "$L")
	e1r=$(shift_coord "$e1" "$SHIFT" "$L")
	s2r=$(shift_coord "$s2" "$SHIFT" "$L")
	e2r=$(shift_coord "$e2" "$SHIFT" "$L")
	# Normalize to increasing ranges on rotated axis
	if ((e1r < s1r)); then
		t=$s1r
		s1r=$e1r
		e1r=$t
	fi
	if ((e2r < s2r)); then
		t=$s2r
		s2r=$e2r
		e2r=$t
	fi

	# Build single-record doubled from rotated
	note "Build doubled (single-record) from rotated ref"
	DBL="${outdir}/pt.double.fa"
	python - "$ROT_FA" "$DBL" <<'PY'
import sys
inp,outp=sys.argv[1],sys.argv[2]
name=None; seq=[]
with open(inp) as f:
    for ln in f:
        if ln.startswith('>'):
            if name is None: name=ln[1:].strip().split()[0]
            else: break
        else: seq.append(ln.strip())
S=''.join(seq)
with open(outp,'w') as w:
    w.write('>pt.double\n'); w.write(S+S+'\n')
PY

	# Now segments are contiguous (no wrap) by design
	S_IRA="$s1r"
	L_IRA="$L_IRA"
	S_IRB="$s2r"
	L_IRB="$L_IRB"
	# LSC starts at 1 on rotated ref
	S_LSC=1
	L_LSC=$((s1r - 1)) # from 1 up to just before IRa
	# SSC lies between IRa_end and IRb_start
	S_SSC=$((e1r + 1))
	L_SSC=$((s2r - S_SSC))

	# If any length is negative/zero due to edge cases, recompute robustly:
	fix_len() {
		local x="$1"
		((x < 0)) && echo 0 || echo "$x"
	}
	L_LSC=$(fix_len "$L_LSC")
	L_SSC=$(fix_len "$L_SSC")

	# Slice
	note "Slice IRa (start=$S_IRA len=$L_IRA)"
	slice_by_start_len "$DBL" "$L" "$S_IRA" "$L_IRA" "IRa" "${outdir}/IRa.fa"
	note "Slice IRb (start=$S_IRB len=$L_IRB)"
	slice_by_start_len "$DBL" "$L" "$S_IRB" "$L_IRB" "IRb" "${outdir}/IRb.fa"
	note "Slice LSC (start=$S_LSC len=$L_LSC)"
	slice_by_start_len "$DBL" "$L" "$S_LSC" "$L_LSC" "LSC" "${outdir}/LSC.fa"
	note "Slice SSC (start=$S_SSC len=$L_SSC)"
	slice_by_start_len "$DBL" "$L" "$S_SSC" "$L_SSC" "SSC" "${outdir}/SSC.fa"

	# Build isomers as single-record FASTAs
	note "Build isomer A: LSC – IRb – SSC – IRa"
	{
		echo ">pt_isomerA"
		seqkit seq -w 0 "${outdir}/LSC.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/IRb.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/SSC.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/IRa.fa" | awk 'NR>1{print}'
	} >"${outdir}/pt_isomerA.fa"

	note "Build isomer B: LSC – IRb – rev(SSC) – IRa"
	seqkit seq -r -p -w 0 "${outdir}/SSC.fa" >"${outdir}/SSC_rc.fa"
	{
		echo ">pt_isomerB"
		seqkit seq -w 0 "${outdir}/LSC.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/IRb.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/SSC_rc.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/IRa.fa" | awk 'NR>1{print}'
	} >"${outdir}/pt_isomerB.fa"

	if [[ $VERBOSE -gt 0 ]]; then
		echo "[SANITY] lengths:" \
			"IRa=$L_IRA IRb=$L_IRB LSC=$L_LSC SSC=$L_SSC  (L=$L)" >&2
		seqkit stats -T "${outdir}/pt_isomerA.fa" "${outdir}/pt_isomerB.fa" || true
	fi

else
	# Fallback
	note "Fallback: single form A from input"
	cp -f "${refbase}.fa" "${outdir}/pt_isomerA.fa"
	echo -e "s1\te1\ts2\te2\tlen" >"${outdir}/ir.pick.tsv"
	# Build rotated = original (no IR known)
	cp -f "${refbase}.fa" "${outdir}/pt.rot.fa"
	# Also build doubled from rotated
	{
		echo ">pt.double"
		seqkit seq -w 0 "${outdir}/pt.rot.fa" | awk 'NR>1{print}'
		seqkit seq -w 0 "${outdir}/pt.rot.fa" | awk 'NR>1{print}'
	} \
		>"${outdir}/pt.double.fa"
fi
