#!/usr/bin/env bash
# polap-bash-polish.sh
# POLAP "polish" subcommand: ONT-only (Racon) or Hybrid (FMLRC2/Polypolish),
# optionally preserve Flye GFA topology by extract → polish → reinject.
# Version: v0.5.1
# License: GPL-3.0+
set -euo pipefail
IFS=$'\n\t'

_show_usage() {
	cat <<'USAGE'
Usage:
  polap-bash-polish.sh polish <assembly.fa|assembly.gfa>
                         [--ont ONT.fq[.gz]] [--sr1 R1.fq[.gz] --sr2 R2.fq[.gz]]
                         [--mode auto|ont|hybrid]
                         [--gfa-preserve]
                         [--polisher auto|fmlrc2|polypolish]
                         [--racon-rounds N] [--threads N]
                         [--minimap-preset map-ont|map-hifi]
                         [--speed normal|fast|turbo]
                         [--recruit auto|draft|bait|off] [--bait BAITS.fa]
                         [--ont-cap X] [--ont-cap-method rasusa|seqtk]
                         [--eval-merqury] [--kmer K] [--prefix NAME]
                         [--ab-test] [--qc-report]

Notes:
  • DEFAULT prioritizes effectiveness: uses ALL reads (recruit=off, ont-cap=0) unless you opt in.
  • If you enable capping, selection for capping ignores MAPQ and keeps reads with:
      aligned_len > 3000 AND identity > 0.75  (identity from PAF: nmatch/alen), then caps to X× (default 200× in A/B).
  • Presets (affect ONT path only; explicit flags still win):
      normal: ont-cap unchanged, rounds=3, racon -w 1000, minimap2 -K 2g
      fast:   ont-cap=80× (if unset), rounds=2, racon -w 1200, minimap2 -K 2g
      turbo:  ont-cap=60× (if unset), rounds=1, racon -w 1500, minimap2 -K 4g

Examples:
  ONT-only (all reads):
    polap-bash-polish.sh polish pt.fa --ont ont.fq.gz --speed normal --prefix pt

  ONT-only (fast path, if you allow capping):
    polap-bash-polish.sh polish pt.fa --ont ont.fq.gz --speed fast --recruit draft --prefix pt_fast

  Hybrid (FMLRC2 default):
    polap-bash-polish.sh polish mt.fa --sr1 R1.fq.gz --sr2 R2.fq.gz --prefix mt

  Hybrid with Polypolish + Merqury:
    polap-bash-polish.sh polish mt.fa --sr1 R1.fq.gz --sr2 R2.fq.gz --polisher polypolish --eval-merqury --kmer 31 --prefix mt_pp

  Preserve Flye GFA topology:
    polap-bash-polish.sh polish flye/assembly_graph.gfa --ont ont.fq.gz --gfa-preserve --prefix mito_graph
USAGE
}

# ----------------- subcommand gate -----------------
_subcmd="${1:-}"
shift || true
[[ -z "${_subcmd:-}" || "${_subcmd}" != "polish" ]] && {
	_show_usage
	exit 1
}

# ----------------- defaults -----------------
asm=""
ont=""
sr1=""
sr2=""
mode="auto" # auto|ont|hybrid
gfa_preserve=0
polisher="auto" # auto|fmlrc2|polypolish
racon_rounds=3
threads=16
minimap_preset="map-ont" # map-ont|map-hifi
speed="normal"           # normal|fast|turbo
minimap_K="2g"

# recruitment/cap (defaults → ALL reads, NO cap)
recruit="off" # off|draft|bait|auto
bait_fa=""
ont_cap=0               # 0 = do not cap unless user/preset sets it
ont_cap_method="rasusa" # rasusa|seqtk

# QC / meta
eval_merqury=0
kmer=31
prefix="polap"
ab_test=0
qc_report=0

# ----------------- parse args -----------------
asm="${1:-}"
shift || true
[[ -z "${asm}" ]] && {
	_show_usage
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--ont)
		ont="$2"
		shift 2
		;;
	--sr1)
		sr1="$2"
		shift 2
		;;
	--sr2)
		sr2="$2"
		shift 2
		;;
	--mode)
		mode="$2"
		shift 2
		;;
	--gfa-preserve)
		gfa_preserve=1
		shift 1
		;;
	--polisher)
		polisher="$2"
		shift 2
		;;
	--racon-rounds)
		racon_rounds="$2"
		shift 2
		;;
	--threads)
		threads="$2"
		shift 2
		;;
	--minimap-preset)
		minimap_preset="$2"
		shift 2
		;;
	--speed)
		speed="$2"
		shift 2
		;;
	--recruit)
		recruit="$2"
		shift 2
		;;
	--bait)
		bait_fa="$2"
		shift 2
		;;
	--ont-cap)
		ont_cap="$2"
		shift 2
		;;
	--ont-cap-method)
		ont_cap_method="$2"
		shift 2
		;;
	--eval-merqury)
		eval_merqury=1
		shift 1
		;;
	--kmer)
		kmer="$2"
		shift 2
		;;
	--prefix)
		prefix="$2"
		shift 2
		;;
	--ab-test)
		ab_test=1
		shift 1
		;;
	--qc-report)
		qc_report=1
		shift 1
		;;
	-h | --help)
		_show_usage
		exit 0
		;;
	*)
		echo "Unknown argument: $1" >&2
		_show_usage
		exit 1
		;;
	esac
done

# ----------------- apply speed preset (only sets defaults; explicit flags win) -----------------
case "$speed" in
normal)
	:
	;;
fast)
	[[ "$ont_cap" -eq 0 ]] && ont_cap=80
	[[ "$racon_rounds" == "3" ]] && racon_rounds=2
	minimap_K="2g"
	;;
turbo)
	[[ "$ont_cap" -eq 0 ]] && ont_cap=60
	[[ "$racon_rounds" == "3" ]] && racon_rounds=1
	minimap_K="4g"
	;;
*) echo "[WARN] Unknown --speed '$speed' (using 'normal')" ;;
esac

# ----------------- setup -----------------
run="$(date +%Y%m%d_%H%M%S)"
odir="out/${prefix}.polish_${run}"
mkdir -p "$odir"
exec > >(tee -i "$odir/run.log") 2>&1

echo "[POLAP] assembly: $asm"
echo "[POLAP] ONT:      ${ont:-'(none)'}"
echo "[POLAP] SR1:      ${sr1:-'(none)'}"
echo "[POLAP] SR2:      ${sr2:-'(none)'}"
echo "[POLAP] mode:     $mode"
echo "[POLAP] polisher: $polisher"
echo "[POLAP] gfa-preserve: $gfa_preserve"
echo "[POLAP] speed:    $speed"
echo "[POLAP] racon-rounds: $racon_rounds"
echo "[POLAP] recruit:  $recruit   bait=${bait_fa:-'(none)'}"
echo "[POLAP] ont-cap:  $ont_cap ($ont_cap_method)   minimap: -x $minimap_preset -K $minimap_K"

[[ -f "$asm" ]] || {
	echo "ERROR: assembly not found: $asm" >&2
	exit 1
}
ext="${asm##*.}"
is_gfa=0
[[ "$ext" =~ ^([Gg][Ff][Aa])$ ]] && is_gfa=1

# Resolve mode
if [[ "$mode" == "auto" ]]; then
	if [[ -n "$sr1" && -n "$sr2" ]]; then
		mode="hybrid"
	elif [[ -n "$ont" ]]; then
		mode="ont"
	else
		echo "ERROR: auto mode found no usable reads; provide --ont or --sr1/--sr2" >&2
		exit 1
	fi
fi
[[ "$mode" == "ont" && -z "$ont" ]] && {
	echo "ERROR: --mode ont requires --ont" >&2
	exit 1
}
[[ "$mode" == "hybrid" && (-z "$sr1" || -z "$sr2") ]] && {
	echo "ERROR: --mode hybrid requires --sr1/--sr2" >&2
	exit 1
}

# ----------------- paths & helpers -----------------
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scripts_dir="$_here/scripts"
mkdir -p "$scripts_dir"

_quiet_zcat() {
	set +e
	zcat -f "$1" 2>/dev/null || cat "$1"
	set -e
}
_fa_len_total() { awk '/^>/{next}{L+=length($0)}END{print (L?L:1)}' "$1"; }
_sh() { printf "'%s'" "${1//\'/\'\\\'\'}"; }

# --- helper files (create if missing) ---
if [[ ! -s "$scripts_dir/gfa_inject_polished.py" ]]; then
	cat >"$scripts_dir/gfa_inject_polished.py" <<'PY'
#!/usr/bin/env python3
# Version: v0.5.1
import sys,gzip
gfa, polished_fa = sys.argv[1], sys.argv[2]
def open_auto(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
seq={}
with open_auto(polished_fa) as f:
    cur=None; buf=[]
    for line in f:
        if line.startswith('>'):
            if cur: seq[cur]=''.join(buf); buf=[]
            cur=line[1:].strip().split()[0]
        else:
            buf.append(line.strip())
    if cur: seq[cur]=''.join(buf)
with open_auto(gfa) as f:
    for line in f:
        if not line.strip(): continue
        if line[0]=='S':
            parts=line.rstrip('\n').split('\t')
            sid=parts[1]
            if sid in seq: parts[2]=seq[sid]
            print('\t'.join(parts))
        else:
            print(line, end='')
PY
	chmod +x "$scripts_dir/gfa_inject_polished.py"
fi

if [[ ! -s "$scripts_dir/fastq_n50.py" ]]; then
	cat >"$scripts_dir/fastq_n50.py" <<'PY'
#!/usr/bin/env python3
# Version: v0.5.1
import sys,gzip
from statistics import median
def openx(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
L=[]
with openx(sys.argv[1]) as f:
    i=0
    for line in f:
        i+=1
        if i%4==2: L.append(len(line.strip()))
L.sort()
tot=sum(L); half=tot/2; acc=0; n50=0
for l in L:
    acc+=l
    if acc>=half:
        n50=l; break
print(f"n={len(L)} total={tot} n50={n50} median={median(L) if L else 0}")
PY
	chmod +x "$scripts_dir/fastq_n50.py"
fi

if [[ ! -s "$scripts_dir/plot_merqury_qv.R" ]]; then
	cat >"$scripts_dir/plot_merqury_qv.R" <<'RS'
# Version: v0.5.1
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_merqury_qv.R <out.png> <A.qv> [B.qv...]")
out <- args[1]; qv_files <- args[-1]
labs <- sub("\\.qv$","",basename(qv_files))
qvs <- sapply(qv_files, function(f) as.numeric(readLines(f)[1]))
png(out, 800, 400, res=120)
barplot(qvs, names.arg=labs, ylab="QV (Phred)", main="Merqury QV")
dev.off()
RS
fi

if [[ ! -s "$scripts_dir/paf_stats.awk" ]]; then
	cat >"$scripts_dir/paf_stats.awk" <<'AWK'
# Version: v0.5.1
# Usage: awk -f paf_stats.awk in.paf > out.tsv
BEGIN{ FS=OFS="\t"; print "n_align","mean_identity","total_aligned" }
{
  nmatch=$10+0; alen=$11+0;
  if (alen>0){ sum_id+=(nmatch/alen) }
  sum_al+=alen; n++
}
END{
  mi = (n>0)? sum_id/n : 0.0;
  printf "%d\t%.6f\t%d\n", n, mi, sum_al
}
AWK
fi

# ----------------- stats -----------------
if [[ -n "${ont:-}" && -f "$ont" ]]; then
	echo "[STATS] ONT read length summary:"
	python3 "$scripts_dir/fastq_n50.py" "$ont" || true
fi
if [[ -n "${sr1:-}" && -f "$sr1" ]]; then
	echo "[STATS] Short reads R1 length summary:"
	python3 "$scripts_dir/fastq_n50.py" "$sr1" || true
fi

# ----------------- extract sequences from GFA (if requested) -----------------
target_fa="$asm"
if ((is_gfa == 1)); then
	if ((gfa_preserve == 1)); then
		echo "[POLAP] GFA preserve ON → extracting S-sequences"
		gfatools gfa2fa "$asm" >"$odir/segments.fa"
		target_fa="$odir/segments.fa"
	else
		echo "[WARN] Input is GFA without --gfa-preserve; polishing exported FASTA only."
		gfatools gfa2fa "$asm" >"$odir/export.fa"
		target_fa="$odir/export.fa"
	fi
fi

# ----------------- recruitment & capping helpers -----------------
_recruit_for_capping() {
	local target="$1" reads="$2" outfq="$3" preset="$4" threads="$5"
	local paf="$odir/recruit.paf" ids="$odir/recruit.ids"
	echo "[RECRUIT] target=$(basename "$target") preset=$preset (len>3000 & id>0.75; MAPQ ignored)"
	minimap2 -x "$preset" --secondary=no -K "$minimap_K" -t "$threads" "$target" "$reads" >"$paf"
	awk '$11+0>3000 {pid=$10/$11; if(pid>0.75) print $1}' "$paf" | sort -u >"$ids"
	local n
	n=$(wc -l <"$ids" || echo 0)
	if [[ "$n" -eq 0 ]]; then
		echo "[RECRUIT] no candidates; using full read set"
		cp -f "$reads" "$outfq"
	else
		echo "[RECRUIT] kept read ids: $n"
		seqkit grep -f "$ids" "$reads" | gzip -dc >"$outfq"
	fi
}

_cap_ont_to() {
	local draft="$1" reads="$2" targetx="$3" outfq="$4"
	if [[ "$targetx" -le 0 ]]; then
		cp -f "$reads" "$outfq"
		return
	fi
	local asm_len
	asm_len=$(_fa_len_total "$draft")
	if command -v rasusa >/dev/null 2>&1; then
		local bases=$((targetx * asm_len))
		echo "[CAP] rasusa to ~${targetx}x (bases=$bases)"
		rasusa -s 42 -i "$reads" -o "$outfq" -b "$bases"
	else
		local total
		total=$(_quiet_zcat "$reads" | awk 'NR%4==2{n+=length($0)}END{print (n?n:1)}')
		local p
		p=$(
			python3 - <<PY
asm_len=$asm_len; total=$total; targetx=$targetx
p = min(1.0, (targetx*asm_len)/total) if total>0 else 1.0
print(f"{p:.6f}")
PY
		)
		echo "[CAP] seqtk sample p=${p} to ~${targetx}x"
		seqtk sample -s42 "$reads" "$p" >"$outfq"
	fi
}

# ----------------- polishing -----------------
polished="$target_fa"
final_out=""

if [[ "$mode" == "ont" ]]; then
	echo "[POLAP] ONT-only polishing with Racon (rounds=$racon_rounds)"
	cur="$target_fa"

	if ((ab_test == 1)); then
		echo "[AB] Running A (ALL reads) and B (recruit+cap) for comparison"
		# A: ALL reads, no cap
		bash "$0" polish "$asm" --ont "$ont" --mode ont --threads "$threads" \
			--speed "$speed" --recruit off --ont-cap 0 --prefix "${prefix}_A" \
			$([[ $gfa_preserve -eq 1 ]] && echo --gfa-preserve) >/dev/null
		# B: recruit on draft, then cap to 200× (or --ont-cap if set)
		local capX=$([ "$ont_cap" -gt 0 ] && echo "$ont_cap" || echo 200)
		bash "$0" polish "$asm" --ont "$ont" --mode ont --threads "$threads" \
			--speed "$speed" --recruit draft --ont-cap "$capX" --prefix "${prefix}_B" \
			$([[ $gfa_preserve -eq 1 ]] && echo --gfa-preserve) >/dev/null
		qc_report=1
		# Prefer A’s result as default output
		A_DIR=$(ls -1dt out/"${prefix}_A".polish_* 2>/dev/null | head -n1)
		[[ -n "$A_DIR" ]] && polished="$(ls -1 "$A_DIR"/*polished.* 2>/dev/null | head -n1 || true)"
		[[ -z "${polished:-}" ]] && polished="$target_fa"
	else
		# Single run
		use_reads="$ont"
		label="ALL"
		if [[ "$recruit" != "off" ]]; then
			recruited="$odir/ont.recruit.fq"
			case "$recruit" in
			draft) _recruit_for_capping "$cur" "$ont" "$recruited" "$minimap_preset" "$threads" ;;
			bait)
				if [[ -n "$bait_fa" && -s "$bait_fa" ]]; then
					_recruit_for_capping "$bait_fa" "$ont" "$recruited" "$minimap_preset" "$threads"
				else
					echo "[RECRUIT] --bait not provided; using all reads"
					cp -f "$ont" "$recruited"
				fi
				;;
			auto) _recruit_for_capping "$cur" "$ont" "$recruited" "$minimap_preset" "$threads" ;;
			off | *) cp -f "$ont" "$recruited" ;;
			esac
			capped="$odir/ont.capped.fq"
			_cap_ont_to "$cur" "$recruited" "$ont_cap" "$capped"
			use_reads="$capped"
			label="RECRUIT+CAP"
		fi

		# Racon rounds
		racon_w=1000
		[[ "$speed" == "fast" ]] && racon_w=1200
		[[ "$speed" == "turbo" ]] && racon_w=1500
		for i in $(seq 1 "$racon_rounds"); do
			echo "[minimap2][$label] round $i  --secondary=no  -K $minimap_K"
			minimap2 -x "$minimap_preset" --secondary=no -K "$minimap_K" -t "$threads" "$cur" "$use_reads" >"$odir/r${i}.paf"
			echo "[racon] round $i (-w $racon_w)"
			racon -t "$threads" -w "$racon_w" "$use_reads" "$odir/r${i}.paf" "$cur" >"$odir/segments.racon${i}.fa"
			cur="$odir/segments.racon${i}.fa"
		done
		polished="$cur"
	fi
elif [[ "$mode" == "hybrid" ]]; then
	[[ "$polisher" == "auto" ]] && polisher="fmlrc2"
	echo "[POLAP] Hybrid polishing using $polisher"
	if [[ "$polisher" == "fmlrc2" ]]; then
		echo "[FMLRC2] building BWT (ropebwt2 → fmlrc2-convert)"
		(_quiet_zcat "$sr1" "$sr2" || cat "$sr1" "$sr2") |
			awk 'NR%4==2' | ropebwt2 -LR | fmlrc2-convert "$odir/comp_msbwt.npy"
		echo "[FMLRC2] polishing"
		fmlrc2 -t "$threads" "$odir/comp_msbwt.npy" "$target_fa" "$odir/segments.fmlrc2.fa"
		polished="$odir/segments.fmlrc2.fa"
	elif [[ "$polisher" == "polypolish" ]]; then
		echo "[Polypolish] mapping ALL placements with BWA-MEM"
		bwa index "$target_fa"
		bwa mem -a -t "$threads" "$target_fa" "$sr1" | samtools sort -@ "$threads" -o "$odir/r1.all.bam"
		bwa mem -a -t "$threads" "$target_fa" "$sr2" | samtools sort -@ "$threads" -o "$odir/r2.all.bam"
		samtools index -@ "$threads" "$odir/r1.all.bam"
		samtools index -@ "$threads" "$odir/r2.all.bam"
		polypolish "$target_fa" "$odir/r1.all.bam" "$odir/r2.all.bam" >"$odir/segments.polypolish.fa"
		polished="$odir/segments.polypolish.fa"
	else
		echo "ERROR: unsupported --polisher '$polisher' (use auto|fmlrc2|polypolish)" >&2
		exit 1
	fi
fi

# ----------------- reinject into GFA or finalize FASTA -----------------
if ((is_gfa == 1)) && ((gfa_preserve == 1)); then
	echo "[POLAP] Reinserting polished S-sequences into GFA"
	python3 "$scripts_dir/gfa_inject_polished.py" "$asm" "$polished" >"$odir/${prefix}.polished.gfa"
	final_out="$odir/${prefix}.polished.gfa"
else
	cp -f "$polished" "$odir/${prefix}.polished.fa"
	final_out="$odir/${prefix}.polished.fa"
fi

# ----------------- QC REPORT -----------------
if ((qc_report == 1)); then
	mkdir -p "$odir/qc"
	_ont_qc() {
		local ref="$1" tag="$2"
		local paf="$odir/qc/${tag}.paf" tsv="$odir/qc/${tag}.pafstats.tsv"
		minimap2 -x "$minimap_preset" --secondary=no -K "$minimap_K" -t "$threads" "$ref" "$ont" >"$paf"
		awk -f "$scripts_dir/paf_stats.awk" "$paf" >"$tsv" || true
	}
	if ((ab_test == 1)); then
		A_DIR=$(ls -1dt out/"${prefix}_A".polish_* 2>/dev/null | head -n1)
		B_DIR=$(ls -1dt out/"${prefix}_B".polish_* 2>/dev/null | head -n1)
		A_OUT=$(ls -1 "$A_DIR"/*polished.* 2>/dev/null | head -n1 || true)
		B_OUT=$(ls -1 "$B_DIR"/*polished.* 2>/dev/null | head -n1 || true)
		[[ -n "${A_OUT:-}" ]] && _ont_qc "$A_OUT" "A_all_reads"
		[[ -n "${B_OUT:-}" ]] && _ont_qc "$B_OUT" "B_recruit_cap"
	else
		_ont_qc "$final_out" "single_run"
	fi
	# Merqury QV if SR present
	if ((eval_merqury == 1)) && [[ -n "${sr1:-}" && -n "${sr2:-}" ]]; then
		if ((ab_test == 1)); then
			meryl k="$kmer" count "$sr1" "$sr2" output "$odir/qc/reads.meryl"
			[[ -n "${A_OUT:-}" ]] && merqury.sh "$odir/qc/reads.meryl" "$A_OUT" "$odir/qc/A" || true
			[[ -n "${B_OUT:-}" ]] && merqury.sh "$odir/qc/reads.meryl" "$B_OUT" "$odir/qc/B" || true
		else
			meryl k="$kmer" count "$sr1" "$sr2" output "$odir/qc/reads.meryl"
			merqury.sh "$odir/qc/reads.meryl" "$final_out" "$odir/qc/SINGLE" || true
		fi
	fi
	echo "[QC] Files:"
	ls -1 "$odir/qc" | sed "s|^|[QC]   $odir/qc/|"
fi

echo "[DONE] Output: $final_out"
