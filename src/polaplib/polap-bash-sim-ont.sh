#!/usr/bin/env bash
# polap-bash-sim-ont.sh v0.1.0
# Simulate ONT reads for organelles + nuclear background, with optional NUMT/NUPT spiking.
# Modes:
#   --mode badread       : Badread (no spiking)
#   --mode badread-nu    : Badread with NUMT/NUPT spiking
#   --mode nanosim       : NanoSim (train from --in-fastq), no spiking   [ELIDED]
#   --mode nanosim-nu    : NanoSim (train from --in-fastq), with spiking [ELIDED]
set -euo pipefail

############################################
# helpers
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
SIM_LOG_LEVEL="${SIM_LOG_LEVEL:-1}" # 0=quiet,1=info,2=debug
vlog() {
	local l=$1
	shift
	[[ $SIM_LOG_LEVEL -ge $l ]] && echo "[INFO]" "$@" >&2
}
dlog() {
	local l=$1
	shift
	[[ $SIM_LOG_LEVEL -ge $l ]] && echo "[DBG]" "$@" >&2
}
die() {
	echo "[ERR]" "$@" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

############################################
# defaults
############################################
SIM_MT_REF=""
SIM_PT_REF=""
SIM_NUC_REF=""
SIM_NUC_SIZE=1000000

SIM_DEPTH_NUC="10x"
SIM_DEPTH_MT="50x"
SIM_DEPTH_PT="500x"

SIM_OUT="sim_ont"
SIM_THREADS=8
SIM_IN_FASTQ=""
SIM_MODE_CLI=""

# ONT tech preset
SIM_TECH="r10.4" # r9,r9.4,r9.4.1,r10,r10.4,kit14,sup,sup_r10,dorado,duplex

# Badread knobs
SIM_LEN_MEAN=15000
SIM_LEN_SD=7000
SIM_CHIMERAS=0.01
SIM_JUNK=0.002
SIM_RANDOM=0.002
SIM_ERR_MODEL="nanopore2023"
SIM_QS_MODEL="nanopore2023"
SIM_BADREAD_IDENT="" # e.g., "96,99.5,1.2" for SUP

SIM_SEED_NUC=101
SIM_SEED_MT=102
SIM_SEED_PT=103

# NanoSim (train/sim) — elided paths remain for future wiring
SIM_SRA=""
SIM_TRAIN_READS=""
SIM_TRAIN_REF=""
SIM_MODEL_DIR=""
SIM_MODEL_PREFIX=""

# NUMT/NUPT spiking (used only in *-nu modes)
SIM_SPIKE_NUMT=2
SIM_SPIKE_NUPT=2
SIM_SPIKE_NUMT_LEN="1000,3000"
SIM_SPIKE_NUPT_LEN="500,1500"
SIM_SPIKE_LEN_MODE="uniform" # or gaussian
SIM_SPIKE_DIV_SUB=0.02
SIM_SPIKE_DIV_INDEL=0.002
SIM_SPIKE_INV_FRAC=0.2
SIM_SPIKE_MODE="insert" # or replace
SIM_SPIKE_CIRC_MT=0     # spiker setting (independent from read circular sampling)
SIM_SPIKE_CIRC_PT=0
SIM_SPIKE_SEED=42
SIM_SPIKER_BIN="${SCRIPT_DIR}/polap-py-spike-numt-nupt.py"

# Circular read sampling defaults
# - Nuclear ALWAYS linear (no option to circularize)
# - Organelles are circular by default; --no-circular flips them to linear
SIM_ORG_CIRC=1

# Rotation of circular refs to avoid fixed breakpoint (default ON)
SIM_JUNC_MARGIN=1000 # min bp away from biological origin
SIM_CIRC_SEED=41     # RNG seed for rotation offset

# Optional sanity map
SIM_SANITY_MAP=1

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --mt mt.fa --pt pt.fa [--nuc nuc.fa|--nuc-size INT] --out DIR [--mode MODE] [options]

Required:
  --mt FILE                 mitochondrion reference (FASTA)
  --pt FILE                 plastid/chloroplast reference (FASTA)

General:
  --nuc FILE                nuclear reference (FASTA). If omitted, synthesize of --nuc-size
  --nuc-size INT            synthetic nuclear size (default ${SIM_NUC_SIZE})
  --depth-nuc STR           nuclear depth (default ${SIM_DEPTH_NUC})
  --depth-mt STR            mito depth   (default ${SIM_DEPTH_MT})
  --depth-pt STR            plastid depth(default ${SIM_DEPTH_PT})
  --out DIR                 output directory (default ${SIM_OUT})
  --threads INT             threads (default ${SIM_THREADS}); also controls Badread parallel chunks
  -v, --verbose             verbose logs
  -q, --quiet               quiet logs
  --no-sanity-map           skip minimap2 sanity mapping at the end

Circular sampling (read simulation):
  --no-circular             make organelles (mt/pt) linear; nuclear is always linear
  --junction-margin INT     min bp away from origin when rotating circular refs (default ${SIM_JUNC_MARGIN})
  --circ-seed INT           RNG seed for rotation (default ${SIM_CIRC_SEED})

Mode (choose one; default auto):
  --mode {badread,badread-nu,nanosim,nanosim-nu}

Badread tuning:
  --tech STR                {r9,r9.4,r9.4.1,r10,r10.4,kit14,sup,sup_r10,dorado,duplex} (default ${SIM_TECH})
  --len-mean INT            override mean read length
  --len-sd   INT            override length SD
  --identity M,N,S          Badread identity triplet (override)
  --err-model STR           Badread --error_model
  --qs-model STR            Badread --qscore_model

NanoSim (train/sim options) [ELIDED]:
  --in-fastq FILE           single ONT FASTQ for training (required for nanosim/-nu unless --model-prefix)
  --train-reads FILE        (alt) training reads
  --train-ref FILE          training reference (required if training a model)
  --model-dir DIR           save new model here
  --model-prefix DIR        use a pretrained NanoSim model (skip training)

Spiking (used only in *-nu modes):
  --spike-numt INT          number of NUMT inserts
  --spike-nupt INT          number of NUPT inserts
  --spike-numt-len A,B      NUMT length spec (uniform min,max or gaussian mean,sd)
  --spike-nupt-len A,B      NUPT length spec
  --spike-len-mode MODE     {uniform|gaussian}
  --spike-div-sub FLOAT     substitution divergence
  --spike-div-indel FLOAT   single-base indel rate
  --spike-inv-frac FLOAT    fraction inverted
  --spike-mode MODE         {insert|replace}
  --spike-circular-mt       mt circular sampling inside spiker (independent from read simulator)
  --spike-circular-pt       pt circular sampling inside spiker
  --spike-seed INT          RNG seed

EOF
}

############################################
# parse args
############################################
while [[ $# -gt 0 ]]; do
	case "$1" in
	--mt)
		SIM_MT_REF="$2"
		shift 2
		;;
	--pt)
		SIM_PT_REF="$2"
		shift 2
		;;
	--nuc)
		SIM_NUC_REF="$2"
		shift 2
		;;
	--nuc-size)
		SIM_NUC_SIZE="$2"
		shift 2
		;;
	--depth-nuc)
		SIM_DEPTH_NUC="$2"
		shift 2
		;;
	--depth-mt)
		SIM_DEPTH_MT="$2"
		shift 2
		;;
	--depth-pt)
		SIM_DEPTH_PT="$2"
		shift 2
		;;
	--out)
		SIM_OUT="$2"
		shift 2
		;;
	--threads)
		SIM_THREADS="$2"
		shift 2
		;;
	--no-sanity-map)
		SIM_SANITY_MAP=0
		shift
		;;
	--no-circular)
		SIM_ORG_CIRC=0
		shift
		;;
	--junction-margin)
		SIM_JUNC_MARGIN="$2"
		shift 2
		;;
	--circ-seed)
		SIM_CIRC_SEED="$2"
		shift 2
		;;
	--mode)
		SIM_MODE_CLI="$2"
		shift 2
		;;
	--tech)
		SIM_TECH="$2"
		shift 2
		;;
	--len-mean)
		SIM_LEN_MEAN="$2"
		shift 2
		;;
	--len-sd)
		SIM_LEN_SD="$2"
		shift 2
		;;
	--identity)
		SIM_BADREAD_IDENT="$2"
		shift 2
		;;
	--err-model)
		SIM_ERR_MODEL="$2"
		shift 2
		;;
	--qs-model)
		SIM_QS_MODEL="$2"
		shift 2
		;;
	--in-fastq)
		SIM_IN_FASTQ="$2"
		shift 2
		;;
	--sra)
		SIM_SRA="$2"
		shift 2
		;;
	--train-reads)
		SIM_TRAIN_READS="$2"
		shift 2
		;;
	--train-ref)
		SIM_TRAIN_REF="$2"
		shift 2
		;;
	--model-dir)
		SIM_MODEL_DIR="$2"
		shift 2
		;;
	--model-prefix)
		SIM_MODEL_PREFIX="$2"
		shift 2
		;;
	--spike-numt)
		SIM_SPIKE_NUMT="$2"
		shift 2
		;;
	--spike-nupt)
		SIM_SPIKE_NUPT="$2"
		shift 2
		;;
	--spike-numt-len)
		SIM_SPIKE_NUMT_LEN="$2"
		shift 2
		;;
	--spike-nupt-len)
		SIM_SPIKE_NUPT_LEN="$2"
		shift 2
		;;
	--spike-len-mode)
		SIM_SPIKE_LEN_MODE="$2"
		shift 2
		;;
	--spike-div-sub)
		SIM_SPIKE_DIV_SUB="$2"
		shift 2
		;;
	--spike-div-indel)
		SIM_SPIKE_DIV_INDEL="$2"
		shift 2
		;;
	--spike-inv-frac)
		SIM_SPIKE_INV_FRAC="$2"
		shift 2
		;;
	--spike-mode)
		SIM_SPIKE_MODE="$2"
		shift 2
		;;
	--spike-circular-mt)
		SIM_SPIKE_CIRC_MT=1
		shift
		;;
	--spike-circular-pt)
		SIM_SPIKE_CIRC_PT=1
		shift
		;;
	--spike-seed)
		SIM_SPIKE_SEED="$2"
		shift 2
		;;
	-v | --verbose)
		SIM_LOG_LEVEL=2
		shift
		;;
	-q | --quiet)
		SIM_LOG_LEVEL=0
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*) die "unknown arg: $1" ;;
	esac
done
[[ -z "$SIM_MT_REF" || -z "$SIM_PT_REF" ]] && {
	usage
	die "--mt and --pt are required"
}

############################################
# prep
############################################
need seqkit
SIM_OUT="$(readlink -f "$SIM_OUT")"
SIM_MT_REF="$(readlink -f "$SIM_MT_REF")"
SIM_PT_REF="$(readlink -f "$SIM_PT_REF")"
[[ -n "$SIM_NUC_REF" ]] && SIM_NUC_REF="$(readlink -f "$SIM_NUC_REF")"
[[ -n "$SIM_IN_FASTQ" ]] && SIM_IN_FASTQ="$(readlink -f "$SIM_IN_FASTQ")"
mkdir -p "$SIM_OUT"
cd "$SIM_OUT"
vlog 1 "polap-bash-sim-ont.sh v0.1.0"
vlog 1 "Output dir: $SIM_OUT"
vlog 2 "Args: mt=$SIM_MT_REF pt=$SIM_PT_REF nuc=${SIM_NUC_REF:-<synthetic>} tech=$SIM_TECH org_circular=$SIM_ORG_CIRC"

make_nuclear() {
	local outfa="$1" size="$2"
	vlog 1 "[nuclear] generating synthetic ${size} bp → ${outfa}"
	python3 - "$size" >"$outfa" <<'PY'
import sys, random
n=int(sys.argv[1]); rng=random.Random(1); b='ACGT'
print(">chrSynthetic")
for i in range(0,n,80):
    print(''.join(rng.choice(b) for _ in range(min(80,n-i))))
PY
}

############################################
# Badread tech preset
############################################
set_ont_tech() {
	case "$SIM_TECH" in
	r9 | r9.4)
		: "${SIM_LEN_MEAN:=10000}"
		: "${SIM_LEN_SD:=6000}"
		: "${SIM_ERR_MODEL:=nanopore2020}"
		: "${SIM_QS_MODEL:=nanopore2020}"
		: "${SIM_BADREAD_IDENT:=87,95,3}"
		;;
	r9.4.1)
		: "${SIM_LEN_MEAN:=12000}"
		: "${SIM_LEN_SD:=6500}"
		: "${SIM_ERR_MODEL:=nanopore2020}"
		: "${SIM_QS_MODEL:=nanopore2020}"
		: "${SIM_BADREAD_IDENT:=90,97,2.5}"
		;;
	r10 | r10.4)
		: "${SIM_LEN_MEAN:=15000}"
		: "${SIM_LEN_SD:=7000}"
		: "${SIM_ERR_MODEL:=nanopore}"
		: "${SIM_QS_MODEL:=nanopore}"
		: "${SIM_BADREAD_IDENT:=94,99,1.6}"
		;;
	kit14 | sup | sup_r10 | dorado)
		: "${SIM_LEN_MEAN:=18000}"
		: "${SIM_LEN_SD:=8000}"
		: "${SIM_ERR_MODEL:=nanopore}"
		: "${SIM_QS_MODEL:=nanopore}"
		: "${SIM_BADREAD_IDENT:=96,99.5,1.2}"
		;;
	duplex)
		: "${SIM_LEN_MEAN:=12000}"
		: "${SIM_LEN_SD:=5000}"
		: "${SIM_ERR_MODEL:=nanopore}"
		: "${SIM_QS_MODEL:=nanopore}"
		: "${SIM_BADREAD_IDENT:=98.5,99.7,0.4}"
		;;
	*)
		vlog 1 "[warn] unknown --tech '$SIM_TECH'; using r10.4-like defaults"
		: "${SIM_LEN_MEAN:=15000}"
		: "${SIM_LEN_SD:=7000}"
		: "${SIM_ERR_MODEL:=nanopore}"
		: "${SIM_QS_MODEL:=nanopore}"
		: "${SIM_BADREAD_IDENT:=94,99,1.6}"
		;;
	esac
	vlog 1 "[tech] $SIM_TECH len≈${SIM_LEN_MEAN}±${SIM_LEN_SD}, identity=${SIM_BADREAD_IDENT}, error_model=${SIM_ERR_MODEL}, qscore_model=${SIM_QS_MODEL}"
}
set_ont_tech

############################################
# Nuclear reference (always linear) + optional spiking
############################################
prepare_nuclear_base() {
	if [[ -n "${SIM_NUC_REF}" ]]; then
		SIM_NUC_REF="$(readlink -f "$SIM_NUC_REF")"
		cp -f "$SIM_NUC_REF" nuclear.base.fa
		vlog 1 "[nuclear] using provided reference: $SIM_NUC_REF"
	else
		make_nuclear nuclear.base.fa "$SIM_NUC_SIZE"
	fi
	cp -f nuclear.base.fa nuclear.fa
}

prepare_nuclear_with_optional_spike() {
	local want_spike="$1" # 0/1
	cp -f nuclear.base.fa nuclear.fa
	if ((want_spike)); then
		[[ -s "$SIM_SPIKER_BIN" ]] || die "spiker not found: $SIM_SPIKER_BIN"
		((SIM_SPIKE_NUMT > 0 || SIM_SPIKE_NUPT > 0)) || die "--mode *-nu requires --spike-* > 0"
		vlog 1 "[spike] NUMT=${SIM_SPIKE_NUMT} NUPT=${SIM_SPIKE_NUPT} (seed=${SIM_SPIKE_SEED})"
		python3 "$SIM_SPIKER_BIN" \
			--nuc nuclear.base.fa --mt "$SIM_MT_REF" --pt "$SIM_PT_REF" \
			--numt "$SIM_SPIKE_NUMT" --nupt "$SIM_SPIKE_NUPT" \
			--numt-len "$SIM_SPIKE_NUMT_LEN" --nupt-len "$SIM_SPIKE_NUPT_LEN" \
			--len-mode "$SIM_SPIKE_LEN_MODE" \
			--div-sub "$SIM_SPIKE_DIV_SUB" --div-indel "$SIM_SPIKE_DIV_INDEL" \
			--inv-frac "$SIM_SPIKE_INV_FRAC" \
			--mode "$SIM_SPIKE_MODE" \
			--seed "$SIM_SPIKE_SEED" \
			${SIM_SPIKE_CIRC_MT:+--circular-mt} \
			${SIM_SPIKE_CIRC_PT:+--circular-pt} \
			--out nuclear.spiked
		cp -f nuclear.spiked.augmented.fa nuclear.fa
		vlog 1 "[spike] truth: nuclear.spiked.inserts.{bed,tsv}"
	fi
}

############################################
# Circular helpers (rotate & double) + coverage scaling
############################################
make_circular_ref_rot() {
	local infa="$1" outfa="$2" margin="$3" seed="$4"
	python3 - "$infa" "$outfa" "$margin" "$seed" <<'PY'
import sys,random
inf, outf, margin, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
name, seq = None, []
with open(inf) as fh:
    for ln in fh:
        ln = ln.strip()
        if not ln: continue
        if ln.startswith('>'):
            if name is None: name = ln[1:].split()[0]
        else:
            seq.append(ln.upper())
s=''.join(seq)
if not s: raise SystemExit("no sequence in "+inf)
L=len(s)
margin=min(max(0,margin), max(0,L//2))
rng=random.Random(seed)
if L>2*margin:
    r=rng.randrange(margin, L-margin)
else:
    r=0
rot = s[r:]+s[:r]
with open(outf,'w') as fo:
    fo.write(f">{name}_circ2_rot{r}\n")
    ss = rot+rot
    for i in range(0, 2*L, 80):
        fo.write(ss[i:i+80]+"\n")
PY
}

scale_quantity_x() {
	local qty="$1" factor="$2"
	if [[ "$qty" =~ ^([0-9]+([.][0-9]+)?)x$ ]]; then
		awk -v q="${BASH_REMATCH[1]}" -v f="$factor" 'BEGIN{printf("%.6fx", q*f)}'
	else
		echo "$qty" # absolute read/bp quantities left as-is
	fi
}

############################################
# Badread parallel helpers (split quantity)
############################################
split_quantity_even() {
	local qty="$1" jobs="$2"
	((jobs <= 1)) && {
		echo "$qty"
		return
	}
	if [[ "$qty" =~ ^([0-9]+([.][0-9]+)?)x$ ]]; then
		local total="${BASH_REMATCH[1]}"
		local each=$(awk -v n="$total" -v j="$jobs" 'BEGIN{printf("%.6f", n/j)}')
		local parts=() i
		for ((i = 1; i < jobs; i++)); do parts+=("${each}x"); done
		local last=$(awk -v n="$total" -v j="$jobs" -v e="$each" 'BEGIN{printf("%.6f", n-(j-1)*e)}')
		parts+=("${last}x")
		echo "${parts[*]}"
	else
		local total=$(awk -v q="$qty" 'BEGIN{printf("%.0f", q)}')
		local chunk=$((total / jobs))
		local parts=() i
		for ((i = 1; i < jobs; i++)); do parts+=("$chunk"); done
		parts+=($((total - chunk * (jobs - 1))))
		echo "${parts[*]}"
	fi
}

run_badread_parallel_one() {
	local ref="$1" qty="$2" mean="$3" sd="$4" em="$5" qm="$6" seed="$7" outf="$8"
	shift 8
	local -a ident_opts=("$@")
	need badread
	local jobs="$SIM_THREADS"
	((jobs > 6)) && jobs=6
	((jobs < 1)) && jobs=1

	if ((jobs == 1)); then
		badread simulate --reference "$ref" --quantity "$qty" \
			--length "$mean","$sd" --error_model "$em" --qscore_model "$qm" \
			"${ident_opts[@]}" --chimeras "$SIM_CHIMERAS" --junk_reads "$SIM_JUNK" --random_reads "$SIM_RANDOM" \
			--seed "$seed" >"$outf"
		return
	fi

	local chunks i=0 pids=() tmpfiles=()
	chunks=$(split_quantity_even "$qty" "$jobs")
	for q in $chunks; do
		i=$((i + 1))
		local tf
		tf="$(mktemp "${outf%.fq}.part${i}.XXXXXX.fq")"
		tmpfiles+=("$tf")
		local s=$((seed + i))
		(
			badread simulate --reference "$ref" --quantity "$q" \
				--length "$mean","$sd" --error_model "$em" --qscore_model "$qm" \
				"${ident_opts[@]}" --chimeras "$SIM_CHIMERAS" --junk_reads "$SIM_JUNK" --random_reads "$SIM_RANDOM" \
				--seed "$s" >"$tf"
		) &
		pids+=($!)
	done
	local ok=1 pid
	for pid in "${pids[@]}"; do wait "$pid" || ok=0; done
	((ok == 1)) || die "badread parallel chunk failed for $ref"
	: >"$outf"
	for tf in "${tmpfiles[@]}"; do cat "$tf" >>"$outf"; done
	for tf in "${tmpfiles[@]}"; do rm -f "$tf"; done
}

############################################
# Simulators
############################################
sim_badread() {
	local ident_opts=()
	[[ -n "$SIM_BADREAD_IDENT" ]] && ident_opts=(--identity "$SIM_BADREAD_IDENT")

	# Decide circular vs linear refs for read simulator (nuclear always linear)
	local REF_NUC="nuclear.fa"
	local REF_MT="$SIM_MT_REF"
	local REF_PT="$SIM_PT_REF"

	local QTY_NUC="$SIM_DEPTH_NUC"
	local QTY_MT="$SIM_DEPTH_MT"
	local QTY_PT="$SIM_DEPTH_PT"

	if ((SIM_ORG_CIRC)); then
		make_circular_ref_rot "$REF_MT" "mt.circ2.fa" "$SIM_JUNC_MARGIN" "$SIM_CIRC_SEED"
		make_circular_ref_rot "$REF_PT" "pt.circ2.fa" "$SIM_JUNC_MARGIN" "$SIM_CIRC_SEED"
		REF_MT="mt.circ2.fa"
		REF_PT="pt.circ2.fa"
		# doubled ref → double the x to preserve effective depth on the original circle
		QTY_MT="$(scale_quantity_x "$QTY_MT" 2.0)"
		QTY_PT="$(scale_quantity_x "$QTY_PT" 2.0)"
		vlog 1 "[circular] mt/pt rotated with margin=${SIM_JUNC_MARGIN}, quantities: mt=${QTY_MT} pt=${QTY_PT}"
	else
		vlog 1 "[circular] disabled for organelles; mt/pt linear"
	fi

	# nuclear
	vlog 1 "[sim] Badread → nuclear ${QTY_NUC}"
	run_badread_parallel_one "$REF_NUC" "${QTY_NUC}" "${SIM_LEN_MEAN}" "${SIM_LEN_SD}" \
		"${SIM_ERR_MODEL}" "${SIM_QS_MODEL}" "${SIM_SEED_NUC}" "nuc.fq" \
		"${ident_opts[@]}"

	# mito
	vlog 1 "[sim] Badread → mito ${QTY_MT}"
	run_badread_parallel_one "$REF_MT" "${QTY_MT}" "${SIM_LEN_MEAN}" "${SIM_LEN_SD}" \
		"${SIM_ERR_MODEL}" "${SIM_QS_MODEL}" "${SIM_SEED_MT}" "mt.fq" \
		"${ident_opts[@]}"

	# plastid
	vlog 1 "[sim] Badread → plastid ${QTY_PT}"
	run_badread_parallel_one "$REF_PT" "${QTY_PT}" "${SIM_LEN_MEAN}" "${SIM_LEN_SD}" \
		"${SIM_ERR_MODEL}" "${SIM_QS_MODEL}" "${SIM_SEED_PT}" "pt.fq" \
		"${ident_opts[@]}"

	# Add prefixes to read names
	for tag in nuc mt pt; do
		tmp="${tag}.tmp"
		awk -v PFX="$tag" '{if(NR%4==1){print "@"PFX"_"substr($0,2)} else print}' ${tag}.fq >$tmp
		mv $tmp ${tag}.fq
	done

	cat nuc.fq mt.fq pt.fq >all.sim.fq
	vlog 1 "[sim] wrote all.sim.fq"
	seqkit stats -Ta all.sim.fq | awk 'NR==2{print "[stats] combined:", $0}' >&2 || true
}

sim_nanosim() {
	die "NanoSim path elided to keep this version focused on Badread + circular sampling. Wire as needed."
}

############################################
# Decide mode
############################################
decide_mode() {
	if [[ -n "$SIM_MODE_CLI" ]]; then
		case "$SIM_MODE_CLI" in
		badread | badread-nu | nanosim | nanosim-nu)
			echo "$SIM_MODE_CLI"
			return
			;;
		*) die "--mode must be one of {badread,badread-nu,nanosim,nanosim-nu}" ;;
		esac
	fi
	if [[ -n "${SIM_SRA}" || -n "${SIM_TRAIN_READS}" || -n "${SIM_MODEL_PREFIX}" || -n "${SIM_IN_FASTQ}" ]]; then
		echo "nanosim"
	else
		echo "badread"
	fi
}

############################################
# run
############################################
MODE="$(decide_mode)"
vlog 1 "[mode] ${MODE}"

vlog 1 "1"
prepare_nuclear_base
vlog 1 "2"
case "$MODE" in
badread)
	prepare_nuclear_with_optional_spike 0
	sim_badread
	;;
badread-nu)
	prepare_nuclear_with_optional_spike 1
	sim_badread
	;;
nanosim)
	prepare_nuclear_with_optional_spike 0
	sim_nanosim
	;;
nanosim-nu)
	prepare_nuclear_with_optional_spike 1
	sim_nanosim
	;;
esac

############################################
# sanity map (optional)
############################################
if ((SIM_SANITY_MAP)); then
	if command -v minimap2 >/dev/null 2>&1; then
		vlog 1 "[sanity] minimap2 all.sim.fq to combined refs"
		# For sanity, index what was fed to the simulator:
		mt_idx_src=$([[ $SIM_ORG_CIRC -eq 1 ]] && echo "mt.circ2.fa" || echo "$SIM_MT_REF")
		pt_idx_src=$([[ $SIM_ORG_CIRC -eq 1 ]] && echo "pt.circ2.fa" || echo "$SIM_PT_REF")
		nuc_idx_src="nuclear.fa" # nuclear always linear
		cat "$mt_idx_src" "$pt_idx_src" "$nuc_idx_src" >combined_refs.fa
		minimap2 -x map-ont -t "$SIM_THREADS" combined_refs.fa all.sim.fq >all_vs_refs.paf
		vlog 1 "[sanity] wrote all_vs_refs.paf"
	else
		vlog 1 "[sanity] minimap2 not found; skip"
	fi
fi

vlog 1 "[done] Outputs in $SIM_OUT/: all.sim.fq, nuc.fq, mt.fq, pt.fq, nuclear.fa"
