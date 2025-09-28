#!/usr/bin/env bash
# polap-bash-mt-bridge.sh
# Iterative mtDNA corridor growth with ONE all-vs-all read mapping.
# Inputs:  R (long reads), P-ref (plastid refs), N-ref (BUSCO proteins), M-seed_0 (initial mt seeds)
# Output:  M-seed_k (final seeds after stopping by delta threshold)
#
# Core ideas:
#   • Do ONE all-vs-all overlap (minimap2 → PAF) on filtered reads R1.
#   • Per round: bin reads to current seeds, sample bridge paths between seed bins, pick high-frequency “middle anchors”, assemble them, update seeds.
#   • Stop when |M| grows by < --delta-stop (e.g., 5%).
#
# Refs (for methods used):
#   - Minimap2 & PAF fields: Li 2018, Bioinformatics 34(18):3094–3100.
#   - Miniprot (protein→genome): Li 2023 (preprint/documentation).
#   - BUSCO lineage protein sets: Manni et al., Mol. Biol. Evol. 38(10):4647–4654.
#
# NOTE:
#   - This script *creates* helper Python/R files when needed (longer than ~10 lines),
#     then invokes them with `python` / `Rscript --vanilla`.
#   - Requires tools in $PATH: minimap2, miniprot (or blastx alternative), samtools, seqkit, awk, gzip, python, Rscript.
#   - You can swap miniprot with BLASTX-ish workflow if preferred.

set -euo pipefail

# -------------------------
# CLI
# -------------------------
_usage() {
	cat <<'EOF'
Usage:
  polap-bash-mt-bridge.sh -r reads.fq[.gz] -p plastid.fa -n busco_proteins.fa -m mt_seeds_0.fa -o outdir
Options:
  -r, --reads         FASTQ(.gz) long reads (R0)
  -p, --plastid-ref   plastid reference FASTA (P-ref)
  -n, --nuc-prot      BUSCO protein FASTA for lineage (N-ref)
  -m, --mt-seeds      initial mt seed contigs FASTA (M-seed_0)
  -o, --outdir        output directory (default: mt_bridge_out)
  -t, --tech          read tech: hifi|ont (affects minimap2 presets; default: hifi)
  --rounds-max        max iterations (default: 6)
  --delta-stop        stop when growth < this fraction (default: 0.05)
  --seed              RNG seed for sampling (default: 13)
  --assembler         assembler command template for middle anchors (default: flye minimal example)
  -v, --verbose       increase verbosity
  --quiet             quiet mode
  --dry               print commands without executing
  -h, --help          show this help
EOF
}

# Defaults
_outdir="mt_bridge_out"
_tech="hifi"
_rounds_max=6
_delta_stop="0.05"
_seed=13
_verbose=0
_quiet=0
_dry=0
_reads=""
_pref=""
_nprot=""
_mtseed0=""
_assembler="flye --threads 8 --out-dir {OUT} --genome-size 500k --${_tech} {READS}"

# Parse
ARGS=()
while [[ $# -gt 0 ]]; do
	case "$1" in
	-r | --reads)
		_reads="$2"
		shift 2
		;;
	-p | --plastid-ref)
		_pref="$2"
		shift 2
		;;
	-n | --nuc-prot)
		_nprot="$2"
		shift 2
		;;
	-m | --mt-seeds)
		_mtseed0="$2"
		shift 2
		;;
	-o | --outdir)
		_outdir="$2"
		shift 2
		;;
	-t | --tech)
		_tech="$2"
		shift 2
		;;
	--rounds-max)
		_rounds_max="$2"
		shift 2
		;;
	--delta-stop)
		_delta_stop="$2"
		shift 2
		;;
	--seed)
		_seed="$2"
		shift 2
		;;
	--assembler)
		_assembler="$2"
		shift 2
		;;
	-v | --verbose)
		_verbose=$((_verbose + 1))
		shift
		;;
	--quiet)
		_quiet=1
		shift
		;;
	--dry)
		_dry=1
		shift
		;;
	-h | --help)
		_usage
		exit 0
		;;
	*)
		ARGS+=("$1")
		shift
		;;
	esac
done

# -------------------------
# Logging (you may replace with your library later)
# -------------------------
_polap_log_ts() { date +"%F %T"; }
_polap_log0() { echo "[$(_polap_log_ts)][ERR] $*" >&2; }
_polap_log1() { [[ ${_quiet:-0} -eq 0 ]] && echo "[$(_polap_log_ts)][INFO] $*" >&2 || true; }

_run() {
	local cmd="$*"
	if [[ "$_dry" -eq 1 ]]; then
		echo "[DRY] $cmd" >&2
	else
		[[ $_verbose -gt 0 ]] && echo "[CMD] $cmd" >&2
		eval "$cmd"
	fi
}

# -------------------------
# Sanity
# -------------------------
[[ -z "$_reads" ]] && _polap_log0 "Missing --reads"
[[ -z "$_reads" ]] && exit 1
[[ -z "$_pref" ]] && _polap_log0 "Missing --plastid-ref"
[[ -z "$_pref" ]] && exit 1
[[ -z "$_nprot" ]] && _polap_log0 "Missing --nuc-prot"
[[ -z "$_nprot" ]] && exit 1
[[ -z "$_mtseed0" ]] && _polap_log0 "Missing --mt-seeds"
[[ -z "$_mtseed0" ]] && exit 1

mkdir -p "$_outdir"
LOG="${_outdir}/pipeline.log"
if [[ "$_dry" -eq 0 ]]; then
	exec > >(tee -a "$LOG") 2>&1
fi
_polap_log1 "Output dir: $_outdir"
_polap_log1 "Tech: $_tech ; Rounds max: $_rounds_max ; Delta stop: $_delta_stop ; Seed: $_seed"
_polap_log1 "Assembler template: $_assembler"

# -------------------------
# Presets
# -------------------------
if [[ "$_tech" == "hifi" ]]; then
	PRESET_AVA="ava-pb"
	PRESET_MAP_SEED="map-hifi"
	PRESET_MAP_PT="map-hifi"
elif [[ "$_tech" == "ont" ]]; then
	PRESET_AVA="ava-ont"
	PRESET_MAP_SEED="map-ont"
	PRESET_MAP_PT="map-ont"
else
	_polap_log0 "Unknown --tech=$_tech (use hifi|ont)"
	exit 1
fi

# -------------------------
# 0) Conservative subtraction: P-ref (plastid) + N-ref (protein)
# -------------------------
STAGE0="${_outdir}/00_subtract"
mkdir -p "$STAGE0"

# Plastid removal
_polap_log1 "Plastid subtraction with minimap2 → remove confident plastid-origin reads"
PT_PAF="${STAGE0}/plastid.paf"
PT_IDS="${STAGE0}/plastid.ids"
_run "minimap2 -x ${PRESET_MAP_PT} --secondary=yes -N 5 -t 8 '${_pref}' '${_reads}' > '${PT_PAF}'"
# Keep confident hits (ID>=0.9, qcov>=0.3, MAPQ>=20). Adjust if needed.
_run "awk 'BEGIN{FS=OFS=\"\\t\"} { if(NF>=12){ident=\$10>0?\$9/\$10:0; qcov=(\$4-\$3)/\$2; if(ident>=0.9 && qcov>=0.3 && \$12>=20) print \$1} }' '${PT_PAF}' | sort -u > '${PT_IDS}'"

# Nuclear removal via miniprot (protein → reads)
_polap_log1 "Nuclear (BUSCO) subtraction with miniprot"
NUC_PAF="${STAGE0}/nuc_miniprot.paf"
NUC_IDS="${STAGE0}/nuc.ids"
_run "miniprot -t 8 '${_nprot}' '${_reads}' > '${NUC_PAF}'"
# Conservative: accept matches with min protein alignment length & MAPQ-like tag 'mp:i:' not standardized; use aln length from PAF[10]
_run "awk 'BEGIN{FS=OFS=\"\\t\"} { if(NF>=12){alen=\$10; if(alen>=150) print \$1} }' '${NUC_PAF}' | sort -u > '${NUC_IDS}'"

# Build R1 (reads after removing confident PT/NUC)
R1="${STAGE0}/R1.fq.gz"
TMP_KEEP="${STAGE0}/keep.ids"
_run "comm -23 <(seqkit fq2fa -n '${_reads}' | sort -u) <(cat '${PT_IDS}' '${NUC_IDS}' | sort -u) > '${TMP_KEEP}'"
_run "seqkit grep -f '${TMP_KEEP}' '${_reads}' -o '${R1}'"
_polap_log1 "R1 built: $(seqkit stats -T '${R1}' | awk 'NR==2{print \$4\" reads, \"\$6\" bp\"}')"

# -------------------------
# 1) ONE all-vs-all on R1 → PAF (used for ALL rounds)
# -------------------------
STAGE1="${_outdir}/01_allvsall"
mkdir -p "$STAGE1"
PAF_ALL="${STAGE1}/allvsall.paf"
_polap_log1 "All-vs-all on R1 (one time only)"
_run "minimap2 -x ${PRESET_AVA} -t 16 --secondary=yes -N 50 --mask-level 0.50 '${R1}' '${R1}' > '${PAF_ALL}'"

# -------------------------
# 2) Gauge coverage from M-seed_0 (map R1→M-seed_0; samtools depth)
# -------------------------
STAGE2="${_outdir}/02_gauge"
mkdir -p "$STAGE2"
BAM0="${STAGE2}/r1_vs_mseed0.bam"
DEP0="${STAGE2}/depth.tsv"
GAUGE="${STAGE2}/gauge.txt" # contains Cstar (seed coverage gauge)

_polap_log1 "Gauge coverage using R1 → M-seed_0"
_run "minimap2 -x ${PRESET_MAP_SEED} -a -t 8 '${_mtseed0}' '${R1}' | samtools sort -@4 -o '${BAM0}'"
_run "samtools index '${BAM0}'"
_run "samtools depth -a '${BAM0}' > '${DEP0}'"

# Short R snippet (<=10 lines): compute median of per-contig medians
R -q --vanilla <<'RSNIP'
dep <- read.table(Sys.getenv("DEP0"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(dep) <- c("ctg","pos","cov")
library(dplyr)
Cstar <- dep %>% group_by(ctg) %>% summarize(m=median(cov, na.rm=TRUE)) %>% summarize(med=median(m, na.rm=TRUE)) %>% pull(med)
cat(Cstar, "\n")
RSNIP
# Capture R output
CSTAR=$(Rscript --vanilla -e 'dep<-read.table("'"$DEP0"'",header=FALSE);colnames(dep)<-c("ctg","pos","cov");suppressPackageStartupMessages(library(dplyr));Cstar<-dep|>group_by(ctg)|>summarize(m=median(cov,na.rm=TRUE))|>summarize(med=median(m,na.rm=TRUE))|>pull(med);cat(Cstar,"\n")' 2>/dev/null | tail -n1 || echo "50")
echo "$CSTAR" >"$GAUGE"
_polap_log1 "Coverage gauge C* ≈ $(cat "$GAUGE")"

# -------------------------
# 3) Helper PY scripts (created if not present)
#    (A) seed-bridge-sampler.py  (B) summarize-walks.py
# -------------------------
HELPERS="${_outdir}/helpers"
mkdir -p "$HELPERS"

if [[ ! -s "${HELPERS}/polap-py-seed-bridge-sampler.py" ]]; then
	_polap_log1 "Writing helper: polap-py-seed-bridge-sampler.py"
	cat >"${HELPERS}/polap-py-seed-bridge-sampler.py" <<'PYEOF'
#!/usr/bin/env python3
import sys, os, argparse, gzip, math, random, itertools
from collections import defaultdict, Counter
import networkx as nx
def _open(p,m="rt"): return gzip.open(p,m) if p.endswith(".gz") else open(p,m)
def parse_paf(l):
  p=l.rstrip("\n").split("\t"); 
  if len(p)<12: return None
  q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=p[0],int(p[1]),int(p[2]),int(p[3]),p[4],p[5],int(p[6]),int(p[7]),int(p[8]),int(p[9]),int(p[10]),int(p[11])
  ident=nm/float(al) if al>0 else 0.0
  return dict(q=q,ql=ql,t=t,tl=tl,al=al,ident=ident,mapq=int(mq))
def wscore(r,formula="ident_x_norm"):
  fr=r["al"]/min(r["ql"],r["tl"]) if min(r["ql"],r["tl"])>0 else 0.0
  if formula=="ident_x_norm": return r["ident"]*fr
  if formula=="ident": return r["ident"]
  if formula=="alen_norm": return fr
  if formula=="mapq_x_ident": return (r["mapq"]/60.0)*r["ident"]
  return r["ident"]
def add(G,a,b,w,d):
  u,v=(a,b) if a<=b else (b,a)
  if u==v: return
  if G.has_edge(u,v):
    if w>G[u][v]["weight"]: G[u][v]["weight"]=w; G[u][v]["data"]=d
  else:
    G.add_edge(u,v,weight=w,data=d)
def paf_to_graph(paf,min_al,min_id,min_w,formula):
  G=nx.Graph()
  with _open(paf,"rt") as f:
    for l in f:
      if not l.strip() or l[0]=="#": continue
      r=parse_paf(l)
      if r is None: continue
      if r["al"]<min_al or r["ident"]<min_id: continue
      w=wscore(r,formula)
      if w<min_w: continue
      add(G,r["q"],r["t"],w,r)
  return G
def load_groups(p):
  D=defaultdict(list)
  with _open(p,"rt") as f:
    for l in f:
      if not l.strip() or l[0]=="#": continue
      s,r=l.rstrip("\n").split("\t")[:2]; D[s].append(r)
  return D
def softmax(ws,T):
  if T<=0: m=max(ws); return [1.0 if w==m else 0.0 for w in ws]
  xs=[w/T for w in ws]; m=max(xs); ex=[math.exp(x-m) for x in xs]; s=sum(ex)
  return [e/s for e in ex] if s>0 else [1.0/len(ws)]*len(ws)
def walk(G,src,dst,rng,T,min_w,max_hops,simple,pr):
  if src not in G or dst not in G: return (False,[src],0.0,0)
  path=[src]; vis={src} if simple else set(); cur=src; logsum=0.0; h=0
  while True:
    if cur==dst: return (True,path,logsum,h)
    if h>=max_hops: return (False,path,logsum,h)
    if rng.random()<pr: path=[src]; vis={src} if simple else set(); cur=src; logsum=0.0; h=0; continue
    nbrs=[]; ws=[]
    for nb in G.neighbors(cur):
      if simple and nb in vis: continue
      w=G[cur][nb]["weight"]
      if w>=min_w: nbrs.append(nb); ws.append(w)
    if not nbrs: return (False,path,logsum,h)
    for i,nb in enumerate(nbrs):
      if nb==dst: ws[i]=max(ws[i],min(1.0,ws[i]*1.05))
    ps=softmax(ws,T); r=rng.random(); a=0.0; idx=len(ps)-1
    for i,p in enumerate(ps):
      a+=p
      if r<=a: idx=i; break
    nxt=nbrs[idx]; w=ws[idx]
    path.append(nxt); 
    if simple: vis.add(nxt)
    logsum+=math.log(max(w,1e-300)); cur=nxt; h+=1
def main():
  ap=argparse.ArgumentParser()
  ap.add_argument("paf"); ap.add_argument("seed_groups"); ap.add_argument("-o","--outdir",default="seed_bridge")
  ap.add_argument("--formula",default="ident_x_norm",choices=["ident_x_norm","ident","alen_norm","mapq_x_ident"])
  ap.add_argument("--min-alen",type=int,default=600); ap.add_argument("--min-ident",type=float,default=0.86)
  ap.add_argument("--min-weight",type=float,default=0.12)
  ap.add_argument("--walks-per-pair",type=int,default=100); ap.add_argument("--max-hops",type=int,default=40)
  ap.add_argument("--temperature",type=float,default=0.25); ap.add_argument("--min-edge-w-walk",type=float,default=0.12)
  ap.add_argument("--simple",action="store_true"); ap.add_argument("--restart",type=float,default=0.02)
  ap.add_argument("--seed",type=int,default=13); ap.add_argument("--only-seeds")
  args=ap.parse_args()
  os.makedirs(args.outdir,exist_ok=True)
  G=paf_to_graph(args.paf,args.min_alen,args.min_ident,args.min_weight,args.formula)
  groups=load_groups(args.seed_groups)
  if args.only_seeds: 
    keep=set([x.strip() for x in args.only_seeds.split(",") if x.strip()])
    groups={k:v for k,v in groups.items() if k in keep}
  nodes=set(G.nodes()); groups={k:[r for r in v if r in nodes] for k,v in groups.items()}
  seeds=sorted([k for k,v in groups.items() if v])
  pairs_dir=os.path.join(args.outdir,"pairs"); os.makedirs(pairs_dir,exist_ok=True)
  rng=random.Random(args.seed); freq=Counter()
  import itertools
  for a,b in itertools.combinations(seeds,2):
    A,B=groups[a],groups[b]; 
    if not A or not B: continue
    outp=os.path.join(pairs_dir,f"{a}__{b}"); os.makedirs(outp,exist_ok=True)
    succ=[]
    for i in range(1,args.walks_per_pair+1):
      src=rng.choice(A); dst=rng.choice(B)
      ok, path, logsum, hops = walk(G,src,dst,rng,args.temperature,args.min_edge_w_walk,args.max_hops,args.simple,args.restart)
      if ok:
        succ.append((logsum,hops,path)); freq.update(path)
        with open(os.path.join(outp,f"walk_{i:05d}.tsv"),"w") as fo:
          print("#idx\tu\tv\tweight",file=fo)
          for j in range(len(path)-1):
            u,v=path[j],path[j+1]; w=G[u][v]["weight"]; print(f"{j+1}\t{u}\t{v}\t{w:.6f}",file=fo)
    succ.sort(key=lambda x:x[0],reverse=True)
    with open(os.path.join(outp,"summary.tsv"),"w") as so:
      print("rank\tlogsum\tedges\tnodes_csv",file=so)
      for rk,(ls,h,pa) in enumerate(succ,1):
        print(f"{rk}\t{ls:.6f}\t{h}\t{','.join(pa)}",file=so)
  # global summary
  rows=sorted(freq.items(), key=lambda x:x[1], reverse=True)
  with open(os.path.join(args.outdir,"read_frequencies.tsv"),"w") as g:
    print("read_id\tcount",file=g)
    for r,c in rows: print(f"{r}\t{c}",file=g)
  with open(os.path.join(args.outdir,"freq_hist.tsv"),"w") as gh:
    from collections import Counter
    H=Counter(freq.values()); print("count\tn_reads",file=gh)
    for k in sorted(H): print(f"{k}\t{H[k]}",file=gh)
PYEOF
	chmod +x "${HELPERS}/polap-py-seed-bridge-sampler.py"
fi

if [[ ! -s "${HELPERS}/polap-py-summarize-walks.py" ]]; then
	_polap_log1 "Writing helper: polap-py-summarize-walks.py"
	cat >"${HELPERS}/polap-py-summarize-walks.py" <<'PYEOF'
#!/usr/bin/env python3
import sys, argparse
from collections import Counter
def parse_summary(p):
  nodes_all=[]
  with open(p) as f:
    for l in f:
      if not l.strip() or l[0]=="#" or l.startswith("rank"): continue
      parts=l.rstrip("\n").split("\t")
      if len(parts)>=4:
        nodes_all.append(parts[3].split(","))
  return nodes_all
def main():
  ap=argparse.ArgumentParser()
  ap.add_argument("summary"); ap.add_argument("-o","--outdir",default="walk_summary")
  ap.add_argument("--min-count",type=int,default=5); ap.add_argument("--top-frac",type=float)
  args=ap.parse_args()
  import os; os.makedirs(args.outdir,exist_ok=True)
  walks=parse_summary(args.summary)
  C=Counter(); [C.update(ns) for ns in walks]
  rows=sorted(C.items(), key=lambda x:x[1], reverse=True)
  with open(os.path.join(args.outdir,"read_frequencies.tsv"),"w") as g:
    print("read_id\tcount",file=g)
    for r,c in rows: print(f"{r}\t{c}",file=g)
  if args.top_fract is not None: pass
PYEOF
	chmod +x "${HELPERS}/polap-py-summarize-walks.py"
fi

# -------------------------
# 4) Iterative rounds
# -------------------------
# Seed FASTA path per round:
SEEDS_CUR="${_outdir}/round_0/m_seeds.fa"
mkdir -p "${_outdir}/round_0"
_run "cp '${_mtseed0}' '${SEEDS_CUR}'"

# Keep track of M size (unique reads attached to any seed)
M_PREV_SIZE=0

for round in $(seq 1 "${_rounds_max}"); do
	RDIR="${_outdir}/round_${round}"
	mkdir -p "$RDIR"
	_polap_log1 "== Round ${round} =="

	# 4.1) Bin reads to current seeds (conservative)
	_polap_log1 "Map R1 → current seeds to form bins"
	PAF_SEED="${RDIR}/r1_vs_seeds.paf"
	_run "minimap2 -x ${PRESET_MAP_SEED} -t 8 --secondary=yes -N 10 '${SEEDS_CUR}' '${R1}' > '${PAF_SEED}'"

	# Build seed_groups.tsv (seed_id \t read_id) with conservative filters
	# Keep hits with identity>=0.88 (HiFi) or 0.85 (ONT), qcov>=0.2, MAPQ>=10
	MIN_IDENT=$([[ "$_tech" == "hifi" ]] && echo "0.88" || echo "0.85")
	SEED_GROUPS="${RDIR}/seed_groups.tsv"
	_run "awk -v minid=${MIN_IDENT} 'BEGIN{FS=OFS=\"\\t\"} { if(NF>=12){ident=\$10>0?\$9/\$10:0; qcov=(\$4-\$3)/\$2; if(ident>=minid && qcov>=0.2 && \$12>=10) print \$6, \$1} }' '${PAF_SEED}' | sort -u > '${SEED_GROUPS}'"

	# 4.2) Compute |M| (union of reads across bins) and delta
	M_READS="${RDIR}/M_reads.ids"
	_run "cut -f2 '${SEED_GROUPS}' | sort -u > '${M_READS}'"
	M_CURR_SIZE=$(wc -l <"${M_READS}" | tr -d '[:space:]')
	if [[ "$M_PREV_SIZE" -gt 0 ]]; then
		num=$(echo "$M_CURR_SIZE - $M_PREV_SIZE" | bc)
		den="$M_PREV_SIZE"
		DELTA=$(
			python - <<PY
num=$num; den=$den
print(0.0 if den==0 else num/den)
PY
		)
	else
		DELTA="1.0"
	fi
	_polap_log1 "|M|_prev=${M_PREV_SIZE} |M|_curr=${M_CURR_SIZE}  Δ=${DELTA}"

	# 4.3) Determine sampling budget and frequency cutoff using C*
	CSTAR_VAL=$(cat "$GAUGE")
	# Simple scaling: N_pair_base=100 scaled by C*/50, clamped [0.5, 3]
	NP_BASE=100
	SCALE=$(
		python - <<PY
c=float("${CSTAR_VAL}"); Cref=50.0
s=max(0.5, min(3.0, (c/Cref if Cref>0 else 1.0)))
print(s)
PY
	)
	WALKS_PER_PAIR=$(
		python - <<PY
print(int(round(${NP_BASE}*${SCALE})))
PY
	)
	FREQ_MIN=$(
		python - <<PY
c=float("${CSTAR_VAL}"); Cref=50.0
# depth-anchored absolute rule: f_cut_abs = 5 + 3*(C*/Cref)
fc=5.0 + 3.0*(c/Cref if Cref>0 else 1.0)
print(int(round(fc)))
PY
	)
	_polap_log1 "Sampling: walks/pair=${WALKS_PER_PAIR} ; frequency cutoff (min-count)=${FREQ_MIN} ; C*=${CSTAR_VAL}"

	# 4.4) Run bridge sampler (one allvsall PAF reused every round)
	BRIDGE_OUT="${RDIR}/bridge"
	mkdir -p "$BRIDGE_OUT"
	_run "python '${HELPERS}/polap-py-seed-bridge-sampler.py' '${PAF_ALL}' '${SEED_GROUPS}' -o '${BRIDGE_OUT}' --walks-per-pair ${WALKS_PER_PAIR} --max-hops 30 --temperature 0.25 --min-edge-w-walk 0.12 --simple --restart 0.02 --seed ${_seed}"

	# 4.5) Pick high-frequency reads
	HF_READS="${RDIR}/highfreq_reads.txt"
	_run "awk -v C=${FREQ_MIN} 'BEGIN{FS=OFS=\"\\t\"} NR>1 && \$2>=C {print \$1}' '${BRIDGE_OUT}/read_frequencies.tsv' | sort -u > '${HF_READS}'"
	_polap_log1 "High-frequency reads: $(wc -l <"${HF_READS}")"

	# 4.6) Assemble middle anchors (extract reads → run assembler)
	MID_FQ="${RDIR}/middle_anchors.fq.gz"
	_run "seqkit grep -f '${HF_READS}' '${R1}' -o '${MID_FQ}'"
	ASM_DIR="${RDIR}/asm"
	mkdir -p "$ASM_DIR"
	# Fill in assembler template
	ASM_CMD="${_assembler}"
	ASM_CMD="${ASM_CMD//\{OUT\}/$ASM_DIR}"
	ASM_CMD="${ASM_CMD//\{READS\}/$MID_FQ}"
	_polap_log1 "Assembler cmd: ${ASM_CMD}"
	_run "${ASM_CMD}"

	# Detect produced contigs (common patterns: assembly.fasta, contigs.fasta, assembly.fna)
	NEW_CONTIGS=""
	for x in "${ASM_DIR}/assembly.fasta" "${ASM_DIR}/contigs.fasta" "${ASM_DIR}/assembly.fna" "${ASM_DIR}/consensus.fasta"; do
		if [[ -s "$x" ]]; then
			NEW_CONTIGS="$x"
			break
		fi
	done
	if [[ -z "$NEW_CONTIGS" ]]; then
		_polap_log1 "WARNING: cannot locate assembler contigs; using middle_anchors.fq.gz (skip update)"
		NEW_SEEDS="${RDIR}/m_seeds.fa"
		_run "seqkit fq2fa '${MID_FQ}' > '${NEW_SEEDS}'"
	else
		NEW_SEEDS="${RDIR}/m_seeds.fa"
		_run "cp '${NEW_CONTIGS}' '${NEW_SEEDS}'"
	fi

	# 4.7) Merge with current seeds → next seeds
	SEEDS_NEXT="${_outdir}/round_${round}/m_seeds_merged.fa"
	_run "cat '${SEEDS_CUR}' '${NEW_SEEDS}' | seqkit seq -u > '${SEEDS_NEXT}'"

	# Prepare next round
	SEEDS_CUR="${SEEDS_NEXT}"
	M_PREV_SIZE="${M_CURR_SIZE}"

	# 4.8) Stopping rule
	# Stop if Δ < delta-stop (unless round==1 where Δ=1.0)
	STOP=$(
		python - <<PY
delta=float("${DELTA}")
thr=float("${_delta_stop}")
print(1 if (delta < thr) else 0)
PY
	)
	if [[ "$STOP" -eq 1 ]]; then
		_polap_log1 "Stop: growth Δ=${DELTA} < threshold ${_delta_stop}"
		break
	fi
done

# -------------------------
# Finalize
# -------------------------
FINAL_ROUND_DIR="$(dirname "${SEEDS_CUR}")"
FINAL_SEEDS="${FINAL_ROUND_DIR}/m_seeds_final.fa"
_run "cp '${SEEDS_CUR}' '${FINAL_SEEDS}'"
_polap_log1 "Done. Final seeds: ${FINAL_SEEDS}"
