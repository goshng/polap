You’re not imagining it—while **Oatk** (and SyncAsm) was _built_ around HiFi, you **can** use its components on **ONT** successfully if you adapt a few things. The “HiFi-only” idea comes from early assumptions about error rates and graph heuristics; once you account for ONT’s profile, the same closed-syncmer/sparse-DBG machinery is perfectly usable.

Below is a practical “what changes” + a ready-to-run **ONT recipe** and guardrails.

---

## Why it works on ONT (principle)

- **Closed syncmers still sample k-mer space unbiasedly.** Errors just reduce matches.
- **Homopolymer compression (HPC)** largely neutralizes ONT’s indel-heavy error mode, making syncmer matching act more “HiFi-like”.
- Adapt **density** (smaller `k/s`) and **graph filters** (looser overlap/Jaccard, edge normalization), add a **coverage prior**, and iterate—then the same graph logic works.

---

## ONT vs HiFi: what you change

**HiFi baseline**

- `k/s`: `121/27` (no HPC)
- `--min-shared`: 4–6
- `--jaccard-min`: 0.01–0.02
- `--edge-norm`: optional
- `--steps`: 2
- SyncAsm seeds: conservative cleanup; no HPC
- Final assembly: allow small bubbles/crossings
- Polishing: often minimal

**ONT-adapted**

- `k/s`: **`41–61 / 17–21`**, and **HPC on** (`--hpc`)
- `--min-shared`: **5–6** (raise if repeats inflate noise)
- `--jaccard-min`: **0.005–0.01** (fuzzier overlaps)
- `--edge-norm`: **on** (de-bias repeat/long-read edges)
- `--steps`: **1–2** (avoid long-range percolation through repeats)
- SyncAsm seeds: conservative cleanup **with HPC**
- Final assembly: allow `--weak-cross 1 --max-bubble 1 --max-tip 1`
- Polishing: **Racon+Medaka** (or your preferred ONT polisher); if available, hybrid polishing or a small HiFi set helps a lot

---

## Minimal ONT recipe (end-to-end)

1. **Pre-filter & coverage prior** (if starting from WGS):

```bash
# quickview (X) on ONT with HPC
syncfilter --mode quickview --hpc -k 41 -s 21 -t 16 -o out.qv reads.ont.fq.gz
# anchors → X-band, PPR (+ X-prior), ensemble
polap select-mt --preset ont --reads reads.ont.fq.gz --mt anchors/mt.id.all.txt \
  -o out/sel --threads 16 --method hybrid --tail 0.05 --tail-mt 0.30 \
  --ppr-k 41 --ppr-s 21 --max-occ 200 --min-shared 5 --jaccard-min 0.0075 \
  --edge-norm --steps 1 --topk-nei 40 --nuc-cut-log10 1.30 --x-slope 0.20 \
  --ensemble majority --emit-fastq
```

2. **Seed assembly (conservative)**:

```bash
# mt-only FASTQ from selection
gzip -dc out/sel.mt.fastq.gz > mt.fastq
# SyncAsm seeds
syncasm -k 41 -s 21 --hpc -t 16 \
  -o mt.seed \
  -c 1 -a 0.05 --unzip-round 0 --max-bubble 0 --max-tip 0 --weak-cross 0 --no-read-ec \
  mt.fastq
```

3. **Score seeds & pick best k** (use your scoring script):

```bash
bash ${_POLAPLIB_DIR}/polap-bash-oatk-score-seeds.sh \
  --seeds mt.seed.asm.utg.gfa --reads mt.fastq --hmm-db db/mt_genes.hmm \
  --out mt.seed --preset ont --threads 16
```

4. **Final (relaxed) pass** (optional if seed already excellent):

```bash
syncasm -k 41 -s 21 --hpc -t 16 \
  -o mt.final \
  -c <picked_c> -a 0.05 --unzip-round 1 --max-bubble 1 --max-tip 1 --weak-cross 1 --no-read-ec \
  mt.fastq
```

5. **Polish** (ONT):

```bash
# e.g., racon (N rounds) + medaka
minimap2 -x map-ont -t 16 mt.final.utg.gfa mt.fastq | samtools sort -@16 -o aln.bam
samtools index aln.bam
# racon …; medaka …  (models depend on your basecaller)
```

> Or let the **autotuner** drive it:
> `polap-bash-autotune-oatk.sh --wgs-mode --reads reads.ont.fq.gz --label mt --anchors anchors/mt.id.all.txt --hmm-db db/mt_genes.hmm --preset ont --threads 16 --final-if "genes_score < 0.90 || breadth < 0.95"`

---

## When ONT works well vs. when to be cautious

**Works well when:**

- ONT basecalling ≥ Q20 (SUP/duplex), coverage per mt ≥ \~40–60×, you use **HPC**, and keep **k/s small**.
- NUMTs exist but you enforce **X-prior** and **edge-norm**, and keep **steps small**.

**Be cautious when:**

- Coverage is very low, or **read N50** is way shorter than repeat spans → fragmentation likely.
- Extremely NUMT-rich genomes (recent transfers at mt-like depth) → rely heavily on **gene anchors**, **X**, and **post-assembly validation**.
- Very old ONT runs (pre-Q20) without re-basecalling → polishing becomes essential; consider hybrid.

---

## Validation checklist (mt/pt)

- **Gene completeness**: `hmmannot` coverage of canonical mt/pt genes (your `genes_score`).
- **Mapping breadth**: fraction of mt/pt reads covering seeds (`breadth`).
- **Coverage band**: unitig X consistent with organelle (not nuclear).
- **Topology sanity**: Bandage; repeats/isomers look expected; no giant nuclear branches.
- **Circularization support** (if expected); **read tiling** continuity.

---

### Bottom line

- The “HiFi-only” myth comes from default assumptions, not core method limits.
- With **HPC**, **smaller k/s**, **looser overlap (Jaccard)**, **edge normalization**, and a **coverage prior**, Oatk/SyncAsm is **perfectly viable on ONT**—especially for **seed-first** workflows followed by polishing.
- Your autotune tools (`polap-bash-autotune-oatk.sh`, `polap-bash-oatk-tune-c.sh`, `polap-bash-oatk-score-seeds.sh`) make the parameter search tractable and repeatable.

If you’d like, I can add an **ONT-focused preset block** to your recipes (including recommended `--final-if`, `k/s`, Jaccard, min-shared, steps) and make the launcher consume it automatically.
