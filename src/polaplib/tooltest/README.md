# polap tool tests

Lightweight smoke/IO tests for CLI tools used by polap.

- **Default tests are fast**: only version/help and tiny IO checks.
- **SLOW tests** (e.g., a tiny flye run) are gated behind `POLAP_TOOLTEST_SLOW=1`.

## Run

```bash
make -C polaplib/tooltest test
# or
bash polaplib/tooltest/Makefile test


Great question — and you’re right to be skeptical.

## 1) Does `tt::mkout` (double colon) work in Bash?

**No.** In Bash, a function name must be a valid shell identifier: letters, digits, and underscores, not starting with a digit. The `:` character isn’t allowed in function names. That `tt::mkout` style works in some other shells (and in languages like Python for “namespacing”), but not in Bash.

**Fix:** rename functions to use underscores, e.g. `tt_mkout`, `tt_banner`, etc.

Below is a corrected `common.sh` you can drop in, followed by a one‑liner to mass‑rename any earlier `tt::name` calls to `tt_name` in your `tooltest` scripts.

---

### ✅ Drop‑in replacement: `polaplib/tooltest/_lib/common.sh` (Bash‑legal names)

```bash
#!/usr/bin/env bash
# polaplib/tooltest/_lib/common.sh
# Minimal test helpers (no external deps). Bash-safe names (underscores).

set -Eeuo pipefail

# Optional: enable your failsafe crash tracer if present
if [[ "${POLAP_FAILSAFE:-0}" == "1" ]]; then
  __TT_THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  __TT_POLAPLIB="$(cd "${__TT_THIS_DIR}/../.." && pwd)"
  if [[ -r "${__TT_POLAPLIB}/polap-lib-failsafe.sh" ]]; then
    # shellcheck disable=SC1091
    source "${__TT_POLAPLIB}/polap-lib-failsafe.sh"
    polap_enable_failsafe
  fi
fi

: "${POLAP_TOOLTEST_OUTDIR:=o/tooltest}"

tt_conda_tag() {
  local tag=""
  if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then tag="(${CONDA_DEFAULT_ENV}) "; fi
  printf '%s' "$tag"
}

tt_ts() { date '+%Y-%m-%d %H:%M:%S'; }

tt_banner() {
  printf '[%s %s%s] %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$*"
}

tt_esc() {
  local s=""
  local a
  for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
  printf '%s' "${s# }"
}

tt_echo_cmd() {
  local esc; esc="$(tt_esc "$@")"
  printf '[%s %s%s] %s ; %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$esc" "$*"
}

tt_require_tool() {
  local t="$1"
  if ! command -v "$t" >/dev/null 2>&1; then
    printf '[%s %s%s] SKIP: %s not found in PATH\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$t" >&2
    return 99
  fi
}

tt_assert_rc() {
  local have="$1" want="${2:-0}" msg="${3:-}"
  if [[ "$have" -ne "$want" ]]; then
    printf '[%s %s%s] ASSERT FAIL rc=%d want=%d %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$have" "$want" "$msg" >&2
    return 1
  fi
}

tt_assert_file() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    printf '[%s %s%s] ASSERT FAIL missing/empty: %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$f" >&2
    return 1
  fi
}

tt_assert_contains() {
  local pattern="$1" file="$2"
  if ! grep -Fq -- "$pattern" "$file"; then
    printf '[%s %s%s] ASSERT FAIL: %s not found in %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$pattern" "$file" >&2
    return 1
  fi
}

tt_mkout() { mkdir -p -- "$1"; }

tt_run() {
  tt_echo_cmd "$@"
  "$@"
}

tt_run_to() {
  local out="$1" err="$2"; shift 2
  tt_echo_cmd "$@"
  "$@" >"$out" 2>"$err"
}

# Fixtures
tt_ensure_mini_fa() {
  local dst="$1"
  [[ -s "$dst" ]] && return 0
  mkdir -p -- "$(dirname "$dst")"
  cat >"$dst" <<'FA'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
>chr2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAA
FA
}

tt_ensure_mini_sam() {
  local dst="$1"
  [[ -s "$dst" ]] && return 0
  mkdir -p -- "$(dirname "$dst")"
  cat >"$dst" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:44
@SQ	SN:chr2	LN:44
r001	0	chr1	1	60	10M	*	0	0	ACGTACGTAA	*	NM:i:0
SAM
}

tt_ensure_mini_aln_fa() {
  local dst="$1"
  [[ -s "$dst" ]] && return 0
  mkdir -p -- "$(dirname "$dst")"
  cat >"$dst" <<'ALN'
>sp1
ACGTACGTACGTACGTACGT
>sp2
ACGTACGTACGTACGTACGA
>sp3
ACGTACGTACGTACGTACGG
>sp4
ACGTACGTACGTACGTACGC
ALN
}

tt_ensure_mini_fq() {
  local dst="$1"
  [[ -s "$dst" ]] && return 0
  mkdir -p -- "$(dirname "$dst")"
  cat >"$dst" <<'FQ'
@r1
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
#####################################
@r2
ACGTACGTACGTACGTACGTACGTACGTACGTACGA
+
######################################
FQ
}
```

**Bulk-rename calls in your existing `run.sh` scripts:**

```bash
# From tooltest root
sed -i 's/tt::/tt_/g' */run.sh ../_lib/common.sh 2>/dev/null || true
```

---

## 2) How to document bioinformatics CLI tools (templates + examples)

A consistent, **lightweight template** pays off fast: it saves your future self and your collaborators from guessing flags, validating inputs, or remembering which preset to use.

Below is a **Markdown template** you can copy into each tool’s folder as `USAGE.md`, followed by two filled examples for **minimap2** and **flye** that fit the way you’re testing tools under `polaplib/tooltest`.

### 📄 `USAGE_TEMPLATE.md`

````markdown
# <TOOL NAME> — Usage & Notes

**Purpose:** What problem does this tool solve in polap?

**Version & Source**
- Expected version (min / pinned): `<version or >=…>`
- How to check: `tool --version`
- Package/channel: `<conda channel or system package>`
- License: `<license or link>`

**Inputs**
- Required:
  - `<file(s)>` — format, constraints, typical size
- Optional / common flags:
  - `--flag`: what it does, defaults, pitfalls

**Outputs**
- Primary:
  - `<file>` — brief description
- Secondary / logs:
  - `<file>` — what to look for

**Quickstart (copy/paste)**
```bash
# minimal, single-line working example
````

**Common Tasks**

1. <Task name>

```bash
# exact command with inline comments for important flags
```

2. <Task name>

```bash
# …
```

**Polap Integration**

* Where used: `<subcommand or pipeline step>`
* Data contracts:

  * Input produced by: `<step>`
  * Output consumed by: `<step>`

**Resource Notes**

* CPU/threads: `<guidance>`
* Memory: `<guidance>`
* Disk/IO: `<guidance>`
* Runtime: `<rough range with small/medium/large>`

**Validation & QA**

* Smoke test: `tool --version` / `tool --help`
* Sanity checks: `<commands to verify outputs>`
* Known failure patterns: `<messages and how to address>`

**Reproducibility**

* Capture exact version(s): `tool --version > o/meta/<tool>.version.txt`
* Deterministic flags: `<flags>`

**Troubleshooting**

* Error: `<message>` → Fix: `<action>`
* Error: `<message>` → Fix: `<action>`

**References**

* Docs: `<link>`
* Paper(s): `<citation>`

````

> Keep these docs **near the tests**:  
> `polaplib/tooltest/<tool>/USAGE.md`  
> so they evolve with your smoke/integration tests.

---

### 🧭 Example: `polaplib/tooltest/minimap2/USAGE.md`

```markdown
# minimap2 — Usage & Notes

**Purpose**  
Long/short read mapper used in polap for read selection, depth stats, and seed mapping.

**Version & Source**
- Expected: ≥ 2.22 (prefer the latest stable)
- Check: `minimap2 --version`
- Conda: `bioconda/minimap2`

**Inputs**
- Reference FASTA (can be indexed to `.mmi` for speed)
- Reads:
  - ONT: FASTQ (gz OK) → preset `-x map-ont`
  - PacBio HiFi: `-x map-hifi`
  - Illumina PE/SE: `-x sr`

**Outputs**
- SAM to stdout by default (pipe to samtools for BAM)
- Optional index: `.mmi` (when using `-d`)

**Quickstart**
```bash
# Index (optional, improves speed on repeated runs)
minimap2 -d ref.mmi ref.fa

# Map ONT reads → BAM
minimap2 -x map-ont -a ref.mmi ont.fq.gz \
  | samtools sort -o aln.sorted.bam
samtools index aln.sorted.bam
````

**Common Tasks**

1. **PacBio HiFi mapping to reference**

```bash
minimap2 -x map-hifi -a ref.fa hifi.fq.gz \
  | samtools sort -o hifi.sorted.bam
samtools index hifi.sorted.bam
```

2. **Illumina short reads**

```bash
# PE reads, output BAM
minimap2 -x sr -a ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o sr.sorted.bam
samtools index sr.sorted.bam
```

3. **Generate `.paf` for graph / chaining stats**

```bash
minimap2 -x map-ont ref.fa ont.fq.gz > aln.paf
```

**Polap Integration**

* Used in: read selection & mapping steps (e.g., `_run_polap_*` mapping menus)
* Contract:

  * Input: FASTA from assembly or seed set; FASTQ reads
  * Output: sorted BAM + index for downstream coverage/annotation

**Resource Notes**

* Threads: `-t N`
* Memory: scales with reference length and index building
* Disk: BAM can be large; ensure enough space

**Validation & QA**

```bash
minimap2 --version > o/meta/minimap2.version.txt
samtools flagstat aln.sorted.bam
samtools idxstats  aln.sorted.bam | head
```

**Troubleshooting**

* *“insufficiently long reads”* → wrong preset (`-x sr` vs `-x map-ont/map-hifi`)
* *“fail to open file”* → path or gzip mismatch
* BAM empty → check read/reference coordinate compatibility and preset

````

---

### 🧭 Example: `polaplib/tooltest/flye/USAGE.md`

```markdown
# Flye — Usage & Notes

**Purpose**  
Long-read de novo assembler used for organelle/whole-genome assembly in polap.

**Version & Source**
- Expected: a recent stable (e.g. ≥ 2.9)
- Check: `flye --version`
- Conda: `bioconda/flye`

**Inputs**
- Reads: choose one of:
  - ONT raw: `--nano-raw reads.fq.gz`
  - PacBio HiFi: `--pacbio-hifi reads.fq.gz`
- Required:
  - `--genome-size <e.g., 500m or 1.6g>`
  - `--out-dir <dir>`

**Outputs**
- Assembly in `<out-dir>/assembly.fasta`
- Logs and intermediate files in `<out-dir>`

**Quickstart**
```bash
# ONT (raw) example
flye --nano-raw ont.fq.gz --genome-size 500m --threads 16 --out-dir asm_ont
# HiFi example
flye --pacbio-hifi hifi.fq.gz --genome-size 500m --threads 16 --out-dir asm_hifi
````

**Common Tasks**

1. **Resume failed/partial run**

```bash
flye --resume --out-dir asm_ont --threads 16
```

2. **Set smaller memory footprint (slower)**

```bash
flye --pacbio-hifi hifi.fq.gz --genome-size 500m --threads 8 --out-dir asm_hifi --meta
```

**Polap Integration**

* Used in: assembly steps (`flye1`, `flye2`, etc.)
* Contract:

  * Input: long reads (ONT/PB) with sufficient coverage (≥20–30× as a rough rule)
  * Output: assembly FASTA for downstream polishing/mapping

**Resource Notes**

* Threads: `--threads N`
* Memory: depends strongly on genome size and coverage
* Runtime: increases with genome size and read count

**Validation & QA**

```bash
flye --version > o/meta/flye.version.txt
ls -lh asm_hifi/assembly.fasta
```

**Troubleshooting**

* *No `assembly.fasta` produced* → check `flye.log` for early termination or not enough reads
* *Very small genome size* → Flye may error or over-prune; consider larger test set for meaningful results

```

---

## 3) Organizing the `tooltest/` docs + tests

A good pattern is **co‑locating** tests and usage docs:

```

polaplib/tooltest/
_lib/common.sh
minimap2/
run.sh
USAGE.md   # <- this file
fixtures/  # (optional; small test data)
flye/
run.sh
USAGE.md
fixtures/
samtools/
run.sh
USAGE.md   # (add later if useful)

````

- **Docs live next to the tests** that validate them.
- Tests write to `o/tooltest/<tool>/…` so they never clash with your pipeline runs.
- Add a tiny **fixture** set per tool (as we already did) to keep tests hermetic.

---

### Quick tip: make your existing test scripts use the fixed names

If any `run.sh` still uses `tt::something`, blitz them to underscores:

```bash
# From polaplib/tooltest
sed -i 's/tt::/tt_/g' */run.sh 2>/dev/null || true
````

That’s it—your helpers will now work in Bash, and your tools will come with small, reproducible docs that grow with your test coverage.

If you want, I can also emit the updated `run.sh` files (samtools/iqtree2/flye) with the underscore functions already wired in.

