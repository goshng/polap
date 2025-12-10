## OATK mtDNA graph workflow

This section documents a typical end‑to‑end workflow for building **Oatk mtDNA assembly manifests**, **QC metrics plots**, and **PT/MT graph sheets** using the Oatk scripts in `polaplib/`.

### 0. Dependencies

You’ll need:

- **Python 3** (for the manifest + metrics table)
- **R** with these packages:

  - `optparse`, `readr`, `dplyr`, `ggplot2`, `patchwork`, `png`, `grid`, `grDevices`

- POLAP scripts in your `$PATH` or referenced directly, e.g.:

  - `polaplib/Makefile.oatk`
  - `polaplib/polap-data-oatk.sh`
  - `polaplib/polap-bash-oatk-make-manifest.sh`
  - `polaplib/polap-bash-oatk-graph.sh`
  - `polaplib/polap-bash-oatk-ptmt-sheet.sh`

---

### 1. Prepare inputs

1. **Facts table** (TSV) harvested from the pipelines, e.g.:

   ```text
   species  tier  inum  organelle  kind  key   value
   Arabidopsis_thaliana  v6  0  pt   file  png  /path/to/pmat-graph.png
   Arabidopsis_thaliana  v6  0  oatk attr  c30_mem_gb  32
   ...
   ```

   Default path (can be overridden):

   ```bash
   oatk-facts.tsv
   ```

2. **Optional species code mapping** (`species-codes.txt`), space‑delimited:

   ```text
   # code species
   AT  Arabidopsis_thaliana
   OS  Oryza_sativa
   ```

3. **PNG list CSV** for the PT/MT sheet (one row per species):

   ```text
   species,code2,pmat_png,tippo_png,himt_png,oatk_png
   Arabidopsis_thaliana,AT,/path/pmat.png,/path/tippo.png,/path/himt.png,/path/oatk.png
   ...
   ```

   Default path:

   ```bash
   oatk-ptmt-sheet-list.csv
   ```

---

### 2. Quick setup with `polap-data-oatk.sh`

From your project root:

```bash
source polaplib/polap-data-oatk.sh
```

This sets handy defaults like:

```bash
$POLAP_OATK_FACTS        # e.g. $POLAP_OATK_ROOT/oatk-facts.tsv
$POLAP_OATK_MANIFEST     # e.g. $POLAP_OATK_ROOT/oatk-manifest.json
$POLAP_OATK_GRAPH_CSV
$POLAP_OATK_GRAPH_PDF
$POLAP_OATK_PTMTSHEET_LIST
$POLAP_OATK_PTMTSHEET_PDF
```

You can override any of these before or after sourcing:

```bash
export POLAP_OATK_ROOT=/data/oatk-run
export POLAP_OATK_FACTS=/data/oatk-run/my-facts.tsv
source polaplib/polap-data-oatk.sh
```

---

### 3. One‑shot build with `Makefile.oatk`

If you like `make`:

```bash
make -f polaplib/Makefile.oatk all
```

This will:

1. Build the manifest JSON
   (`$POLAP_OATK_MANIFEST`, default: `oatk-manifest.json`)
2. Build the long‑format metrics CSV
   (`$POLAP_OATK_GRAPH_CSV`, default: `oatk-graph-metrics.csv`)
3. Build the QC boxplot PDF
   (`$POLAP_OATK_GRAPH_PDF`, default: `oatk-graph-metrics.pdf`)
4. Build the PT/MT/Oatk sheet PDF
   (`$POLAP_OATK_PTMTSHEET_PDF`, default: `oatk-ptmt-sheet.pdf`)

You can also run pieces:

```bash
# just the manifest
make -f polaplib/Makefile.oatk manifest

# just the metrics plots
make -f polaplib/Makefile.oatk graphs

# just the PT/MT/Oatk sheet
make -f polaplib/Makefile.oatk ptmt-sheet
```

Override individual paths on the command line as needed:

```bash
make -f polaplib/Makefile.oatk \
  OATK_FACTS=/data/oatk-facts.tsv \
  OATK_MANIFEST=/results/oatk-manifest.json \
  OATK_GRAPH_PDF=/results/oatk-metrics.pdf
```

---

### 4. Using the bash wrappers directly

You can also drive each step explicitly.

#### 4.1 Manifest assembly

```bash
polaplib/polap-bash-oatk-make-manifest.sh \
  --facts oatk-facts.tsv \
  --set v5-0-auto \
  --tier v6 \
  --inum 0 \
  --out oatk-manifest.json \
  --species-codes species-codes.txt \
  --pretty
```

This wraps `scripts/oatk_manifest_assemble.py` and writes a structured JSON manifest with:

- per‑species blocks (`pt`, `mt`, `ptpt`, `pmat2`, `tippo`, `himt`, `oatk`, `oatk-*`)
- numeric coercion for data/short1data/short2data
- string‑only metrics for pipeline blocks (`pmat2`, `tippo`, `himt`, `oatk`, `oatk-*`)
- optional `code2` (2‑letter species code)
- optional GFA stats if `scripts/gfa_stats_json.py` is available

#### 4.2 QC metrics plots

```bash
polaplib/polap-bash-oatk-graph.sh \
  --manifest oatk-manifest.json \
  --csv oatk-graph-metrics.csv \
  --pdf oatk-graph-metrics.pdf \
  --title "Mitogenome assembly QC metrics (PMAT, TIPPo, HiMT, Oatk)" \
  --ncol-per-row 3
```

This does:

1. `polap-py-oatk-graph-table.py`
   manifest → long‑format CSV with:

   - `species`, `pipeline` ∈ {`pmat`, `tippo`, `himt`, `oatk`}
   - `mem_gb`, `time_hours`, `disk_gb`,
   - `geneset_completeness_prop`,
   - `num_contigs`, `N50`,
   - `fragmentation_index`, `contig_length_cv`, `max_contig_prop`, `total_length`

2. `polap-r-oatk-graph.R`
   CSV → boxplots for the **main** metrics:

   - memory, runtime, disk
   - geneset completeness
   - number of contigs
   - N50

You can pass additional layout options to the R script (after the usual args), e.g.:

```bash
polaplib/polap-bash-oatk-graph.sh \
  --manifest oatk-manifest.json \
  --pdf oatk-graph-metrics.pdf \
  --ncol-per-row 2 \
  --page-width-in 8 \
  --page-height-in 6
```

#### 4.3 PT/MT/Oatk sheet

Assuming you have `oatk-ptmt-sheet-list.csv` prepared:

```bash
polaplib/polap-bash-oatk-ptmt-sheet.sh \
  --list oatk-ptmt-sheet-list.csv \
  --out oatk-ptmt-sheet.pdf \
  --rows-per-page 6 \
  --start-page 1 \
  --title "Mitogenome assembly graphs" \
  --subtitle "PMAT / TIPPo / HiMT / Oatk"
```

This wraps `polap-r-oatk-ptmt-sheet.R` and produces a multi‑page PDF with:

- one row per species
- columns: PMAT, TIPPo, HiMT, Oatk
- species label (code2 + name) on each row
- configurable rows per page, page size, margins, and gap between images

---

### 5. Minimal “just run it” sequence

Assuming you’ve filled in `oatk-facts.tsv`, `species-codes.txt`, and `oatk-ptmt-sheet-list.csv`:

```bash
# optional, set base directory
export POLAP_OATK_ROOT=/data/oatk-demo
mkdir -p "$POLAP_OATK_ROOT"
cd "$POLAP_OATK_ROOT"

# set defaults
source /path/to/polaplib/polap-data-oatk.sh

# (a) build manifest
/path/to/polaplib/polap-bash-oatk-make-manifest.sh

# (b) build metrics plots
/path/to/polaplib/polap-bash-oatk-graph.sh

# (c) build PT/MT/Oatk sheet
/path/to/polaplib/polap-bash-oatk-ptmt-sheet.sh
```

You should end up with:

- `oatk-manifest.json`
- `oatk-graph-metrics.csv`
- `oatk-graph-metrics.pdf`
- `oatk-ptmt-sheet.pdf`

# Makefile.oatk

with polap source folder at the present directory:
make -f polap/src/polaplib/Makefile.oatk
make -f polap/src/polaplib/Makefile.oatk print-vars
make -f polap/src/polaplib/Makefile.oatk manifest
make -f polap/src/polaplib/Makefile.oatk graph

Draw four different kinds of plots.
First, boxplots of one of the following values: total_mem_gb, total_time_hours, disk_gb, completeness_prop.
X-axis can be ordered by pattern or cov.
pattern is t-h, h-t, t-x, x-h, x-x.
cov is 2,5,8,11,20,30,40,50.
for only x-x, cov is 11,20,30,40,50.
One way plotting boxplots of X-axis can be grouped by pattern, and cov is listed within each pattern.
e.g., t-h-02, t-h-05, t-h-08, t-h-11, t-h-20, t-h-30, t-h-40, t-h-50,
h-t-02, h-t-05, h-t-08, h-t-11, h-t-20, h-t-30, h-t-40, h-t-50,
t-x-02, t-x-05, t-x-08, t-x-11, t-x-20, t-x-30, t-x-40, t-x-50,
x-h-02, x-h-05, x-h-08, x-h-11, x-h-20, x-h-30, x-h-40, x-h-50,
x-x-11, x-x-20, x-x-30, x-x-40, x-x-50,
and tippo, himt, and oatk-30 (or x-30).

Another way of plotting boxplots of X-axis is grouped by cov, and patterns are listed with each cov.
e.g.,
t-h-02, h-t-02, t-x-02, x-h-02, 
t-h-05, h-t-05, t-x-05, x-h-05, 
t-h-08, h-t-08, t-x-08, x-h-08, 
t-h-11, h-t-11, t-x-11, x-h-11, x-x-11,
t-h-20, h-t-20, t-x-20, x-h-20, x-x-20,
t-h-30, h-t-30, t-x-30, x-h-30, x-x-30,
t-h-40, h-t-40, t-x-40, x-h-40, x-x-40,
t-h-50, h-t-50, t-x-50, x-h-50, x-x-50,
and tippo, himt, and oatk-30 (or x-30).

Second, aggregate each group of the boxplots in the first.
Plotting boxplots can be listed by only patterns.
e.g.,
t-h, h-t, t-x, x-h, x-x, and tippo, himt, and oatk-30 (or x-30).
each group aggregates all in that group; 
i.e., t-h is all t-h-02, t-h-05, t-h-08, t-h-11, t-h-20, t-h-30, t-h-40, t-h-50,

Plotting can be listed by only cov.
e.g.,
02, 05, 08, 11, 20, 30, 40, 50, and tippo, himt, and oatk-30 (or x-30).
11 group includes t-h-11, h-t-11, t-x-11, x-h-11, x-x-11,

So, up until now, 4 different kinds of boxplots for each of values: total_mem_gb, total_time_hours, disk_gb, completeness_prop.

Third, this is not about boxplots. 
Y-axis is over all the 81 datasets.
X-axis is yes or no for completeness_prop is equal to 1 or at least 0.95 for one particular case. Each case is listed as follows:
t-h-02, t-h-05, t-h-08, t-h-11, t-h-20, t-h-30, t-h-40, t-h-50,
h-t-02, h-t-05, h-t-08, h-t-11, h-t-20, h-t-30, h-t-40, h-t-50,
t-x-02, t-x-05, t-x-08, t-x-11, t-x-20, t-x-30, t-x-40, t-x-50,
x-h-02, x-h-05, x-h-08, x-h-11, x-h-20, x-h-30, x-h-40, x-h-50,
x-x-11, x-x-20, x-x-30, x-x-40, x-x-50,
and tippo, himt, and oatk-30 (or x-30).

Then, a graph of grid can be 81-by-(the number of cases) or 81-by-40.
Mark one entry for yes, and no mark for no.
We want to see marks near 0 (eg, 2, 5, 8, 11) for the filtering patterns such as t-h, h-t, t-x, x-h.
Marks are some dispersed for x-x patterns and tippo, himt, and oatk-30.

Fourth, for the X-axis of the grid plot in the third, sum the number of yes cases, and draw barplots along with all pattern-cov:
t-h-02, t-h-05, t-h-08, t-h-11, t-h-20, t-h-30, t-h-40, t-h-50,
h-t-02, h-t-05, h-t-08, h-t-11, h-t-20, h-t-30, h-t-40, h-t-50,
t-x-02, t-x-05, t-x-08, t-x-11, t-x-20, t-x-30, t-x-40, t-x-50,
x-h-02, x-h-05, x-h-08, x-h-11, x-h-20, x-h-30, x-h-40, x-h-50,
x-x-11, x-x-20, x-x-30, x-x-40, x-x-50,
and tippo, himt, and oatk-30 (or x-30).
