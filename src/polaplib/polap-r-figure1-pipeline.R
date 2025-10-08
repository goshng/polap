#!/usr/bin/env Rscript
# polap-r-figure1-pipeline.R
# Version : v1.0.5  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build Figure 1 (vertical layout) with unicode icons & labeled arrows.
# Outputs:
#   man/figures/figure1-pipeline.svg
#   man/figures/figure1-pipeline.pdf (if rsvg present)

suppressPackageStartupMessages({
  library(DiagrammeR)
  library(DiagrammeRsvg)
})

# ---- Output paths -----------------------------------------------------------
out_dir <- file.path(getwd(), "man", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_svg <- file.path(out_dir, "figure1-pipeline.svg")
out_pdf <- file.path(out_dir, "figure1-pipeline.pdf")

# ---- Colors -----------------------------------------------------------------
C1 <- "#5DADE2" # Data
C2 <- "#48C9B0" # Processing
C3 <- "#58D68D" # Assembly
C4 <- "#F5B041" # Annotation
C5 <- "#AF7AC5" # Reporting

# ---- DOT graph (vertical) ---------------------------------------------------
# Use HTML-like labels so we can bold and add line breaks; include emoji icons.
dot <- sprintf('
digraph polap {
  graph [rankdir=TB, nodesep=0.6, ranksep=0.8, fontname="Helvetica"];
  node  [shape=box, style="rounded,filled", color="#555555",
         fontname="Helvetica", fontsize=12, width=3.8, penwidth=1.2];
  edge  [color="#555555", penwidth=1.2, arrowsize=0.8,
         fontname="Helvetica", fontsize=11];

  N1 [label=< <b>📥 1. Data Input</b><br/>FASTQ + manifest<br/>polap get, download<br/>SRRxxxx.fastq.gz, manifest.json >,
      fillcolor="%s"];
  N2 [label=< <b>🧹 2. Read Processing</b><br/>SeqKit, minimap2<br/>filtered-reads.fq.gz<br/>tableS1-dataset-summary.tsv >,
      fillcolor="%s"];
  N3 [label=< <b>🧬 3. Assembly</b><br/>Flye + Bandage<br/>pt.1.gfa / mt.1.gfa<br/>polap-readassemble-1-miniasm >,
      fillcolor="%s"];
  N4 [label=< <b>📑 4. Annotation</b><br/>anno-scan (AWK sums)<br/>anno-pt.csv / anno-mt.csv<br/>contig-annotation-depth >,
      fillcolor="%s"];
  N5 [label=< <b>📊 5. Reporting</b><br/>sheet-ptmt (PT/MT pages)<br/>table-s1, benchmarks >,
      fillcolor="%s"];

  N1 -> N2 [label="FASTQ + manifest"];
  N2 -> N3 [label="filtered reads"];
  N3 -> N4 [label="GFA + depth"];
  N4 -> N5 [label="CSV + manifest"];

}
', C1, C2, C3, C4, C5)

g <- grViz(dot)

# ---- Export SVG -------------------------------------------------------------
svg_txt <- DiagrammeRsvg::export_svg(g)
writeLines(svg_txt, out_svg)

# ---- Export PDF (optional) --------------------------------------------------
ok <- FALSE
if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_pdf(charToRaw(svg_txt), file = out_pdf)
  ok <- TRUE
}

cat("[OK] Wrote:", out_svg, if (ok) paste("and", out_pdf) else "(install rsvg for PDF)", "\n")
