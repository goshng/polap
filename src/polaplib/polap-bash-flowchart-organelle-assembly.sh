#!/usr/bin/env bash
# polap-bash-organelle-assembly-flow.sh
# Render a flowchart describing the four organelle assembly cases in polap/bolap.
# Requirements: graphviz (dot)
# Usage:
#   polap-bash-organelle-assembly-flow.sh -o outdir [-f png|pdf|svg] [--theme light|dark] [--title "Custom Title"]
#   polap-bash-organelle-assembly-flow.sh --help

set -euo pipefail

# Defaults
OUTDIR="."
FMT="png"     # png|pdf|svg
THEME="light" # light|dark
TITLE="Plant Organelle Assembly in polap (with bolap benchmarking)"
ENGINE="dot" # dot|neato|fdp (dot recommended)
BASENAME="polap_organelle_assembly"
GV_FILE=""
OPEN_AFTER=0

print_help() {
	cat <<'EOF'
Render the polap organelle assembly flowchart.

Options:
  -o, --outdir DIR        Output directory (default: .)
  -f, --format FMT        Output format: png|pdf|svg (default: png)
  -t, --theme THEME       Theme: light|dark (default: light)
  --title STR             Title text (default shown)
  --engine NAME           Graphviz engine: dot|neato|fdp (default: dot)
  -n, --name BASENAME     Base name for outputs (default: polap_organelle_assembly)
  --open                  Try to open the rendered file after build (macOS: open; Linux: xdg-open)
  -h, --help              Show this help
EOF
}

# Parse CLI
while [[ $# -gt 0 ]]; do
	case "$1" in
	-o | --outdir)
		OUTDIR="$2"
		shift 2
		;;
	-f | --format)
		FMT="$2"
		shift 2
		;;
	-t | --theme)
		THEME="$2"
		shift 2
		;;
	--title)
		TITLE="$2"
		shift 2
		;;
	--engine)
		ENGINE="$2"
		shift 2
		;;
	-n | --name)
		BASENAME="$2"
		shift 2
		;;
	--open)
		OPEN_AFTER=1
		shift
		;;
	-h | --help)
		print_help
		exit 0
		;;
	*)
		echo "Unknown argument: $1" >&2
		exit 1
		;;
	esac
done

mkdir -p "$OUTDIR"
GV_FILE="${OUTDIR}/${BASENAME}.gv"
OUT_FILE="${OUTDIR}/${BASENAME}.${FMT}"

# Theme
if [[ "$THEME" == "dark" ]]; then
	BG="#0e1116"
	FG="#e6edf3"
	EDGE="#9da7b3"
	ACC="#6cb6ff"
	ACC2="#7ee787"
	ACC3="#ffab70"
else
	BG="white"
	FG="#0b1221"
	EDGE="#5b6b7a"
	ACC="#2563eb"
	ACC2="#059669"
	ACC3="#b45309"
fi

# Check graphviz
if ! command -v dot >/dev/null 2>&1; then
	echo "Error: graphviz 'dot' is required on PATH." >&2
	exit 2
fi

cat >"$GV_FILE" <<EOF
digraph polap_assembly {
  graph [label="${TITLE}\\n(polap wrapper; bolap benchmarking)",
         labelloc="t", fontsize=20, fontname="Helvetica",
         bgcolor="${BG}", pad="0.2", ranksep="0.5", nodesep="0.35"];
  node  [shape=rectangle, style="rounded,filled", color="${EDGE}", fontcolor="${FG}",
         fontname="Helvetica", fillcolor="${BG}", penwidth=1.5];
  edge  [color="${EDGE}", penwidth=1.4, arrowsize=0.9];

  # Legend
  subgraph cluster_legend {
    label="Legend"; fontsize=12; color="${EDGE}";
    l1 [label="Flye v2.9.6", shape=box, fillcolor="${ACC}", fontcolor="white", style="rounded,filled"];
    l2 [label="Miniasm (seed contigs)", shape=box, fillcolor="${ACC2}", fontcolor="white", style="rounded,filled"];
    l3 [label="Oatk (HiFi-centric)", shape=box, fillcolor="${ACC3}", fontcolor="white", style="rounded,filled"];
  }

  # Inputs
  start [label="Start", shape=circle, fillcolor="${ACC}", fontcolor="white"];
  ont   [label="ONT reads", shape=folder];
  hifi  [label="HiFi reads", shape=folder];
  ptref [label="Protein homology panel\\n(e.g., miniprot/minimap2)"];
  nukes [label="Remove nuclear reads\\n(e.g., BUSCO/miniprot/minimap2)"];

  # ONT -> ptDNA
  subgraph cluster_pt_ont {
    label="Case 1: ptDNA from ONT"; color="${EDGE}";
    s1   [label="Select pt-origin reads\\nvia protein homology", shape=box, fillcolor="${BG}"];
    s2   [label="Assemble ptDNA\\nFlye v2.9.6", shape=box, fillcolor="${ACC}", fontcolor="white"];
  }

  # ONT -> mtDNA
  subgraph cluster_mt_ont {
    label="Case 2: mtDNA from ONT"; color="${EDGE}";
    m1   [label="Filter pt & nuclear reads\\n(minimize contamination)", shape=box];
    m2   [label="Seed contigs\\nMiniasm", shape=box, fillcolor="${ACC2}", fontcolor="white"];
    m3   [label="Assemble mtDNA\\nFlye v2.9.6", shape=box, fillcolor="${ACC}", fontcolor="white"];
  }

  # HiFi -> ptDNA
  subgraph cluster_pt_hifi {
    label="Case 3: ptDNA from HiFi"; color="${EDGE}";
    h1a  [label="Option A: protein-homology read selection\\nthen Flye v2.9.6", shape=box, fillcolor="${ACC}", fontcolor="white"];
    h1b  [label="Option B: Oatk (syncasm + pathfinder)", shape=box, fillcolor="${ACC3}", fontcolor="white"];
  }

  # HiFi -> mtDNA
  subgraph cluster_mt_hifi {
    label="Case 4: mtDNA from HiFi"; color="${EDGE}";
    h2   [label="Oatk (syncasm + pathfinder)", shape=box, fillcolor="${ACC3}", fontcolor="white"];
  }

  # Edges
  start -> ont; start -> hifi;

  # ONT splits
  ont -> s1;
  s1 -> s2;

  ont -> m1;
  ptref -> s1;
  nukes -> m1;

  m1 -> m2 -> m3;

  # HiFi splits
  hifi -> h1a;
  hifi -> h1b;
  hifi -> h2;

  # Rank hints
  {rank=same; ont; hifi;}
}
EOF

# Render
case "$ENGINE" in
dot | neato | fdp) : ;;
*)
	echo "Unsupported engine: $ENGINE" >&2
	exit 3
	;;
esac

CMD=("$ENGINE" "-T${FMT}" "$GV_FILE" "-o" "$OUT_FILE")
echo "[INFO] Rendering: ${CMD[*]}" >&2
"${CMD[@]}"

echo "[OK] Wrote: $GV_FILE"
echo "[OK] Wrote: $OUT_FILE"

if ((OPEN_AFTER)); then
	if command -v xdg-open >/dev/null 2>&1; then xdg-open "$OUT_FILE" >/dev/null 2>&1 || true; fi
	if command -v open >/dev/null 2>&1; then open "$OUT_FILE" || true; fi
fi
