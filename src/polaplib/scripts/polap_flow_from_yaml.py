#!/usr/bin/env python3
"""
Read a minimal YAML and emit a Graphviz .gv with the same topology used above.
Usage:
  python scripts/polap_flow_from_yaml.py config.yaml out.gv
YAML keys (optional):
  title: str
  theme: light|dark
  flye_version: "2.9.6"
"""
import sys, yaml


def theme_colors(theme):
    if theme == "dark":
        return dict(
            BG="#0e1116",
            FG="#e6edf3",
            EDGE="#9da7b3",
            ACC="#6cb6ff",
            ACC2="#7ee787",
            ACC3="#ffab70",
        )
    return dict(
        BG="white",
        FG="#0b1221",
        EDGE="#5b6b7a",
        ACC="#2563eb",
        ACC2="#059669",
        ACC3="#b45309",
    )


def main():
    if len(sys.argv) != 3:
        print("Usage: polap_flow_from_yaml.py config.yaml out.gv", file=sys.stderr)
        sys.exit(1)
    cfg = (
        yaml.safe_load(open(sys.argv[1]))
        if sys.argv[1] != "-"
        else yaml.safe_load(sys.stdin)
    )
    title = cfg.get(
        "title", "Plant Organelle Assembly in polap (with bolap benchmarking)"
    )
    flye_ver = cfg.get("flye_version", "2.9.6")
    colors = theme_colors(cfg.get("theme", "light"))
    gv = f"""digraph polap_assembly {{
  graph [label="{title}\\n(polap wrapper; bolap benchmarking)",
         labelloc="t", fontsize=20, fontname="Helvetica",
         bgcolor="{colors['BG']}", pad="0.2", ranksep="0.5", nodesep="0.35"];
  node  [shape=rectangle, style="rounded,filled", color="{colors['EDGE']}", fontcolor="{colors['FG']}",
         fontname="Helvetica", fillcolor="{colors['BG']}", penwidth=1.5];
  edge  [color="{colors['EDGE']}", penwidth=1.4, arrowsize=0.9];

  subgraph cluster_legend {{
    label="Legend";
    l1 [label="Flye v{flye_ver}", shape=box, fillcolor="{colors['ACC']}", fontcolor="white", style="rounded,filled"];
    l2 [label="Miniasm (seed contigs)", shape=box, fillcolor="{colors['ACC2']}", fontcolor="white", style="rounded,filled"];
    l3 [label="Oatk (HiFi-centric)", shape=box, fillcolor="{colors['ACC3']}", fontcolor="white", style="rounded,filled"];
  }}

  start [label="Start", shape=circle, fillcolor="{colors['ACC']}", fontcolor="white"];
  ont   [label="ONT reads", shape=folder];
  hifi  [label="HiFi reads", shape=folder];
  ptref [label="Protein homology panel\\n(e.g., miniprot/minimap2)"];
  nukes [label="Remove nuclear reads\\n(e.g., BUSCO/miniprot/minimap2)"];

  subgraph cluster_pt_ont {{
    label="Case 1: ptDNA from ONT";
    s1   [label="Select pt-origin reads\\nvia protein homology"];
    s2   [label="Assemble ptDNA\\nFlye v{flye_ver}", shape=box, fillcolor="{colors['ACC']}", fontcolor="white"];
  }}

  subgraph cluster_mt_ont {{
    label="Case 2: mtDNA from ONT";
    m1   [label="Filter pt & nuclear reads\\n(minimize contamination)"];
    m2   [label="Seed contigs\\nMiniasm", shape=box, fillcolor="{colors['ACC2']}", fontcolor="white"];
    m3   [label="Assemble mtDNA\\nFlye v{flye_ver}", shape=box, fillcolor="{colors['ACC']}", fontcolor="white"];
  }}

  subgraph cluster_pt_hifi {{
    label="Case 3: ptDNA from HiFi";
    h1a  [label="Option A: protein-homology read selection\\nthen Flye v{flye_ver}", shape=box, fillcolor="{colors['ACC']}", fontcolor="white"];
    h1b  [label="Option B: Oatk (syncasm + pathfinder)", shape=box, fillcolor="{colors['ACC3']}", fontcolor="white"];
  }}

  subgraph cluster_mt_hifi {{
    label="Case 4: mtDNA from HiFi";
    h2   [label="Oatk (syncasm + pathfinder)", shape=box, fillcolor="{colors['ACC3']}", fontcolor="white"];
  }}

  start -> ont; start -> hifi;
  ont -> s1; s1 -> s2;
  ont -> m1; ptref -> s1; nukes -> m1; m1 -> m2 -> m3;
  hifi -> h1a; hifi -> h1b; hifi -> h2;
}}
"""
    open(sys.argv[2], "w").write(gv)


if __name__ == "__main__":
    main()
