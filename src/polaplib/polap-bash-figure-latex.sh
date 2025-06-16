#!/bin/bash

# Convert a CSV file and a text file to a latex file for creating a PDF file of
# assembly graphs. Assembly graph png files are necessary before executing this
# command. Such assembly graph png files can be generated using Bandage software.
# Adjust the number of columns for displaying figures using <num_columns>.
# The start page numebr is updated using <page>.
# Use <groupby> to group assembly graphs: we tested using species.
# Use <percent_of_textwidth> to adjust the size of assembly graph figure.
#
# Example:
# input.csv (group, taxon, tool name, assembly graph png filename)
# ---------
# getorganelle,Anthoceros agrestis,GetOrganelle,./Test_species/t1/0/getorganelle/embplant_pt.K115.complete.graph1.selected_graph.png
# ptgaul,Anthoceros agrestis,ptGAUL,./Test_species/t1/0/ptgaul/flye_cpONT/assembly_graph.png
#
# text.txt
# --------
# Any text that can be used as a latex document.

csv_file="$1"
text_file="$2"
ncol="${3:-2}"
startpage="${4:-1}"
groupby="${5:-none}" # tool | species | none
percent="${6:-15}"   # figure width % of \textwidth

if [[ -z "$csv_file" || -z "$ncol" ]]; then
  echo "Usage: $0 <input.csv> <text.txt> [num_columns:2] [page:1] [groupby: none|tool|species] [percent_of_textwidth:15]" >&2
  exit 1
fi

width_frac=$(awk -v p="$percent" 'BEGIN { printf "%.4f", p/100 }')

echo "\\documentclass{article}"
echo "\\usepackage[bottom=2.5cm, top=2cm, left=2cm, right=2cm]{geometry}"
echo "\\usepackage{graphicx}"
echo "\\usepackage{caption}"
echo "\\usepackage{float}"
echo "\\usepackage{fancyhdr}"
echo "\\pagestyle{fancy}"
echo "\\fancyhf{}"
echo "\\cfoot{\\thepage}"
echo "\\renewcommand{\\headrulewidth}{0pt}"
echo "\\setcounter{page}{${startpage}}"
echo "\\begin{document}"

cat "${text_file}"

groupkey=""
count=0
group_count=0
first=1

print_group_header() {
  echo "\\begin{figure}[H]"
  echo "\\centering"
  echo "\\mbox{}" # prevents LaTeX from dropping single figure
  if [[ "$groupby" != "none" ]]; then
    echo "\\noindent{\\Large\\textit{$1}}\\\\[0.8em]"
  fi
}

end_group() {
  if ((count % ncol != 0)); then
    echo "\\\\[1em]"
  fi
  if ((group_count > 1)); then
    echo "\\vspace{1em}\\hrule\\vspace{1em}"
  fi
  echo "\\end{figure}"
  count=0
}

# No sorting: process in CSV order
while IFS=',' read -r tool species caption path; do
  [[ -z "$path" ]] && continue

  # Determine group key
  case "$groupby" in
  tool) current_group="$tool" ;;
  species) current_group="$species" ;;
  *) current_group="" ;;
  esac

  # 2025-06-11: single-group problem: single group has errors.
  # New group detected
  # if [[ "$current_group" != "$groupkey" ]]; then
  #   [[ "$groupkey" != "" ]] && end_group
  #   groupkey="$current_group"
  #   print_group_header "$groupkey"
  # fi
  # New group detected
  if ((first)); then
    groupkey="$current_group"
    print_group_header "$groupkey"
    ((group_count++))
    first=0
  elif [[ "$current_group" != "$groupkey" ]]; then
    end_group
    groupkey="$current_group"
    print_group_header "$groupkey"
    ((group_count++))
  fi

  echo "\\begin{minipage}[b]{${width_frac}\\textwidth}"
  echo "\\centering"
  echo "\\includegraphics[width=\\linewidth]{$path}"
  echo "\\captionsetup{labelformat=empty}"
  echo "\\caption*{$caption}"
  echo "\\end{minipage}%"

  ((count++))
  if ((count % ncol == 0)); then
    echo "\\\\[1em]"
  else
    echo -n " "
  fi
done <"$csv_file"

end_group

echo "\\end{document}"
