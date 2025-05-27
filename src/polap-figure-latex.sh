#!/bin/bash

csv_file="$1"
ncol="$2"
groupby="${3:-none}" # tool | species | none
percent="${4:-22}"   # figure width % of \textwidth

if [[ -z "$csv_file" || -z "$ncol" ]]; then
  echo "Usage: $0 <input.csv> <num_columns> [groupby: tool|species|none] [percent_of_textwidth]" >&2
  exit 1
fi

width_frac=$(awk -v p="$percent" 'BEGIN { printf "%.4f", p/100 }')

echo "\\documentclass{article}"
echo "\\usepackage[margin=1cm]{geometry}"
echo "\\usepackage{graphicx}"
echo "\\usepackage{caption}"
echo "\\usepackage{float}"
echo "\\begin{document}"

groupkey=""
count=0

print_group_header() {
  if [[ "$groupby" != "none" ]]; then
    echo "\\section*{$1}"
  fi
  echo "\\begin{figure}[H]"
  echo "\\centering"
}

end_group() {
  if ((count % ncol != 0)); then
    echo "\\\\[1em]"
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

  # New group detected
  if [[ "$current_group" != "$groupkey" ]]; then
    [[ "$groupkey" != "" ]] && end_group
    groupkey="$current_group"
    print_group_header "$groupkey"
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

[[ "$groupkey" != "" ]] && end_group

echo "\\end{document}"
