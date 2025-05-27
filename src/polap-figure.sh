#!/bin/bash

# Usage: ./generate_md_grid.sh input.csv 4 > output.md

csv_file="$1"
ncol="$2"

if [[ -z "$csv_file" || -z "$ncol" ]]; then
	echo "Usage: $0 <input.csv> <num_columns>" >&2
	exit 1
fi

echo
echo "# Grid of Figures"
echo

rowcount=0

# Print header row (empty to suppress Pandoc's table headers)
for ((i = 0; i < ncol; i++)); do
	printf "|           "
done
echo "|"

for ((i = 0; i < ncol; i++)); do
	printf "|:----------:"
done
echo "|"

# Read CSV and output rows of images and captions
while IFS=',' read -r tool species caption path; do
	[[ -z "$path" ]] && continue # skip empty lines

	# Generate image with optional sizing
	printf "| ![](%s){width=100%%} " "$path"

	((rowcount++))
	if ((rowcount % ncol == 0)); then
		echo "|"

		# Now print caption row
		for ((i = rowcount - ncol; i < rowcount; i++)); do
			# Reread the same block from the CSV
			line=$(sed -n "$((i + 1))p" "$csv_file")
			cap=$(echo "$line" | awk -F',' '{print $3}')
			printf "| **%s** " "$cap"
		done
		echo "|"

		# Print alignment row again (for Pandoc markdown)
		for ((i = 0; i < ncol; i++)); do
			printf "|:----------:"
		done
		echo "|"
	fi
done <"$csv_file"

# Print trailing incomplete row if needed
remainder=$((rowcount % ncol))
if ((remainder > 0)); then
	for ((i = remainder; i < ncol; i++)); do
		printf "|           "
	done
	echo "|"

	# Caption row for last incomplete line
	for ((i = rowcount - remainder; i < rowcount; i++)); do
		line=$(sed -n "$((i + 1))p" "$csv_file")
		cap=$(echo "$line" | awk -F',' '{print $3}')
		printf "| **%s** " "$cap"
	done
	for ((i = remainder; i < ncol; i++)); do
		printf "|           "
	done
	echo "|"
fi
