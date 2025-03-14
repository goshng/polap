#!/bin/bash

# Function to process the input file and write results to output file
extract_edges() {
	local input_file="$1"
	local output_file="$2"

	# Ensure output file is empty before writing
	>"$output_file"

	# Process lines starting with 'P'
	awk '$1 == "P" {print $2, $3}' "$input_file" | while read -r contig edges; do
		# Extract unique edge numbers without '+' or '-'
		unique_edges=$(echo "$edges" | tr ',' '\n' | sed 's/[+-]//' | sort -u)

		# Count unique edges
		edge_count=$(echo "$unique_edges" | wc -l)

		# If exactly 3 unique edges, write to output file
		if [ "$edge_count" -eq 3 ]; then
			echo "$unique_edges" >>"$output_file"
		fi
	done
}

# Example usage:
extract_edges "$1" "$2"
