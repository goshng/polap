#!/bin/bash

# Script to generate a Table of Contents for a README.md file

INPUT_FILE="README.md"
OUTPUT_FILE="README.md"
TOC_START="<!-- TOC START -->"
TOC_END="<!-- TOC END -->"

# Check if README.md exists
if [[ ! -f $INPUT_FILE ]]; then
	echo "Error: $INPUT_FILE not found."
	exit 1
fi

# Function to create an anchor from a heading
generate_anchor() {
	local heading="$1"
	# Convert to lowercase, replace spaces with hyphens, remove special characters
	echo "$heading" | tr '[:upper:]' '[:lower:]' | sed -E 's/[[:space:]]+/-/g; s/[^a-z0-9\-]//g'
}

# Extract headings and build TOC
TOC=""
while IFS= read -r line; do
	if [[ $line =~ ^(#{1,6})[[:space:]](.+) ]]; then
		heading_level="${#BASH_REMATCH[1]}" # Number of # characters
		heading_text="${BASH_REMATCH[2]}"   # Heading text
		anchor=$(generate_anchor "$heading_text")
		indent=$(printf '  %.0s' $(seq 1 $((heading_level - 1)))) # Indentation based on level
		TOC+="$indent- [${heading_text}](#$anchor)\n"
	fi
done <"$INPUT_FILE"

# Add TOC markers
TOC_HEADER="## Table of Contents"
TOC_FULL="$TOC_HEADER\n\n$TOC"

# Replace or insert TOC in README.md
if grep -q "$TOC_START" "$INPUT_FILE"; then
	# Replace existing TOC
	sed -i "/$TOC_START/,/$TOC_END/c\\$TOC_START\n$TOC_FULL\n$TOC_END" "$OUTPUT_FILE"
else
	# Prepend TOC at the top
	echo -e "$TOC_START\n$TOC_FULL\n$TOC_END\n" | cat - "$INPUT_FILE" >temp && mv temp "$OUTPUT_FILE"
fi

echo "Table of Contents has been generated and added to $OUTPUT_FILE."
