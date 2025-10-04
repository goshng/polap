#!/bin/bash

# Usage: ./replace_include.sh input.md > output.md
# This script replaces `!include filename.md` with the content of the file.

while IFS= read -r line; do
	if [[ "$line" =~ ^!include[[:space:]]+(.+)$ ]]; then
		file="${BASH_REMATCH[1]}"
		if [[ -f "$file" ]]; then
			cat "$file"
		else
			echo "Warning: Could not open file: $file" >&2
		fi
	else
		echo "$line"
	fi
done <"$1"
