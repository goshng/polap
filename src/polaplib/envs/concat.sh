#!/usr/bin/env bash
# combine-environments.sh
# Version: v0.1.0
# Combine all subdirectory environment.yml files into one file with headers.

set -euo pipefail

output_file="combined-environment.yml"
>"$output_file"

# Find all environment.yml files (depth = 2)
find . -mindepth 2 -maxdepth 2 -type f -name "environment.yml" | sort | while read -r file; do
	dir_name=$(basename "$(dirname "$file")")
	echo "# ───────────────────────────────────────────────" >>"$output_file"
	echo "# From: $dir_name/environment.yml" >>"$output_file"
	echo "# ───────────────────────────────────────────────" >>"$output_file"
	cat "$file" >>"$output_file"
	echo -e "\n" >>"$output_file"
done

echo "✅ Combined file created: $output_file"
