#!/bin/bash

log_file="completed_bioprojects.txt"

l=0
for csv_file in runinfo/*.csv; do
	# l=$((l + 1))
	# if ((l > 10)); then
	# 	exit
	# fi
	base_name=$(basename "$csv_file" .csv)

	# Skip if file is empty
	if [[ ! -s "$csv_file" ]]; then
		continue
	fi

	# Create folder and copy CSV as bioproject.runinfo
	mkdir -p "$base_name"
	cp "$csv_file" "$base_name/bioproject.runinfo"

	# Run command
	bash polap/src/polap.sh get-bioproject pacbio --bioproject "$base_name" -o "$base_name"

	# Check if both long-read and short-read files exist and are non-empty
	long_read="$base_name/00-bioproject/1-sra-long-read.tsv"
	# short_read="$base_name/00-bioproject/1-sra-short-read.tsv"

	if [[ -s "$long_read" ]]; then
		echo "$base_name" >>"$log_file"
	else
		rm -rf "$base_name"
	fi
done
