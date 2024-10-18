#!/bin/bash

# Define your files
file1="$1" # File to prepend

# Create a temporary file that will hold the combined content
temp_file=$(mktemp)

# Combine the content of file1 and file2 into the temporary file
#

cat <(head -16 run-polap-function-template.sh) "$file1" >"$temp_file"

# Replace file2 with the new combined content
mv "$temp_file" "$file1"

echo "$file1" now has a header of the POLAP GPLv3 copyright.
