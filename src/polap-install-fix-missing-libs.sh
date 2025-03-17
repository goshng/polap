#!/bin/bash

set -e # Exit on error

if [[ $# -ne 1 ]]; then
	echo "Usage: $0 <executable>"
	exit 1
fi

EXECUTABLE=$1

if [[ ! -f "$EXECUTABLE" || ! -x "$EXECUTABLE" ]]; then
	echo "Error: File '$EXECUTABLE' not found or not executable."
	exit 1
fi

echo "Checking missing libraries for: $EXECUTABLE"

# Find missing libraries
MISSING_LIBS=$(ldd "$EXECUTABLE" | awk '/not found/ {print $1}')

if [[ -z "$MISSING_LIBS" ]]; then
	echo "✅ No missing libraries found!"
	exit 0
fi

echo "⚠️  Missing libraries detected:"
echo "$MISSING_LIBS"

LIB_PATHS=()
for LIB in $MISSING_LIBS; do
	FOUND_PATH=$(ldconfig -p | grep "$LIB" | awk '{print $NF}' | head -n 1)
	if [[ -n "$FOUND_PATH" ]]; then
		echo "✔ Found $LIB at $FOUND_PATH"
		LIB_PATHS+=("$(dirname "$FOUND_PATH")")
	else
		echo "❌ Library $LIB not found in system paths!"
	fi
done

# Remove duplicates and construct new RPATH
NEW_RPATH=$(echo "${LIB_PATHS[@]}" | tr ' ' ':')

if [[ -n "$NEW_RPATH" ]]; then
	echo "🔗 Setting new RPATH: $NEW_RPATH"
	patchelf --set-rpath "$NEW_RPATH" "$EXECUTABLE"
	echo "✅ RPATH updated successfully!"
else
	echo "⚠️  No valid library paths found. Cannot update RPATH."
fi
