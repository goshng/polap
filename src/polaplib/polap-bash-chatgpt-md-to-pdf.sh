#!/usr/bin/env bash
# File: chatgpt-md-to-pdf.sh
# Version: v0.1.0
#
# Convert a ChatGPT-style Markdown file to PDF using pandoc.
#
# Usage:
#   chatgpt-md-to-pdf.sh INPUT.md [OUTPUT.pdf]
#
# If OUTPUT.pdf is omitted, it will be derived from INPUT.md
# (e.g., chat.md -> chat.pdf).

set -Eeuo pipefail

usage() {
	cat <<EOF
Usage: $(basename "$0") INPUT.md [OUTPUT.pdf]

Convert a ChatGPT Markdown export to PDF using pandoc.

Examples:
  $(basename "$0") chat.md
  $(basename "$0") chat.md chat-conversation.pdf

Environment variables:
  PDF_ENGINE   Pandoc PDF engine (default: xelatex)
  PANDOC_OPTS  Extra options passed to pandoc (default: basic layout)
EOF
}

# --- Check arguments --------------------------------------------------------

if [[ $# -lt 1 || $# -gt 2 ]]; then
	usage
	exit 1
fi

in_md=$1
if [[ ! -f "$in_md" ]]; then
	echo "Error: input file not found: $in_md" >&2
	exit 1
fi

# Derive output PDF if not given
if [[ $# -eq 2 ]]; then
	out_pdf=$2
else
	base=${in_md##*/} # strip path
	base=${base%.*}   # strip extension
	out_pdf="${base}.pdf"
fi

# --- Check dependencies -----------------------------------------------------

if ! command -v pandoc >/dev/null 2>&1; then
	echo "Error: pandoc not found in PATH." >&2
	exit 1
fi

# PDF engine can be overridden from the environment
PDF_ENGINE=${PDF_ENGINE:-xelatex}

# Reasonable defaults; can be overridden with PANDOC_OPTS
: "${PANDOC_OPTS:=\
--from markdown --pdf-engine=${PDF_ENGINE} -V geometry:margin=2cm -V fontsize=11pt -V papersize=a4paper }"

# --- Run pandoc -------------------------------------------------------------

echo "Converting:"
echo "  Input : $in_md"
echo "  Output: $out_pdf"
echo "  Engine: $PDF_ENGINE"
echo

# shellcheck disable=SC2086
pandoc $PANDOC_OPTS "$in_md" -o "$out_pdf"

echo "Done: $out_pdf"
