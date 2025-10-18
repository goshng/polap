#!/usr/bin/env awk -f
# scripts/extract_genus_species_prefixes.awk
# Version: v0.2.1
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Purpose:
#   Extract function name prefixes before `_genus_species()` from Bash files.
#   Handles:
#     - optional "function" keyword
#     - hyphens (-) and underscores (_) in names
#     - regex filter via -v re='REGEX' (empty = match all)
#
# Usage:
#   awk -f scripts/extract_genus_species_prefixes.awk [-v re='regex'] file...
#
# Example:
#   awk -v re="^run" -f scripts/extract_genus_species_prefixes.awk polap-lib-data.sh

BEGIN {
  OFS = "\n"
  # If 're' not provided, default to empty → match all
  if (re == "") re = ""
}

# Return 1 if name matches re (or re is empty).
function match_name(name,   norm) {
  if (re == "") return 1
  # also try underscore-normalized variant so queries match both hyphen/underscore forms
  norm = name
  gsub(/-/, "_", norm)
  return (name ~ re) || (norm ~ re)
}

# Match function headers:
#   NAME_genus_species() {
#   function NAME_genus_species() {
# NAME may contain letters, digits, underscore, or hyphen.
match($0, /^[[:space:]]*(function[[:space:]]+)?([[:alnum:]_-]+)_genus_species[[:space:]]*\(\)[[:space:]]*\{/, m) {
  name = m[2]
  if (match_name(name)) print name
}
