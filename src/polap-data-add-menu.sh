#!/bin/bash

cmd="${1:-xyz}"
file="${2:-polap-data-v2.sh}"
cmd_snake="${cmd//-/_}"

# Backup original
cp "$file" "$file.bak"

awk -v cmd="$cmd" -v cmd_snake="$cmd_snake" '
function print_help_block() {
  print "help_message_" cmd_snake "=$(";
  print "  cat <<HEREDOC";
  print "";
  print "  get ${_local_host}:$PWD/<outdir>-a.tar.gz";
  print "HEREDOC";
  print ")";
  print "";  # <- blank line after help block
}
function print_function_block() {
  print cmd "_genus_species() {";
  print "  local _brg_outdir=\"${1}\"";
  print "  echo \"Preparing archive for ${_brg_outdir} ...\"";
  print "  # tar zcf \"${_brg_outdir}-a.tar.gz\" \"${_brg_outdir}\"";
  print "}";
  print "";  # <- blank line after function block
}
function print_case_block() {
  print "    " cmd ")"
  print "      if [[ -z \"${_arg2}\" || \"${_arg2}\" == arg2 ]]; then"
  print "        echo \"Help: ${subcmd1} <outdir>\""
  print "        echo \"  ${0##*/} ${subcmd1} Arabidopsis_thaliana\""
  print "        _subcmd1_clean=\"${subcmd1//-/_}\""
  print "        declare -n ref=\"help_message_${_subcmd1_clean}\""
  print "        echo \"$ref\""
  print "        exit 0"
  print "      fi"
  print "      ${subcmd1}_genus_species \"${_arg2}\""
  print "      ;;"
}

{
  print  # print the original line

  if ($0 ~ /##### INSERT_HELP_HERE #####/) {
    print_help_block()
  }

  if ($0 ~ /##### INSERT_FUNCTION_HERE #####/) {
    print_function_block()
  }

  if ($0 ~ /##### INSERT_CASE_HERE #####/) {
    print_case_block()
  }
}
' "$file.bak" >"$file"

if [[ -s "$file" ]]; then
  echo "✅ Inserted blocks for '${cmd}' into ${file}. Backup saved as ${file}.bak"
else
  echo "❌ Output file is empty! Restoring backup..."
  mv "$file.bak" "$file"
fi
