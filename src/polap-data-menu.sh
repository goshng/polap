#!/bin/bash

print_usage() {
  cat <<EOF
Usage: $0 [-c cmd] [-f file] [file] [cmd]

Options:
  -c <cmd>    Command name (default: xyz)
  -f <file>   Target shell script file (default: polap-data-v2.sh)
  -h          Show this help message

You can also provide 'cmd' and 'file' as positional arguments:
  $0 my-command target.sh

Required placeholders must exist in the target file:
  ##### INSERT_HELP_HERE #####
  ##### INSERT_FUNCTION_HERE #####
  ##### INSERT_CASE_HERE #####
  ##### INSERT_COMMAND_HERE #####
EOF
}

if [[ $# -eq 0 ]]; then
  print_usage
  exit 0
fi

# Defaults
cmd=""
file=""

# Parse options
while [[ $# -gt 0 ]]; do
  case "$1" in
  -c)
    cmd="$2"
    shift 2
    ;;
  -f)
    file="$2"
    shift 2
    ;;
  -h | --help)
    print_usage
    exit 0
    ;;
  --)
    shift
    break
    ;;
  -*)
    echo "❌ Unknown option: $1" >&2
    print_usage
    exit 1
    ;;
  *) break ;;
  esac
done

# Fallback to positional arguments
[[ -z "$file" && -n "$1" ]] && file="$1"
[[ -z "$cmd" && -n "$2" ]] && cmd="$2"

# Set defaults if still unset
cmd="${cmd:-xyz}"
file="${file:-polaplib/polap-lib-data.sh}"

# Allow the first argument to be a subcommand
if [[ ! -r "$file" ]]; then
  cmd="$file"
  file=polaplib/polap-lib-data.sh
fi

# Validate file exists and readable
if [[ ! -r "$file" ]]; then
  echo "❌ Error: File '$file' does not exist or is not readable." >&2
  exit 1
fi

# echo "cmd: $cmd"
# echo "file: $file"
# exit

# Validate required placeholder lines
for placeholder in "##### INSERT_HELP_HERE #####" \
  "##### INSERT_FUNCTION_HERE #####" \
  "##### INSERT_CASE_HERE #####"; do
  if ! grep -qF "$placeholder" "$file"; then
    echo "❌ Error: Missing required placeholder in $file: $placeholder" >&2
    exit 1
  fi
done

cmd_snake="${cmd//-/_}"

# Backup original
cp "$file" "$file.bak"

awk -v cmd="$cmd" -v cmd_snake="$cmd_snake" '
function print_help_block() {
  print "help_message_" cmd_snake "=$(";
  print "  cat <<HEREDOC";
  print "";
  print "  menu title";
  print "HEREDOC";
  print ")";
  print "";
}
function print_function_block() {
  print cmd "_genus_species_for() {";
  print "  local _brg_outdir=\"${1}\"";
  print "  echo \"Preparing archive for ${_brg_outdir} ...\"";
  print "  # tar zcf \"${_brg_outdir}-a.tar.gz\" \"${_brg_outdir}\"";
  print "}";
  print "";

  print cmd "_genus_species() {";
  print "  local _brg_outdir=\"${1:-all}\"";
	print "  local _brg_inum=\"${2:-0}\"";
	print "  local _brg_polished=\"${3:-hifiasm}\"";
	print "  local _brg_fc=\"${4:-30}\"";
  print "";
	print "  if [[ \"${_brg_outdir}\" == \"all\" ]]; then";
	print "    for _v1 in \"${Sall[@]}\"; do";
	print "      " cmd "_genus_species_for \"${_v1}\" \"${@:2}\"";
	print "    done";
	print "  elif [[ \"${_brg_outdir}\" == \"each\" ]]; then";
	print "    for _v1 in \"${Sall[@]}\"; do";
	print "      " cmd "_genus_species_for \"${_v1}\" \"${@:2}\"";
	print "    done";
	print "  else";
	print "    " cmd "_genus_species_for \"$@\"";
	print "  fi";
  print "}";
  print "";
}
function print_case_block() {
  print "    " cmd ")"
  print "      if [[ -z \"${_arg2}\" || \"${_arg2}\" == arg2 || \"${_arg2}\" == \"-h\" || \"${_arg2}\" == \"--help\" ]]; then"
  print "        echo \"Help: ${subcmd1} <outdir> [inum:0|N]\""
  print "        echo \"  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana\""
  print "        _subcmd1_clean=\"${subcmd1//-/_}\""
  print "        declare -n ref=\"help_message_${_subcmd1_clean}\""
  print "        echo \"$ref\""
  print "        exit 0"
  print "      fi"
	print "      [[ \"${_arg3}\" == arg3 ]] && _arg3=\"\""
  print "      ${subcmd1}_genus_species \"${_arg2}\""
  print "      ${subcmd1}_genus_species \"${cmd_args[@]}\""
  print "      ${subcmd1}_genus_species \"${cmd_args_ref[@]}\""
  print "      ;;"
}
function print_command_block() {
  print "  " cmd ")"
  print "    handled=1"
  print "    ;;"
}

{
  print

  if ($0 ~ /##### INSERT_HELP_HERE #####/) {
    print_help_block()
  }

  if ($0 ~ /##### INSERT_FUNCTION_HERE #####/) {
    print_function_block()
  }

  if ($0 ~ /##### INSERT_CASE_HERE #####/) {
    print_case_block()
  }

  if ($0 ~ /##### INSERT_COMMAND_HERE #####/) {
    print_command_block()
  }
}
' "$file.bak" >"$file"

if [[ -s "$file" ]]; then
  echo "✅ Inserted blocks for '${cmd}' into ${file}. Backup saved as ${file}.bak"
else
  echo "❌ Output file is empty! Restoring backup..."
  mv "$file.bak" "$file"
fi
