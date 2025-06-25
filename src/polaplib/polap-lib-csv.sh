################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
  set -u
  return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
  return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

# Learn Bash: CSV config parsing refactor
# source polap/src/read_csv_config_dynamic.sh
# source polap/src/print_species_field_summary.sh
#
# _polap_data_csv=config.csv
# read_csv_config_dynamic
# print_species_field_summary --values --markdown --count-missing

read_csv_config_dynamic() {
  local csv_file="${csv_file:-${PWD}/${_polap_data_csv}}"
  [[ -s "$csv_file" ]] || csv_file="${_POLAPLIB_DIR}/${_polap_data_csv}"

  if [[ ! -f "$csv_file" ]]; then
    echo "[ERROR] CSV file not found: $csv_file" >&2
    return 1
  fi

  # Read and clean header
  local -a headers cleaned_headers
  IFS=',' read -r -a headers < <(head -n1 "$csv_file")
  for col in "${headers[@]}"; do
    col=$(echo "$col" | tr -d '[:space:]' | sed 's/^_//' | tr -cd '[:alnum:]_')
    cleaned_headers+=("$col")
  done

  # Declare associative arrays
  for col in "${cleaned_headers[@]}"; do
    declare -gA "_$col"
  done
  declare -gA _species

  # Parse and populate arrays
  while IFS=',' read -r -a fields; do
    local -A row
    for i in "${!cleaned_headers[@]}"; do
      row["${cleaned_headers[$i]}"]="${fields[$i]}"
    done

    local species="${row[species]}"
    [[ -z "$species" || "$species" == \#* ]] && continue

    _species["$species"]="$species"
    for col in "${cleaned_headers[@]}"; do
      eval "_${col}[\"$species\"]=\"\${row[$col]}\""
    done
  done < <(tail -n +2 "$csv_file")

  # Create Sall with all species folder names
  mapfile -t Sall < <(
    for key in "${!_long[@]}"; do
      echo "${key%%-*}"
    done | sort -u
  )
}

print_species_field_summary() {
  local output_format="csv"
  local missing_only=false
  local show_values=false
  local na_value="NA"
  local fields_file=""
  local out_file=""
  local match_pattern=""
  local sort_field=""
  local markdown=false
  local count_missing=false
  local -a selected_fields=()
  local add_field=""

  # Parse arguments
  for arg in "$@"; do
    case "$arg" in
    --tsv) output_format="tsv" ;;
    --pretty) output_format="pretty" ;;
    --csv) output_format="csv" ;;
    --markdown)
      output_format="markdown"
      markdown=true
      ;;
    --missing-only) missing_only=true ;;
    --values) show_values=true ;;
    --count-missing) count_missing=true ;;
    --na=*) na_value="${arg#--na=}" ;;
    --fields=*) IFS=',' read -r -a selected_fields <<<"${arg#--fields=}" ;;
    --fields-file=*) fields_file="${arg#--fields-file=}" ;;
    --out=*) out_file="${arg#--out=}" ;;
    --match=*) match_pattern="${arg#--match=}" ;;
    --sort-by=*) sort_field="${arg#--sort-by=}" ;;
    --add-field=*) add_field="${arg#--add-field=}" ;;
    *)
      echo "[ERROR] Unknown option: $arg" >&2
      return 1
      ;;
    esac
  done

  # Load fields from file
  if [[ -n "$fields_file" && -f "$fields_file" ]]; then
    mapfile -t selected_fields < <(grep -v '^\s*#' "$fields_file" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
  fi

  # If --add-field is used, dynamically create the field
  if [[ -n "$add_field" ]]; then
    local new_field="${add_field%%=*}"
    local default_value="${add_field#*=}"
    [[ "$new_field" != _* ]] && new_field="_$new_field"

    # Declare new associative array
    declare -gA "$new_field"

    # Populate with default value if not already set
    for s in "${!_species[@]}"; do
      if ! eval "[[ -v ${new_field}[\$s] ]]"; then
        eval "${new_field}[\$s]=\"\$default_value\""
      fi
    done
  fi

  # All associative arrays (excluding _species)
  local -a all_fields=()
  for var in $(compgen -A variable | grep '^_' | grep -v '^_species$'); do
    if declare -p "$var" 2>/dev/null | grep -q 'declare \-A'; then
      all_fields+=("$var")
    fi
  done

  # Determine fields
  local -a fields=()
  if [[ ${#selected_fields[@]} -eq 0 ]]; then
    fields=("${all_fields[@]}")
  else
    for f in "${selected_fields[@]}"; do
      [[ "$f" != _* ]] && f="_$f"
      if declare -p "$f" 2>/dev/null | grep -q 'declare \-A'; then
        fields+=("$f")
      else
        echo "[WARNING] Field $f not found." >&2
      fi
    done
  fi

  [[ ${#_species[@]} -eq 0 ]] && {
    echo "[WARNING] No species entries." >&2
    return
  }

  local sep=","
  [[ "$output_format" == "tsv" ]] && sep=$'\t'
  [[ "$markdown" == true ]] && sep='|'

  # Filter species
  local -a species_list=()
  for s in "${!_species[@]}"; do
    [[ -n "$match_pattern" && ! "$s" =~ $match_pattern ]] && continue
    species_list+=("$s")
  done

  # Sort if needed
  if [[ -n "$sort_field" ]]; then
    [[ "$sort_field" != _* ]] && sort_field="_$sort_field"
    if declare -p "$sort_field" 2>/dev/null | grep -q 'declare \-A'; then
      mapfile -t species_list < <(
        for s in "${species_list[@]}"; do
          val=$(eval "printf '%s' \"\${${sort_field}[\$s]:-$na_value}\"")
          printf "%s\t%s\n" "$val" "$s"
        done | sort -k1,1 | cut -f2
      )
    else
      echo "[WARNING] Unknown sort field: $sort_field" >&2
    fi
  fi

  # Header
  local -a lines=()
  local header="species"
  for f in "${fields[@]}"; do header+="${sep}${f#_}"; done
  lines+=("$header")

  [[ "$markdown" == true ]] && {
    local mdline="---"
    for _ in "${fields[@]}"; do mdline+="${sep}---"; done
    lines+=("$mdline")
  }

  # Missing counter
  declare -A missing_count_per_field
  for f in "${fields[@]}"; do missing_count_per_field["$f"]=0; done

  # Rows
  for s in "${species_list[@]}"; do
    local -a values=()
    local missing=0

    for f in "${fields[@]}"; do
      if eval "[[ -v ${f}[\$s] ]]"; then
        if [[ "$show_values" == true ]]; then
          val=$(eval "printf '%s' \"\${${f}[\$s]}\"")
          values+=("$val")
        else
          values+=("1")
        fi
      else
        values+=("$na_value")
        ((missing++))
        ((missing_count_per_field["$f"]++))
      fi
    done

    [[ "$missing_only" == true && $missing -eq 0 ]] && continue

    local line="$s"
    for v in "${values[@]}"; do line+="${sep}${v}"; done
    lines+=("$line")
  done

  # Count-missing summary
  if [[ "$count_missing" == true ]]; then
    if [[ "$output_format" == "pretty" ]]; then
      echo "Missing fields per column:"
      for f in "${fields[@]}"; do
        printf "%-12s %3d\n" "${f#_}" "${missing_count_per_field[$f]}"
      done
    else
      local miss_line="MissingCount"
      for f in "${fields[@]}"; do
        miss_line+="${sep}${missing_count_per_field[$f]}"
      done
      lines+=("$miss_line")
    fi
  fi

  # Output
  if [[ -n "$out_file" ]]; then
    printf "%s\n" "${lines[@]}" >"$out_file"
    echo "[INFO] Summary written to $out_file"
  else
    printf "%s\n" "${lines[@]}"
  fi
}

dump_parsed_config() {
  local output_file="${1:-config_dump.sh}"
  local arr_name s value
  local -a arrays=()
  local -a output_lines=()

  # Collect only associative arrays (excluding _species)
  for var in $(compgen -A variable | grep '^_' | grep -v '^_species$'); do
    if declare -p "$var" 2>/dev/null | grep -q 'declare \-A'; then
      arrays+=("$var")
    fi
  done

  # Add array declarations
  output_lines+=("# Auto-generated config dump")
  output_lines+=("# Source this file with: source $output_file")
  output_lines+=("")
  for arr_name in "${arrays[@]}"; do
    output_lines+=("declare -gA $arr_name")
  done
  output_lines+=("")

  # Write to file
  printf "%s\n" "${output_lines[@]}" >"$output_file"

  # Add grouped entries by species
  for s in $(printf "%s\n" "${!_species[@]}" | sort); do
    echo "# Species: $s"
    for arr_name in "${arrays[@]}"; do
      if eval "[[ -v ${arr_name}[\$s] ]]"; then
        value=$(eval "printf '%s' \"\${${arr_name}[\$s]}\"")
        printf '%s[%q]=%q\n' "$arr_name" "$s" "$value"
      fi
    done
    echo
  done >>"$output_file"

  echo "[INFO] Configuration grouped and dumped to $output_file"
}
