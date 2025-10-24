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
# Convert numbers between different units.
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

#!/usr/bin/env bash
# polap-lib-dataset.sh
# JSON-based dataset manager for polap-dataset-read.json

DATASET_FILE="${_brg_dataset:-dataset.json}"
[[ -f "$DATASET_FILE" ]] || {
	echo "Error: $DATASET_FILE not found" >&2
	return 1
}

########################################
# Load dataset.json into variable
########################################
_dataset_load() {
	if [[ ! -f "$DATASET_FILE" ]]; then
		echo "{}" >"$DATASET_FILE"
	fi
}

########################################
# Save updated JSON to file
########################################
_dataset_save() {
	local tmp_file="${DATASET_FILE}.tmp"
	cat >"$tmp_file"
	mv "$tmp_file" "$DATASET_FILE"
}

########################################
# Backup
########################################
_dataset_backup() {
	if [[ -s "$DATASET_FILE" ]]; then
		cp -p "$DATASET_FILE" "${DATASET_FILE}.bak.$(date +%Y%m%d%H%M%S)"
	else
		echo "No file to backup"
	fi
}

########################################
# Add new record
########################################
dataset-add-record() {
	local species="$1"
	shift
	_dataset_load

	if jq -e --arg s "$species" 'has($s)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Species $species already exists." >&2
		return 1
	fi

	local json_obj="{}"
	for kv in "$@"; do
		local k="${kv%%=*}"
		local v="${kv#*=}"
		json_obj=$(jq --arg k "$k" --arg v "$v" '. + {($k): $v}' <<<"$json_obj")
	done

	jq --arg s "$species" --argjson obj "$json_obj" '. + {($s): $obj}' "$DATASET_FILE" | _dataset_save
	echo "[INFO] Added $species"
}

########################################
# Delete a record
########################################
dataset-delete-record() {
	local species="$1"
	_dataset_load
	_dataset_backup

	jq --arg s "$species" 'del(.[$s])' "$DATASET_FILE" | _dataset_save
	echo "[INFO] Deleted $species"
}

########################################
# Search
########################################
dataset-search() {
	local match="$1"
	jq --arg m "$match" 'to_entries | map(select(.key | test($m))) | from_entries' "$DATASET_FILE"
}

########################################
# Update a field in one or more species
########################################
dataset-update-field() {
	local field="$1"
	local value="$2"
	local match="${3:-.*}"
	_dataset_load
	_dataset_backup

	jq --arg f "$field" --arg v "$value" --arg m "$match" '
    with_entries(if (.key | test($m)) then .value[$f] = $v | . else . end)
  ' "$DATASET_FILE" | _dataset_save
	echo "[INFO] Updated field $field for species matching $match"
}

########################################
# Add a new field to all records
########################################
dataset-add-field() {
	local field="$1"
	local default="${2:-null}"
	_dataset_load
	_dataset_backup

	jq --arg f "$field" --arg v "$default" '
    with_entries(.value[$f] = $v)
  ' "$DATASET_FILE" | _dataset_save
	echo "[INFO] Added field $field to all records"
}

########################################
# Export JSON dataset back to CSV
########################################
dataset-export-csv() {
	_dataset_load
	local fields=$(jq -r 'to_entries | map(.value | keys) | flatten | unique | .[]' "$DATASET_FILE")
	local header="species"
	for f in $fields; do header+="|$f"; done
	echo "$header" | tr '|' ','

	jq -r --arg fields "$fields" '
    to_entries[] | [.key] + [.value[$fields | split(" ")][]] | @csv
  ' "$DATASET_FILE"
}

########################################
# REDEFINED: Add a new record using a template (corrected grammar)
# Usage:
#  A) dataset-add-template NEW_ID TEMPLATE_ID [field=value ...]
#  B) dataset-add-template --auto-increment BASE_GENUS_SPECIES TEMPLATE_ID [field=value ...]
########################################
dataset-add-template() {
	local new_id=""
	local template_id=""

	_dataset_load

	if [[ "$1" == "--auto-increment" ]]; then
		shift
		local base="$1"
		shift || true
		template_id="$1"
		shift || true
		if [[ -z "$base" || -z "$template_id" ]]; then
			echo "Usage: dataset-add-template --auto-increment BASE_GENUS_SPECIES TEMPLATE_ID [field=value ...]" >&2
			return 1
		fi
		local max_index
		max_index=$(jq -r 'keys[]' "$DATASET_FILE" | grep -E "^${base}-[0-9]+$" | sed 's/.*-//' | sort -nr | head -n1)
		[[ -z "$max_index" ]] && max_index=-1
		local next_index=$((max_index + 1))
		new_id="${base}-${next_index}"
		echo "[INFO] Auto-generated species ID: $new_id"
	else
		new_id="$1"
		shift || true
		template_id="$1"
		shift || true
		if [[ -z "$new_id" || -z "$template_id" ]]; then
			echo "Usage: dataset-add-template NEW_ID TEMPLATE_ID [field=value ...]" >&2
			return 1
		fi
	fi

	if jq -e --arg n "$new_id" 'has($n)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Target $new_id already exists." >&2
		return 1
	fi
	if ! jq -e --arg t "$template_id" 'has($t)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Template $template_id not found." >&2
		return 1
	fi

	local clone
	clone=$(jq --arg t "$template_id" '.[$t]' "$DATASET_FILE")

	local kv k v
	for kv in "$@"; do
		k="${kv%%=*}"
		v="${kv#*=}"
		clone=$(jq --arg k "$k" --arg v "$v" '. + {($k): $v}' <<<"$clone")
	done

	jq --arg n "$new_id" --argjson obj "$clone" '. + {($n): $obj}' "$DATASET_FILE" | _dataset_save
	echo "[INFO] Created $new_id from $template_id"
}

########################################
# REDEFINED: Rename a species record with auto-updates
########################################
dataset-rename-species() {
	local old_id="$1"
	local new_id="$2"
	_dataset_load

	if [[ -z "$old_id" || -z "$new_id" ]]; then
		echo "Usage: dataset-rename-species OLD_ID NEW_ID" >&2
		return 1
	fi
	if ! jq -e --arg o "$old_id" 'has($o)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Species $old_id not found." >&2
		return 1
	fi
	if jq -e --arg n "$new_id" 'has($n)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Target $new_id already exists." >&2
		return 1
	fi

	local base species2
	base="${new_id%-*}"
	species2="${base//_/ }"

	jq --arg o "$old_id" --arg n "$new_id" --arg base "$base" --arg species2 "$species2" '
    . + {($n): .[$o]} | del(.[$o]) |
    .[$n].species = $n |
    .[$n].folder = $base |
    .[$n].species2 = $species2
  ' "$DATASET_FILE" | _dataset_save

	echo "[INFO] Renamed $old_id -> $new_id (updated fields: species, folder, species2)"
}

########################################
# Convert CSV to JSON dataset (with --infer-types support)
########################################
dataset-import-csv() {
	# Usage:
	#   dataset-import-csv [--infer-types] [--key-field=species] path/to/file.csv
	# Notes:
	#   - --infer-types: convert numeric, boolean, and NA/null blanks to proper JSON types
	#   - --key-field: which header column is used as the JSON key (default: species)
	#   - Existing JSON will be backed up before import
	local infer=false
	local key_field="species"
	while [[ "$1" == --* ]]; do
		case "$1" in
		--infer-types)
			infer=true
			shift
			;;
		--key-field=*)
			key_field="${1#--key-field=}"
			shift
			;;
		*)
			echo "[ERROR] Unknown option: $1" >&2
			return 1
			;;
		esac
	done

	local csv_file="$1"
	[[ -f "$csv_file" ]] || {
		echo "[ERROR] File not found: $csv_file" >&2
		return 1
	}

	_dataset_backup
	echo "[INFO] Importing CSV: $csv_file (key=$key_field infer=$infer)"
	jq -n '{}' >"$DATASET_FILE"

	# Read header
	IFS=',' read -r -a headers < <(head -n1 "$csv_file")

	# Locate key column index
	local key_idx=-1
	for i in "${!headers[@]}"; do
		h="${headers[$i]//[$'\r\n']/}"
		if [[ "$h" == "$key_field" ]]; then
			key_idx=$i
			break
		fi
	done
	if ((key_idx < 0)); then
		echo "[ERROR] Key field '$key_field' not found in header." >&2
		return 1
	fi

	tail -n +2 "$csv_file" | while IFS=',' read -r -a fields; do
		local rec_key="${fields[$key_idx]}"
		[[ -z "$rec_key" ]] && continue

		# Build JSON object for this row
		local obj='{}'
		for i in "${!headers[@]}"; do
			local k="${headers[$i]//[$'\r\n']/}"
			local v="${fields[$i]}"

			# Choose set expression based on infer flag
			if $infer; then
				obj=$(jq -c --arg k "$k" --arg v "$v" '
          def parse($s):
            if ($s|length)==0 or ($s|ascii_downcase) == "na" or ($s|ascii_downcase)=="null" then null
            elif ($s|ascii_downcase)=="true" then true
            elif ($s|ascii_downcase)=="false" then false
            elif ($s|test("^-?[0-9]+$")) then ($s|tonumber)
            elif ($s|test("^-?[0-9]*\\.[0-9]+$")) then ($s|tonumber)
            else $s end;
          . + {($k): (parse($v))}
        ' <<<"$obj")
			else
				obj=$(jq -c --arg k "$k" --arg v "$v" '. + {($k): $v}' <<<"$obj")
			fi
		done

		# Append to dataset
		jq --arg s "$rec_key" --argjson o "$obj" '. + {($s): $o}' "$DATASET_FILE" >"${DATASET_FILE}.tmp"
		mv "${DATASET_FILE}.tmp" "$DATASET_FILE"
	done

	echo "[INFO] Successfully imported CSV -> $DATASET_FILE"
}

########################################
# View dataset
########################################
dataset-view() {
	# Usage: dataset-view [fields_csv] [match_regex]
	# If fields_csv is omitted or empty, prints full records.
	local fields="${1:-}"
	local match="${2:-.*}"
	_dataset_load

	local jq_filter
	if [[ -n "$fields" ]]; then
		local jq_fields=""
		IFS=',' read -ra cols <<<"$fields"
		for col in "${cols[@]}"; do
			jq_fields+="\"$col\": .\"$col\","
		done
		jq_fields="{${jq_fields%,}}"
		jq_filter="map({key: .key, value: (${jq_fields})})"
	else
		jq_filter="map({key: .key, value: .value})"
	fi

	jq --arg m "$match" 'to_entries | map(select(.key | test($m))) | '"$jq_filter" "$DATASET_FILE"
}

########################################
# Tabular view (column-aligned / TSV/CSV)
########################################
dataset-view-table() {
	# Usage: dataset-view-table [--fields=a,b,c] [--match=regex] [--na=NA] [--sep=$'\t'] [--no-header]
	local fields_csv="" match=".*" na="NA" sep=$'\t' header=true

	# Safe arg parsing with set -u
	while (($#)) && [[ $1 == --* ]]; do
		case "$1" in
		--fields=*) fields_csv="${1#--fields=}" ;;
		--match=*) match="${1#--match=}" ;;
		--na=*) na="${1#--na=}" ;;
		--sep=*) sep="${1#--sep=}" ;;
		--no-header) header=false ;;
		*)
			echo "[ERROR] Unknown option: $1" >&2
			return 1
			;;
		esac
		shift
	done

	_dataset_load

	# Build fields list (JSON array)
	local fields_json
	if [[ -n "$fields_csv" ]]; then
		fields_json="[\"${fields_csv//,/\",\"}\"]"
	else
		fields_json=$(jq -c '
      to_entries
      | map(select(.value|type=="object") | .value | keys)
      | flatten | unique
    ' "$DATASET_FILE")
	fi

	# Header
	if $header; then
		if [[ -n "$fields_csv" ]]; then
			printf "species%s%s\n" "$sep" "${fields_csv//,/$sep}"
		else
			if [[ "$sep" == $'\t' ]]; then
				jq -r '. as $F | "species\t" + ($F | join("\t"))' <<<"$fields_json"
			else
				printf "species"
				jq -r '.[]' <<<"$fields_json" | while IFS= read -r col; do
					printf "%s%s" "$sep" "$col"
				done
				printf "\n"
			fi
		fi
	fi

	# Rows (keep '.' = dataset; guard per-record)
	local tsv
	tsv=$(jq -r --arg m "$match" --argjson F "$fields_json" --arg na "$na" '
    def safe($x; $na): if $x==null or ($x|tostring|length)==0 then $na else $x end;
    . as $D
    | (keys | map(select(test($m)))) as $S
    | $S[] as $id
    | ($D[$id] as $rec | if ($rec|type)!="object" then {} else $rec end) as $R
    | [$id] + ($F | map(. as $k | safe($R[$k]; $na)))
    | @tsv
  ' "$DATASET_FILE")

	if [[ "$sep" == $'\t' ]]; then
		printf "%s\n" "$tsv"
	else
		printf "%s\n" "$tsv" | sed $'s/\t/'"$sep"'/g'
	fi
}

########################################
# Get a field value for a given species key
########################################
# Get a field value for a given species key (with aliasing + default)
dataset-get-field() {
	# Usage: dataset-get-field SPECIES_KEY FIELD_NAME [DEFAULT]
	local species="$1"
	local field="$2"
	local defval="${3-}" # optional default
	if [[ -z "$species" || -z "$field" ]]; then
		echo "Usage: dataset-get-field SPECIES FIELD [DEFAULT]" >&2
		return 1
	fi

	# Field aliases
	# case "$field" in
	# long) field="long_sra" ;;
	# short) field="short_sra" ;;
	# esac

	_dataset_load
	jq -r --arg s "$species" --arg f "$field" --arg d "$defval" '
    ( .[$s][$f] // $d // "" )
  ' "$DATASET_FILE"
}

# Global composite-key cache: _dataset["<species>|<field>"]=value
declare -gA _dataset

dataset-cache-fields() {
	# Usage: dataset-cache-fields field1 field2 ...
	# Example: dataset-cache-fields long short bioproject
	_dataset=() # reset cache
	_dataset_load

	# require at least one field
	if (($# == 0)); then
		echo "[ERROR] dataset-cache-fields: need at least one field (e.g., long short)" >&2
		return 2
	fi

	# fields exactly as given (no aliasing)
	local -a fields=("$@")

	# all species keys
	local -a _keys
	mapfile -t _keys < <(jq -r 'keys[]' "$DATASET_FILE")

	local k f v
	for k in "${_keys[@]}"; do
		for f in "${fields[@]}"; do
			v="$(jq -r --arg s "$k" --arg f "$f" '.[$s][$f] // ""' "$DATASET_FILE")"
			_dataset["$k|$f"]="$v"
		done
	done
}

########################################
# CLI dispatcher: polap dataset ...
########################################
dataset-cli() {
	local subcommand="$1"
	shift
	case "$subcommand" in
	add) dataset-add-record "$@" ;;
	del | delete) dataset-delete-record "$@" ;;
	search) dataset-search "$@" ;;
	update) dataset-update-field "$@" ;;
	add-field) dataset-add-field "$@" ;;
	view) dataset-view "$@" ;;
	view-grouped) dataset-view-grouped "$@" ;;
	report) dataset-report-species "$@" ;;
	validate) dataset-validate-platforms ;;
	tui) dataset-tui ;;
	import) dataset-import-csv "$@" ;;
	export) dataset-export-csv ;;
	*)
		echo "Usage: polap dataset <subcommand> [...]"
		echo "  subcommands: add, delete, update, add-field, view, view-grouped, report, search, validate, tui, import, export"
		return 1
		;;
	esac
}

########################################
# Generate all Markdown reports as HTML
########################################
dataset-generate-html-reports() {
	_dataset_load
	mkdir -p dataset-reports
	for s in $(jq -r 'keys[]' "$DATASET_FILE"); do
		local mdfile="dataset-reports/${s}.md"
		local htmlfile="dataset-reports/${s}.html"
		dataset-report-species "$s" >"$mdfile"
		pandoc "$mdfile" -o "$htmlfile" --standalone --toc --metadata title="$s"
		echo "[INFO] Generated $htmlfile"
	done
}

########################################
# REDEFINED: Add a new record using a template (corrected grammar)
# Usage:
#  A) dataset-add-template NEW_ID TEMPLATE_ID [field=value ...]
#  B) dataset-add-template --auto-increment BASE_GENUS_SPECIES TEMPLATE_ID [field=value ...]
########################################
dataset-add-template() {
	local new_id=""
	local template_id=""

	_dataset_load

	if [[ "$1" == "--auto-increment" ]]; then
		shift
		local base="$1"
		shift || true
		template_id="$1"
		shift || true
		if [[ -z "$base" || -z "$template_id" ]]; then
			echo "Usage: dataset-add-template --auto-increment BASE_GENUS_SPECIES TEMPLATE_ID [field=value ...]" >&2
			return 1
		fi
		local max_index
		max_index=$(jq -r 'keys[]' "$DATASET_FILE" | grep -E "^${base}-[0-9]+$" | sed 's/.*-//' | sort -nr | head -n1)
		[[ -z "$max_index" ]] && max_index=-1
		local next_index=$((max_index + 1))
		new_id="${base}-${next_index}"
		echo "[INFO] Auto-generated species ID: $new_id"
	else
		new_id="$1"
		shift || true
		template_id="$1"
		shift || true
		if [[ -z "$new_id" || -z "$template_id" ]]; then
			echo "Usage: dataset-add-template NEW_ID TEMPLATE_ID [field=value ...]" >&2
			return 1
		fi
	fi

	if jq -e --arg n "$new_id" 'has($n)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Target $new_id already exists." >&2
		return 1
	fi
	if ! jq -e --arg t "$template_id" 'has($t)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Template $template_id not found." >&2
		return 1
	fi

	local clone
	clone=$(jq --arg t "$template_id" '.[$t]' "$DATASET_FILE")

	local kv k v
	for kv in "$@"; do
		k="${kv%%=*}"
		v="${kv#*=}"
		clone=$(jq --arg k "$k" --arg v "$v" '. + {($k): $v}' <<<"$clone")
	done

	jq --arg n "$new_id" --argjson obj "$clone" '. + {($n): $obj}' "$DATASET_FILE" | _dataset_save
	echo "[INFO] Created $new_id from $template_id"
}

########################################
# REDEFINED: Rename a species record with auto-updates
########################################
dataset-rename-species() {
	local old_id="$1"
	local new_id="$2"
	_dataset_load

	if [[ -z "$old_id" || -z "$new_id" ]]; then
		echo "Usage: dataset-rename-species OLD_ID NEW_ID" >&2
		return 1
	fi
	if ! jq -e --arg o "$old_id" 'has($o)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Species $old_id not found." >&2
		return 1
	fi
	if jq -e --arg n "$new_id" 'has($n)' "$DATASET_FILE" >/dev/null; then
		echo "[ERROR] Target $new_id already exists." >&2
		return 1
	fi

	local base species2
	base="${new_id%-*}"
	species2="${base//_/ }"

	jq --arg o "$old_id" --arg n "$new_id" --arg base "$base" --arg species2 "$species2" '
    . + {($n): .[$o]} | del(.[$o]) |
    .[$n].species = $n |
    .[$n].folder = $base |
    .[$n].species2 = $species2
  ' "$DATASET_FILE" | _dataset_save

	echo "[INFO] Renamed $old_id -> $new_id (updated fields: species, folder, species2)"
}
