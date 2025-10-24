#!/usr/bin/env bash
set -euo pipefail

# polap-bash-config-new.sh
# Create a new flat YAML config by enforcing ordered leading flags, then
# forwarding remaining options to polap-py-config-new.py.
#
# Ordered flags (prefix of argv):
#   --preset NAME                 (REQUIRED)
#   [--template TEMP]             (OPTIONAL; name or path)
#   [--config-dir DIR]            (OPTIONAL; default: ~/.polap/profiles)
#
# Remaining options are interpreted as:
#   --wgs-mode | --no-wgs-mode | --no-wgs_mode   # forwarded
#   --hpc | --no-hpc                              # forwarded
#   --key value                                  # saved as top-level via --kv key=value
#   --bool key true|false                        # saved as boolean
#
# Examples:
#   source polap-bash-config-new.sh
#   polap_config_new --preset mt_hifi --config-dir ~/.polap/profiles \
#     --wgs-mode \
#     --reads /data/all.fq.gz \
#     --anchors /data/anchors/mt.id.all.txt \
#     --hmm-db /data/hmm/mt_genes.hmm \
#     --outdir /work/autotune \
#     --threads 16 \
#     --k 121 --s 27 --no-hpc \
#     --final-if "genes_score < 0.95 || breadth < 0.97" \
#     --c-radius "0,10,20" \
#     --min-shared 5 --jaccard-min 0.015 --topk-nei 50 --steps 2

polap_config_new() {
	# ordered header
	local PRESET=""
	local TEMPLATE=""
	local CONFIG_DIR=""

	# tool
	local PY="${_POLAPLIB_DIR:-.}/polap-py-config-new.py"

	# forwarding accumulators
	local FWD_ARGS=()
	local KV_ARGS=()
	local BOOL_ARGS=()

	# parse ordered leading flags
	local ARGS=("$@")
	local idx=0

	if [[ ${#ARGS[@]} -lt 2 || "${ARGS[0]}" != "--preset" ]]; then
		echo "[error] expected ordered prefix: --preset NAME [--template TEMP] [--config-dir DIR]" >&2
		return 2
	fi
	PRESET="${ARGS[1]}"
	idx=2

	if [[ $idx -lt ${#ARGS[@]} && "${ARGS[$idx]}" == "--template" ]]; then
		TEMPLATE="${ARGS[$((idx + 1))]}"
		idx=$((idx + 2))
	fi
	if [[ $idx -lt ${#ARGS[@]} && "${ARGS[$idx]}" == "--config-dir" ]]; then
		CONFIG_DIR="${ARGS[$((idx + 1))]}"
		idx=$((idx + 2))
	fi
	[[ -z "$CONFIG_DIR" ]] && CONFIG_DIR="${HOME}/.polap/profiles"

	# process the rest
	local i
	i=$idx
	while [[ $i -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$i]}"

		case "$tok" in
		--wgs-mode | --no-wgs-mode | --no-wgs_mode | --hpc | --no-hpc)
			FWD_ARGS+=("$tok")
			i=$((i + 1))
			;;

		--bool)
			# --bool key true|false
			if [[ $((i + 2)) -ge ${#ARGS[@]} ]]; then
				echo "[error] --bool expects: --bool key true|false" >&2
				return 2
			fi
			local bkey="${ARGS[$((i + 1))]}"
			local bval="${ARGS[$((i + 2))]}"
			BOOL_ARGS+=("--bool" "$bkey" "$bval")
			i=$((i + 3))
			;;

		--*)
			# treat as key value -> --kv key=value
			local key="${tok#--}"
			# normalize hyphens to underscores for top-level keys
			key="${key//-/_}"
			if [[ $((i + 1)) -ge ${#ARGS[@]} || "${ARGS[$((i + 1))]}" == --* ]]; then
				echo "[error] missing value for $tok" >&2
				return 2
			fi
			local val="${ARGS[$((i + 1))]}"
			KV_ARGS+=("--kv" "${key}=${val}")
			i=$((i + 2))
			;;

		*)
			echo "[error] unknown token: $tok" >&2
			return 2
			;;
		esac
	done

	# sanity
	[[ -f "$PY" || -x "$PY" ]] || {
		echo "[error] not found: $PY" >&2
		return 127
	}

	# build call
	local CALL=(python3 "$PY" --config-dir "$CONFIG_DIR" --preset "$PRESET")
	[[ -n "$TEMPLATE" ]] && CALL+=(--template "$TEMPLATE")
	CALL+=("${FWD_ARGS[@]}" "${KV_ARGS[@]}" "${BOOL_ARGS[@]}")

	# run
	"${CALL[@]}"
}

# Edit (overlay) a flat YAML config, resolved by --preset + --config-dir,
# or explicitly by --path. Supports in-place (default), --backup, --output,
# and --dry-run. Writes only top-level keys (flat YAML).
#
# Examples:
#   polap_config_edit \
#     --preset mt_hifi --config-dir ~/.polap/profiles \
#     --set threads 32 \
#     --set outdir /scratch/polap \
#     --bool hpc false \
#     --backup
#
#   polap_config_edit \
#     --path ~/.polap/profiles/mt_ont.yaml \
#     --wgs-mode \
#     --hpc \
#     --set k 41 --set s 21 \
#     --output /tmp/mt_ont.edited.yaml --dry-run

polap_config_edit() {
	# -------- ordered addressing & options --------
	local PRESET=""
	local CONFIG_DIR=""
	local PATH_YAML=""
	local OUTPUT_YAML=""
	local BACKUP=0
	local DRY_RUN=0

	# engine
	local PY_EDIT="${_POLAPLIB_DIR:-.}/polap-py-config-new.py"

	# accumulators
	local FWD_ARGS=()  # boolean flags forwarded as-is to python: --wgs-mode/--no-*, --hpc/--no-*
	local KV_ARGS=()   # --kv key=value
	local BOOL_ARGS=() # --bool key true|false

	# -------- parse argv --------
	local ARGS=("$@")
	local i=0
	while [[ $i -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$i]}"
		case "$tok" in
		# addressing
		--preset)
			PRESET="${ARGS[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--config-dir)
			CONFIG_DIR="${ARGS[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--path)
			PATH_YAML="${ARGS[$((i + 1))]:-}"
			i=$((i + 2))
			;;

		# destination & safety
		--output)
			OUTPUT_YAML="${ARGS[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--backup)
			BACKUP=1
			i=$((i + 1))
			;;
		--dry-run)
			DRY_RUN=1
			i=$((i + 1))
			;;

		# boolean flags (forwarded)
		--wgs-mode | --no-wgs-mode | --no-wgs_mode | --hpc | --no-hpc)
			FWD_ARGS+=("$tok")
			i=$((i + 1))
			;;

		# explicit boolean k/v
		--bool)
			if [[ -z "${ARGS[$((i + 1))]:-}" || -z "${ARGS[$((i + 2))]:-}" ]]; then
				echo "[error] --bool expects: --bool key true|false" >&2
				return 2
			fi
			local bkey="${ARGS[$((i + 1))]}"
			local bval="${ARGS[$((i + 2))]}"
			BOOL_ARGS+=("--bool" "$bkey" "$bval")
			i=$((i + 3))
			;;

		# generic key value -> --kv key=value
		--*)
			local key="${tok#--}"
			key="${key//-/_}"
			if [[ -z "${ARGS[$((i + 1))]:-}" || "${ARGS[$((i + 1))]}" == --* ]]; then
				echo "[error] missing value for $tok" >&2
				return 2
			fi
			local val="${ARGS[$((i + 1))]}"
			KV_ARGS+=("--kv" "${key}=${val}")
			i=$((i + 2))
			;;
		-h | --help)
			cat <<'EOF'
Usage:
  polap_config_edit [--preset NAME --config-dir DIR | --path FILE.yaml] \
    [--output FILE.yaml] [--backup] [--dry-run] \
    [--wgs-mode|--no-wgs-mode|--no-wgs_mode] [--hpc|--no-hpc] \
    [--bool key true|false] [--set key value]...

Examples:
  polap_config_edit --preset mt_hifi --config-dir ~/.polap/profiles \
    --set threads 32 --set outdir /scratch/polap --bool hpc false --backup

  polap_config_edit --path ~/.polap/profiles/mt_ont.yaml \
    --wgs-mode --hpc --set k 41 --set s 21 --output /tmp/mt_ont.edited.yaml
EOF
			return 0
			;;
		*)
			echo "[error] unknown token: $tok" >&2
			return 2
			;;
		esac
	done

	# -------- resolve target path --------
	if [[ -z "$PATH_YAML" ]]; then
		[[ -n "$PRESET" ]] || {
			echo "[error] require --preset or --path" >&2
			return 2
		}
		[[ -n "$CONFIG_DIR" ]] || CONFIG_DIR="${HOME}/.polap/profiles"
		PATH_YAML="${CONFIG_DIR%/}/${PRESET}.yaml"
	fi
	local DEST
	if [[ -n "$OUTPUT_YAML" ]]; then
		DEST="$OUTPUT_YAML"
	else
		DEST="$PATH_YAML"
	fi

	# -------- sanity --------
	[[ -f "$PY_EDIT" || -x "$PY_EDIT" ]] || {
		echo "[error] not found: $PY_EDIT" >&2
		return 127
	}

	# -------- dry-run / backup / call plan --------
	if [[ $DRY_RUN -eq 1 ]]; then
		echo "[dry-run] template: $PATH_YAML" >&2
		echo "[dry-run] dest    : $DEST" >&2
		if [[ -z "$OUTPUT_YAML" ]]; then
			[[ $BACKUP -eq 1 ]] && echo "[dry-run] would backup -> ${PATH_YAML}.bak" >&2
			echo "[dry-run] would edit in-place" >&2
		else
			echo "[dry-run] would write new file (no in-place edit)" >&2
		fi
		echo "[dry-run] forwarded flags: ${FWD_ARGS[*]:-(none)}" >&2
		echo "[dry-run] bool args      : ${BOOL_ARGS[*]:-(none)}" >&2
		echo "[dry-run] kv args        : ${KV_ARGS[*]:-(none)}" >&2
		return 0
	fi

	# backup if in-place and file exists
	if [[ -z "$OUTPUT_YAML" && -s "$PATH_YAML" && $BACKUP -eq 1 ]]; then
		local bak="${PATH_YAML}.bak"
		if [[ -e "$bak" ]]; then
			local n=1
			while [[ -e "${bak}.${n}" ]]; do n=$((n + 1)); done
			bak="${bak}.${n}"
		fi
		cp -p "$PATH_YAML" "$bak"
		echo "[backup] saved original -> $bak" >&2
	fi

	# ensure dest parent directory
	mkdir -p "$(dirname -- "$DEST")"

	# call python writer using current YAML as template
	if [[ -s "$PATH_YAML" ]]; then
		python3 "$PY_EDIT" \
			--config-dir "$(dirname -- "$DEST")" \
			--preset "$(basename -- "$DEST" .yaml)" \
			--path "$DEST" \
			--template "$PATH_YAML" \
			"${FWD_ARGS[@]}" "${KV_ARGS[@]}" "${BOOL_ARGS[@]}"
	else
		# no template available; write overlay directly
		python3 "$PY_EDIT" \
			--config-dir "$(dirname -- "$DEST")" \
			--preset "$(basename -- "$DEST" .yaml)" \
			--path "$DEST" \
			"${FWD_ARGS[@]}" "${KV_ARGS[@]}" "${BOOL_ARGS[@]}"
	fi
}
