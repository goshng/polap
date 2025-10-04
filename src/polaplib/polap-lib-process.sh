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
# Bash functions for monitoring system usage such as memory, CPU, storage
# TEST-SCC: not yet
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

_polap_lib_process-run_script_with_args_in_bash_c() {
	local script="$1"
	shift

	# Quote each argument safely
	local quoted_args=()
	for arg in "$@"; do
		quoted_args+=("\"$arg\"")
	done

	# Build the final command string
	local command="bash \"$script\" ${quoted_args[*]}"

	# Run the command in a new bash -c subshell and return the PID
	bash -c "$command"
}

_polap_lib_process-start_memtracker() {
	local log_file="$1"
	local interval="${2:-60}"

	# This function runs only if opt_f_flag is defined to be other than false.
	if [[ "${opt_f_flag:-false}" == "false" ]]; then
		return
	fi

	# input: "/a/b/c/memlog-msbwt-node.csv"
	# output: msbwt-node
	local run_title="${log_file##*/}"
	run_title="${run_title#*-}"
	run_title="${run_title%.csv}"

	echo "timestamp,total_used_kb,cpu_load_1min,disk_free_gb,max_rss_cmd" >"$log_file"
	# echo "timestamp,total_used_kb,cpu_load_1min,disk_free_gb" >"$log_file"
	# local start_ts
	# start_ts=$(date +%s)
	# echo "$start_ts" >"${log_file}.start"

	(
		while true; do
			local ts
			ts=$(date +%s)

			# Memory
			local total_kb available_kb used_kb
			total_kb=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
			# echo "${total_kb}"
			available_kb=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
			# echo "${available_kb}"
			used_kb=$((total_kb - available_kb))
			# echo "${used_kb}"

			# CPU load (1 min avg)
			local cpu_load
			cpu_load=$(awk '{print $1}' /proc/loadavg)
			# echo "${cpu_load}"

			# Disk space (GB, root mount)
			local disk_free_gb
			disk_free_gb=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
			# echo "${disk_free_gb}"

			# Most memory-hungry command (RSS)
			# local max_cmd
			# max_cmd=$(ps -eo rss,cmd --no-headers | sort -nr | head -n1 | sed 's/^[0-9]* //')
			# echo "${max_cmd}"
			# if ps -eo rss,cmd --no-headers >/dev/null 2>&1; then
			# 	max_cmd=$(ps -eo rss,cmd --no-headers | sort -nr | head -n1 | sed 's/^[0-9]* //')
			# else
			# 	max_cmd="NA"
			# fi

			# Requires: bash with set -euo pipefail
			# Goal: set max_cmd to the command line of the largest-RSS process, or "NA".

			# Safest defaults
			max_cmd="NA"

			# Check once whether ps supports the query
			if ps -eo rss,cmd --no-headers >/dev/null 2>&1; then
				# Capture ps output without letting failures trip -e
				ps_out="$(ps -eo rss,cmd --no-headers 2>/dev/null || true)"

				if [[ -n "$ps_out" ]]; then
					# Pick the largest-RSS line; guard each stage against -e
					top_line="$(printf '%s\n' "$ps_out" | LC_ALL=C sort -nr 2>/dev/null | head -n1 || true)"

					if [[ -n "$top_line" ]]; then
						# Strip the leading RSS field (+spaces) → leave only the command
						# no-awk: use sed with character classes
						max_cmd="$(sed -e 's/^[[:space:]]*[0-9][0-9]*[[:space:]]*//' <<<"$top_line" || true)"
						[[ -z "$max_cmd" ]] && max_cmd="NA"
					fi
				fi
			else
				max_cmd="NA"
			fi

			echo "$ts,$used_kb,$cpu_load,$disk_free_gb,\"$max_cmd\""
			# echo "$ts,$used_kb,$cpu_load,$disk_free_gb"

			sleep "$interval"
		done
	) >>"$log_file" &
	echo $! >"${log_file}.pid"
	if [[ "${_POLAP_DEBUG}" == "1" ]]; then
		echo "[INFO] Job: [$run_title] started with PID $(<"${log_file}.pid")"
	fi
}

_polap_lib_process-end_memtracker() {
	local log_file="$1"
	local summary_base="${2:-memtrack-summary.txt}"
	local _arg_verbose="${3:-no_verbose}"

	# This function runs only if opt_f_flag is defined to be other than false.
	if [[ "${opt_f_flag:-false}" == "false" ]]; then
		return
	fi

	local timestamp
	timestamp=$(date +"%Y%m%d_%H%M")
	local summary_file="${summary_base}.${timestamp}.txt"

	if [[ ! -f "$log_file" ]]; then
		echo "[ERROR] Log file not found: $log_file" >&2
		return 1
	fi

	# Stop logger
	if [[ -f "${log_file}.pid" ]]; then
		kill "$(cat "${log_file}.pid")" 2>/dev/null || true
		rm -f "${log_file}.pid"
	fi

	local line_count
	line_count=$(wc -l <"$log_file" | xargs)
	if [[ "$line_count" -le 1 ]]; then
		echo "[ERROR] Log file has no log: $log_file" >&2
		return 1
	fi

	local start_ts end_ts elapsed
	start_ts=$(awk -F',' 'NR==2 {print $1; exit}' "$log_file")
	end_ts=$(awk -F',' 'END {print $1}' "$log_file")
	elapsed=$((end_ts - start_ts))

	awk -F',' -v start_ts="$start_ts" -v end_ts="$end_ts" -v elapsed="$elapsed" \
		-v summary_file="$summary_file" '
    BEGIN {
      min_disk = 999999
      peak_cpu = 0
      peak_mem = 0
      start_disk_free = -1
    }
    NR == 1 { next }
    NR == 2 {
      start_used = $2
      start_disk_free = $4
    }
    {
      if ($2 > peak_mem) {
        peak_mem = $2
        peak_cmd = $5
      }
      if ($3 > peak_cpu) peak_cpu = $3
      if ($4 < min_disk) min_disk = $4
    }
    END {
      delta_mem = peak_mem - start_used
      delta_disk = start_disk_free - min_disk
      h = int(elapsed / 3600)
      m = int((elapsed % 3600) / 60)
      s = elapsed % 60
      elapsed_str = sprintf("%02d:%02d:%02d (%.2f h)", h, m, s, elapsed / 3600)

      print "Physical Memory Usage Summary:"                       > summary_file
      printf "Start used:      %d KB (%.2f GB)\n", start_used, start_used/1048576  >> summary_file
      printf "Peak used:       %d KB (%.2f GB)\n", peak_mem, peak_mem/1048576      >> summary_file
      printf "Net increase:    %d KB (%.2f GB)\n", delta_mem, delta_mem/1048576    >> summary_file
      printf "Elapsed time:    %s\n", elapsed_str                                  >> summary_file
      printf "Peak CPU load:   %.2f\n", peak_cpu                                   >> summary_file
      printf "Start disk free: %d GB\n", start_disk_free                           >> summary_file
      printf "Min disk free:   %d GB\n", min_disk                                   >> summary_file
      printf "Disk used:       %d GB\n", delta_disk                                 >> summary_file
      printf "Top memory cmd:  %s\n", peak_cmd                                     >> summary_file
    }
  ' "$log_file"

	# Show summary to terminal
	if [[ "${_arg_verbose}" == "verbose" ]]; then
		cat "$summary_file"
	fi

	# Update symlink to point to latest
	ln -sf "$(basename "$summary_file")" "$summary_base"
}

_polap_lib_process-analyze_memtracker_log() {
	local log_file="$1"
	local summary_base="${2:-memtrack-summary.txt}"
	local timestamp
	timestamp=$(date +"%Y%m%d_%H%M")
	local summary_file="${summary_base}.${timestamp}.txt"

	if [[ ! -f "$log_file" ]]; then
		echo "[ERROR] Log file not found: $log_file" >&2
		return 1
	fi

	local line_count
	line_count=$(wc -l <"$log_file" | xargs)
	if [[ "$line_count" -le 1 ]]; then
		echo "[ERROR] Log file has no log: $log_file" >&2
		return 1
	fi

	local start_ts end_ts elapsed
	start_ts=$(awk -F',' 'NR==2 {print $1; exit}' "$log_file")
	# end_ts=$(awk -F',' 'END {print $1}' "$log_file")
	end_ts=$(tail -n2 "$log_file" | head -n1 | cut -d',' -f1)
	elapsed=$((end_ts - start_ts))

	awk -F',' -v start_ts="$start_ts" -v end_ts="$end_ts" -v elapsed="$elapsed" \
		-v summary_file="$summary_file" '
    BEGIN {
      min_disk = 999999
      peak_cpu = 0
      peak_mem = 0
      start_disk_free = -1
    }
    NR == 1 { next }
    NR == 2 {
      start_used = $2
      start_disk_free = $4
    }
    {
      if ($2 > peak_mem) {
        peak_mem = $2
        peak_cmd = $5
      }
      if ($3 > peak_cpu) peak_cpu = $3
      if ($4 < min_disk) min_disk = $4
    }
    END {
      delta_mem = peak_mem - start_used
      delta_disk = start_disk_free - min_disk
      h = int(elapsed / 3600)
      m = int((elapsed % 3600) / 60)
      s = elapsed % 60
      elapsed_str = sprintf("%02d:%02d:%02d (%.2f h)", h, m, s, elapsed / 3600)

      print "Physical Memory Usage Summary:"                       > summary_file
      printf "Start used:      %d KB (%.2f GB)\n", start_used, start_used/1048576  >> summary_file
      printf "Peak used:       %d KB (%.2f GB)\n", peak_mem, peak_mem/1048576      >> summary_file
      printf "Net increase:    %d KB (%.2f GB)\n", delta_mem, delta_mem/1048576    >> summary_file
      printf "Elapsed time:    %s\n", elapsed_str                                  >> summary_file
      printf "Peak CPU load:   %.2f\n", peak_cpu                                   >> summary_file
      printf "Start disk free: %d GB\n", start_disk_free                           >> summary_file
      printf "Min disk free:   %d GB\n", min_disk                                   >> summary_file
      printf "Disk used:       %d GB\n", delta_disk                                 >> summary_file
      printf "Top memory cmd:  %s\n", peak_cmd                                     >> summary_file
    }
  ' "$log_file"

	cat "$summary_file"
	ln -sf "$(basename "$summary_file")" "$summary_base"
}
