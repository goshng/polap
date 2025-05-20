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

	echo "timestamp,total_used_kb,cpu_load_1min,disk_free_gb,max_rss_cmd" >"$log_file"
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
			available_kb=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
			used_kb=$((total_kb - available_kb))

			# CPU load (1 min avg)
			local cpu_load
			cpu_load=$(awk '{print $1}' /proc/loadavg)

			# Disk space (GB, root mount)
			local disk_free_gb
			disk_free_gb=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')

			# Most memory-hungry command (RSS)
			local max_cmd
			max_cmd=$(ps -eo rss,cmd --no-headers | sort -nr | head -n1 | sed 's/^[0-9]* //')

			echo "$ts,$used_kb,$cpu_load,$disk_free_gb,\"$max_cmd\"" >>"$log_file"

			sleep "$interval"
		done
	) &
	echo $! >"${log_file}.pid"
	echo "[INFO] Logging started with PID $(<"${log_file}.pid")"
}

_polap_lib_process-end_memtracker() {
	local log_file="$1"
	local summary_base="${2:-memtrack-summary.txt}"
	local timestamp
	timestamp=$(date +"%Y%m%d_%H%M")
	local summary_file="${summary_base}.${timestamp}.txt"
	local md_summary_file="${summary_base}.${timestamp}.md"

	if [[ ! -f "$log_file" ]]; then
		echo "[ERROR] Log file not found: $log_file" >&2
		return 1
	fi

	# Stop logger
	if [[ -f "${log_file}.pid" ]]; then
		kill "$(cat "${log_file}.pid")" 2>/dev/null || true
		rm -f "${log_file}.pid"
	fi

	local line_count=$(wc -l <"$log_file" | xargs)
	if [[ "$line_count" -le 1 ]]; then
		echo "[ERROR] Log file has no log: $log_file" >&2
		return 1
	fi

	local start_ts end_ts elapsed
	start_ts=$(awk -F',' 'NR==2 {print $1; exit}' "$log_file")
	end_ts=$(awk -F',' 'END {print $1}' "$log_file")
	elapsed=$((end_ts - start_ts))

	awk -F',' -v start_ts="$start_ts" -v end_ts="$end_ts" -v elapsed="$elapsed" \
		-v summary_file="$summary_file" -v md_file="$md_summary_file" '
    BEGIN {
      min_disk = 999999
      peak_cpu = 0
      peak_mem = 0
    }
    NR == 1 { next }
    NR == 2 { start_used = $2 }
    {
      if ($2 > peak_mem) {
        peak_mem = $2
        peak_cmd = $5
      }
      if ($3 > peak_cpu) peak_cpu = $3
      if ($4 < min_disk) min_disk = $4
    }
    END {
      delta = peak_mem - start_used
      h = int(elapsed / 3600)
      m = int((elapsed % 3600) / 60)
      s = elapsed % 60
      elapsed_str = sprintf("%02d:%02d:%02d (%.2f h)", h, m, s, elapsed / 3600)

      # Plain text
      printf "🧠 Physical Memory Usage Summary:\n"                       > summary_file
      printf "Start used:      %d KB (%.2f GB)\n", start_used, start_used/1048576  >> summary_file
      printf "Peak used:       %d KB (%.2f GB)\n", peak_mem, peak_mem/1048576      >> summary_file
      printf "Net increase:    %d KB (%.2f GB)\n", delta, delta/1048576            >> summary_file
      printf "Elapsed time:    %s\n", elapsed_str                                  >> summary_file
      printf "Peak CPU load:   %.2f\n", peak_cpu                                   >> summary_file
      printf "Min disk space:  %d GB\n", min_disk                                  >> summary_file
      printf "Top memory cmd:  %s\n", peak_cmd                                     >> summary_file

      # Markdown
      print "### 🧠 Physical Memory Usage Summary" > md_file
      print "| Metric           | Value                         |" >> md_file
      print "|------------------|-------------------------------|" >> md_file
      printf "| Start used       | %d KB (%.2f GB)              |\n", start_used, start_used/1048576 >> md_file
      printf "| Peak used        | %d KB (%.2f GB)              |\n", peak_mem, peak_mem/1048576 >> md_file
      printf "| Net increase     | %d KB (%.2f GB)              |\n", delta, delta/1048576 >> md_file
      printf "| Elapsed time     | %s                          |\n", elapsed_str >> md_file
      printf "| Peak CPU load    | %.2f                          |\n", peak_cpu >> md_file
      printf "| Min disk space   | %d GB                         |\n", min_disk >> md_file
      printf "| Top memory cmd   | `%s`                         |\n", peak_cmd >> md_file
    }
  ' "$log_file"

	# Show plain-text summary to terminal
	cat "$summary_file"

	# Update symlinks
	ln -sf "$(basename "$summary_file")" "$summary_base"
	ln -sf "$(basename "$md_summary_file")" "${summary_base%.txt}.md"

	echo "[INFO] Summary written to: $summary_file"
	echo "[INFO] Markdown written to: $md_summary_file"
	echo "[INFO] Symlinks updated: $summary_base, ${summary_base%.txt}.md"
}
