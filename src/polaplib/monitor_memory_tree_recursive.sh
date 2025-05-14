#!/bin/bash

TARGET_PID="$1"
LOGFILE=${2:-memlog_pss_nosudo.csv}
INTERVAL=${3:-60}

echo "timestamp,total_pss_kb,tree_cpu_percent,cpu_load_1min,disk_free_gb,max_pss_cmd" >"$LOGFILE"

get_descendants() {
  local pid=$1
  echo "$pid"
  for child in $(pgrep -P "$pid"); do
    get_descendants "$child"
  done
}

while kill -0 "$TARGET_PID" 2>/dev/null; do
  all_pids=$(get_descendants "$TARGET_PID")

  total_pss=0
  max_pss=0
  max_pid=""
  for pid in $all_pids; do
    if [[ -r "/proc/$pid/smaps_rollup" ]]; then
      pss=$(grep ^Pss: /proc/$pid/smaps_rollup 2>/dev/null | awk '{print $2}')
      total_pss=$((total_pss + pss))
      if [[ "$pss" -gt "$max_pss" ]]; then
        max_pss=$pss
        max_pid=$pid
      fi
    fi
  done

  # Get process name of max PSS PID
  if [[ -n "$max_pid" && -r "/proc/$max_pid/cmdline" ]]; then
    max_cmd=$(tr '\0' ' ' </proc/$max_pid/cmdline | cut -c1-100)
  else
    max_cmd="N/A"
  fi

  total_cpu=$(ps -o %cpu= -p $(echo "$all_pids" | tr '\n' ',' | sed 's/,$//') | awk '{sum+=$1} END{print sum}')
  load_avg_1min=$(awk '{print $1}' /proc/loadavg)
  disk_free_gb=$(df -BG . | awk 'NR==2 {gsub("G","",$4); print $4}')
  timestamp=$(date +%s)

  printf "%s,%d,%.2f,%.2f,%s,\"%s\"\n" "$timestamp" "$total_pss" "$total_cpu" "$load_avg_1min" "$disk_free_gb" "$max_cmd" >>"$LOGFILE"

  sleep "$INTERVAL"
done
