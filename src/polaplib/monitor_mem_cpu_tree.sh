#!/bin/bash

TARGET_PID=$1
LOGFILE=$2
INTERVAL=${3:-1} # seconds

echo "timestamp,total_rss_kb,total_cpu_percent" >"$LOGFILE"

get_descendants() {
	local pid=$1
	echo $pid
	for child in $(ps --no-headers -o pid --ppid "$pid"); do
		get_descendants "$child"
	done
}

while kill -0 "$TARGET_PID" 2>/dev/null; do
	pids=$(get_descendants "$TARGET_PID")
	pid_list=$(echo $pids | tr ' ' ',')

	total_rss=$(ps --no-headers -o rss -p "$pid_list" | awk '{sum+=$1} END{print sum}')
	total_cpu=$(ps --no-headers -o %cpu -p "$pid_list" | awk '{sum+=$1} END{print sum}')
	timestamp=$(date +%s)

	echo "$timestamp,$total_rss,$total_cpu" >>"$LOGFILE"
	sleep "$INTERVAL"
done
