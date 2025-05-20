#!/bin/bash
set -euo pipefail
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

short_read1="$1"
short_read2="$2"
output_dest="$3"
outdir="$4"

script_name="$(basename "$0" .sh)"
logfile="${outdir}/log/${script_name}.log"

echo "[INFO] Polishing pipeline started at $(date)" >"$logfile"
start_time=$(date +%s)

source "$(conda info --base)/etc/profile.d/conda.sh"

# Function to verify dependencies
check_command() {
  if ! command -v "$1" &>/dev/null; then
    echo "[ERROR] Required command '$1' not found in PATH." >>"$logfile"
    exit 1
  fi
}

{
  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap-fmlrc" ]]; then
    echo "[INFO] Activating conda environment 'polap-fmlrc'..."
    conda activate polap-fmlrc
  fi

  if [[ "$CONDA_DEFAULT_ENV" != "polap-fmlrc" ]]; then
    echo "[ERROR] Failed to enter conda environment 'polap-fmlrc'"
    exit 1
  fi

  echo "[INFO] Running polishing pipeline in environment: $CONDA_DEFAULT_ENV"

  # Verify required tools
  for cmd in ropebwt2 msbwt awk sort tr zcat; do
    check_command "$cmd"
  done

  echo "[INFO] Removing previous msbwt directory: ${outdir}/msbwt"
  rm -rf "${outdir}/msbwt"

  echo "[INFO] Creating ${outdir}/msbwt ..."
  sleep 5

  if [[ "$short_read1" = *.fastq || "$short_read1" = *.fq ]]; then
    cat "$short_read1" "${short_read2:-/dev/null}" |
      awk 'NR % 4 == 2' | sort | tr NT TN |
      ropebwt2 -LR 2>"$output_dest" |
      tr NT TN |
      msbwt convert "${outdir}/msbwt" >/dev/null 2>&1

  elif [[ "$short_read1" = *.fq.gz || "$short_read1" = *.fastq.gz ]]; then
    zcat "$short_read1" "${short_read2:-/dev/null}" |
      awk 'NR % 4 == 2' | sort | tr NT TN |
      ropebwt2 -LR 2>"$output_dest" |
      tr NT TN |
      msbwt convert "${outdir}/msbwt" >/dev/null 2>&1

  else
    echo "[ERROR] Unrecognized input file format:"
    echo "[ERROR] short_read1: $short_read1"
    echo "[ERROR] short_read2: $short_read2"
    exit 1
  fi

  pipeline_status=$?
  conda deactivate
} >>"$logfile" 2>&1

end_time=$(date +%s)
elapsed=$((end_time - start_time))

if [[ "${pipeline_status:-1}" -eq 0 ]]; then
  printf "[INFO] Polishing pipeline finished at %s (elapsed: %d sec)\n" "$(date)" "$elapsed" >>"$logfile"
else
  printf "[ERROR] Pipeline failed at %s (after %d sec, exit code %d)\n" "$(date)" "$elapsed" "${pipeline_status:-1}" >>"$logfile"
  exit "${pipeline_status:-1}"
fi
