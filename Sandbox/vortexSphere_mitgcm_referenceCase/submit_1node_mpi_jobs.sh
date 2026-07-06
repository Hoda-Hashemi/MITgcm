#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
JOBS_DIR="$SCRIPT_DIR/jobs"
LOG_DIR="$SCRIPT_DIR/logs"

mkdir -p "$LOG_DIR"

mapfile -t JOBS < <(find "$JOBS_DIR" -maxdepth 1 -type f -name "job_*mpi_1node.slurm" | sort -V)

if [ "${#JOBS[@]}" -eq 0 ]; then
    echo "No one-node MPI Slurm jobs found in $JOBS_DIR" >&2
    exit 1
fi

echo "Submitting ${#JOBS[@]} one-node MPI jobs"
echo "Logs will be written under: $LOG_DIR"

for job in "${JOBS[@]}"; do
    echo "Submitting $(basename "$job")"
    sbatch --chdir="$LOG_DIR" "$job"
done
