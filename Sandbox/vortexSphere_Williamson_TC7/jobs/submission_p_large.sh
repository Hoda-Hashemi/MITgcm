#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
MITGCM_DIR=$(cd "$CASE_DIR/../.." && pwd)
PYTHON="$MITGCM_DIR/.venv/bin/python"
VALIDATE_SCRIPT="$MITGCM_DIR/Sandbox/Scripts/tc7_validate_initial_conditions.py"
RAW_DIR="$CASE_DIR/input/raw"

if [ ! -x "$PYTHON" ]; then
    PYTHON=python3
fi

missing=0
for label in 19781221_0000 19790116_0000 19790109_0000; do
    raw_file="$RAW_DIR/tc7_${label}_initial_conditions.npz"
    if [ ! -f "$raw_file" ]; then
        echo "Missing $raw_file" >&2
        missing=1
    fi
done

if [ "$missing" -ne 0 ]; then
    cat >&2 <<EOF
TC7 is not submit-ready: missing one or more prepared case NPZ files.

Prepare them with:
  $CASE_DIR/tools/prepare_tc7_from_ncep.py --case all
EOF
    exit 2
fi

for label in 19781221_0000 19790116_0000 19790109_0000; do
    "$PYTHON" "$VALIDATE_SCRIPT" "$RAW_DIR/tc7_${label}_initial_conditions.npz"
done

cd "$SCRIPT_DIR/large" || exit 1

sbatch job_tc7_c1.slurm
sbatch job_tc7_c2.slurm
sbatch job_tc7_c3.slurm
