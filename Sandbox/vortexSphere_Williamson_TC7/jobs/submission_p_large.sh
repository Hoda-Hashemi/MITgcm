#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
MITGCM_DIR=$(cd "$CASE_DIR/../.." && pwd)
PYTHON="$MITGCM_DIR/.venv/bin/python"
RAW_FILE="$CASE_DIR/input/raw/tc7_initial_conditions.npz"
VALIDATE_SCRIPT="$MITGCM_DIR/Sandbox/Scripts/tc7_validate_initial_conditions.py"

if [ ! -x "$PYTHON" ]; then
    PYTHON=python3
fi

if [ ! -f "$RAW_FILE" ]; then
    cat >&2 <<EOF
TC7 is not submit-ready: missing $RAW_FILE

Download the 500 hPa analysis file, convert it to tc7_initial_conditions.npz,
and rerun this wrapper. See:
  $CASE_DIR/input/raw/README_TC7.md
EOF
    exit 2
fi

"$PYTHON" "$VALIDATE_SCRIPT" "$RAW_FILE"

cd "$SCRIPT_DIR/large" || exit 1

sbatch job_tc7_c1.slurm
