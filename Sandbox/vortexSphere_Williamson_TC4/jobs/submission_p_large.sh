#!/bin/bash

set -euo pipefail

CASE_DIR=$(cd "$(dirname "$0")/.." && pwd)
BUILD_DIR="$CASE_DIR/build/c1_large48"

for required in \
    "$CASE_DIR/code/tc4_forcing.F" \
    "$CASE_DIR/code/external_forcing_surf.F" \
    "$CASE_DIR/input/gendata_ref.py"
do
    if [ ! -s "$required" ]; then
        echo "Missing required TC4 file: $required" >&2
        exit 2
    fi
done

if [ ! -x "$BUILD_DIR/mitgcmuv" ]; then
    echo "TC4 executable not built yet: $BUILD_DIR/mitgcmuv" >&2
    echo "Run the job once or build with jobs/large/job_tc4_c1.slurm's genmake settings." >&2
    exit 2
fi

cd large || exit 1

sbatch job_tc4_c1.slurm
