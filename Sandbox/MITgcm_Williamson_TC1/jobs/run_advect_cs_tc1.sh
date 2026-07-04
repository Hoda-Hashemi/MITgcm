#!/bin/bash

set -euo pipefail

: "${MITGCM_DIR:?MITGCM_DIR is required}"
: "${SETUP_DIR:?SETUP_DIR is required}"
: "${CASE_NAME:?CASE_NAME is required}"
: "${ALPHA_LABEL:?ALPHA_LABEL is required}"
: "${ALPHA_VALUE:?ALPHA_VALUE is required}"

export EXP_DIR=$SETUP_DIR/existingTutorials/advect_cs
export OUTPUT_ROOT=$MITGCM_DIR/Sandbox/output/existingTutorials/test1/MITGCM_Williamson_TC1/advect_cs
export BUILD_DIR=$EXP_DIR/build/$CASE_NAME
export RUN_DIR=$OUTPUT_ROOT/alpha_${ALPHA_LABEL}/run_alpha_${ALPHA_LABEL}
export LOG_DIR=$SETUP_DIR/logs
export MPI_TASKS=${SLURM_NTASKS:-6}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-6}

START_TIME=$(date +%s)

finish() {
    local status=$1
    local end_time elapsed
    end_time=$(date +%s)
    elapsed=$((end_time-START_TIME))
    echo "======================================"
    echo "End: $(date)"
    echo "Status: $status"
    echo "Elapsed: $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"
    echo "======================================"
}

require_safe_dir() {
    local dir=$1
    case "$dir" in
        "$EXP_DIR"/build/*|"$OUTPUT_ROOT"/alpha_*/run_alpha_*) ;;
        *)
            echo "ERROR: refusing unsafe cleanup target: $dir" >&2
            exit 1
            ;;
    esac
}

empty_dir_contents() {
    local dir=$1
    require_safe_dir "$dir"
    mkdir -p "$dir"
    find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
}

stage_inputs() {
    cp "$EXP_DIR"/input/data "$RUN_DIR"/
    cp "$EXP_DIR"/input/data.pkg "$RUN_DIR"/
    cp "$EXP_DIR"/input/data.diagnostics "$RUN_DIR"/
    cp "$EXP_DIR"/input/eedata "$RUN_DIR"/
    cp "$EXP_DIR"/input/S.init "$RUN_DIR"/
    cp "$EXP_DIR"/input/T.init "$RUN_DIR"/
    printf "%s\n" "$ALPHA_VALUE" > "$RUN_DIR/tc1_alpha.txt"

    for face in 001 002 003 004 005 006; do
        ln -sf "$MITGCM_DIR/verification/tutorial_held_suarez_cs/input/grid_cs32.face${face}.bin" \
            "$RUN_DIR/grid_cs32.face${face}.bin"
    done
}

clean_run_after_plot() {
    "$MITGCM_DIR/.venv/bin/python" "$SETUP_DIR/tools/postprocess_advect_cs.py" --clean-run-products
}

echo "======================================"
echo "MITgcm advect_cs Williamson TC1"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Nodes:  ${SLURM_JOB_NODELIST:-unknown}"
echo "Case:   $CASE_NAME"
echo "Alpha:  $ALPHA_VALUE"
echo "Start:  $(date)"
echo "======================================"

mkdir -p "$BUILD_DIR" "$RUN_DIR" "$LOG_DIR"

echo "Cleaning build directory: $BUILD_DIR"
empty_dir_contents "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

"$MITGCM_DIR/tools/genmake2" \
    -mpi \
    -rootdir="$MITGCM_DIR" \
    -mods="$EXP_DIR/code" \
    -of="$MITGCM_DIR/tools/build_options/linux_amd64_gfortran"
make depend
make -j "$MAKE_JOBS"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created" >&2
    exit 1
fi

echo "Cleaning run directory: $RUN_DIR"
empty_dir_contents "$RUN_DIR"
cp mitgcmuv "$RUN_DIR/"
stage_inputs

cd "$RUN_DIR" || exit 1

echo "======================================"
echo "RUNNING MITgcm"
echo "MPI tasks = $MPI_TASKS"
echo "======================================"

set +e
mpirun -np "$MPI_TASKS" ./mitgcmuv
STATUS=$?
set -e

if [ "$STATUS" -ne 0 ]; then
    finish "$STATUS"
    exit "$STATUS"
fi

if [ "${RUN_POSTPROCESS:-0}" = "1" ]; then
    clean_run_after_plot
fi

finish 0
exit 0
