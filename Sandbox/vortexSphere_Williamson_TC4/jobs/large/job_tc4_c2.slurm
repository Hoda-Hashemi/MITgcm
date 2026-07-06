#!/bin/bash
#SBATCH --job-name=TC4_c2_u0_40_large48
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC4/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC4/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12
source /home/hmh85/scratch/MITgcm/.venv/bin/activate

export MITGCM_DIR=/home/hmh85/scratch/MITgcm
export CASE_DIR=$MITGCM_DIR/Sandbox/vortexSphere_Williamson_TC4
export CASE_NAME=c2
export TC4_U0_VALUE=40.0
export DELTA_T=60.0
export TOTAL_SECONDS=432000
export BUILD_DIR=$CASE_DIR/build/c2_large48
export BUILD_MODS_DIR=$BUILD_DIR/code
export RUN_DIR=$CASE_DIR/run_u0_40
export LOG_DIR=$CASE_DIR/logs
export MPI_TASKS=${SLURM_NTASKS:-48}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-48}
export N_TIME_STEPS=$(python3 -c "print(int(round(${TOTAL_SECONDS}/${DELTA_T})))")

START_TIME=$(date +%s)
TOOLS_DIR=$MITGCM_DIR/tools
BUILD_OPTIONS=$TOOLS_DIR/build_options/linux_amd64_gfortran
PLOT_SCRIPT=$MITGCM_DIR/Sandbox/Scripts/TC4.py
CONSERVATION_SCRIPT=$MITGCM_DIR/Sandbox/Scripts/conservation_diagnostics.py
SITE_SCRIPT=$MITGCM_DIR/Sandbox/Scripts/github_pages.py

finish() {
    local status=$1 end_time elapsed
    end_time=$(date +%s)
    elapsed=$((end_time-START_TIME))
    echo "End: $(date) | Status: $status | Elapsed: $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"
}

require_safe_dir() {
    local dir=$1
    case "$dir" in
        "$CASE_DIR"/build/*|"$CASE_DIR"/run_u0_40) ;;
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

patch_size_for_48_ranks() {
    python3 - "$BUILD_MODS_DIR/SIZE.h" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text()
for key, value in {
    "sNx": 120,
    "sNy": 180,
    "nPx": 12,
    "nPy": 4,
}.items():
    text = re.sub(
        rf"(PARAMETER\s*\(\s*{key}\s*=\s*)\d+(\s*\))",
        rf"\g<1>{value:5d}\2",
        text,
    )
path.write_text(text)
PY
}

stage_case_inputs() {
    find "$CASE_DIR/input" -maxdepth 1 -type f \
        ! -name "eta_init*.bin" \
        ! -name "u_init*.bin" \
        ! -name "v_init*.bin" \
        ! -name "bathymetry*.bin" \
        -exec cp {} "$RUN_DIR/" \;
    sed -i "s/^ nTimeSteps=.*/ nTimeSteps=${N_TIME_STEPS},/" "$RUN_DIR/data"
    sed -i "s/^ deltaT=.*/ deltaT=${DELTA_T},/" "$RUN_DIR/data"
    (cd "$RUN_DIR" && python3 gendata_ref.py)
}

echo "======================================"
echo "MITgcm TC4 c2 compile/run/postprocess"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Nodes:  ${SLURM_JOB_NODELIST:-unknown}"
echo "MPI:    ${MPI_TASKS}"
echo "U0:     ${TC4_U0_VALUE} m/s"
echo "Start:  $(date)"
echo "======================================"

mkdir -p "$LOG_DIR"
empty_dir_contents "$BUILD_DIR"
empty_dir_contents "$RUN_DIR"

cd "$BUILD_DIR" || exit 1
mkdir -p "$BUILD_MODS_DIR"
find "$CASE_DIR/code" -maxdepth 1 -type f -exec cp {} "$BUILD_MODS_DIR/" \;
patch_size_for_48_ranks

"$TOOLS_DIR/genmake2" -mpi -mods="$BUILD_MODS_DIR" -of="$BUILD_OPTIONS"
make depend
make -j "${MAKE_JOBS:-8}"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created" >&2
    exit 1
fi

cp mitgcmuv "$RUN_DIR/"
stage_case_inputs

cd "$RUN_DIR" || exit 1
set +e
mpirun -np "${MPI_TASKS}" ./mitgcmuv
STATUS=$?
set -e
if [ "$STATUS" -ne 0 ]; then
    finish "$STATUS"
    exit "$STATUS"
fi

set +e
TC4_RUN_DIRS="$RUN_DIR" python3 "$PLOT_SCRIPT"
STATUS=$?
set -e
if [ "$STATUS" -ne 0 ]; then
    finish "$STATUS"
    exit "$STATUS"
fi

python3 "$CONSERVATION_SCRIPT" "$RUN_DIR" --cases TC4
python3 "$SITE_SCRIPT"

finish 0
exit 0
