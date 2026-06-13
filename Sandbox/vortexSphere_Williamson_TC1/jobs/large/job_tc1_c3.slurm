#!/bin/bash
#SBATCH --job-name=TC1_c3_large
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodelist=onode16
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail
module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12
source /home/hmh85/scratch/MITgcm/.venv/bin/activate

export MITGCM_DIR=/home/hmh85/scratch/MITgcm
export CASE_DIR=$MITGCM_DIR/Sandbox/vortexSphere_Williamson_TC1
export CASE_NAME=c3
export ALPHA_VALUE=1.5207963267948966
export RUN_DIR=$CASE_DIR/run_alpha_1.52
export BUILD_DIR=$CASE_DIR/build/$CASE_NAME
export MPI_TASKS=${SLURM_NTASKS:-48}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-48}

START_TIME=$(date +%s)

TOOLS_DIR=$MITGCM_DIR/tools
BUILD_OPTIONS=$TOOLS_DIR/build_options/linux_amd64_gfortran
BUILD_MODS_DIR=$BUILD_DIR/code
PLOT_SCRIPT=$MITGCM_DIR/Sandbox/Scripts/TC1.py

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
        "$CASE_DIR"/build/*|"$CASE_DIR"/run_alpha_*) ;;
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

stage_case_inputs() {
    find "$CASE_DIR/input" -maxdepth 1 -type f -exec cp {} "$RUN_DIR/" \;
    sed -i "s/^ALPHA_RAD = .*/ALPHA_RAD = ${ALPHA_VALUE}/" "$RUN_DIR/gendata_ref.py"
    (cd "$RUN_DIR" && python3 gendata_ref.py)
}

clean_run_after_plot() {
    local path name
    shopt -s nullglob dotglob
    for path in "$RUN_DIR"/*; do
        name=${path##*/}
        case "$name" in
            data|data.pkg|data.diagnostics|eedata|eedata.mth|gendata_ref.py)
                continue
                ;;
            S_init.bin|u_init.bin|v_init.bin|bathymetry_flat4000.bin|mitgcmuv)
                continue
                ;;
            available_diagnostics.log|run_inventory.log|STDOUT*|STDERR*)
                continue
                ;;
        esac
        rm -rf -- "$path"
    done
    shopt -u nullglob dotglob
}

echo "======================================"
echo "MITgcm COMPILE AND RUN"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Nodes:  ${SLURM_JOB_NODELIST:-unknown}"
echo "Start:  $(date)"
echo "======================================"

mkdir -p "$BUILD_DIR" "$RUN_DIR"

echo "Cleaning build directory: $BUILD_DIR"
empty_dir_contents "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

mkdir -p "$BUILD_MODS_DIR"
find "$CASE_DIR/code" -maxdepth 1 -type f -exec cp {} "$BUILD_MODS_DIR/" \;

"$TOOLS_DIR/genmake2" -mpi -mods="$BUILD_MODS_DIR" -of="$BUILD_OPTIONS"
make depend
make -j "${MAKE_JOBS}"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created" >&2
    exit 1
fi

echo "Cleaning run directory: $RUN_DIR"
empty_dir_contents "$RUN_DIR"
cp mitgcmuv "$RUN_DIR/"
stage_case_inputs

cd "$RUN_DIR" || exit 1

echo "======================================"
echo "RUNNING MITGCM"
echo "MPI tasks = ${MPI_TASKS}"
echo "======================================"

set +e
mpirun -np "${MPI_TASKS}" ./mitgcmuv
STATUS=$?
set -e

if [ "$STATUS" -ne 0 ]; then
    finish "$STATUS"
    exit "$STATUS"
fi

echo "======================================"
echo "PLOTTING TC1"
echo "RUN_DIR = $RUN_DIR"
echo "======================================"

set +e
TC1_RUN_DIRS="$RUN_DIR" python3 "$PLOT_SCRIPT"
STATUS=$?
set -e

if [ "$STATUS" -ne 0 ]; then
    finish "$STATUS"
    exit "$STATUS"
fi

clean_run_after_plot
finish 0
exit 0

