#!/bin/bash
#SBATCH --job-name=TC4_c1_normal16
#SBATCH --account=hmh85
#SBATCH --partition=normal
#SBATCH --nodelist=onode08
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
##SBATCH --exclusive
#SBATCH --mem=0
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
export CASE_NAME=c1
export TC4_U0_VALUE=20.0
export DELTA_T=60.0
export TOTAL_SECONDS=432000
export BUILD_DIR=$CASE_DIR/build/c1_normal16
export RUN_DIR=$CASE_DIR/run_u0_20
export LOG_DIR=$CASE_DIR/logs
export MPI_TASKS=${SLURM_NTASKS:-48}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-48}
export N_TIME_STEPS=$(python3 -c "print(int(round(${TOTAL_SECONDS}/${DELTA_T})))")

START_TIME=$(date +%s)
TOOLS_DIR=$MITGCM_DIR/tools
BUILD_OPTIONS=$TOOLS_DIR/build_options/linux_amd64_gfortran
MODS_DIR=$CASE_DIR/code
PLOT_SCRIPT=$MITGCM_DIR/Sandbox/Scripts/TC4.py

finish() {
    local status=$1 end_time elapsed
    end_time=$(date +%s)
    elapsed=$((end_time-START_TIME))
    echo "End: $(date) | Status: $status | Elapsed: $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"
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

mkdir -p "$LOG_DIR" "$BUILD_DIR" "$RUN_DIR"
cd "$BUILD_DIR" || exit 1
"$TOOLS_DIR/genmake2" -mpi -mods="$MODS_DIR" -of="$BUILD_OPTIONS"
make depend
make -j "${MAKE_JOBS:-8}"
cp mitgcmuv "$RUN_DIR/"
stage_case_inputs
cd "$RUN_DIR" || exit 1
set +e
mpirun -np "${MPI_TASKS:-$SLURM_NTASKS}" ./mitgcmuv
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
finish "$STATUS"
exit "$STATUS"
