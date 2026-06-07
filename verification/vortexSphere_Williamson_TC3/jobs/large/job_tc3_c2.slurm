#!/bin/bash
#SBATCH --job-name=TC3_c2_large
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode16
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00

#SBATCH --output=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC3/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC3/logs/slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12

VENV_DIR=/home/hmh85/scratch/MITgcm/.venv
source "$VENV_DIR/bin/activate"

CASE_NAME=c2
ALPHA=1.0471975511965976
MITGCM_DIR=/home/hmh85/scratch/MITgcm
CASE_DIR=$MITGCM_DIR/verification/vortexSphere_Williamson_TC3
BUILD_DIR=$CASE_DIR/build/$CASE_NAME
MAKE_JOBS=${SLURM_CPUS_ON_NODE:-$(nproc)}
MPI_TASKS=${SLURM_NTASKS:-16}
NODE_TAG=${SLURM_JOB_NODELIST:-onode16}
NODE_TAG=${NODE_TAG//[^A-Za-z0-9_.-]/_}
RUN_DIR=$CASE_DIR/run_${MPI_TASKS}mpi_${NODE_TAG}_alpha_$CASE_NAME
RUN_INPUT_FILES=(data data.diagnostics data.pkg eedata gendata.py)

echo "======================================"
echo "MITgcm compile + run"
echo "Case:   $CASE_NAME"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Nodes:  ${SLURM_JOB_NODELIST:-local}"
echo "Start:  $(date)"
echo "======================================"

START_TIME=$(date +%s)

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

make CLEAN >/dev/null 2>&1 || true
rm -f Makefile

"$MITGCM_DIR/tools/genmake2" \
    -rootdir="$MITGCM_DIR" \
    -mpi \
    -mods="$CASE_DIR/code" \
    -of="$MITGCM_DIR/tools/build_options/linux_amd64_gfortran"

make depend
make -j "$MAKE_JOBS"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created"
    exit 1
fi

mkdir -p "$RUN_DIR"

find "$RUN_DIR" -maxdepth 1 -type f \
    \( -name "*.data" -o -name "*.meta" -o -name "*.bin" -o -name "STDOUT*" -o -name "STDERR*" \) \
    -delete

for file in "${RUN_INPUT_FILES[@]}"; do
    cp "$CASE_DIR/input/$file" "$RUN_DIR/"
done

cd "$RUN_DIR"
python3 gendata.py --alpha "$ALPHA"

cp "$BUILD_DIR/mitgcmuv" "$RUN_DIR/"

echo "======================================"
echo "Running MITgcm"
echo "MPI tasks = $MPI_TASKS"
echo "alpha = $ALPHA"
echo "======================================"

if mpirun -np "$MPI_TASKS" ./mitgcmuv; then
    STATUS=0
else
    STATUS=$?
fi

if grep -q "ABNORMAL END\\|\\*\\*\\* ERROR \\*\\*\\*" \
    STDERR.* STDOUT.* 2>/dev/null; then
    echo "ERROR: MITgcm reported an abnormal termination"
    STATUS=1
fi

if [ ! -f Eta.0000000000.data ]; then
    echo "ERROR: MITgcm did not write Eta.0000000000.data"
    STATUS=1
fi

if [ "$STATUS" -eq 0 ]; then
    RUN_NAME=$(basename "$RUN_DIR")
    PLOT_DIR=$CASE_DIR/Scripts/output/$RUN_NAME

    echo "======================================"
    echo "Plotting snapshots"
    echo "======================================"

    if [ ! -f "$RUN_DIR/Eta.0000000000.data" ]; then
        echo "ERROR: missing $RUN_DIR/Eta.0000000000.data" >&2
        exit 1
    fi

    mkdir -p "$PLOT_DIR"

    export RUN_DIR PLOT_DIR

    echo "python = $(which python)"
    python --version
    echo "RUN_DIR = $RUN_DIR"
    echo "PLOT_DIR = $PLOT_DIR"

    cd "$CASE_DIR"
    python "$CASE_DIR/Scripts/HPC/snapshots.py" "$RUN_DIR"

    echo "Plots saved in: $PLOT_DIR"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME-START_TIME))

echo "======================================"
echo "End: $(date)"
echo "Status: $STATUS"
echo "Elapsed: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
echo "======================================"

exit "$STATUS"
