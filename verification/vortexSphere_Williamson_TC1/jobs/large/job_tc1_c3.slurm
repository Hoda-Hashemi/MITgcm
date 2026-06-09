#!/bin/bash
#SBATCH --job-name=TC1_alpha_1.52
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodelist=onode14
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC1/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC1/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12

VENV_DIR=/home/hmh85/scratch/MITgcm/.venv
source "$VENV_DIR/bin/activate"

ALPHA=1.5207963267948966
MITGCM_DIR=/home/hmh85/scratch/MITgcm
CASE_DIR=$MITGCM_DIR/verification/vortexSphere_Williamson_TC1
INPUT_DIR=$CASE_DIR/input
BUILD_DIR=$CASE_DIR/build/alpha_1.52_large_48
MODS_DIR=$BUILD_DIR/mods
RUN_DIR=$CASE_DIR/run_alpha_1.52
MPI_TASKS=${SLURM_NTASKS:-48}
MAKE_JOBS=${SLURM_CPUS_ON_NODE:-48}
RUN_INPUT_FILES=(data data.diagnostics data.pkg eedata gendata_ref.py)

echo "======================================"
echo "Williamson TC1 alpha=$ALPHA"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Node:   ${SLURM_JOB_NODELIST:-local}"
echo "Tasks:  $MPI_TASKS"
echo "Start:  $(date)"
echo "======================================"

if [ "$MPI_TASKS" -ne 48 ]; then
    echo "ERROR: this job builds a 48-rank decomposition"
    exit 1
fi

START_TIME=$(date +%s)

# Build a private 48-rank decomposition: 12 x 4 tiles.
mkdir -p "$BUILD_DIR" "$MODS_DIR"
cp "$CASE_DIR/code/"* "$MODS_DIR/"
sed -E -i.bak \
    -e 's/sNx =[[:space:]]*[0-9]+[[:space:]]*[)]/sNx =  120 )/' \
    -e 's/sNy =[[:space:]]*[0-9]+[[:space:]]*[)]/sNy =  180 )/' \
    -e 's/nPx =[[:space:]]*[0-9]+[[:space:]]*[)]/nPx =   12 )/' \
    -e 's/nPy =[[:space:]]*[0-9]+[[:space:]]*[)]/nPy =    4 )/' \
    "$MODS_DIR/SIZE.h"
rm -f "$MODS_DIR/SIZE.h.bak"

cd "$BUILD_DIR"
make CLEAN >/dev/null 2>&1 || true
rm -f Makefile

"$MITGCM_DIR/tools/genmake2" \
    -rootdir="$MITGCM_DIR" \
    -mpi \
    -mods="$MODS_DIR" \
    -of="$MITGCM_DIR/tools/build_options/linux_amd64_gfortran"

make depend
make -j "$MAKE_JOBS"
test -x mitgcmuv

# Remove old model outputs, then install fresh inputs.
mkdir -p "$RUN_DIR"
find "$RUN_DIR" -maxdepth 1 -type f \
    \( -name "*.data" -o -name "*.meta" -o -name "*.bin" \
       -o -name "STDOUT*" -o -name "STDERR*" \
       -o -name "pickup*" -o -name "available_diagnostics.log" \) \
    -delete

for file in "${RUN_INPUT_FILES[@]}"; do
    cp "$INPUT_DIR/$file" "$RUN_DIR/"
done

cp "$BUILD_DIR/mitgcmuv" "$RUN_DIR/"

cd "$RUN_DIR"
sed -E -i.bak \
    -e 's/^[[:space:]]*nTimeSteps=.*/ nTimeSteps= 1036800 ,/' \
    -e 's/^[[:space:]]*deltaT=.*/ deltaT=1.0,/' \
    -e 's/^[[:space:]]*dumpFreq=.*/ dumpFreq=86400.,/' \
    -e 's/^[[:space:]]*monitorFreq=.*/ monitorFreq=86400.,/' \
    data
rm -f data.bak
sed -E -i.bak \
    "s/^ALPHA_RAD[[:space:]]*=.*/ALPHA_RAD = $ALPHA/" \
    gendata_ref.py
rm -f gendata_ref.py.bak
python3 gendata_ref.py

echo "Run configuration:"
grep -E "nTimeSteps|deltaT=|dumpFreq|monitorFreq" data
grep -E "^ALPHA_RAD" gendata_ref.py

set +e
mpirun -np "$MPI_TASKS" ./mitgcmuv
STATUS=$?
set -e

if grep -q "ABNORMAL END\\|\\*\\*\\* ERROR \\*\\*\\*" \
    STDERR.* STDOUT.* 2>/dev/null; then
    echo "ERROR: MITgcm reported an abnormal termination"
    STATUS=1
fi

if ! grep -q "PROGRAM MAIN: Execution ended Normally" STDOUT.0000; then
    echo "ERROR: normal MITgcm completion marker was not found"
    STATUS=1
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "======================================"
echo "End:     $(date)"
echo "Status:  $STATUS"
echo "Elapsed: $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"
echo "Run dir: $RUN_DIR"
echo "======================================"

exit "$STATUS"
