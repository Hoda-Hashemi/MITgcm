#!/bin/bash

#SBATCH --job-name=mitgcm_ref_compile_run
#SBATCH --account=hmh85

#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=0-06:00:00

#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

echo "======================================"
echo "MITgcm compile + run"
echo "Job ID: $SLURM_JOB_ID"
echo "Start:  $(date)"
echo "Nodes:  $SLURM_JOB_NODELIST"
echo "======================================"

START_TIME=$(date +%s)

CASE_DIR=/home/hmh85/scratch/MITgcm/verification/vortexSphere_mitgcm_referenceCase
BUILD_DIR=$CASE_DIR/build
RUN_DIR=$CASE_DIR/run

# -----------------------------
# Compile
# -----------------------------
cd $BUILD_DIR || exit 1

echo "Cleaning build directory..."
rm -f Makefile
make CLEAN 2>/dev/null || true

echo "Generating Makefile..."
../../../tools/genmake2 \
  -mpi \
  -mods=../code \
  -of=../../../tools/build_options/linux_amd64_gfortran

echo "Making dependencies..."
make depend

echo "Compiling..."
make -j 8

if [ ! -f mitgcmuv ]; then
    echo "ERROR: mitgcmuv was not created."
    exit 1
fi

# -----------------------------
# Run
# -----------------------------
echo "Copying mitgcmuv to run directory..."
cp mitgcmuv $RUN_DIR/

cd $RUN_DIR || exit 1

echo "Running MITgcm..."
mpirun -np 8 ./mitgcmuv

STATUS=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "======================================"
echo "End:     $(date)"
echo "Status:  $STATUS"
echo "Elapsed: $((ELAPSED/3600)) h $(((ELAPSED%3600)/60)) min $((ELAPSED%60)) sec"
echo "======================================"

exit $STATUS

