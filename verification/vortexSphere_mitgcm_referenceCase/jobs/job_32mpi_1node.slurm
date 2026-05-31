#!/bin/bash

#SBATCH --job-name=mitgcm_32mpi
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode14
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6

echo "======================================"
echo "MITgcm compile + run"
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes:  $SLURM_JOB_NODELIST"
echo "Start:  $(date)"
echo "======================================"

START_TIME=$(date +%s)

CASE_DIR=/scratch/p8190783-hmh85/MITgcm/verification/vortexSphere_mitgcm_referenceCase
BUILD_DIR=$CASE_DIR/build
RUN_DIR=$CASE_DIR/run_32mpi_1node

# -----------------------------
# Compile
# -----------------------------
cd $BUILD_DIR || exit 1

rm -f Makefile
make CLEAN >/dev/null 2>&1 || true

../../../tools/genmake2 \
    -mpi \
    -mods=../code \
    -of=../../../tools/build_options/linux_amd64_gfortran

make depend
make -j 8

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created"
    exit 1
fi

# -----------------------------
# Prepare run directory
# -----------------------------
mkdir -p $RUN_DIR

find $RUN_DIR -maxdepth 1 -type f \
    \( -name "*.data" -o -name "*.meta" -o -name "STDOUT*" -o -name "STDERR*" \) \
    -delete

cp mitgcmuv $RUN_DIR/
cp $CASE_DIR/input/* $RUN_DIR/

cd $RUN_DIR || exit 1

echo "======================================"
echo "Running MITgcm"
echo "MPI tasks = $SLURM_NTASKS"
echo "======================================"

mpirun -np $SLURM_NTASKS ./mitgcmuv

STATUS=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME-START_TIME))

echo "======================================"
echo "End: $(date)"
echo "Status: $STATUS"
echo "Elapsed: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
echo "======================================"

exit $STATUS
