#!/bin/bash
#SBATCH --job-name=mitgcm
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6

CASE_DIR=/home/hmh85/scratch/MITgcm/verification/vortexSphere_mitgcm_referenceCase
BUILD_DIR=${CASE_DIR}/build
RUN_DIR=${CASE_DIR}/run

echo "======================================"
echo "MITgcm multi-node compile and run"
echo "Job ID:   ${SLURM_JOB_ID}"
echo "Start:    $(date)"
echo "Nodes:    ${SLURM_JOB_NODELIST}"
echo "Tasks:    ${SLURM_NTASKS}"
echo "======================================"

START_TIME=$(date +%s)

cd "${BUILD_DIR}"

rm -f Makefile
make CLEAN >/dev/null 2>&1 || true

../../../tools/genmake2 \
  -mpi \
  -mods=../code \
  -of=../../../tools/build_options/linux_amd64_gfortran

make depend
make -j 8

test -f mitgcmuv

rm -rf "${RUN_DIR:?}"/*
cp mitgcmuv "${RUN_DIR}/"
cp "${CASE_DIR}/input/"* "${RUN_DIR}/"

cd "${RUN_DIR}"

mpirun -np 8 ./mitgcmuv
STATUS=$?

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "======================================"
echo "End:      $(date)"
echo "Status:   ${STATUS}"
echo "Elapsed:  $((ELAPSED/3600)) h $(((ELAPSED%3600)/60)) min $((ELAPSED%60)) sec"
echo "======================================"

exit "${STATUS}"
