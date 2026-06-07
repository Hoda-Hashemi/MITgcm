#!/bin/bash
#SBATCH --job-name=williamson_TC1
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode15
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6

TC1_DIR=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC1
REF_DIR=/home/hmh85/scratch/MITgcm/verification/vortexSphere_mitgcm_referenceCase

BUILD_DIR=$TC1_DIR/build
INPUT_DIR=$TC1_DIR/input
RUN_DIR=$TC1_DIR/run_48mpi_alpha_PIOVER2

RUN_NAME=$(basename "$RUN_DIR")
PLOT_DIR=$TC1_DIR/Scripts/HPC/$RUN_NAME

mkdir -p "$RUN_DIR" "$PLOT_DIR"

echo "Compile..."
cd "$BUILD_DIR" || exit 1
rm -f Makefile
make CLEAN 2>/dev/null || true

../../../tools/genmake2 \
  -mpi \
  -mods=../code \
  -of=../../../tools/build_options/linux_amd64_gfortran

make depend
make -j 8

[ -f mitgcmuv ] || { echo "ERROR: mitgcmuv not created"; exit 1; }

echo "Prepare run..."
rm -f "$RUN_DIR"/*
cp mitgcmuv "$RUN_DIR/"
cp "$INPUT_DIR"/* "$RUN_DIR/"

echo "Run MITgcm..."
cd "$RUN_DIR" || exit 1
mpirun -np 48 ./mitgcmuv
STATUS=$?

if [ "$STATUS" -ne 0 ]; then
  echo "MITgcm failed. Skipping plots."
  exit "$STATUS"
fi

echo "Plot..."
export RUN_DIR
export PLOT_DIR

cd "$TC1_DIR" || exit 1
python3 "$REF_DIR/Scripts/HPC/snapshots.py"

echo "Plots saved in: $PLOT_DIR"

