#!/usr/bin/env bash

#SBATCH --job-name= TC1
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

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO_DIR=$(cd "$CASE_DIR/../.." && pwd)
BUILD_DIR="$CASE_DIR/build"
RUN_DIR="$CASE_DIR/run"
GENMAKE=${GENMAKE:-"$REPO_DIR/tools/genmake2"}
BUILD_OPTIONS=${BUILD_OPTIONS:-"$REPO_DIR/tools/build_options/linux_amd64_gfortran"}

mkdir -p "$BUILD_DIR" "$RUN_DIR"

python3 "$CASE_DIR/input/gendata_tc1.py"

cd "$BUILD_DIR"
rm -f Makefile

"$GENMAKE" -mpi -mods=../code -of="$BUILD_OPTIONS"
make depend
make -j "${SLURM_CPUS_PER_TASK:-1}"

find "$RUN_DIR" -maxdepth 1 -type f \( -name "*.data" -o -name "*.meta" -o -name "STDOUT*" -o -name "STDERR*" \) -delete
cp mitgcmuv "$RUN_DIR/"
cp "$CASE_DIR"/input/* "$RUN_DIR/"

cd "$RUN_DIR"

if command -v srun >/dev/null 2>&1; then
  srun -n "${SLURM_NTASKS:-1}" ./mitgcmuv
else
  mpirun -np "${SLURM_NTASKS:-1}" ./mitgcmuv
fi

