#!/bin/bash
#SBATCH --job-name=sim_mitgcm
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode13
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -e

module purge
module load python/3.12

REFERENCE_CASE="/home/hmh85/scratch/MITgcm/verification/vortexSphere_mitgcm_referenceCase"
PLOT_SCRIPT="$REFERENCE_CASE/Scripts/HPC/snapshots.py"

source /home/hmh85/scratch/MITgcm/.venv/bin/activate

cd "$REFERENCE_CASE"

# for RUN_DIR in run_*mpi*; do
#     if [ -d "$RUN_DIR" ]; then
#         python "$PLOT_SCRIPT" "$REFERENCE_CASE/$RUN_DIR"
#     fi
# done

RUN_DIR="run_48mpi_1node"
python "$PLOT_SCRIPT" "$REFERENCE_CASE/$RUN_DIR"
