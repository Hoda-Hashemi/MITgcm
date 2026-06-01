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
module load python

REFERENCE_CASE="/path/to/referenceCase"
PLOT_SCRIPT="$REFERENCE_CASE/Scripts/HPC/save_snapshots.py"

cd "$REFERENCE_CASE"

for RUN_DIR in run_*mpi*; do
    if [ -d "$RUN_DIR" ]; then
        echo "Plotting $RUN_DIR"
        python "$PLOT_SCRIPT" "$REFERENCE_CASE/$RUN_DIR"
    fi
done

