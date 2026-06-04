#!/bin/bash
#SBATCH --job-name=plot_tc1_c4_large
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

#SBATCH --output=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC1/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC1/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -e

module purge
module load python/3.12

VENV_DIR=/home/hmh85/scratch/MITgcm/.venv
source "$VENV_DIR/bin/activate"

TC1_DIR=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC1

RUN_DIR=$TC1_DIR/run_48mpi_alpha_c4
RUN_NAME=$(basename "$RUN_DIR")
PLOT_DIR=$TC1_DIR/Scripts/output/$RUN_NAME

test -f "$RUN_DIR/Eta.0000000000.data"
mkdir -p "$PLOT_DIR"

export RUN_DIR
export PLOT_DIR

echo "python = $(which python)"
python --version
echo "RUN_DIR = $RUN_DIR"
echo "PLOT_DIR = $PLOT_DIR"

cd "$TC1_DIR" || exit 1
python3 "$TC1_DIR/Scripts/HPC/snapshots.py" "$RUN_DIR"

echo "Plots saved in: $PLOT_DIR"
