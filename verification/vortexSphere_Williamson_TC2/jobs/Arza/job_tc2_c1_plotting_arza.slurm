#!/bin/bash
#SBATCH --job-name=plot_tc2_c1_arza
#SBATCH --account=hmh85
#SBATCH --partition=arza
#SBATCH --nodelist=anode02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60GB
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

module purge
module load python/3.12

TC2_DIR=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2
RUN_DIR=$TC2_DIR/run_48mpi_alpha_c1
RUN_NAME=$(basename "$RUN_DIR")
PLOT_DIR=$TC2_DIR/Scripts/output/$RUN_NAME
VENV_DIR=/home/hmh85/scratch/MITgcm/.venv

source "$VENV_DIR/bin/activate"

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

cd "$TC2_DIR"
python "$TC2_DIR/Scripts/HPC/snapshots.py" "$RUN_DIR"

echo "Plots saved in: $PLOT_DIR"
