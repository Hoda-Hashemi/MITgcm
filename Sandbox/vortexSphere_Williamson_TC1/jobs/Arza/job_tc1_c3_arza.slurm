#!/bin/bash
#SBATCH --job-name=TC1_c3_arza
#SBATCH --account=hmh85
#SBATCH --partition=arza
#SBATCH --nodelist=anode14
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=60GB
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail
module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12
source /home/hmh85/scratch/MITgcm/.venv/bin/activate

export MITGCM_DIR=/home/hmh85/scratch/MITgcm
export CASE_DIR=$MITGCM_DIR/Sandbox/vortexSphere_Williamson_TC1
export CASE_NAME=c3
export PREPARE_MODE=tc1
export ALPHA_VALUE=1.5207963267948966
export BUILD_DIR=$CASE_DIR/build/c3_arza16
export RUN_DIR=$CASE_DIR/run_alpha_1.52
export MPI_TASKS=${SLURM_NTASKS:-16}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-16}
export FINAL_FIELD=S
export FINAL_ITER=0001036800

bash "$MITGCM_DIR/Sandbox/Scripts/run_williamson_job.sh"

