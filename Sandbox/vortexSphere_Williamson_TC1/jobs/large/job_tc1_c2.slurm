#!/bin/bash
#SBATCH --job-name=TC1_c2_large
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodelist=onode14
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
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
export CASE_NAME=c2
export PREPARE_MODE=tc1
export ALPHA_VALUE=0.05
export BUILD_DIR=$CASE_DIR/build/c2_large16
export RUN_DIR=$CASE_DIR/run_alpha_0.05
export MPI_TASKS=${SLURM_NTASKS:-16}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-16}
export FINAL_FIELD=S
export FINAL_ITER=0001036800

bash "$MITGCM_DIR/Sandbox/Scripts/run_williamson_job.sh"

