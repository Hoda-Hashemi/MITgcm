#!/bin/bash
#SBATCH --job-name=TC1_advect_c1
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/p8190783-hmh85/MITgcm/Sandbox/MITgcm_Williamson_TC1/logs/slurm-%j.out
#SBATCH --error=/scratch/p8190783-hmh85/MITgcm/Sandbox/MITgcm_Williamson_TC1/logs/slurm-%j.err

set -euo pipefail
module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12

export MITGCM_DIR=/scratch/p8190783-hmh85/MITgcm
export SETUP_DIR=$MITGCM_DIR/Sandbox/MITgcm_Williamson_TC1
export CASE_NAME=c1
export ALPHA_LABEL=0
export ALPHA_VALUE=0.0

"$SETUP_DIR/jobs/run_advect_cs_tc1.sh"
