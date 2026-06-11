#!/bin/bash
#SBATCH --job-name=TC2_c4_large
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodelist=onode15
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC2/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC2/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail
module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12
source /home/hmh85/scratch/MITgcm/.venv/bin/activate

export MITGCM_DIR=/home/hmh85/scratch/MITgcm
export CASE_DIR=$MITGCM_DIR/Sandbox/vortexSphere_Williamson_TC2
export CASE_NAME=c4
export PREPARE_MODE=tc2
export ALPHA_VALUE=1.5707963267948966
export BUILD_DIR=$CASE_DIR/build/c4_large48
export RUN_DIR=$CASE_DIR/run_alpha_1.57
export MPI_TASKS=${SLURM_NTASKS:-48}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-48}
export FINAL_FIELD=Eta
export FINAL_ITER=0000017280

bash "$MITGCM_DIR/Sandbox/Scripts/run_williamson_job.sh"

