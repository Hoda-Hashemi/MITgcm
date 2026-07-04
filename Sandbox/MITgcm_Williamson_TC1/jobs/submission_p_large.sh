#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/large" || exit 1

sbatch job_tc1_c1.slurm
sbatch job_tc1_c2.slurm
sbatch job_tc1_c3.slurm
sbatch job_tc1_c4.slurm
