#!/bin/bash

set -euo pipefail

cd jobs/large || exit 1

runid=$(sbatch job_tc1_c1.slurm | awk '{print $4}')
runid=$(sbatch --dependency=afterok:$runid job_tc1_c2.slurm | awk '{print $4}')
runid=$(sbatch --dependency=afterok:$runid job_tc1_c3.slurm | awk '{print $4}')
sbatch --dependency=afterok:$runid job_tc1_c4.slurm
