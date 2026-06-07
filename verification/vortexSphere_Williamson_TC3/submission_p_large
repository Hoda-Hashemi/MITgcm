#!/bin/bash

set -euo pipefail

cd jobs/large || exit 1

runid=$(sbatch job_tc3_c1.slurm | awk '{print $4}')
sbatch --dependency=afterok:$runid job_tc3_c2.slurm
