#!/bin/bash

set -euo pipefail

cd large || exit 1

sbatch job_tc3_c1.slurm
sbatch job_tc3_c2.slurm
