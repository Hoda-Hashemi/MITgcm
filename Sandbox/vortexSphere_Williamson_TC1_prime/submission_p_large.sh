#!/bin/bash

set -euo pipefail

cd jobs/large || exit 1

sbatch job_tc1_prime_a0.slurm
sbatch job_tc1_prime_a157.slurm
