#!/bin/bash

set -euo pipefail

cd large || exit 1

sbatch job_tc6_c1.slurm
