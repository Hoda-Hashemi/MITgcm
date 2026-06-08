#!/bin/bash

cd jobs/Arza || exit 1

sbatch job_tc3_c1_arza.slurm
sbatch job_tc3_c2_arza.slurm
