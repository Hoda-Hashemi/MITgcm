#!/bin/bash

cd jobs/Arza || exit 1

runid=$(sbatch job_tc2_c1_arza.slurm | awk '{print $4}')
sbatch --dependency=afterok:$runid job_tc2_c1_plotting_arza.slurm

runid=$(sbatch job_tc2_c2_arza.slurm | awk '{print $4}')
sbatch --dependency=afterok:$runid job_tc2_c2_plotting_arza.slurm

runid=$(sbatch job_tc2_c3_arza.slurm | awk '{print $4}')
sbatch --dependency=afterok:$runid job_tc2_c3_plotting_arza.slurm

runid=$(sbatch job_tc2_c4_arza.slurm | awk '{print $4}')
sbatch --dependency=afterok:$runid job_tc2_c4_plotting_arza.slurm