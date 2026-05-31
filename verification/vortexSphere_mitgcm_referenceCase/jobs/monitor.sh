BSCRIPT=${1:-job_gaussian_3.slurm}

JOBID=$(sbatch "$JOBSCRIPT" | awk '{print $4}')

echo "Submitted job: $JOBID"
echo

while true; do
  clear
  echo "JOBID: $JOBID"
  echo

  echo "========== QUEUE =========="
  squeue -j "$JOBID" -o "%.10i %.12P %.20j %.8u %.2t %.12M %.6D %R"

  echo
  echo "========== TIME / STATE =========="
  sacct -j "$JOBID" --format=JobID,Partition,State,Elapsed,ExitCode

  sleep 10
done

