#!/bin/bash
#SBATCH --job-name=TC2_c3_large
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode14
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00

#SBATCH --output=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/jobs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/jobs/slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6

CASE_NAME=c3
ALPHA=1.5207963267948966
EXPECTED_INPUT_BYTES=$((1440 * 720 * 4))
FINAL_ITER=0000017280
MITGCM_DIR=/home/hmh85/scratch/MITgcm
CASE_DIR=$MITGCM_DIR/verification/vortexSphere_Williamson_TC2
BUILD_DIR=$CASE_DIR/build/$CASE_NAME
MAKE_JOBS=${SLURM_CPUS_ON_NODE:-$(nproc)}
MPI_TASKS=${SLURM_NTASKS:-48}
NODE_TAG=${SLURM_JOB_NODELIST:-onode14}
NODE_TAG=${NODE_TAG//[^A-Za-z0-9_.-]/_}
RUN_DIR=$CASE_DIR/run_${MPI_TASKS}mpi_${NODE_TAG}_alpha_$CASE_NAME
RUN_INPUT_FILES=(
    data
    data.diagnostics
    data.pkg
    eedata
    "gendata_${CASE_NAME}.py"
    bathymetry_flat4000.bin
)

echo "======================================"
echo "MITgcm compile + run"
echo "Case:   $CASE_NAME"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Nodes:  ${SLURM_JOB_NODELIST:-local}"
echo "Start:  $(date)"
echo "======================================"

START_TIME=$(date +%s)

if [ "$MPI_TASKS" -ne 48 ]; then
    echo "ERROR: this job builds a 48-rank decomposition"
    exit 1
fi

for file in "${RUN_INPUT_FILES[@]}"; do
    if [ ! -r "$CASE_DIR/input/$file" ]; then
        echo "ERROR: missing input file: $CASE_DIR/input/$file"
        exit 1
    fi
done

for file in bathymetry_flat4000.bin "eta_init_${CASE_NAME}.bin" "u_init_${CASE_NAME}.bin" "v_init_${CASE_NAME}.bin"; do
    if [ ! -r "$CASE_DIR/input/$file" ]; then
        echo "ERROR: missing input file: $CASE_DIR/input/$file"
        exit 1
    fi
    input_bytes=$(wc -c < "$CASE_DIR/input/$file")
    if [ "$input_bytes" -ne "$EXPECTED_INPUT_BYTES" ]; then
        echo "ERROR: $file has $input_bytes bytes; expected $EXPECTED_INPUT_BYTES"
        exit 1
    fi
done

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

make CLEAN >/dev/null 2>&1 || true
rm -f Makefile

"$MITGCM_DIR/tools/genmake2" \
    -rootdir="$MITGCM_DIR" \
    -mpi \
    -mods="$CASE_DIR/code" \
    -of="$MITGCM_DIR/tools/build_options/linux_amd64_gfortran"

make depend
make -j "$MAKE_JOBS"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created"
    exit 1
fi

mkdir -p "$RUN_DIR"

find "$RUN_DIR" -maxdepth 1 -type f \
    \( -name "*.data" -o -name "*.meta" -o -name "*.bin" -o -name "gendata*.py" -o -name "STDOUT*" -o -name "STDERR*" \) \
    -delete

for file in "${RUN_INPUT_FILES[@]}"; do
    cp "$CASE_DIR/input/$file" "$RUN_DIR/"
done

cp "$CASE_DIR/input/eta_init_${CASE_NAME}.bin" "$RUN_DIR/eta_init.bin"
cp "$CASE_DIR/input/u_init_${CASE_NAME}.bin" "$RUN_DIR/u_init.bin"
cp "$CASE_DIR/input/v_init_${CASE_NAME}.bin" "$RUN_DIR/v_init.bin"

echo "======================================"
echo "Running MITgcm"
echo "MPI tasks = $MPI_TASKS"
echo "alpha = $ALPHA"
echo "======================================"

cd "$RUN_DIR"
cp "$BUILD_DIR/mitgcmuv" "$RUN_DIR/"

if mpirun -np "$MPI_TASKS" ./mitgcmuv; then
    STATUS=0
else
    STATUS=$?
fi

if grep -q "ABNORMAL END\\|\\*\\*\\* ERROR \\*\\*\\*" \
    STDERR.* STDOUT.* 2>/dev/null; then
    echo "ERROR: MITgcm reported an abnormal termination"
    STATUS=1
fi

if [ ! -f "Eta.${FINAL_ITER}.data" ]; then
    echo "ERROR: MITgcm did not reach final iteration $FINAL_ITER"
    STATUS=1
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME-START_TIME))

echo "======================================"
echo "End: $(date)"
echo "Status: $STATUS"
echo "Elapsed: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
echo "======================================"

exit "$STATUS"
