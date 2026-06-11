#!/bin/bash
#SBATCH --job-name=TC1p_a157_large
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodelist=onode15
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1_prime/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1_prime/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail
module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12
source /home/hmh85/scratch/MITgcm/.venv/bin/activate

export MITGCM_DIR=/home/hmh85/scratch/MITgcm
export CASE_DIR=$MITGCM_DIR/Sandbox/vortexSphere_Williamson_TC1_prime
export CASE_NAME=a157
export PREPARE_MODE=tc1
export ALPHA_VALUE=1.5707963267948966
export BUILD_DIR=$CASE_DIR/build/a157_large48
export RUN_DIR=$CASE_DIR/run_alpha_1.57
export MPI_TASKS=${SLURM_NTASKS:-48}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-48}
export PATCH_SIZE_H=1
export SIZE_SNX=120
export SIZE_SNY=180
export SIZE_NPX=12
export SIZE_NPY=4
export FINAL_FIELD=Eta
export FINAL_ITER=0000017280


START_TIME=$(date +%s)

echo "======================================"
echo "MITgcm compile + run"
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes:  $SLURM_JOB_NODELIST"
echo "Start:  $(date)"
echo "======================================"

TOOLS_DIR=$MITGCM_DIR/tools
BUILD_OPTIONS=$TOOLS_DIR/build_options/linux_amd64_gfortran
BUILD_MODS_DIR=$BUILD_DIR/code

prepare_size_h() {
    if [ -n "${SIZE_SNX:-}" ] && [ -n "${SIZE_SNY:-}" ] && [ -n "${SIZE_NPX:-}" ] && [ -n "${SIZE_NPY:-}" ]; then
        python3 - "$BUILD_MODS_DIR/SIZE.h" "$SIZE_SNX" "$SIZE_SNY" "$SIZE_NPX" "$SIZE_NPY" <<'PY2'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
snx, sny, npx, npy = sys.argv[2:6]
text = path.read_text()
patterns = {
    r'^(\s*PARAMETER \( sNx =\s*)\d+(\s*\).*)$': snx,
    r'^(\s*PARAMETER \( sNy =\s*)\d+(\s*\).*)$': sny,
    r'^(\s*PARAMETER \( nPx =\s*)\d+(\s*\).*)$': npx,
    r'^(\s*PARAMETER \( nPy =\s*)\d+(\s*\).*)$': npy,
}
for pattern, value in patterns.items():
    text, count = re.subn(pattern, rf'\g<1>{value}\g<2>', text, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f'Failed to patch SIZE.h with pattern: {pattern}')
path.write_text(text)
PY2
    fi
}

stage_case_inputs() {
    find "$CASE_DIR/input" -maxdepth 1 -type f -exec cp {} "$RUN_DIR/" \;

    case "$CASE_DIR" in
        */vortexSphere_Williamson_TC1_prime)
            sed -i "s/^ALPHA_RAD = .*/ALPHA_RAD = ${ALPHA_VALUE}/" "$RUN_DIR/gendata_ref.py"
            (cd "$RUN_DIR" && python3 gendata_ref.py)
            ;;
        */vortexSphere_Williamson_TC1)
            sed -i "s/^ALPHA_RAD = .*/ALPHA_RAD = ${ALPHA_VALUE}/" "$RUN_DIR/gendata_ref.py"
            (cd "$RUN_DIR" && python3 gendata_ref.py)
            ;;
        */vortexSphere_Williamson_TC2)
            case "$CASE_NAME" in
                c1) suffix="" ;;
                c2) suffix="_c2" ;;
                c3) suffix="_c3" ;;
                c4) suffix="_c4" ;;
                *)
                    echo "ERROR: unsupported CASE_NAME=$CASE_NAME" >&2
                    exit 1
                    ;;
            esac
            if [ -n "$suffix" ]; then
                cp "$RUN_DIR/eta_init${suffix}.bin" "$RUN_DIR/eta_init.bin"
                cp "$RUN_DIR/u_init${suffix}.bin" "$RUN_DIR/u_init.bin"
                cp "$RUN_DIR/v_init${suffix}.bin" "$RUN_DIR/v_init.bin"
            fi
            ;;
        *)
            echo "ERROR: unsupported CASE_DIR=$CASE_DIR" >&2
            exit 1
            ;;
    esac
}

mkdir -p "$BUILD_DIR" "$RUN_DIR"

# -----------------------------
# Compile
# -----------------------------
cd "$BUILD_DIR" || exit 1

rm -f Makefile
make CLEAN >/dev/null 2>&1 || true
rm -rf "$BUILD_MODS_DIR"
mkdir -p "$BUILD_MODS_DIR"
find "$CASE_DIR/code" -maxdepth 1 -type f -exec cp {} "$BUILD_MODS_DIR/" \;
prepare_size_h

"$TOOLS_DIR/genmake2"     -mpi     -mods="$BUILD_MODS_DIR"     -of="$BUILD_OPTIONS"

make depend
make -j "${MAKE_JOBS:-8}"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created"
    exit 1
fi

# -----------------------------
# Prepare run directory
# -----------------------------
find "$RUN_DIR" -maxdepth 1 -type f     \( -name "*.data" -o -name "*.meta" -o -name "*.bin" -o -name "*.py"     -o -name "STDOUT*" -o -name "STDERR*" -o -name "dyn*" -o -name "pickup*"     -o -name "slurm-*" -o -name "available_diagnostics.log" -o -name "run.log"     -o -name "mitgcmuv" -o -name "data" -o -name "data.pkg"     -o -name "data.diagnostics" -o -name "eedata" -o -name "eedata.mth" \)     -delete

cp mitgcmuv "$RUN_DIR/"
stage_case_inputs

cd "$RUN_DIR" || exit 1

echo "======================================"
echo "Running MITgcm"
echo "MPI tasks = ${MPI_TASKS:-$SLURM_NTASKS}"
echo "======================================"

mpirun -np "${MPI_TASKS:-$SLURM_NTASKS}" ./mitgcmuv

STATUS=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME-START_TIME))

echo "======================================"
echo "End: $(date)"
echo "Status: $STATUS"
echo "Elapsed: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
echo "======================================"

exit $STATUS
