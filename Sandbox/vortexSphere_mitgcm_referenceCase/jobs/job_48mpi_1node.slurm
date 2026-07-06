#!/bin/bash

#SBATCH --job-name=mitgcm_48mpi
#SBATCH --account=hmh85

#SBATCH --partition=large
#SBATCH --nodelist=onode15
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6

set -euo pipefail

echo "======================================"
echo "MITgcm compile + run"
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes:  $SLURM_JOB_NODELIST"
echo "Start:  $(date)"
echo "======================================"

START_TIME=$(date +%s)
RANKS=${SLURM_NTASKS:?SLURM_NTASKS is not set}

case "$RANKS" in
    8)  sNx=360; sNy=360; nPx=4; nPy=2 ;;
    12) sNx=360; sNy=240; nPx=4; nPy=3 ;;
    16) sNx=360; sNy=180; nPx=4; nPy=4 ;;
    24) sNx=240; sNy=180; nPx=6; nPy=4 ;;
    32) sNx=180; sNy=180; nPx=8; nPy=4 ;;
    48) sNx=180; sNy=120; nPx=8; nPy=6 ;;
    *)
        echo "ERROR: unsupported one-node MPI rank count: $RANKS" >&2
        exit 1
        ;;
esac

CASE_DIR=/scratch/p8190783-hmh85/MITgcm/verification/vortexSphere_mitgcm_referenceCase
MITGCM_ROOT=$(cd "$CASE_DIR/../.." && pwd)
MODS_DIR=$CASE_DIR/code_${RANKS}mpi
BUILD_DIR=$CASE_DIR/build_${RANKS}mpi
RUN_DIR=$CASE_DIR/run_${RANKS}mpi_1node

echo "Grid layout: sNx=$sNx sNy=$sNy nPx=$nPx nPy=$nPy"
echo "Mods dir:    $MODS_DIR"
echo "Build dir:   $BUILD_DIR"
echo "Run dir:     $RUN_DIR"

# -----------------------------
# Generate rank-specific code dir
# -----------------------------
mkdir -p "$MODS_DIR"
find "$MODS_DIR" -maxdepth 1 -type f -delete
cp "$CASE_DIR"/code/* "$MODS_DIR"/

cat > "$MODS_DIR/SIZE.h" <<EOF
C    !===  GRID SPACING & DOMAIN SIZE =====================================
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr

      PARAMETER ( sNx =  $sNx )        ! grid points in x-dir per tile
      PARAMETER ( sNy =  $sNy )        ! grid points in y-dir per tile
      PARAMETER ( OLx =    2 )        ! overlap in x-dir
      PARAMETER ( OLy =    2 )        ! overlap in y-dir
      PARAMETER ( nSx =    1 )
      PARAMETER ( nSy =    1 )
      PARAMETER ( nPx =    $nPx )
      PARAMETER ( nPy =    $nPy )
      PARAMETER ( Nx  = sNx*nSx*nPx )
      PARAMETER ( Ny  = sNy*nSy*nPy )
      PARAMETER ( Nr  =    1 )        ! vertical levels

      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )
EOF

# -----------------------------
# Compile
# -----------------------------
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

rm -f Makefile
make CLEAN >/dev/null 2>&1 || true

"$MITGCM_ROOT/tools/genmake2" \
    -mpi \
    -mods="$MODS_DIR" \
    -of="$MITGCM_ROOT/tools/build_options/linux_amd64_gfortran"

make depend
make -j "${MAKE_JOBS:-8}"

if [ ! -x mitgcmuv ]; then
    echo "ERROR: mitgcmuv not created"
    exit 1
fi

# -----------------------------
# Prepare run directory
# -----------------------------
mkdir -p "$RUN_DIR"

find "$RUN_DIR" -maxdepth 1 -type f \
    \( -name "*.data" -o -name "*.meta" -o -name "STDOUT*" -o -name "STDERR*" -o -name "pickup*" -o -name "monitor*" \) \
    -delete

cp mitgcmuv "$RUN_DIR"/
cp "$CASE_DIR"/input/* "$RUN_DIR"/

cat > "$RUN_DIR/mpi_layout.txt" <<EOF
ranks=$RANKS
sNx=$sNx
sNy=$sNy
nPx=$nPx
nPy=$nPy
Nx=1440
Ny=720
Nr=1
EOF

cd "$RUN_DIR" || exit 1

echo "======================================"
echo "Running MITgcm"
echo "MPI tasks = $SLURM_NTASKS"
echo "======================================"

set +e
mpirun -np "$SLURM_NTASKS" ./mitgcmuv
STATUS=$?
set -e

END_TIME=$(date +%s)
ELAPSED=$((END_TIME-START_TIME))

echo "======================================"
echo "End: $(date)"
echo "Status: $STATUS"
echo "Elapsed: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
echo "======================================"

exit $STATUS
