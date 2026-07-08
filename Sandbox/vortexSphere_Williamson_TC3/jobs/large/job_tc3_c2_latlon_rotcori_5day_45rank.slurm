#!/bin/bash
#SBATCH --job-name=TC3_c2_rotcor_5d_45r
#SBATCH --account=hmh85
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=45
#SBATCH --ntasks-per-node=45
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC3/logs/slurm-%j.out
#SBATCH --error=/home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC3/logs/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hmh85@aub.edu.lb

set -euo pipefail
module purge
module load gcc/10.1.0
module load mpi/openmpi/4.1.4-slurm-18.08.6
module load python/3.12
source /home/hmh85/scratch/MITgcm/.venv/bin/activate

export MITGCM_DIR=/home/hmh85/scratch/MITgcm
export CASE_DIR=$MITGCM_DIR/Sandbox/vortexSphere_Williamson_TC3
export CASE_NAME=c2
export ALPHA_VALUE=1.0471975511965976
export DELTA_T=0.5
export TOTAL_SECONDS=432000
export DUMP_SECONDS=86400
export BUILD_DIR=$CASE_DIR/build/c2_latlon_rotcori_5day_large45
export RUN_DIR=$CASE_DIR/run_alpha_1.0472_latlon_rotcori_5day_45rank
export LOG_DIR=$CASE_DIR/logs
export MPI_TASKS=${SLURM_NTASKS:-45}
export MAKE_JOBS=${SLURM_CPUS_ON_NODE:-45}
export N_TIME_STEPS=$(python3 -c "print(int(round(${TOTAL_SECONDS}/${DELTA_T})))")

START_TIME=$(date +%s)
TOOLS_DIR=$MITGCM_DIR/tools
BUILD_OPTIONS=$TOOLS_DIR/build_options/linux_amd64_gfortran
MODS_DIR=$CASE_DIR/code_smoke45
PLOT_SCRIPT=$MITGCM_DIR/Sandbox/Scripts/TC3.py
SETUP_CHECK=$CASE_DIR/tools/check_setup_tc3.py

finish() {
    local status=$1 end_time elapsed
    end_time=$(date +%s)
    elapsed=$((end_time-START_TIME))
    echo "======================================"
    echo "End: $(date)"
    echo "Status: $status"
    echo "Elapsed: $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"
    echo "======================================"
}

set_or_append_data_value() {
    local key=$1
    local value=$2
    local file=$3
    if grep -qi "^[[:space:]]*${key}[[:space:]]*=" "$file"; then
        sed -i "s/^[[:space:]]*${key}[[:space:]]*=.*/ ${key}=${value},/" "$file"
    else
        sed -i "/^[[:space:]]*useCoriolis[[:space:]]*=/a\\ ${key}=${value}," "$file"
    fi
}

stage_case_inputs() {
    find "$CASE_DIR/input" -maxdepth 1 -type f \
        ! -name "data.exch2" \
        ! -name "eta_init*.bin" \
        ! -name "u_init*.bin" \
        ! -name "v_init*.bin" \
        ! -name "bathymetry*.bin" \
        ! -name "fCoriC.bin" \
        ! -name "fCoriG.bin" \
        ! -name "fCorCs.bin" \
        -exec cp {} "$RUN_DIR/" \;
    touch "$RUN_DIR/.use_latlon_grid"
    sed -i "s/^ALPHA_RAD = .*/ALPHA_RAD = ${ALPHA_VALUE}/" "$RUN_DIR/gendata_ref.py"
    sed -i "s/^ nTimeSteps=.*/ nTimeSteps=${N_TIME_STEPS},/" "$RUN_DIR/data"
    sed -i "s/^ deltaT=.*/ deltaT=${DELTA_T},/" "$RUN_DIR/data"
    sed -i "s/^ dumpFreq=.*/ dumpFreq=${DUMP_SECONDS}.,/" "$RUN_DIR/data"
    sed -i "s/^ monitorFreq=.*/ monitorFreq=${DUMP_SECONDS}.,/" "$RUN_DIR/data"
    sed -i "s/^  frequency(3) = -.*/  frequency(3) = -${DUMP_SECONDS}.,/" "$RUN_DIR/data.diagnostics"
    sed -i "s/^  frequency(2) = -.*/  frequency(2) = -${DUMP_SECONDS}.,/" "$RUN_DIR/data.diagnostics"
    sed -i "s/^[[:space:]]*usingCurvilinearGrid[[:space:]]*=.*/ usingSphericalPolarGrid = .TRUE.,/" "$RUN_DIR/data"
    sed -i "/^[[:space:]]*radius_fromHorizGrid[[:space:]]*=/d" "$RUN_DIR/data"
    sed -i "/^[[:space:]]*usingSphericalPolarGrid[[:space:]]*=/a\\ delX = 1440*0.25,\\n delY = 720*0.25,\\n xgOrigin = 0.,\\n ygOrigin = -90.," "$RUN_DIR/data"
    set_or_append_data_value "selectCoriMap" "3" "$RUN_DIR/data"
    cat > "$RUN_DIR/eedata" <<'EOF'
# Example "eedata" file
 &EEPARMS
 &
EOF
    rm -f "$RUN_DIR/data.exch2"
    (cd "$RUN_DIR" && python3 gendata_ref.py)
    python3 "$SETUP_CHECK" --run-dir "$RUN_DIR" --alpha "$ALPHA_VALUE"
    grep -q "usingSphericalPolarGrid = .TRUE." "$RUN_DIR/data"
    grep -qi "selectCoriMap[[:space:]]*=[[:space:]]*3" "$RUN_DIR/data"
    if grep -q "usingCurvilinearGrid" "$RUN_DIR/data"; then
        echo "ERROR: lat-lon run contains usingCurvilinearGrid" >&2
        exit 1
    fi
    if grep -q "useCubedSphereExchange" "$RUN_DIR/eedata"; then
        echo "ERROR: lat-lon run contains useCubedSphereExchange" >&2
        exit 1
    fi
    if [ -e "$RUN_DIR/data.exch2" ]; then
        echo "ERROR: lat-lon run contains data.exch2" >&2
        exit 1
    fi
}

clean_run_dir() {
    find "$RUN_DIR" -maxdepth 1 -type f \
        \( -name "*.data" -o -name "*.meta" -o -name "*.bin" \
        -o -name "STDOUT*" -o -name "STDERR*" -o -name "dyn*" \
        -o -name "tile*.mitgrid" -o -name "fCoriG.face*.bin" \
        -o -name "pickup*" -o -name "available_diagnostics.log" \
        -o -name "run.log" -o -name "run_inventory.log" \
        -o -name "mitgcmuv" -o -name "data" -o -name "data.pkg" \
        -o -name "data.diagnostics" -o -name "data.exch2" \
        -o -name "eedata" -o -name "eedata.mth" \
        -o -name "gendata_ref.py" -o -name ".use_latlon_grid" \) \
        -delete
}

mkdir -p "$LOG_DIR" "$BUILD_DIR" "$RUN_DIR"
cd "$BUILD_DIR" || exit 1
"$TOOLS_DIR/genmake2" -mpi -mods="$MODS_DIR" -of="$BUILD_OPTIONS"
make depend
make -j "${MAKE_JOBS:-8}"
clean_run_dir
cp mitgcmuv "$RUN_DIR/"
stage_case_inputs
cd "$RUN_DIR" || exit 1

set +e
mpirun -np "${MPI_TASKS:-$SLURM_NTASKS}" ./mitgcmuv
STATUS=$?
set -e
if [ "$STATUS" -ne 0 ]; then
    finish "$STATUS"
    exit "$STATUS"
fi

set +e
TC3_RUN_DIRS="$RUN_DIR" python3 "$PLOT_SCRIPT"
STATUS=$?
set -e
finish "$STATUS"
exit "$STATUS"
