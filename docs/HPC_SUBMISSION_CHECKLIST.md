# Williamson HPC Submission Runbook

Audit date: 2026-06-16, Asia/Beirut.

Scope: TC1 through TC7 under `Sandbox/vortexSphere_Williamson_TC*`. This is an audit and submission runbook. It does not modify experiment setup, job scripts, or plotting scripts.

## Decision Summary

| Case | Tonight action | Current evidence | Final status |
| --- | --- | --- | --- |
| TC1 | Submit/rerun after optional output cleanup | Current runs end normally; TC1 alpha `1.57` has stale extra iterations/assets | Submit-ready, cleanup recommended |
| TC2 | Submit only after archiving stale run dirs | Input templates are fixed, but current run dirs do not contain `selectCoriMap=3` or `fCori*.bin` | Needs restage, then submit |
| TC3 | Submit only after archiving stale run dirs | Input templates are fixed, but current run dirs do not contain `selectCoriMap=3` or `fCori*.bin` | Needs restage, then submit |
| TC4 | Do not submit | Job has `TC4_FORCE_READY` guard; forcing hook is not validated | Blocked/scaffold |
| TC5 | Debug-only, not final | MITgcm ended normally, but all saved `Eta/U/V/dynDiag/dyn_Aux` records after iter `0000000000` are non-finite | Invalid output |
| TC6 | Submit/rerun | Current run ends normally and final fields are finite through iter `0000040320` | Submit-ready |
| TC7 | Do not submit | Missing `input/raw/tc7_initial_conditions.npz` | Blocked/missing data |

The safest final-submit set tonight is TC1, TC2, TC3, and TC6, with the TC2/TC3 run directories archived before submission so the corrected templates are actually copied into new run directories. TC5 should not be submitted as a final science run until the non-finite output cause is fixed.

## Expected Output Matrix

| Case | Expected output | Initialization | Key flags/files | Validation status | Correction/action |
| --- | --- | --- | --- | --- | --- |
| TC1 | Passive tracer/flow should complete 12 days and preserve tracer mass; rotated alphas should remain finite | Flat 4000 m bathymetry, analytic solid-body velocity/tracer fields from `input/gendata_ref.py`, alpha patched per job | `implicitFreeSurface=.TRUE.`, `exactConserv=.TRUE.`, no explicit viscosity/diffusion | Current runs complete; conservation summary reports mass/tracer preserved to roundoff or well preserved | Archive stale TC1 outputs if a clean publication rerun is desired |
| TC2 | Steady-state analytic `eta/u/v` should remain steady; error norms should be small and not grow by alpha | Flat 4000 m bathymetry, analytic rotated steady state; alpha patched per job | Must use `selectCoriMap=3` plus `fCoriC.bin`, `fCorCs.bin`, `fCoriG.bin` | Templates fixed, current run dirs stale and invalid for rotated-alpha validation | Archive `run_alpha_*`, submit fresh, regenerate plots/errors/conservation |
| TC3 | Geostrophic adjustment/rotated benchmark should remain finite and match expected alpha behavior after corrected Coriolis | Flat TC3 bathymetry, analytic `eta/u/v`; alpha patched per job | Must use `selectCoriMap=3` plus `fCoriC.bin`, `fCorCs.bin`, `fCoriG.bin` | Templates fixed, current run dirs stale and invalid for rotated-alpha validation | Archive `run_alpha_*`, submit fresh, regenerate plots/errors/conservation |
| TC4 | Forced flow should balance analytic forcing terms and remain steady over 5 days | TC2-like initial state is present, but the required TC4 forcing terms are not implemented/validated | `TC4_FORCE_READY` job guard; `code/README_TC4_FORCING.md` documents scaffold state | Blocked | Implement/validate forcing hook before setting `TC4_FORCE_READY=1` |
| TC5 | Mountain-flow response should remain finite through 15 days with interpretable day 5/10/15 fields | Mountain bathymetry from `wf.tc5_bathymetry`, analytic initial zonal flow/free surface | `implicitFreeSurface=.TRUE.`, `exactConserv=.TRUE.`, no explicit viscosity/diffusion | Invalid: current saved fields are finite at iter `0000000000`, then non-finite from iter `0000001440` onward | Debug bathymetry/sign/free-surface/wet-cell setup before final submission |
| TC6 | Rossby-Haurwitz wave should remain finite through 14 days with finite final fields and diagnostics | Flat bathymetry, analytic TC6 wave `eta/u/v` | `deltaT=30`, `implicitFreeSurface=.TRUE.`, `exactConserv=.TRUE.` | Submit-ready; fields finite through iter `0000040320` | Rerun/regenerate products if final clean assets are needed |
| TC7 | Observed-analysis case should initialize from gridded geopotential height and winds | Requires `eta_m`, `u_m_s`, `v_m_s` arrays in `input/raw/tc7_initial_conditions.npz` | Validator `Sandbox/Scripts/tc7_validate_initial_conditions.py` runs before compile | Blocked | Download/prepare NPZ, run validator, then submit |

## Critical Findings

TC2 and TC3 are not safe to resubmit into the current run directories. The input templates now contain the corrected rotated-Coriolis setup, but the current run directories are stale:

```bash
grep -n "selectCoriMap" Sandbox/vortexSphere_Williamson_TC2/run_alpha_*/data
grep -n "selectCoriMap" Sandbox/vortexSphere_Williamson_TC3/run_alpha_*/data
find Sandbox/vortexSphere_Williamson_TC2/run_alpha_* Sandbox/vortexSphere_Williamson_TC3/run_alpha_* -name 'fCori*.bin'
```

At audit time these checks found no `selectCoriMap`, no rotated-Coriolis writer, and no `fCori*.bin` files in the current TC2/TC3 run dirs. The job scripts call `stage_case_inputs`, but if stale run output remains, old fields can remain mixed into the run directory and the existing output is not valid validation evidence.

MITgcm itself expects the corrected filenames when `selectCoriMap=3`:

```bash
rg -n "selectCoriMap|fCoriC.bin|fCorCs.bin|fCoriG.bin" model/src/ini_cori.F model/inc/PARAMS.h
```

TC5 is not a plotting-only problem. A metadata-aware finite scan of the MDS files showed:

```text
TC5 Eta/U/V/dynDiag/dyn_Aux:
  iter 0000000000: all entries finite
  iter 0000001440 and later: 0 finite entries
```

This means the model output becomes non-finite after the initial record. The plotting scripts are correctly exposing a bad output state.

The CFL audit also flags the old TC2 alpha `0.05` run as above the conservative 0.5 observed-advective margin:

```text
TC2 alpha=0.05: initial advective CFL 0.332, observed old-run advective CFL 0.560
```

Because the current TC2 output is stale and not a valid rotated-Coriolis run, this is a watch item rather than a final blocker. If the fresh TC2 alpha `0.05` rerun still exceeds that margin, lower `DELTA_T` from `10.0` s to about `8.9` s or smaller before using it as final evidence.

## Preflight Commands

Run from the login node:

```bash
cd /home/hmh85/scratch/MITgcm
pwd -P
squeue -u "$USER"
sinfo -p large -N -o "%N %T %C %m" | grep -E "NODELIST|onode13|onode14|onode15|onode16"
```

Compile/check Python entrypoints:

```bash
.venv/bin/python -m py_compile \
  Sandbox/Scripts/TC1.py \
  Sandbox/Scripts/TC2.py \
  Sandbox/Scripts/TC3.py \
  Sandbox/Scripts/TC5.py \
  Sandbox/Scripts/TC6.py \
  Sandbox/Scripts/conservation_diagnostics.py \
  Sandbox/Scripts/cfl_deltaT_audit.py \
  Sandbox/Scripts/tc7_validate_initial_conditions.py
```

CFL and input checks:

```bash
.venv/bin/python Sandbox/vortexSphere_Williamson_TC1/tools/check_cfl_tc1.py --all
.venv/bin/python Sandbox/vortexSphere_Williamson_TC2/tools/check_cfl_tc2.py --all
.venv/bin/python Sandbox/vortexSphere_Williamson_TC3/tools/check_cfl_tc3.py --all
.venv/bin/python Sandbox/vortexSphere_Williamson_TC4/tools/check_cfl_tc4.py --all
.venv/bin/python Sandbox/vortexSphere_Williamson_TC5/tools/check_cfl_tc5.py --all
.venv/bin/python Sandbox/vortexSphere_Williamson_TC6/tools/check_cfl_tc6.py --all
.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments
```

Expected CFL result: TC2, TC3, TC4, TC5, and TC6 should pass their initial-condition checks; TC7 should fail until the raw NPZ exists. TC1's standalone checker uses a template `deltaT=1` for the high-alpha cases and can warn there, while the high-alpha TC1 jobs use `DELTA_T=0.75`.

Blocker checks:

```bash
grep -n "TC4_FORCE_READY" Sandbox/vortexSphere_Williamson_TC4/jobs/large/job_tc4_c1.slurm
test -f Sandbox/vortexSphere_Williamson_TC7/input/raw/tc7_initial_conditions.npz
.venv/bin/python Sandbox/Scripts/tc7_validate_initial_conditions.py
```

The TC7 validator is expected to fail until the NPZ is provided.

## Archive Before Fresh TC2/TC3 Submission

Use archive moves, not deletion, so the current invalid/stale outputs remain available for comparison.

```bash
cd /home/hmh85/scratch/MITgcm
STAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p "Sandbox/run_archive/$STAMP" \
         "Sandbox/output/archive/$STAMP" \
         "docs/assets/williamson/archive/$STAMP" \
         "docs/fragments/archive/$STAMP"

for d in \
  Sandbox/vortexSphere_Williamson_TC2/run_alpha_0 \
  Sandbox/vortexSphere_Williamson_TC2/run_alpha_0.05 \
  Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.52 \
  Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.57 \
  Sandbox/vortexSphere_Williamson_TC3/run_alpha_0 \
  Sandbox/vortexSphere_Williamson_TC3/run_alpha_1.0472
do
  [ -d "$d" ] && mv "$d" "Sandbox/run_archive/$STAMP/"
done

for d in Sandbox/output/TestCase2 Sandbox/output/TestCase3; do
  [ -e "$d" ] && mv "$d" "Sandbox/output/archive/$STAMP/"
done

for d in docs/assets/williamson/TestCase2 docs/assets/williamson/TestCase3; do
  [ -e "$d" ] && mv "$d" "docs/assets/williamson/archive/$STAMP/"
done

for f in docs/fragments/testcase2.html docs/fragments/testcase3.html; do
  [ -e "$f" ] && mv "$f" "docs/fragments/archive/$STAMP/"
done
```

Optional clean TC1 publication archive:

```bash
STAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p "Sandbox/run_archive/$STAMP" \
         "Sandbox/output/archive/$STAMP" \
         "docs/assets/williamson/archive/$STAMP" \
         "docs/fragments/archive/$STAMP"

for d in Sandbox/vortexSphere_Williamson_TC1/run_alpha_0 \
         Sandbox/vortexSphere_Williamson_TC1/run_alpha_0.05 \
         Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.52 \
         Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.57
do
  [ -d "$d" ] && mv "$d" "Sandbox/run_archive/$STAMP/"
done

[ -e Sandbox/output/TestCase1 ] && mv Sandbox/output/TestCase1 "Sandbox/output/archive/$STAMP/"
[ -e docs/assets/williamson/TestCase1 ] && mv docs/assets/williamson/TestCase1 "docs/assets/williamson/archive/$STAMP/"
[ -e docs/fragments/testcase1.html ] && mv docs/fragments/testcase1.html "docs/fragments/archive/$STAMP/"
```

Do not use `Sandbox/vortexSphere_Williamson_TC{2,3,4,5,6,7}/jobs/clean_run.sh` or `clean_build.sh` yet; those wrappers delegate to repo-level `jobs/clean_run.sh` and `jobs/clean_build.sh`, which are not present in this checkout.

## Slurm Dry Run

```bash
sbatch --test-only Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c1.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c2.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c3.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c4.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c1.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c2.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c3.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c4.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC3/jobs/large/job_tc3_c1.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC3/jobs/large/job_tc3_c2.slurm
sbatch --test-only Sandbox/vortexSphere_Williamson_TC6/jobs/large/job_tc6_c1.slurm
```

Run TC5 with `--test-only` only if you are preparing a debug rerun:

```bash
sbatch --test-only Sandbox/vortexSphere_Williamson_TC5/jobs/large/job_tc5_c1.slurm
```

## Final Submit Commands

After archiving stale TC2/TC3 run dirs, submit:

```bash
cd /home/hmh85/scratch/MITgcm
sbatch Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c1.slurm
sbatch Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c2.slurm
sbatch Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c3.slurm
sbatch Sandbox/vortexSphere_Williamson_TC1/jobs/large/job_tc1_c4.slurm
sbatch Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c1.slurm
sbatch Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c2.slurm
sbatch Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c3.slurm
sbatch Sandbox/vortexSphere_Williamson_TC2/jobs/large/job_tc2_c4.slurm
sbatch Sandbox/vortexSphere_Williamson_TC3/jobs/large/job_tc3_c1.slurm
sbatch Sandbox/vortexSphere_Williamson_TC3/jobs/large/job_tc3_c2.slurm
sbatch Sandbox/vortexSphere_Williamson_TC6/jobs/large/job_tc6_c1.slurm
```

Wrapper alternative, only from each case's `jobs` directory:

```bash
cd /home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC1/jobs && bash submission_p_large.sh
cd /home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC2/jobs && bash submission_p_large.sh
cd /home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC3/jobs && bash submission_p_large.sh
cd /home/hmh85/scratch/MITgcm/Sandbox/vortexSphere_Williamson_TC6/jobs && bash submission_p_large.sh
```

## No-submit List

```bash
# TC4: blocked/scaffold. The forcing hook is not implemented/validated.
# sbatch Sandbox/vortexSphere_Williamson_TC4/jobs/large/job_tc4_c1.slurm

# TC5: debug-only. Current final output is non-finite after day 0.
# sbatch Sandbox/vortexSphere_Williamson_TC5/jobs/large/job_tc5_c1.slurm

# TC7: blocked. Missing input/raw/tc7_initial_conditions.npz.
# sbatch Sandbox/vortexSphere_Williamson_TC7/jobs/large/job_tc7_c1.slurm
```

## Monitoring

```bash
squeue -u "$USER"
tail -f Sandbox/vortexSphere_Williamson_TC1/logs/slurm-*.out
tail -f Sandbox/vortexSphere_Williamson_TC2/logs/slurm-*.out
tail -f Sandbox/vortexSphere_Williamson_TC3/logs/slurm-*.out
tail -f Sandbox/vortexSphere_Williamson_TC6/logs/slurm-*.out
```

Normal MITgcm completion check:

```bash
grep -R "PROGRAM MAIN: Execution ended Normally" \
  Sandbox/vortexSphere_Williamson_TC1/run_alpha_*/STDOUT.0000 \
  Sandbox/vortexSphere_Williamson_TC2/run_alpha_*/STDOUT.0000 \
  Sandbox/vortexSphere_Williamson_TC3/run_alpha_*/STDOUT.0000 \
  Sandbox/vortexSphere_Williamson_TC6/run_standard/STDOUT.0000
```

Expected final state files:

```bash
test -s Sandbox/vortexSphere_Williamson_TC1/run_alpha_0/Eta.0000017280.data
test -s Sandbox/vortexSphere_Williamson_TC1/run_alpha_0.05/Eta.0000103680.data
test -s Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.52/Eta.0001382400.data
test -s Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.57/Eta.0001382400.data
test -s Sandbox/vortexSphere_Williamson_TC2/run_alpha_0/Eta.0000017280.data
test -s Sandbox/vortexSphere_Williamson_TC2/run_alpha_0.05/Eta.0000103680.data
test -s Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.52/Eta.0001382400.data
test -s Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.57/Eta.0001382400.data
test -s Sandbox/vortexSphere_Williamson_TC3/run_alpha_0/Eta.0000007200.data
test -s Sandbox/vortexSphere_Williamson_TC3/run_alpha_1.0472/Eta.0000576000.data
test -s Sandbox/vortexSphere_Williamson_TC6/run_standard/Eta.0000040320.data
```

After fresh TC2/TC3 jobs start or finish, confirm that the run directories were restaged with the corrected Coriolis files:

```bash
for d in Sandbox/vortexSphere_Williamson_TC2/run_alpha_* Sandbox/vortexSphere_Williamson_TC3/run_alpha_*; do
  echo "== $d =="
  grep -n "selectCoriMap" "$d/data"
  ls "$d"/fCoriC.bin "$d"/fCorCs.bin "$d"/fCoriG.bin
done
```

## Post-run Regeneration

Each Slurm job runs its own case plotting script at the end. If you need to regenerate after all jobs complete:

```bash
cd /home/hmh85/scratch/MITgcm
TC1_RUN_DIRS=Sandbox/vortexSphere_Williamson_TC1/run_alpha_0:Sandbox/vortexSphere_Williamson_TC1/run_alpha_0.05:Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.52:Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.57 .venv/bin/python Sandbox/Scripts/TC1.py
TC2_RUN_DIRS=Sandbox/vortexSphere_Williamson_TC2/run_alpha_0:Sandbox/vortexSphere_Williamson_TC2/run_alpha_0.05:Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.52:Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.57 .venv/bin/python Sandbox/Scripts/TC2.py
TC3_RUN_DIRS=Sandbox/vortexSphere_Williamson_TC3/run_alpha_0:Sandbox/vortexSphere_Williamson_TC3/run_alpha_1.0472 .venv/bin/python Sandbox/Scripts/TC3.py
TC6_RUN_DIRS=Sandbox/vortexSphere_Williamson_TC6/run_standard .venv/bin/python Sandbox/Scripts/TC6.py
.venv/bin/python Sandbox/Scripts/conservation_diagnostics.py --cases TC1,TC2,TC3,TC6
.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments
.venv/bin/python Sandbox/Scripts/github_pages.py
```

Do not regenerate TC5 final products into the publication site until its output is finite.

## Finite-field Smoke Test

Use this after any rerun before trusting plots:

```bash
.venv/bin/python - <<'PY'
from pathlib import Path
import re
import numpy as np

paths = [
    "Sandbox/vortexSphere_Williamson_TC1/run_alpha_0/Eta.0000017280.data",
    "Sandbox/vortexSphere_Williamson_TC1/run_alpha_0.05/Eta.0000103680.data",
    "Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.52/Eta.0001382400.data",
    "Sandbox/vortexSphere_Williamson_TC1/run_alpha_1.57/Eta.0001382400.data",
    "Sandbox/vortexSphere_Williamson_TC2/run_alpha_0/Eta.0000017280.data",
    "Sandbox/vortexSphere_Williamson_TC2/run_alpha_0.05/Eta.0000103680.data",
    "Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.52/Eta.0001382400.data",
    "Sandbox/vortexSphere_Williamson_TC2/run_alpha_1.57/Eta.0001382400.data",
    "Sandbox/vortexSphere_Williamson_TC3/run_alpha_0/Eta.0000007200.data",
    "Sandbox/vortexSphere_Williamson_TC3/run_alpha_1.0472/Eta.0000576000.data",
    "Sandbox/vortexSphere_Williamson_TC6/run_standard/Eta.0000040320.data",
]

def dtype_for(datafile):
    meta = datafile.with_suffix(".meta")
    text = meta.read_text(errors="ignore")
    match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'", text)
    prec = match.group(1).lower() if match else "float64"
    return ">f8" if "64" in prec else ">f4"

for name in paths:
    path = Path(name)
    arr = np.fromfile(path, dtype=dtype_for(path))
    finite = np.isfinite(arr).sum()
    print(f"{name}: finite {finite}/{arr.size}")
    if finite != arr.size:
        raise SystemExit(f"non-finite values found in {name}")
PY
```

## TC7 Required Data

Before TC7 can run, create:

```text
Sandbox/vortexSphere_Williamson_TC7/input/raw/tc7_initial_conditions.npz
```

Required arrays:

- `eta_m`, shape `(720, 1440)`.
- `u_m_s`, shape `(720, 1440)`.
- `v_m_s`, shape `(720, 1440)`.
- Optional `bathymetry_m`, shape `(720, 1440)`.

Then run:

```bash
.venv/bin/python Sandbox/Scripts/tc7_validate_initial_conditions.py
.venv/bin/python Sandbox/vortexSphere_Williamson_TC7/tools/check_cfl_tc7.py --all
```
