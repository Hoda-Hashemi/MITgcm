## Time-step safety check (CFL and deltaT)

This check asks whether the submitted model time step is small enough for the resolved velocities and grid spacing in each Williamson run.

For each submitted Williamson job, the audited Courant numbers are `C_adv = |u| deltaT / dx + |v| deltaT / dy` and `C_g = sqrt(gH) deltaT / dx`. The spherical-polar grid is 1440 x 720 with 0.25 degree spacing; the smallest center-cell zonal metric is 60.65 m in the polar row and `dy = 27798.73 m`.

The gravity-wave column is an explicit-wave diagnostic only. All active runs set `implicitFreeSurface=.TRUE.`, so the external gravity wave is handled by the implicit free-surface solve rather than by the explicit advective CFL limit.

The template `input/data` files define the grid and default schedule, but submitted jobs are controlled by `jobs/large/job_*.slurm`: each job exports `DELTA_T`, computes `nTimeSteps=round(TOTAL_SECONDS/DELTA_T)`, copies the template input directory, and rewrites the run-local `data`. The table therefore separates template, submitted job, and archived run deltaT.

TC1 has the most visible discrepancy: its legacy `tools/check_cfl_tc1.py` prints the template `deltaT=1 s`, while the submitted jobs use 60 s, 10 s, or 0.75 s and the archived run-local `data` files confirm those values. Use this audit table for submitted-run CFL decisions.

| case | run | alpha(rad) | template dt(s) | job dt(s) | run dt(s) | steps | H(m) | sqrt(gH) | init adv CFL | saved adv CFL | max speed | explicit Cg,x | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TC1 | alpha=0 | 0 | 1.00 | 60.0 | 60.0 | 17280 | 4000 | 198 | 0.083 | n/a | 38.6 | 196 | no advective deltaT change indicated |
| TC1 | alpha=0.05 | 0.05 | 1.00 | 10.0 | 10.0 | 103680 | 4000 | 198 | 0.332 | n/a | 38.6 | 32.7 | no advective deltaT change indicated |
| TC1 | alpha=1.52 | 1.5208 | 1.00 | 0.75 | 0.75 | 1382400 | 4000 | 198 | 0.490 | n/a | 1095 | 2 | no advective deltaT change indicated |
| TC1 | alpha=1.57 | 1.5708 | 1.00 | 0.75 | 0.75 | 1382400 | 4000 | 198 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC2 | alpha=0 | 0 | 60.0 | 60.0 | 60.0 | 17280 | 2997 | 171 | 0.083 | n/a | 38.6 | 170 | no advective deltaT change indicated |
| TC2 | alpha=0.05 | 0.05 | 60.0 | 10.0 | 10.0 | 103680 | 2997 | 171 | 0.332 | n/a | 38.6 | 28.3 | no advective deltaT change indicated |
| TC2 | alpha=1.52 | 1.5208 | 60.0 | 0.75 | 0.75 | 1382400 | 2997 | 171 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC2 | alpha=1.57_cs32 | 1.5708 | 60.0 | 60.0 | 60.0 | 17280 | 2997 | 171 | 38.2 | n/a | 38.7 | 170 | reduce deltaT before rerun |
| TC3 | alpha=0 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 2997 | 171 | 0.124 | n/a | 38.6 | 170 | no advective deltaT change indicated |
| TC3 | alpha=1.0472 | 1.0472 | 60.0 | 0.75 | 0.75 | 576000 | 2997 | 171 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC4 | run_u0_20 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 10194 | 316 | n/a | n/a | n/a | 313 | cannot verify advective CFL until inputs/run output exist |
| TC5 | run_standard_cs32 | 0 | 30.0 | 30.0 | 30.0 | 43200 | 5960 | 242 | 0.022 | n/a | 20.0 | 120 | no advective deltaT change indicated |
| TC6 | standard | 0 | 30.0 | 30.0 | 30.0 | 40320 | 8000 | 280 | 0.126 | n/a | 100.0 | 139 | no advective deltaT change indicated |
| TC7 | run_analysis_cs32 | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.237 | n/a | 55.4 | 115 | no advective deltaT change indicated |

### Decisions

No completed run exceeds advective CFL 1.0. TC2 alpha=0.05 reaches about 0.56 in the saved fields, above the conservative 0.5 margin; deltaT <= 8.93 s would keep the saved-output maximum under 0.5.

TC5 does not look like a simple CFL failure: the initial advective CFL is small, but archived fields become non-finite before the required late-day checks and the CG residuals later print NaN. Treat TC5 as needing a run-health rerun with the corrected H0, smaller timestep, and explicit viscosity before using later-day plots.

TC4 and TC7 should be read from the job schedule, not from their template `nTimeSteps`. TC4 now has completed `run_u0_20` output and its saved advective CFL can be read from archived U/V fields, while TC7 needs a finite rerun with the filtered 25 s setup before saved-output CFL can be reported.

TC7 cannot be accepted from the current archived `run_analysis` output because the saved fields become non-finite after initialization; rerun the filtered 25 s setup.

### Check commands

```bash
cd /home/hmh85/scratch/MITgcm
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments
```

```bash
cd Sandbox/vortexSphere_Williamson_TC1
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc1.py --all

cd Sandbox/vortexSphere_Williamson_TC2
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc2.py --all

cd Sandbox/vortexSphere_Williamson_TC3
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc3.py --all

cd Sandbox/vortexSphere_Williamson_TC4
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc4.py --all

cd Sandbox/vortexSphere_Williamson_TC5
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc5.py --all

cd Sandbox/vortexSphere_Williamson_TC6
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc6.py --all

cd Sandbox/vortexSphere_Williamson_TC7
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_tc7.py --all
```

For TC3-TC7, use the project venv; the system `python3` is too old for the shared `williamson_cfl.py` helper.
