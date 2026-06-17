## Time-step safety check (CFL and deltaT)

This check asks whether the submitted model time step is small enough for the resolved velocities and grid spacing in each Williamson run.

For each submitted Williamson job, the audited Courant numbers are `C_adv = |u| deltaT / dx + |v| deltaT / dy` and `C_g = sqrt(gH) deltaT / dx`. The spherical-polar grid is 1440 x 720 with 0.25 degree spacing; the smallest center-cell zonal metric is 60.65 m in the polar row and `dy = 27798.73 m`.

The gravity-wave column is an explicit-wave diagnostic only. All active runs set `implicitFreeSurface=.TRUE.`, so the external gravity wave is handled by the implicit free-surface solve rather than by the explicit advective CFL limit.

The template `input/data` files define the grid and default schedule, but submitted jobs are controlled by `jobs/large/job_*.slurm`: each job exports `DELTA_T`, computes `nTimeSteps=round(TOTAL_SECONDS/DELTA_T)`, copies the template input directory, and rewrites the run-local `data`. The table therefore separates template, submitted job, and archived run deltaT.

TC1 has the most visible discrepancy: its legacy `tools/check_cfl_tc1.py` prints the template `deltaT=1 s`, while the submitted jobs use 60 s, 10 s, or 0.75 s and the archived run-local `data` files confirm those values. Use this audit table for submitted-run CFL decisions.

| case | run | alpha(rad) | template dt(s) | job dt(s) | run dt(s) | steps | H(m) | sqrt(gH) | init adv CFL | saved adv CFL | max speed | explicit Cg,x | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TC2 | alpha=0 | 0 | 60.0 | 60.0 | 60.0 | 17280 | 2997 | 171 | 0.083 | 0.083 | 38.6 | 170 | no advective deltaT change indicated |
| TC2 | alpha=0.05 | 0.05 | 60.0 | 10.0 | 10.0 | 103680 | 2997 | 171 | 0.332 | 0.561 | 38.6 | 28.3 | above 0.5 margin; use deltaT <= 8.91 s for margin |
| TC2 | alpha=1.52 | 1.5208 | 60.0 | 0.75 | 60.0 | 1382400 | 2997 | 171 | 0.477 | 38.2 | 38.7 | 2 | investigate non-finite run; not explained by initial advective CFL |
| TC2 | alpha=1.57 | 1.5708 | 60.0 | 0.75 | 60.0 | 1382400 | 2997 | 171 | 0.477 | 38.2 | 38.7 | 2 | investigate non-finite run; not explained by initial advective CFL |
| TC3 | alpha=0 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 2997 | 171 | 0.124 | 0.124 | 38.6 | 170 | no advective deltaT change indicated |
| TC3 | alpha=1.0472 | 1.0472 | 60.0 | 0.75 | 0.75 | 576000 | 2997 | 171 | 0.477 | 0.486 | 50.1 | 2 | no advective deltaT change indicated |
| TC4 | run_u0_20 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 10194 | 316 | n/a | 0.139 | 40.0 | 313 | no advective deltaT change indicated |
| TC5 | standard | 0 | 60.0 | 60.0 | 60.0 | 21600 | 5960 | 242 | 0.043 | 0.909 | 40.1 | 239 | above 0.5 margin; use deltaT <= 33.00 s for margin |
| TC6 | standard | 0 | 30.0 | 30.0 | 30.0 | 40320 | 8000 | 280 | 0.126 | 0.131 | 95.2 | 139 | no advective deltaT change indicated |
| TC7 | analysis | 0 | 60.0 | 60.0 | 60.0 | 7200 | 8000 | 280 | 21.2 | 21.2 | 56.4 | 277 | reduce deltaT before rerun |

### Decisions

No completed run exceeds advective CFL 1.0. TC2 alpha=0.05 reaches about 0.56 in the saved fields, above the conservative 0.5 margin; deltaT <= 8.93 s would keep the saved-output maximum under 0.5.

TC5 does not look like a simple CFL failure: the initial advective CFL is 0.043, but archived fields become non-finite after iteration 1440 and the CG residuals later print NaN. Treat TC5 as needing a run-health fix or a targeted shorter-deltaT rerun before using later-day plots.

TC4 and TC7 should be read from the job schedule, not from their template `nTimeSteps`: both submitted jobs target 5 days at 60 s, i.e. 7200 steps. TC4 now has completed `run_u0_20` output and its saved advective CFL can be read from archived U/V fields, while TC7 has input data staged but needs a smaller deltaT before a final rerun.

TC7 cannot be fully audited yet because completed `run_analysis` output is not archived. The staged analyzed input is present, and the preflight advective CFL indicates the current 60 s job schedule is too large for final validation.

### Check commands

```bash
cd /home/hmh85/scratch/MITgcm
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments
```

```bash
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
