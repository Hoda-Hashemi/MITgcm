## Time-step safety check (CFL and deltaT)

This check asks whether the submitted model time step is small enough for the resolved velocities and grid spacing in each Williamson run.

For each submitted Williamson job, the audited Courant numbers are `C_adv = |u| deltaT / dx + |v| deltaT / dy` and `C_g = sqrt(gH) deltaT / dx`. The spherical-polar grid is 1440 x 720 with 0.25 degree spacing; the smallest center-cell zonal metric is 60.65 m in the polar row and `dy = 27798.73 m`.

The gravity-wave column is an explicit-wave diagnostic only. All active runs set `implicitFreeSurface=.TRUE.`, so the external gravity wave is handled by the implicit free-surface solve rather than by the explicit advective CFL limit.

The template `input/data` files define the grid and default schedule, but vortexSphere submitted jobs are controlled by `jobs/large/job_*.slurm`: each job exports `DELTA_T`, computes `nTimeSteps=round(TOTAL_SECONDS/DELTA_T)`, copies the template input directory, and rewrites the run-local `data`. The table therefore separates template, submitted job, and archived run deltaT.

TC1 has two entries here: vortexSphere TC1 uses the submitted lat-lon jobs, while `TC1 cubed` is the MITgcm `advect_cs` tutorial and keeps the tutorial `deltaT=2700 s`. The cubed CFL values come from MITgcm monitor output.

| case | run | alpha(rad) | template dt(s) | job dt(s) | run dt(s) | steps | H(m) | sqrt(gH) | init adv CFL | saved adv CFL | max speed | explicit Cg,x | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TC1 | alpha=0 | 0 | 1.00 | 60.0 | 60.0 | 17280 | 4000 | 198 | 0.083 | n/a | 38.6 | 196 | no advective deltaT change indicated |
| TC1 | alpha=0.05 | 0.05 | 1.00 | 10.0 | 10.0 | 103680 | 4000 | 198 | 0.332 | n/a | 38.6 | 32.7 | no advective deltaT change indicated |
| TC1 | alpha=1.52 | 1.5208 | 1.00 | 0.75 | 0.75 | 1382400 | 4000 | 198 | 0.490 | n/a | 1095 | 2 | no advective deltaT change indicated |
| TC1 | alpha=1.57 | 1.5708 | 1.00 | 0.75 | 0.75 | 1382400 | 4000 | 198 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC2 | alpha=0 | 0 | 60.0 | 60.0 | 60.0 | 17280 | 2997 | 171 | 0.083 | n/a | 38.6 | 170 | no advective deltaT change indicated |
| TC2 | alpha=0.05 | 0.05 | 60.0 | 10.0 | 10.0 | 103680 | 2997 | 171 | 0.332 | n/a | 38.6 | 28.3 | no advective deltaT change indicated |
| TC2 | alpha=1.52 | 1.5208 | 60.0 | 0.75 | 0.75 | 1382400 | 2997 | 171 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC2 | alpha=1.57_cs32 | 1.5708 | 60.0 | 60.0 | 60.0 | 17280 | 2997 | 171 | 0.021 | n/a | 54.6 | 170 | no advective deltaT change indicated |
| TC3 | alpha=0 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 2997 | 171 | 0.124 | n/a | 38.6 | 170 | no advective deltaT change indicated |
| TC3 | alpha=1.0472 | 1.0472 | 60.0 | 0.75 | 0.75 | 576000 | 2997 | 171 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC4 | run_u0_20 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 10194 | 316 | 0.139 | n/a | 40.0 | 313 | no advective deltaT change indicated |
| TC4 | run_u0_40 | 0 | 60.0 | 60.0 | 60.0 | 7200 | 10194 | 316 | 0.192 | 0.193 | 56.4 | 313 | no advective deltaT change indicated |
| TC5 | run_standard_cs32 | 0 | 30.0 | 30.0 | 30.0 | 43200 | 5960 | 242 | 0.005 | n/a | 28.3 | 120 | no advective deltaT change indicated |
| TC6 | standard | 0 | 30.0 | 30.0 | 30.0 | 40320 | 8000 | 280 | 0.126 | n/a | 100.0 | 139 | no advective deltaT change indicated |
| TC7 | run_c1_19781221_0000 | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.009 | n/a | 56.6 | 115 | no advective deltaT change indicated |
| TC7 | run_c2_19790116_0000 | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.012 | n/a | 73.3 | 115 | no advective deltaT change indicated |
| TC7 | run_c3_19790109_0000 | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.011 | n/a | 70.0 | 115 | no advective deltaT change indicated |
| TC1 cubed | advect_cs alpha=0 | 0 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.768 | n/a | n/a | above 0.5 margin; use deltaT <= 1757.71 s for margin |
| TC1 cubed | advect_cs alpha=0.05 | 0.05 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.786 | n/a | n/a | above 0.5 margin; use deltaT <= 1717.65 s for margin |
| TC1 cubed | advect_cs alpha=1.52 | 1.5208 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.786 | n/a | n/a | above 0.5 margin; use deltaT <= 1717.65 s for margin |
| TC1 cubed | advect_cs alpha=1.57 | 1.5708 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.768 | n/a | n/a | above 0.5 margin; use deltaT <= 1757.71 s for margin |

### Decisions

No completed run exceeds advective CFL 1.0. TC2 alpha=0.05 reaches about 0.56 in the saved fields, above the conservative 0.5 margin; deltaT <= 8.93 s would keep the saved-output maximum under 0.5.

TC5 is now a completed CS32 rerun: the initial and monitored advective CFL remain small, the final saved state fields are finite through day 15, and the mountain is verified as static bathymetry rather than an eta bump. TC7 uses cubed-sphere compact initial fields for the submitted three-date suite.

`n/a` means the audit could not read a finite CFL source for that column: missing archived U/V fields, unavailable initial-velocity hook, or a cubed-sphere row where the spherical-polar gravity-wave metric is not used. TC4 now includes both `run_u0_20` and completed `run_u0_40` output.

TC7 has three completed analyzed-state rows: 21 Dec 1978, 16 Jan 1979, and 9 Jan 1979. The table values are CS32 compact-input CFL checks for the completed 25 s, 48-rank runs.

### Check commands

```bash
cd /home/hmh85/scratch/MITgcm
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments
```

```bash
cd Sandbox/MITgcm_Williamson_TC1
/home/hmh85/scratch/MITgcm/.venv/bin/python tools/postprocess_advect_cs.py --check-only

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
