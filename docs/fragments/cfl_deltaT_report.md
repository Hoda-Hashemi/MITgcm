## Time-step safety check (CFL and deltaT)

This check asks whether the submitted model time step is small enough for the resolved velocities and grid spacing in each Williamson run.

For each submitted Williamson job, the audited Courant numbers are `C_adv = |u| deltaT / dx + |v| deltaT / dy` and `C_g = sqrt(gH) deltaT / dx`. The spherical-polar grid is 1440 x 720 with 0.25 degree spacing; the smallest center-cell zonal metric is 60.65 m in the polar row and `dy = 27798.73 m`.

The gravity-wave column is an explicit-wave diagnostic only. All active runs set `implicitFreeSurface=.TRUE.`, so the external gravity wave is handled by the implicit free-surface solve rather than by the explicit advective CFL limit.

The template `input/data` files define the grid and default schedule, but vortexSphere submitted jobs are controlled by `jobs/large/job_*.slurm`: each job exports `DELTA_T`, computes `nTimeSteps=round(TOTAL_SECONDS/DELTA_T)`, copies the template input directory, and rewrites the run-local `data`. The table therefore separates template, submitted job, and archived run deltaT.

TC1 has two entries here: vortexSphere TC1 uses the submitted lat-lon jobs, while `TC1 cubed` is the MITgcm `advect_cs` tutorial and keeps the tutorial `deltaT=2700 s`. The cubed CFL values come from MITgcm monitor output.

The `run dir` and `job file` columns are the exact archived sources used for the row; `run dt(s)` is read from that run-local `data` file when it exists.

| case | run | run dir | job file | alpha(rad) | template dt(s) | job dt(s) | run dt(s) | steps | H(m) | sqrt(gH) | init adv CFL | saved adv CFL | max speed | explicit Cg,x | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TC1 | alpha=0 | run_alpha_0 | job_tc1_c1.slurm | 0 | 1.00 | 60.0 | 60.0 | 17280 | 4000 | 198 | 0.083 | n/a | 38.6 | 196 | no advective deltaT change indicated |
| TC1 | alpha=0.05 | run_alpha_0.05 | job_tc1_c2.slurm | 0.05 | 1.00 | 10.0 | 10.0 | 103680 | 4000 | 198 | 0.332 | n/a | 38.6 | 32.7 | no advective deltaT change indicated |
| TC1 | alpha=1.52 | run_alpha_1.52 | job_tc1_c3.slurm | 1.5208 | 1.00 | 0.75 | 0.75 | 1382400 | 4000 | 198 | 0.490 | n/a | 1095 | 2 | no advective deltaT change indicated |
| TC1 | alpha=1.57 | run_alpha_1.57 | job_tc1_c4.slurm | 1.5708 | 1.00 | 0.75 | 0.75 | 1382400 | 4000 | 198 | 0.477 | n/a | 38.7 | 2 | no advective deltaT change indicated |
| TC2 | alpha=0 | run_alpha_0_latlon | job_tc2_c1.slurm | 0 | 60.0 | 60.0 | 60.0 | 17280 | 2997 | 171 | 0.083 | 0.083 | 38.6 | 170 | no advective deltaT change indicated |
| TC2 | alpha=0.05 | run_alpha_0.05_latlon_rotcori_12day | job_tc2_c2_rotcori_12day.slurm | 0.05 | 60.0 | 10.0 | 10.0 | 103680 | 2997 | 171 | 0.332 | 0.603 | 38.6 | 28.3 | above 0.5 margin; use deltaT <= 8.30 s for margin |
| TC2 | alpha=1.52 | run_alpha_1.52_latlon_rotcori_12day | job_tc2_c3_rotcori_12day.slurm | 1.5208 | 60.0 | 1.00 | 1.00 | 1036800 | 2997 | 171 | 0.636 | 0.728 | 44.2 | 3 | above 0.5 margin; use deltaT <= 0.69 s for margin |
| TC2 | alpha=1.57 | run_alpha_1.57_latlon_rotcori_12day_dt0p5 | job_tc2_c4_rotcori_12day_dt0p5.slurm | 1.5708 | 60.0 | 0.50 | 0.50 | 2073600 | 2997 | 171 | 0.318 | 0.377 | 46.2 | 1 | no advective deltaT change indicated |
| TC3 | alpha=0 | run_alpha_0_cs32 | job_tc3_c1.slurm | 0 | 60.0 | 60.0 | n/a | 7200 | 2997 | 171 | 0.021 | n/a | 54.5 | 170 | no advective deltaT change indicated |
| TC3 | alpha=1.0472 | run_alpha_1.0472_cs32 | job_tc3_c2.slurm | 1.0472 | 60.0 | 0.75 | n/a | 576000 | 2997 | 171 | 2.38e-04 | n/a | 49.9 | 2 | no advective deltaT change indicated |
| TC4 | run_u0_20 | run_u0_20 | job_tc4_c1.slurm | 0 | 60.0 | 60.0 | 60.0 | 7200 | 10194 | 316 | 0.139 | n/a | 40.0 | 313 | no advective deltaT change indicated |
| TC4 | run_u0_40 | run_u0_40 | job_tc4_c2.slurm | 0 | 60.0 | 60.0 | 60.0 | 7200 | 10194 | 316 | 0.192 | 0.193 | 56.4 | 313 | no advective deltaT change indicated |
| TC5 | run_standard_cs32 | run_standard_cs32 | job_tc5_c1.slurm | 0 | 30.0 | 30.0 | 30.0 | 43200 | 5960 | 242 | 0.005 | n/a | 28.3 | 120 | no advective deltaT change indicated |
| TC6 | standard | run_standard | job_tc6_c1.slurm | 0 | 30.0 | 30.0 | 30.0 | 40320 | 8000 | 280 | 0.126 | n/a | 100.0 | 139 | no advective deltaT change indicated |
| TC7 | run_c1_19781221_0000 | run_c1_19781221_0000 | job_tc7_c1.slurm | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.009 | n/a | 56.6 | 115 | no advective deltaT change indicated |
| TC7 | run_c2_19790116_0000 | run_c2_19790116_0000 | job_tc7_c2.slurm | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.012 | n/a | 73.3 | 115 | no advective deltaT change indicated |
| TC7 | run_c3_19790109_0000 | run_c3_19790109_0000 | job_tc7_c3.slurm | 0 | 25.0 | 25.0 | 25.0 | 17280 | 8000 | 280 | 0.011 | n/a | 70.0 | 115 | no advective deltaT change indicated |
| TC1 cubed | advect_cs alpha=0 | alpha_0/run_alpha_0 | advect_cs tutorial | 0 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.768 | n/a | n/a | above 0.5 margin; use deltaT <= 1757.71 s for margin |
| TC1 cubed | advect_cs alpha=0.05 | alpha_0.05/run_alpha_0.05 | advect_cs tutorial | 0.05 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.786 | n/a | n/a | above 0.5 margin; use deltaT <= 1717.65 s for margin |
| TC1 cubed | advect_cs alpha=1.52 | alpha_1.52/run_alpha_1.52 | advect_cs tutorial | 1.5208 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.786 | n/a | n/a | above 0.5 margin; use deltaT <= 1717.65 s for margin |
| TC1 cubed | advect_cs alpha=1.57 | alpha_1.57/run_alpha_1.57 | advect_cs tutorial | 1.5708 | 2700 | 2700 | 2700 | 384 | 100000 | n/a | n/a | 0.768 | n/a | n/a | above 0.5 margin; use deltaT <= 1757.71 s for margin |

### Decisions

No completed run exceeds advective CFL 1.0. In the verified TC2 suite, alpha=0.05 and alpha=1.52 sit above the conservative 0.5 saved-output margin; use deltaT <= 8.30 s and <= 0.69 s, respectively, if that extra margin is required. TC2 alpha=1.57 uses the completed 0.5 s rotated-Coriolis run and remains below the 0.5 margin.

TC5 is now a completed CS32 rerun: the initial and monitored advective CFL remain small, the final saved state fields are finite through day 15, and the mountain is verified as static bathymetry rather than an eta bump. TC6 is the completed Rossby-Haurwitz `run_standard` row at deltaT=30 s for 40320 steps, and TC7 uses cubed-sphere compact initial fields for the submitted three-date suite.

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
