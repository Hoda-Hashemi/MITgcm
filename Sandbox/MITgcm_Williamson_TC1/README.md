# MITgcm Williamson TC1 existing tutorial setup

This folder stages MITgcm's existing `verification/advect_cs` tutorial as the TC1 pure-advection comparison.

Logic:
- `advect_cs` is a passive-tracer cubed-sphere advection tutorial.
- It is not a free-surface test; `Eta` is diagnostic here, not the target.
- The target is passive tracer return after one 12-day solid-body rotation.
- `ini_vel.F` reads `tc1_alpha.txt` and applies the rotated Williamson TC1 velocity.
- Output is written under `Sandbox/output/existingTutorials/test1/MITGCM_Williamson_TC1/advect_cs`.

Large-queue jobs:
- `jobs/large/job_tc1_c1.slurm`: alpha 0
- `jobs/large/job_tc1_c2.slurm`: alpha 0.05
- `jobs/large/job_tc1_c3.slurm`: alpha 1.5207963267948966
- `jobs/large/job_tc1_c4.slurm`: alpha 1.5707963267948966
- `jobs/submission_p_large.sh`: submits the four jobs when run manually

Main paths:
- setup: `Sandbox/MITgcm_Williamson_TC1/existingTutorials/advect_cs`
- jobs: `Sandbox/MITgcm_Williamson_TC1/jobs/large/job_tc1_c*.slurm`
- postprocess helper: `Sandbox/MITgcm_Williamson_TC1/tools/postprocess_advect_cs.py`
