# TC4 Forcing Hook

Williamson TC4 is a forced nonlinear exact-solution test: a translating
low is held on its analytic trajectory by source terms in the shallow-water
momentum and continuity equations.

Current status:

- `tc4_forcing.F` evaluates the analytic TC4 source terms on the sphere.
- `external_forcing_surf.F` injects momentum forcing through
  `surfaceForcingU/V` and the height source through `EmPmR`.
- `input/gendata_ref.py` writes TC4 translating-low initial fields.
- `input/data` uses the TC4 reference depth `gh0/g = 10193.67991845056 m`
  and enables `useRealFreshWaterFlux` so the free-surface source is active.
- `jobs/large/job_tc4_c1.slurm` is the large `u0=20 m/s` standard case and
  stages output in `run_u0_20`.

The standard Williamson TC4 suite also includes `u0=40 m/s`. The forcing code
and generator read `TC4_U0_VALUE`, so a second job can use the same executable
once a `u0=40` Slurm wrapper/run directory is added.

Verified preflight:

- TC4 initial-condition generation with `TC4_U0_VALUE=20.0`.
- MPI build in `build/c1_large48` using the same GCC/OpenMPI modules as the
  Slurm job.
- Submission wrapper syntax and executable-presence checks.

Remaining scientific validation after the first run:

- Compare day-5 fields with the analytic TC4 solution.
- Add/confirm TC4 error-norm postprocessing against the time-dependent exact
  height and velocity fields.
