# Scientific Methodology and MITgcm Implementation Notes

Last updated: 2026-07-07

This note documents the methodology implemented in this branch for the
`Sandbox/vortexSphere_*` MITgcm experiments, especially the Williamson
test-case suite `vortexSphere_Williamson_TC1` through
`vortexSphere_Williamson_TC7`.  It is written as a scientific Markdown file
with LaTeX-style equations so that sections can be copied into a LaTeX report.

The generated website for inspecting figures, movies, CFL tables, and
conservation tables is:

- `docs/index.html`
- `docs/fragments/testcase1.html` through `docs/fragments/testcase7.html`
- `docs/fragments/cfl_deltaT_report.md`
- `docs/fragments/conservation_report.md`
- `docs/assets/williamson/submitted_time_step_cfl_table.tex`

The plotted output artifacts are written under:

- `Sandbox/output/TestCase1` through `Sandbox/output/TestCase7`
- `Sandbox/output/TestCase*/Snapshots*/.../manifest.json`
- `Sandbox/output/TestCase*/Diagnosis/conservation/summary.csv`
- `Sandbox/output/TestCase3/Diagnosis/error/.../TC3_error_norms.csv`

## 1. Governing Model

The active Williamson cases are implemented as one-layer rotating shallow-water
or passive-tracer problems in MITgcm.  The prognostic fields are horizontal
velocity, free-surface anomaly, and, where enabled, a passive tracer:

$$
{\bf u}=(u,v), \qquad \eta, \qquad q .
$$

The spherical longitude-latitude coordinates are denoted by
$(\lambda,\theta)$, with planetary radius $a$.  The active shallow-water form is

$$
\frac{D{\bf u}}{Dt}
 + f \hat{k}\times{\bf u}
 =
 -g\nabla_s \eta
 + {\cal D}_{\bf u}
 + {\bf F},
$$

$$
\frac{\partial \eta}{\partial t}
 + \nabla_s\cdot\left[(H+\eta){\bf u}\right]
 =
 S_\eta ,
$$

and the passive TC1 tracer satisfies

$$
\frac{Dq}{Dt}={\cal D}_q+S_q .
$$

In MITgcm these equations are advanced by the standard momentum, continuity,
and implicit free-surface machinery.  The submitted cases use
`rigidLid=.FALSE.`, `implicitFreeSurface=.TRUE.`, `exactConserv=.TRUE.`, and
the finite-volume spherical metric terms supplied by MITgcm.  The external
gravity-wave CFL is therefore diagnostic only; the free-surface pressure is
handled by the implicit solver.

The linearized pressure variable is

$$
\phi_s = g\eta .
$$

The pressure-correction step can be summarized as

$$
{\bf u}^{n+1} = {\bf u}^* - \Delta t \nabla_s \phi_s^{n+1},
$$

$$
\frac{\eta^{n+1}-\eta^n}{\Delta t}
 + \nabla_s\cdot\left(H{\bf u}^{n+1}\right)=0 .
$$

Substitution gives the elliptic free-surface system solved by MITgcm:

$$
A\phi_s^{n+1}=b .
$$

The relevant MITgcm source landmarks are:

- `model/src/forward_step.F`: top-level time-step driver.
- `model/src/dynamics.F`: momentum tendency assembly.
- `model/src/solve_for_pressure.F`: implicit free-surface solve.
- `model/src/ini_cg2d.F`: two-dimensional elliptic operator assembly.
- `model/src/momentum_correction_step.F`: pressure-gradient correction.
- `model/src/integr_continuity.F`: free-surface continuity update.
- `model/src/ini_cori.F`: Coriolis-map initialization, including
  `selectCoriMap=3` custom Coriolis reads.

## 2. How MITgcm Communicates With Each Experiment

Each experiment directory follows the MITgcm verification/tutorial pattern:

```text
Sandbox/vortexSphere_Williamson_TCX/
  code/                  compile-time options and SIZE.h
  input/                 template data, data.pkg, data.diagnostics, gendata_ref.py
  jobs/large/            Slurm scripts that stage, compile, run, and postprocess
  run_*                  archived run-local copies and MITgcm MDS output
  tools/                 case-local CFL and setup wrappers
```

The communication boundary between Python setup and MITgcm is deliberately
file-based:

1. `jobs/large/job_*.slurm` exports the case, rotation angle, time step,
   duration, build directory, and run directory.
2. The Slurm script copies `input/` into a run-local directory.
3. The script rewrites run-local `data` entries such as `deltaT`,
   `nTimeSteps`, `dumpFreq`, grid type, and, where needed, `selectCoriMap`.
4. `gendata_ref.py` generates big-endian binary initial conditions:
   `bathymetry*.bin`, `eta_init.bin`, `u_init.bin`, `v_init.bin`, and, for
   rotated Coriolis cases, `fCoriC.bin`, `fCoriG.bin`, and `fCorCs.bin`.
5. MITgcm reads these files through `PARM05` in `data`:

$$
bathyFile,\quad pSurfInitFile,\quad uVelInitFile,\quad vVelInitFile .
$$

6. MITgcm writes MDS `.data/.meta` output.
7. `Sandbox/Scripts/TC*.py` and shared diagnostics scripts turn MDS output into
   snapshots, conservation summaries, error norms, and website fragments.

The main shared Python implementation files are:

- `Sandbox/Scripts/williamson_fields.py`: analytic fields and constants.
- `Sandbox/Scripts/williamson_cube.py`: CS32 compact cubed-sphere packing.
- `Sandbox/Scripts/williamson_setup_checks.py`: setup validation.
- `Sandbox/Scripts/williamson_cfl.py`: initial-condition CFL audit.
- `Sandbox/Scripts/cfl_deltaT_audit.py`: website CFL table generation.
- `Sandbox/Scripts/conservation_diagnostics.py`: conservation diagnostics.
- `Sandbox/Scripts/github_pages.py`: website assembly.

## 3. Grid, Decomposition, and Time-Step Summary

The latitude-longitude Williamson cases use a global
$1440\times720$ grid with $0.25^\circ$ spacing unless otherwise noted.  The
CS32 cubed-sphere runs use MITgcm compact cubed-sphere fields with six
$32\times32$ faces.

| Case | Main grid | Compile decomposition | Main duration | Submitted time step | Notes |
| --- | --- | --- | ---: | ---: | --- |
| TC1 | lat-lon 1440 x 720 | `sNx=120`, `sNy=180`, `nPx=12`, `nPy=4` | 12 d | 60, 10, 0.75 s | passive cosine-bell transport |
| TC2 | lat-lon 1440 x 720 | `sNx=180`, `sNy=120`, `nPx=8`, `nPy=6` | 12 d | 60, 10, 1, 0.5 s | rotated Coriolis required for tilted cases |
| TC3 | lat-lon archive and CS32 | lat-lon 45-rank shared-node archive; CS32 requires 48 ranks | 5 d | 60, 0.75, 0.5 s | compact-support jet; rotated lat-lon needs matching Coriolis |
| TC4 | lat-lon 1440 x 720 | `sNx=360`, `sNy=180`, `nPx=4`, `nPy=4` | 5 d | 60 s | forced translating low, `U0=20` and `40` m/s |
| TC5 | CS32 | `sNx=16`, `sNy=8`, `nPx=48`, `nPy=1` | 15 d | 30 s | mountain is bathymetry, not eta |
| TC6 | lat-lon 1440 x 720 | `sNx=120`, `sNy=180`, `nPx=12`, `nPy=4` | 14 d | 30 s | Rossby-Haurwitz wave |
| TC7 | CS32 | `sNx=16`, `sNy=8`, `nPx=48`, `nPy=1` | 5 d | 25 s | analyzed 500 mb states |

The submitted-run CFL audit is stored in
`docs/assets/williamson/submitted_time_step_cfl_table.tex`.  The audit uses

$$
C_{\rm adv} =
\frac{|u|\Delta t}{\Delta x}
+\frac{|v|\Delta t}{\Delta y},
\qquad
C_g =
\frac{\sqrt{gH}\Delta t}{\Delta x}.
$$

The advective CFL is the decision criterion used here.  The gravity-wave CFL is
reported only as an explicit-wave diagnostic because `implicitFreeSurface` is
enabled.

## 4. Initial Conditions Implemented

### 4.1 TC1: Cosine-Bell Passive Advection

TC1 advects a compact cosine bell by a prescribed solid-body velocity field.
The bell is centered at

$$
(\lambda_c,\theta_c)=(3\pi/2,0),
$$

with spherical distance $r$ and compact radius $R=a/3$:

$$
q(\lambda,\theta,0)
=
\frac{q_0}{2}
\left[
  1+\cos\left(\frac{\pi r}{R}\right)
\right],
\qquad r<R,
$$

$$
q(\lambda,\theta,0)=0,\qquad r\ge R.
$$

The velocity is the Williamson solid-body field, tilted by angle $\alpha$:

$$
u = u_0\left(\cos\theta\cos\alpha
 + \cos\lambda\sin\theta\sin\alpha\right),
$$

$$
v = -u_0\sin\lambda\sin\alpha .
$$

The implemented values are $u_0=2\pi a/(12\,{\rm days})$, $q_0=1000$,
$a=6.371\times10^6$ m, and active alphas
$0$, $0.05$, $1.5207963268$, and $\pi/2$.

### 4.2 TC2: Global Steady Nonlinear Zonal Geostrophic Flow

TC2 initializes a nonlinear shallow-water steady state.  Define the rotated
latitude coordinate through

$$
\mu(\lambda,\theta;\alpha)
=
-\cos\lambda\cos\theta\sin\alpha
 + \sin\theta\cos\alpha .
$$

Then

$$
u = u_0\left(\cos\theta\cos\alpha
 + \cos\lambda\sin\theta\sin\alpha\right),
$$

$$
v = -u_0\sin\lambda\sin\alpha ,
$$

and

$$
gh =
gh_0 -
\left(a\Omega u_0+\frac{u_0^2}{2}\right)\mu^2,
\qquad
\eta = \frac{gh}{g}-H_0.
$$

The implemented constants are

$$
u_0=\frac{2\pi a}{12\,{\rm days}},\qquad
gh_0=2.94\times10^4\ {\rm m^2\,s^{-2}},\qquad
H_0=gh_0/g .
$$

The tilted TC2 cases require custom rotated Coriolis.  The run-local scripts set
`selectCoriMap=3`, and MITgcm reads

$$
f_C = 2\Omega\mu_C,\qquad
f_G = 2\Omega\mu_G,\qquad
f_{\cos}=2\Omega\sqrt{1-\mu_C^2}
$$

from `fCoriC.bin`, `fCoriG.bin`, and `fCorCs.bin`.

### 4.3 TC3: Steady Geostrophic Flow With Compact Support

TC3 replaces the global solid-body jet with a smooth compact-support jet.  In
this context, compact support means the analytic wind profile is nonzero only
inside a finite latitude band and is exactly zero outside it.

The implemented compact profile uses

$$
\theta_0=-\frac{\pi}{6},\qquad
\theta_1=\frac{\pi}{2},\qquad
x_e=0.3,
$$

$$
x=x_e\frac{\theta'-\theta_0}{\theta_1-\theta_0},
$$

and

$$
U(\theta') =
u_0
\exp\left(
 -\frac{1}{x}
 -\frac{1}{x_e-x}
 +\frac{4}{x_e}
\right),
\qquad
\theta_0 < \theta' < \theta_1,
$$

with $U=0$ outside the interval.  Here $\theta'$ is the rotated latitude
defined from $\mu=\sin\theta'$.

The free-surface height is obtained by numerically integrating the geostrophic
balance:

$$
\frac{dH}{d\theta'}
=
-\frac{a U(\theta')}{g}
\left[
2\Omega\sin\theta'
+ \frac{U(\theta')\tan\theta'}{a}
\right],
$$

then subtracting $H_0$ to write `eta_init.bin`.

For the rotated C2 run, the same `selectCoriMap=3` mechanism as TC2 is used.
This was a critical implementation fix: rotating only `eta_init/u_init/v_init`
while leaving MITgcm on the unrotated spherical Coriolis map produced a
non-steady state.

### 4.4 TC4: Forced Nonlinear Translating Low

TC4 is an exact-solution forced nonlinear shallow-water test.  The implemented
background flow is

$$
u_b(\theta)=u_0\sin^{14}(2\theta),
$$

and the localized streamfunction disturbance is represented as

$$
\psi =
\psi_0
\exp\left[
-\sigma
\frac{1-c}{1+c}
\right],
$$

where

$$
c = \sin\theta_c\sin\theta
+\cos\theta_c\cos\theta\cos\left(\lambda-\frac{u_0}{a}t-\lambda_c\right).
$$

The implemented constants are

$$
gh_0=10^5\ {\rm m^2\,s^{-2}},\quad
H_0=gh_0/g,\quad
\theta_c=\pi/4,\quad
\lambda_c=-7\pi/12,\quad
\sigma=12.74244^2.
$$

Both required velocities, $u_0=20$ m/s and $u_0=40$ m/s, are represented by
separate jobs and output directories.

### 4.5 TC5: Zonal Flow Over an Isolated Mountain

TC5 uses a zonal flow over a conical mountain.  The corrected implementation
places the mountain in bathymetry rather than in eta:

$$
b(\lambda,\theta) = -\left(H_0-h_s(\lambda,\theta)\right).
$$

The mountain is centered at

$$
\lambda_m=3\pi/2,\qquad \theta_m=\pi/6,
$$

with radius $R_m=\pi/9$ and peak height $h_{s0}=2000$ m:

$$
h_s =
h_{s0}\left(1-\frac{r_m}{R_m}\right),
\qquad r_m < R_m,
$$

and $h_s=0$ outside the mountain.  The flow and eta are

$$
u = u_0\cos\theta,\qquad v=0,\qquad
\eta =
-\frac{\left(a\Omega u_0+\frac{u_0^2}{2}\right)\sin^2\theta}{g}.
$$

The completed verified run uses CS32, 48 MPI ranks, $\Delta t=30$ s, and a
15-day duration.

### 4.6 TC6: Rossby-Haurwitz Wave

TC6 implements the Rossby-Haurwitz wave with wave number $r=4$,
$K=7.848\times10^{-6}$, $\omega=7.848\times10^{-6}$, and $H_0=8000$ m.  The
implemented velocities are

$$
u =
a\omega\cos\theta

+aK\cos^{r-1}\theta
\left(r\sin^2\theta-\cos^2\theta\right)
\cos(r\lambda),
$$

$$
v =
-aKr\cos^{r-1}\theta\sin\theta\sin(r\lambda).
$$

The geopotential height is generated by the analytic expression in
`williamson_fields.py::tc6_height`, then written as

$$
\eta = h-H_0.
$$

The verified run is `run_standard`, $\Delta t=30$ s for 40320 steps, with
maximum initial advective CFL 0.126.

### 4.7 TC7: Analyzed 500 mb Initial State

TC7 uses analyzed atmospheric height and wind fields as shallow-water initial
conditions rather than a closed-form exact solution.  The implementation loads
raw analyzed fields from `TC7_RAW_FILE` or `input/raw/tc7_*_initial_conditions.npz`,
then applies preprocessing:

- reference depth $H_0=8000$ m,
- zonal wavenumber cutoff 42,
- 8 Shapiro smoothing passes,
- cubed-sphere interpolation and vector projection for CS32 runs.

The active verified dates are:

- `run_c1_19781221_0000`
- `run_c2_19790116_0000`
- `run_c3_19790109_0000`

TC7 is judged by finite state health, mass conservation, energy-proxy drift,
and smooth day-0 through day-5 evolution, not by exact analytic error norms.

## 5. Cubed-Sphere Packing and Vector Projection

The CS32 cases use MITgcm compact cubed-sphere I/O.  `williamson_cube.py`
copies `tile001.mitgrid` through `tile006.mitgrid`, evaluates analytic fields
on face centers and face-normal velocity locations, projects east/north
velocity components onto cube-face axes, and writes compact arrays with shape

$$
6\times32\times32.
$$

The compact binary layout is built as a vertical stack of transposed faces:

$$
{\rm compact} =
\begin{bmatrix}
F_1^T\\
F_2^T\\
\vdots\\
F_6^T
\end{bmatrix}.
$$

The preflight script `williamson_cube_preflight.py` verifies the required
`data`, `eedata`, `data.exch2`, tile files, finite initial fields, and an
initial advective CFL below the chosen safety threshold.

## 6. Diagnostics and Verification

### 6.1 Error Norms

For steady cases, the primary numerical diagnostic is the difference between
the saved model state and the analytic initial or target state:

$$
e_\phi(t)=\phi_{\rm MITgcm}(t)-\phi_{\rm ref}(t).
$$

Reported norms include

$$
L_2(e)=
\left(
\frac{\int_\Omega e^2\,dA}{\int_\Omega dA}
\right)^{1/2},
\qquad
L_\infty(e)=\max_\Omega |e|.
$$

TC2 and TC3 verification uses these steady-state drift norms.  TC4 uses the
forced analytic target.  TC5 and TC7 are not exact-solution tests and therefore
use finite-state health and qualitative/diagnostic comparisons.

### 6.2 Conservation Diagnostics

The conservation diagnostics are generated from MITgcm MDS output.  The
implemented quantities are

$$
h = H+\eta,
\qquad
M = \rho_0\int_\Omega h\,dA,
$$

$$
E_m =
\rho_0\int_\Omega hK\,dA
+
\frac{1}{2}\rho_0 g\int_\Omega \eta^2\,dA,
$$

$$
\zeta =
\frac{1}{a\cos\theta}
\left[
\frac{\partial v}{\partial \lambda}
-
\frac{\partial(u\cos\theta)}{\partial \theta}
\right],
$$

$$
q = \frac{\zeta+f}{h},
\qquad
Z = \frac{1}{2}\int_\Omega hq^2\,dA,
$$

and relative drift

$$
\Delta_r X(t)=\frac{X(t)-X(0)}{|X(0)|}.
$$

When `momKE` is available it is used as $K$; otherwise $K=(u^2+v^2)/2$ is
reconstructed from centered velocities.  For rotated runs, $f$ is read from
`fCoriC.bin` when available; otherwise it defaults to
$2\Omega\sin\theta$.

The generated conservation report is:

- `docs/fragments/conservation_report.md`
- `Sandbox/output/TestCase*/Diagnosis/conservation/summary.csv`

## 7. Website and Output Products

The website is generated from local run outputs and fragment files.  The main
viewing entry point is:

```text
docs/index.html
```

Case-specific pages are assembled from:

```text
docs/fragments/testcase1.html
docs/fragments/testcase2.html
docs/fragments/testcase3.html
docs/fragments/testcase4.html
docs/fragments/testcase5.html
docs/fragments/testcase6.html
docs/fragments/testcase7.html
```

The website currently labels the suite as follows:

| Case | Website status | Main statement |
| --- | --- | --- |
| TC1 | Verified | cosine-bell transport and conservation assets published |
| TC2 | Verified | corrected rotated-Coriolis suite used for tilted lat-lon cases |
| TC3 | Validated with caveat | alpha=0 lat-lon is clean; alpha=1.0472 lat-lon is finite but drifty; CS32 alpha=1.0472 rerun queued |
| TC4 | Verified | required $U_0=20$ and $U_0=40$ runs complete |
| TC5 | Verified | corrected CS32 mountain-in-bathymetry run complete |
| TC6 | Verified | Rossby-Haurwitz run complete through day 14 |
| TC7 | Verified | three analyzed-state CS32 runs complete |

Useful website and output regeneration commands are:

```bash
cd /home/hmh85/scratch/MITgcm
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/conservation_diagnostics.py --write-fragments
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/github_pages.py
```

## 8. Development Problems, Bugs, and Fixes

The main scientific and implementation issues found in this branch were:

1. Rotated Coriolis mismatch in TC2 and TC3.

   The tilted analytic fields use a rotated latitude $\theta'$, but MITgcm's
   default spherical Coriolis uses geographic latitude $\theta$.  If only
   `eta_init.bin`, `u_init.bin`, and `v_init.bin` are rotated, the model is not
   initialized in geostrophic balance and emits spurious gravity-wave drift.
   For lat-lon runs, the fix is to set `selectCoriMap=3` and write matching
   `fCoriC.bin`, `fCoriG.bin`, and `fCorCs.bin`.  For CS32 runs, use the
   cubed-sphere grid/exchange setup and do not carry a lat-lon `selectCoriMap`
   override.

2. Latitude-longitude polar CFL sensitivity.

   Near the pole, $\Delta x=a\cos\theta\Delta\lambda$ becomes very small.
   Tilted flows can place nonzero zonal velocity in the polar row.  This forced
   much smaller `deltaT` values for high-tilt TC1, TC2, and TC3 cases.  For the
   TC3 rotated C2 lat-lon run, `deltaT=0.5` s gives
   `max_total_CFL=0.318292`.

3. Invalid MPI decompositions.

   MITgcm's `SIZE.h` decomposition must exactly tile the domain.  A request for
   46 ranks on the 1440 x 720 lat-lon grid is not valid because 46 contains the
   factor 23, which divides neither 1440 nor 720.  The nearest packed valid
   shared-node layout used here is 45 ranks:

   ```text
   nPx=9, nPy=5, sNx=160, sNy=144.
   ```

   The CS32 compact layout is different: the six $32\times32$ faces are packed
   as 48 MITgcm tiles with `sNx=16`, `sNy=8`, `nPx=48`, and `nPy=1`, so a
   45-rank CS32 job is not a valid equal-tile setup.

4. Cubed-sphere packing and vector orientation.

   CS32 runs require correct compact face ordering, face transposition,
   `data.exch2`, `eedata` cubed-sphere exchange settings, and vector projection
   from east/north components to face-local MITgcm velocity axes.  The helper
   `williamson_cube_preflight.py` was added to catch missing tile files, wrong
   `SIZE.h`, missing `data.exch2`, non-finite input fields, and excessive CFL.

5. TC5 mountain placement.

   The mountain must be static bathymetry, not an eta perturbation.  Earlier
   incorrect setups placed topography-like height in `eta_init.bin`.  The fixed
   run writes the mountain through `bathymetry_mountain_tc5.bin` and keeps eta
   as the geopotential drop only.

6. TC7 analyzed-field preprocessing.

   TC7 raw analyzed fields need smoothing, spectral cutoff, finite-field
   validation, and careful interpolation to CS32.  The current pipeline uses a
   zonal wavenumber cutoff of 42 and 8 Shapiro smoothing passes before writing
   compact cubed-sphere fields.

7. Python environment mismatch.

   The system `python3` is Python 3.6 and cannot run the shared helpers that
   use modern type annotations.  Job scripts load `python/3.12` and activate
   `/home/hmh85/scratch/MITgcm/.venv`.  Local checks should use the project
   virtual environment.

8. Slurm resource sharing.

   Several early jobs used `--exclusive` and `--mem=0`, which prevents sharing
   even when only a subset of CPUs is needed.  Later smoke and TC3 lat-lon jobs
   use explicit memory such as `--mem=32G`, no exclusive allocation, and valid
   rank counts to fit on shared `large` nodes.  The current CS32 TC3 C2 rerun
   uses `--mem=64G`, no exclusive allocation, and 48 ranks pinned to `onode13`.

## 9. HPC Execution Metadata

The table below records representative Slurm metadata available from `sacct`.
Failed and cancelled rows are development history, not final scientific
validation.

| Job | Name | State | Node(s) | CPUs | Elapsed | Notes |
| ---: | --- | --- | --- | ---: | ---: | --- |
| 817496 | `TC1_advect_c1` | COMPLETED | onode13 | 64 | 00:03:16 | TC1 tutorial-style run |
| 817497 | `TC1_advect_c2` | COMPLETED | onode14 | 64 | 00:03:16 | TC1 tutorial-style run |
| 817498 | `TC1_advect_c4` | COMPLETED | onode15 | 64 | 00:04:12 | TC1 tutorial-style run |
| 817499 | `TC1_advect_c3` | COMPLETED | onode13 | 64 | 00:03:12 | TC1 tutorial-style run |
| 817507 | `TC4_c2_u0_40_large48` | COMPLETED | onode13 | 48 | 00:44:50 | TC4 $U_0=40$ run |
| 817514 | `TC7_c1_19781221_48` | COMPLETED | onode13 | 64 | 00:10:51 | TC7 analyzed-state run |
| 817515 | `TC7_c2_19790116_48` | COMPLETED | onode14 | 64 | 00:11:27 | TC7 analyzed-state run |
| 817516 | `TC7_c3_19790109_48` | COMPLETED | onode15 | 64 | 00:10:51 | TC7 analyzed-state run |
| 817517 | `TC5_c1_cs32_48` | COMPLETED | onode13 | 64 | 00:15:26 | corrected CS32 TC5 |
| 817591 | `TC2_c1_cs32_48` | COMPLETED | onode13 | 64 | 00:16:56 | TC2 CS32 alpha=0 |
| 817596 | `TC2_c1_latlon_48` | COMPLETED | onode13 | 64 | 00:26:21 | TC2 lat-lon alpha=0 |
| 817609 | `TC2_c2_rotcor_12d` | COMPLETED | onode13 | 64 | 01:20:06 | TC2 alpha=0.05 rotated Coriolis |
| 817611 | `TC2_c3_rotcor_12d` | COMPLETED | onode14 | 64 | 07:46:19 | TC2 alpha=1.52 rotated Coriolis |
| 817643 | `TC2_c4_rotcor_12d05` | COMPLETED | onode13 | 64 | 15:11:20 | TC2 alpha=pi/2, dt=0.5 s |
| 817526 | `TC3_c1_cs32_48` | COMPLETED | onode13 | 64 | 00:08:31 | TC3 CS32 alpha=0 attempt |
| 818978 | `TC3_c2_smoke_30r` | COMPLETED | onode13 | 30 | 00:39:43 | pre-fix smoke: finite but not steady-valid |
| 818979 | `TC3_c1_smoke_30r` | COMPLETED | onode13 | 30 | 00:06:18 | alpha=0 smoke |
| 818981 | `TC3_c1_5d_45r` | COMPLETED | onode13 | 45 | 00:19:42 | full 5-day alpha=0 lat-lon run |
| 818988 | `TC3_c2_rotcor_5d_45r` | COMPLETED | onode13 | 45 | 08:07:39 | full 5-day alpha=1.0472 rotated-Coriolis lat-lon run |
| 819272 | `TC3_c2_cs32_o13` | PENDING (Resources) as of 2026-07-08 12:09 | requested onode13 | 48 | 00:00:00 | CS32 alpha=1.0472 rerun on large partition |

Job 818988 completed the 45-rank lat-lon rotated-Coriolis C2 run on `onode13`
with `deltaT=0.5` s.  The follow-up CS32 C2 submission is job 819272
(`TC3_c2_cs32_o13`), pinned to `onode13` on the `large` partition with 48 MPI
ranks and `deltaT=0.75` s; at submit time it is pending for resources.

## 10. Earlier vortexSphere Reference and Poisson Work

Two earlier directories document the transition from prototype mathematics to
MITgcm-native experiments:

### `vortexSphere_mitgcm_referenceCase`

This directory is a one-layer, global, latitude-longitude MITgcm reference
experiment.  It uses the native MITgcm implicit free-surface solve, a
$1440\times720$ grid, one vertical layer, and generated Gaussian initial
conditions.  The current generator includes:

- cosine bump free-surface adjustment,
- geostrophic Gaussian balanced eddy,
- Gaussian free-surface patch with zero initial velocity.

The active generator selection in `input/gendata_ref.py` writes the Gaussian
patch:

$$
\eta = A\exp\left[-\frac{1}{2}\left(\frac{d}{\sigma}\right)^2\right],
\qquad
A=0.1\ {\rm m},\quad
\sigma=10^\circ.
$$

### `vortexSphere_mitgcm_Poisson`

This directory contains an earlier screened/unscreened spherical Poisson
prototype and bathymetry conversion tools.  It is useful as derivation and
debugging context for the free-surface elliptic solve, but it is not the final
validated experiment pipeline.  The final Williamson suite uses MITgcm's native
pressure/free-surface operator rather than a custom inserted Poisson solver.

## 11. Reproducibility Checklist

Use this sequence for a clean scientific rerun:

1. Activate the MITgcm environment:

   ```bash
   module purge
   module load gcc/10.1.0
   module load mpi/openmpi/4.1.4-slurm-18.08.6
   module load python/3.12
   source /home/hmh85/scratch/MITgcm/.venv/bin/activate
   ```

2. Submit the appropriate `jobs/large/job_*.slurm`.
3. Confirm that the run-local `data` file contains the submitted `deltaT`,
   `nTimeSteps`, grid setting, and Coriolis setting.
4. Run the case-local setup preflight when available:

   ```bash
   python tools/check_setup_tcX.py --run-dir RUN_DIR --alpha ALPHA
   ```

5. Run the case-local CFL audit:

   ```bash
   python tools/check_cfl_tcX.py --all
   ```

6. Postprocess with the case script:

   ```bash
   TCX_RUN_DIRS=RUN_DIR python Sandbox/Scripts/TCX.py
   ```

7. Regenerate website fragments and open `docs/index.html`.

## 12. Interpretation Notes

The website status flags are not identical to mathematical proof of exact
conservation.  They mean that the current run products pass the intended
case-specific health checks:

- exact-solution cases use error norms and visual drift,
- non-exact tests use finite saved fields, coherent dynamics, and diagnostics,
- conservation tables separate mass preservation from energy/PV drift,
- failed/cancelled Slurm jobs remain useful debugging history but are not final
  validation rows.

For TC3 C2 lat-lon specifically, a steady output requires both of the
following:

1. the rotated analytic initial state,
2. the matching rotated Coriolis map through `selectCoriMap=3`.

Without the second condition the run is a smoke/no-crash test only; it is not a
balanced Williamson TC3 validation.

For TC3 C2 CS32, the corresponding setup requirement is native cubed-sphere
input staging: `data.exch2`, `useCubedSphereExchange=.TRUE.`, compact CS32
grid files, and no lat-lon `selectCoriMap` override.
