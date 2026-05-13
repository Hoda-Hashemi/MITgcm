# Vortex Sphere MITgcm Reference Case

## Abstract

`verification/vortexSphere_mitgcm_referenceCase` is a one-layer free-surface primitive-equation reference experiment in MITgcm. The case runs on a global latitude-longitude grid over realistic bathymetry, retains spherical Coriolis and metric terms, and uses MITgcm's standard implicit linear free-surface solver rather than a custom pressure equation inserted into `code/`. The directory structure and build workflow follow the standard MITgcm tutorial pattern used in [`verification/tutorial_global_oce_latlon`](../tutorial_global_oce_latlon), while the scientific configuration is reduced to a barotropic-style, one-layer reference problem.

## Mathematical Formulation

The prognostic fields are

$$
u(\lambda,\theta,t), \qquad
v(\lambda,\theta,t), \qquad
\eta(\lambda,\theta,t),
$$

where $\lambda$ is longitude, $\theta$ is latitude, and $\eta$ is the free-surface anomaly. The clean continuum interpretation of the active dynamics is

$$
\frac{\partial u}{\partial t} - fv
= -\frac{1}{R \cos\theta}\frac{\partial (g\eta)}{\partial \lambda} + \mathcal{M}_u,
$$

$$
\frac{\partial v}{\partial t} + fu
= -\frac{1}{R}\frac{\partial (g\eta)}{\partial \theta} + \mathcal{M}_v,
$$

$$
\frac{\partial \eta}{\partial t} + \nabla_h \cdot (h\mathbf{u}_h) = 0,
\qquad
\mathbf{u}_h = (u,v),
$$

with

$$
f(\theta) = 2 \Omega \sin\theta,
\qquad
\phi_s = g \eta.
$$

The active constants are

$$
R = 6.371 \times 10^6\ \mathrm{m},
\qquad
\Omega = 7.292 \times 10^{-5}\ \mathrm{s^{-1}},
\qquad
g = 9.81\ \mathrm{m\,s^{-2}},
\qquad
\rho_0 = 1000\ \mathrm{kg\,m^{-3}}.
$$

The experiment uses the linear free-surface formulation, so the continuity equation is written with the resting depth $d$ rather than the total thickness $h = d + \eta$. The case is therefore best described as a one-layer, rotating, free-surface ocean reference problem in which bathymetry enters through the resting geometry and the free surface evolves as a perturbation on top of that geometry.

Because the grid is spherical, MITgcm also retains the metric contributions, which in continuum notation may be summarized as

$$
\mathcal{M}_u \approx \frac{u v \tan\theta}{R},
\qquad
\mathcal{M}_v \approx -\frac{u^2 \tan\theta}{R}.
$$

These are carried in MITgcm through staggered finite-volume operators rather than pointwise formulas.

## Experiment Configuration

### Grid and time step

$$
N_x = 1440, \qquad
N_y = 720, \qquad
\Delta \lambda = \Delta \theta = 0.25^\circ,
$$

$$
N_r = 1, \qquad
\Delta r_1 = 11000\ \mathrm{m}, \qquad
\Delta t = 900\ \mathrm{s}.
$$

The compile-time layout in [`code/SIZE.h`](./code/SIZE.h) uses

- `Nx = 1440`, `Ny = 720`, `Nr = 1`
- `sNx = 360`, `sNy = 360`
- `nPx = 4`, `nPy = 2`

so the natural MPI decomposition is `8` ranks.

### Active and inactive physics

The runtime configuration in [`input/data`](./input/data) and the MITgcm defaults reported in `run/STDOUT.0000` imply:

- active linear free surface: `rigidLid = .FALSE.`, `implicitFreeSurface = .TRUE.`, `exactConserv = .TRUE.`
- active Coriolis and spherical geometry: `useCoriolis = .TRUE.`, `selectCoriMap = 2`, `selectMetricTerms = 1`
- implicit free-surface pressure solve: `implicSurfPress = 1`, `implicDiv2DFlow = 1`, `nonlinFreeSurf = 0`
- one active layer and no tracer evolution: `Nr = 1`, `tempStepping = .FALSE.`, `saltStepping = .FALSE.`
- no explicit momentum advection: `momAdvection = .FALSE.`
- no explicit physical viscosity or bottom drag: `viscAh = 0.`, `viscA4 = 0.`, `bottomDragLinear = 0.`, `bottomDragQuadratic = 0.`
- hydrostatic formulation: `nonHydrostatic = .FALSE.`

Scientifically, the experiment is close to a linear barotropic model, but it is not perfectly linear because the spherical metric terms remain active.

## Geometry, Bathymetry, and Partial Cells

The bathymetry file is [`input/bathymetry.bin`](./input/bathymetry.bin), using the standard MITgcm sign convention

$$
b(\lambda,\theta) = 0 \quad \text{on land},
\qquad
b(\lambda,\theta) < 0 \quad \text{in ocean}.
$$

The resting ocean depth is therefore

$$
d(\lambda,\theta) = -b(\lambda,\theta).
$$

With one vertical layer, MITgcm represents the effective local depth through partial cells:

$$
Depth_{ij} = d_{ij} = \Delta r_1 \, hFacC_{ij1}.
$$

For this case,

$$
\Delta r_1 = 11000\ \mathrm{m}.
$$

Hence,

$$
hFacC = 0 \Rightarrow \text{land},
$$

$$
0 < hFacC < 1 \Rightarrow \text{partially wet cell},
$$

$$
hFacC = 1 \Rightarrow \text{full 11000 m wet cell}.
$$

The settings

$$
hFacMin = 10^{-5},
\qquad
hFacMinDr = 0.5\ \mathrm{m}
$$

allow shallow and partially wet cells to preserve realistic bathymetric depth instead of collapsing the ocean to a single uniform layer thickness.

## MITgcm Time-Stepping and Pressure Solve

At the top level, the forward step follows

$$
\texttt{DYNAMICS}
\;\longrightarrow\;
\texttt{SOLVE\_FOR\_PRESSURE}
\;\longrightarrow\;
\texttt{MOMENTUM\_CORRECTION\_STEP}
\;\longrightarrow\;
\texttt{INTEGR\_CONTINUITY}.
$$

The momentum predictor and pressure correction can be summarized as

$$
u^{n+1} = u^* - \Delta t_{\mathrm{mom}} \partial_x \phi_s^{n+1},
$$

$$
v^{n+1} = v^* - \Delta t_{\mathrm{mom}} \partial_y \phi_s^{n+1},
$$

$$
\frac{\eta^{n+1} - \eta^n}{\Delta t_{\mathrm{fs}}} + \nabla_h \cdot \bigl(d \mathbf{u}^{n+1}\bigr) = 0.
$$

Substituting the corrected velocity into continuity yields the elliptic free-surface equation solved by MITgcm:

$$
\nabla_h \cdot \bigl(d \nabla_h \phi_s^{n+1}\bigr) - \frac{1}{g \Delta t_{\mathrm{mom}} \Delta t_{\mathrm{fs}}}\phi_s^{n+1} = \frac{1}{\Delta t_{\mathrm{mom}}}\nabla_h \cdot \bigl(d \mathbf{u}^*\bigr) - \frac{1}{g \Delta t_{\mathrm{mom}} \Delta t_{\mathrm{fs}}}\phi_s^n.
$$

MITgcm discretizes this as

$$
A \phi_s^{n+1} = b,
\qquad
\phi_s^{n+1} = g \eta^{n+1}.
$$

The matrix is assembled in [`model/src/ini_cg2d.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/ini_cg2d.F) through the coefficients `aW2d`, `aS2d`, and `aC2d`, using the face areas, partial-cell factors, and the linear free-surface storage term. The right-hand side is filled in [`model/src/solve_for_pressure.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/solve_for_pressure.F) from the divergence of the provisional transport plus the previous free-surface state.

The active solver controls from [`input/data`](./input/data) are

$$
cg2dMaxIters = 5000,
\qquad
cg2dTargetResidual = 10^{-10},
\qquad
cg2dPreCondFreq = 1.
$$

## Source-Code Landmarks

The main source files for this experiment are:

- [`model/src/forward_step.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/forward_step.F): top-level time-step driver
- [`model/src/dynamics.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/dynamics.F): explicit momentum tendency assembly
- [`model/src/solve_for_pressure.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/solve_for_pressure.F): implicit free-surface solve and right-hand-side assembly
- [`model/src/ini_cg2d.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/ini_cg2d.F): construction of the 2-D elliptic operator
- [`model/src/momentum_correction_step.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/momentum_correction_step.F): pressure-gradient velocity correction
- [`model/src/integr_continuity.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/integr_continuity.F): continuity update
- [`model/src/update_masks_etc.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/update_masks_etc.F): bathymetry to `Depth` and `hFac*`
- [`model/src/ini_linear_phisurf.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/model/src/ini_linear_phisurf.F): initialization of $\phi_s = g \eta$
- [`pkg/mom_fluxform/mom_fluxform.F`](/Users/hodahashemi/Documents/Github/Playground/MITgcm/pkg/mom_fluxform/mom_fluxform.F): active horizontal momentum package

## Directory Layout

| Path | Role |
| --- | --- |
| [`input/`](./input) | runtime namelists and binary inputs |
| [`code/`](./code) | experiment-specific compile-time overrides |
| [`build/`](./build) | generated build directory created by `genmake2` |
| [`run/`](./run) | generated execution directory containing links, executable, and model output |
| [`Scripts/`](./Scripts) | plotting and post-processing utilities |

Important files include:

- [`input/data`](./input/data): main runtime namelist
- [`input/data.pkg`](./input/data.pkg): runtime package selection
- [`input/data.diagnostics`](./input/data.diagnostics): diagnostic stream definitions
- [`input/bathymetry.bin`](./input/bathymetry.bin): resting bathymetry and land mask
- [`input/eta_init.bin`](./input/eta_init.bin): initial free-surface anomaly
- [`input/gendata_ref.py`](./input/gendata_ref.py) and [`input/gendata.m`](./input/gendata.m): input generators
- [`code/SIZE.h`](./code/SIZE.h): grid size and MPI layout
- [`code/DIAGNOSTICS_SIZE.h`](./code/DIAGNOSTICS_SIZE.h): diagnostic array sizing
- [`code/packages.conf`](./code/packages.conf): compiled package list

## Relation to the Standard MITgcm Tutorial Layout

The experiment follows the standard MITgcm verification workflow used by [`verification/tutorial_global_oce_latlon`](../tutorial_global_oce_latlon):

1. configure a local `build/` directory with `genmake2`
2. compile the executable in `build/`
3. populate a local `run/` directory with links or copies of the runtime inputs
4. execute the model from `run/`
5. inspect diagnostics and post-process the outputs

The main scientific differences from the standard tutorial configuration are:

- one active vertical layer instead of a multi-level ocean configuration
- no active temperature or salinity stepping
- no imposed forcing fields beyond bathymetry and the initial free-surface anomaly
- no explicit momentum advection
- focus on the MITgcm barotropic free-surface solve over realistic bathymetry

## Build and Run Workflow

The typical workflow is

```bash
cd verification/vortexSphere_mitgcm_referenceCase/build
../../../tools/genmake2 -mpi -mods=../code -of=../../../tools/build_options/darwin_arm64_gfortran
make depend
make -j8
```

Then populate the run directory:

```bash
cd ../run
ln -sf ../input/data .
ln -sf ../input/data.pkg .
ln -sf ../input/data.diagnostics .
ln -sf ../input/eedata .
ln -sf ../input/bathymetry.bin .
ln -sf ../input/eta_init.bin .
ln -sf ../build/mitgcmuv .
```

The executable may then be launched, for example, with

```bash
mpirun -np 8 ./mitgcmuv
```

The input generators in [`input/gendata_ref.py`](./input/gendata_ref.py) and [`input/gendata.m`](./input/gendata.m) may be used to regenerate the initial anomaly field before building or rerunning the case.

## Outputs and Diagnostics

The run directory typically contains:

- grid and geometry fields such as `Depth`, `hFacC`, `hFacW`, and `hFacS`
- prognostic fields such as `Eta`, `U`, `V`, `W`, `PH`, and `PHL`
- diagnostic bundles such as `dynDiag`, `dyn_Aux`, and `dynStDiag`
- logs such as `STDOUT.*`, `STDERR.*`, and `available_diagnostics.log`

Post-processing figures are written under [`Scripts/results/`](./Scripts/results). The plotting utilities live in [`Scripts/mitgcm_plotter.py`](./Scripts/mitgcm_plotter.py) and [`Scripts/scripts.py`](./Scripts/scripts.py).

## Repository Hygiene

`build/` and `run/` are generated working directories and are intended to remain git-ignored. Plot outputs and Python cache directories are also excluded through experiment-local `.gitignore` files so that the experiment can be documented and rerun without committing generated artifacts.

