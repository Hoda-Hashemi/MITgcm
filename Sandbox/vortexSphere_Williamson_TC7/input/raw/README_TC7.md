# TC7 Static Initial Analysis File

`gendata_ref.py` expects one static initial-condition file:

```text
tc7_initial_conditions.npz
```

Place it in this directory.

Required arrays:

- `eta_m`: shape `(720, 1440)`, free-surface anomaly in meters, cell centers.
- `u_m_s`: shape `(720, 1440)`, zonal velocity in m/s, MITgcm U-face layout.
- `v_m_s`: shape `(720, 1440)`, meridional velocity in m/s, MITgcm V-face layout.

Optional array:

- `bathymetry_m`: shape `(720, 1440)`, MITgcm bathymetry in meters, negative below sea level. If omitted, the generator writes a flat `-8000 m` bathymetry.

Grid convention:

- `eta_m` centers: longitude `0.125, 0.375, ..., 359.875 deg`; latitude `-89.875, -89.625, ..., 89.875 deg`.
- `u_m_s` U faces: longitude `0.0, 0.25, ..., 359.75 deg`; latitude at cell centers.
- `v_m_s` V faces: longitude at cell centers; latitude `-90.0, -89.75, ..., 89.75 deg`. The generator zeros the first and last rows.

This is a static initial condition only. It is not a continuous wind forcing file and is not time dependent.

Validate before staging or submitting:

```bash
cd /home/hmh85/scratch/MITgcm
.venv/bin/python Sandbox/Scripts/tc7_validate_initial_conditions.py
```

The TC7 SLURM job runs the same validation before compiling.

## Provenance and Download Notes

The canonical Williamson TC7 is an analyzed 500 mb height/wind initial state
with reference-solution provenance. A direct ERA5 file is useful only if the
project accepts ERA5 as the analysis source and records the chosen analysis
time and preprocessing. Do not label an ERA5-derived file as the original
Williamson reference dataset unless that provenance has been established.

The current repository has one TC7 job, `job_tc7_c1.slurm`. For that job use
the first Williamson TC7 state:

```text
0000 GMT 21 December 1978
```

The Williamson paper also lists `0000 GMT 16 January 1979` and
`0000 GMT 9 January 1979`. Those should become separate jobs/files if we decide
to run them; do not silently replace `c1` with a different date.

If using Copernicus/ERA5 as the source for `c1`, the browser download is:

1. Log in to the Copernicus Climate Data Store.
2. Open "ERA5 hourly data on pressure levels from 1940 to present"
   (`reanalysis-era5-pressure-levels`).
3. Select product type `Reanalysis`.
4. Select variables `Geopotential`, `U component of wind`, and
   `V component of wind`.
5. Select pressure level `500 hPa`.
6. Select `1978` / `December` / `21` / `00:00`.
7. Select the full globe, 0.25 degree grid if the form exposes grid controls,
   and NetCDF output.
8. Submit the form, wait for the request to complete, and download the file as:

```text
input/raw/era5_tc7_19781221_0000_500hpa.nc
```

To convert that file into `tc7_initial_conditions.npz`, preprocess onto the
MITgcm grid conventions above. ERA5 geopotential is in m2/s2; divide by
`9.80665` to get geopotential height in meters, then set
`eta_m = geopotential_height_m - 8000.0` if using the default flat
`-8000 m` bathymetry. Winds must be in m/s and interpolated/placed on the U
and V face layouts listed above.

The scoped converter does that interpolation and writes provenance metadata:

```bash
cd /home/hmh85/scratch/MITgcm
.venv/bin/python Sandbox/vortexSphere_Williamson_TC7/tools/prepare_tc7_from_era5.py \
  Sandbox/vortexSphere_Williamson_TC7/input/raw/era5_tc7_19781221_0000_500hpa.nc
.venv/bin/python Sandbox/Scripts/tc7_validate_initial_conditions.py
```

The converter requires `xarray` and a NetCDF backend such as `netCDF4`; those
packages are not part of the current project virtual environment.

Programmatic download is blocked until the local environment has:

- a CDS personal access token in `$HOME/.cdsapirc`
- the user has accepted the dataset terms in the CDS browser
- `cdsapi>=0.7.7` installed in the Python environment

## Public NOAA/NCEP Fallback Used Locally

If Copernicus credentials are unavailable, this repository can prepare a
clearly marked NOAA/NCEP-derived TC7 input from public PSL THREDDS subset
downloads:

```bash
cd /home/hmh85/scratch/MITgcm
.venv/bin/python Sandbox/vortexSphere_Williamson_TC7/tools/prepare_tc7_from_ncep.py
```

This downloads 500 hPa `hgt`, `uwnd`, and `vwnd` slices for
`1978-12-21T00:00:00Z`, writes `tc7_initial_conditions.npz`, and stores
provenance in `tc7_initial_conditions.metadata.json`. Do not label this file
as the original Williamson reference archive unless that provenance is later
confirmed.
