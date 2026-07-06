# Williamson TC7 Raw Input Files

TC7 uses analyzed 500 hPa geopotential height and horizontal wind fields.
This folder stages NOAA PSL NCEP/NCAR Reanalysis pressure-level slices as a
practical local source, not a confirmed copy of the original Williamson archive.
Source catalog:
https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/pressure/catalog.html

Prepared experiments:

| Case | Williamson state | Prepared NPZ |
| --- | --- | --- |
| `c1` | `0000 GMT 21 December 1978` | `tc7_19781221_0000_initial_conditions.npz` |
| `c2` | `0000 GMT 16 January 1979` | `tc7_19790116_0000_initial_conditions.npz` |
| `c3` | `0000 GMT 9 January 1979` | `tc7_19790109_0000_initial_conditions.npz` |

Each prepared NPZ contains:

- `eta_m = hgt_m - 8000`
- `u_m_s`
- `v_m_s`
- `bathymetry_m = -8000`

Regenerate all prepared inputs from the available local NetCDF files, downloading
any missing 500 hPa slice from NOAA PSL if needed:

```bash
../../tools/prepare_tc7_from_ncep.py --case all
```

The job scripts select the case explicitly through `TC7_RAW_FILE`.
