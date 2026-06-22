#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

NX = 1440
NY = 720
DEG_PER_CELL = 0.25
DEFAULT_GRAVITY = 9.80665
DEFAULT_REFERENCE_DEPTH = 8000.0

VAR_CANDIDATES = {
    "geopotential": ("z", "geopotential"),
    "u_wind": ("u", "u_component_of_wind"),
    "v_wind": ("v", "v_component_of_wind"),
}
LAT_CANDIDATES = ("latitude", "lat")
LON_CANDIDATES = ("longitude", "lon")
LEVEL_CANDIDATES = ("pressure_level", "level", "isobaricInhPa", "plev")
TIME_CANDIDATES = ("valid_time", "time")


def default_output_path() -> Path:
    return Path(__file__).resolve().parents[1] / "input" / "raw" / "tc7_initial_conditions.npz"


def import_xarray():
    try:
        import xarray as xr  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "ERROR: converting ERA5 NetCDF needs xarray plus a NetCDF backend.\n"
            "Install in the project environment with:\n"
            "  .venv/bin/python -m pip install xarray netCDF4"
        ) from exc
    return xr


def find_dataset_name(names: list[str], candidates: tuple[str, ...], label: str) -> str:
    lowered = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    raise ValueError(f"could not find {label}; looked for {', '.join(candidates)}")


def find_optional_name(names: list[str], candidates: tuple[str, ...]) -> str | None:
    lowered = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def select_single_slice(da, pressure_level: float, time_utc: str | None):
    coord_names = list(da.coords) + list(da.dims)
    level_name = find_optional_name(coord_names, LEVEL_CANDIDATES)
    if level_name and da.sizes.get(level_name, 1) > 1:
        da = da.sel({level_name: pressure_level})

    time_name = find_optional_name(list(da.coords) + list(da.dims), TIME_CANDIDATES)
    if time_name and da.sizes.get(time_name, 1) > 1:
        if not time_utc:
            raise ValueError(
                f"{da.name} contains multiple times; pass --time UTC_ISO, "
                "for example --time 1978-12-21T00:00:00"
            )
        da = da.sel({time_name: np.datetime64(time_utc)})

    for dim in list(da.dims):
        if dim in LAT_CANDIDATES or dim in LON_CANDIDATES:
            continue
        if da.sizes.get(dim, 1) == 1:
            da = da.isel({dim: 0})

    return da.squeeze(drop=True)


def extract_lat_lon_field(ds, variable_name: str, pressure_level: float, time_utc: str | None):
    da = select_single_slice(ds[variable_name], pressure_level=pressure_level, time_utc=time_utc)
    names = list(da.coords) + list(da.dims)
    lat_name = find_dataset_name(names, LAT_CANDIDATES, f"{variable_name} latitude coordinate")
    lon_name = find_dataset_name(names, LON_CANDIDATES, f"{variable_name} longitude coordinate")
    da = da.transpose(lat_name, lon_name)
    values = np.asarray(da.values, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError(f"{variable_name} did not reduce to a 2-D latitude/longitude field")
    return (
        values,
        np.asarray(da[lat_name].values, dtype=np.float64),
        np.asarray(da[lon_name].values, dtype=np.float64),
    )


def sorted_unique_longitudes(lon: np.ndarray, field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon = np.mod(lon, 360.0)
    order = np.argsort(lon)
    lon = lon[order]
    field = field[:, order]
    keep = np.r_[True, np.diff(lon) > 1.0e-10]
    return lon[keep], field[:, keep]


def sorted_latitudes(lat: np.ndarray, field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(lat)
    return lat[order], field[order, :]


def interp_periodic_latlon(
    field: np.ndarray,
    src_lat: np.ndarray,
    src_lon: np.ndarray,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> np.ndarray:
    src_lon, field = sorted_unique_longitudes(src_lon, field)
    src_lat, field = sorted_latitudes(src_lat, field)

    lon_ext = np.concatenate(([src_lon[-1] - 360.0], src_lon, [src_lon[0] + 360.0]))
    along_lon = np.empty((src_lat.size, target_lon.size), dtype=np.float64)
    for j in range(src_lat.size):
        row_ext = np.concatenate(([field[j, -1]], field[j, :], [field[j, 0]]))
        along_lon[j, :] = np.interp(target_lon, lon_ext, row_ext)

    out = np.empty((target_lat.size, target_lon.size), dtype=np.float64)
    for i in range(target_lon.size):
        out[:, i] = np.interp(target_lat, src_lat, along_lon[:, i])
    return out


def target_grids() -> dict[str, np.ndarray]:
    return {
        "lon_center": (np.arange(NX, dtype=np.float64) + 0.5) * DEG_PER_CELL,
        "lon_u": np.arange(NX, dtype=np.float64) * DEG_PER_CELL,
        "lat_center": -90.0 + (np.arange(NY, dtype=np.float64) + 0.5) * DEG_PER_CELL,
        "lat_v": -90.0 + np.arange(NY, dtype=np.float64) * DEG_PER_CELL,
    }


def write_metadata(path: Path, args: argparse.Namespace, source_vars: dict[str, str]) -> None:
    metadata = {
        "source_file": str(args.netcdf.expanduser().resolve()),
        "source_dataset": "ERA5 hourly data on pressure levels from 1940 to present",
        "source_is_original_williamson_reference": False,
        "selected_time_utc": args.time,
        "pressure_level_hpa": args.pressure_level,
        "geopotential_conversion": f"height_m = geopotential / {args.gravity}",
        "eta_definition": f"eta_m = height_m - {args.reference_depth}",
        "target_grid": {
            "shape": [NY, NX],
            "spacing_degrees": DEG_PER_CELL,
            "eta": "cell centers: lon 0.125..359.875, lat -89.875..89.875",
            "u": "U faces: lon 0.0..359.75, lat cell centers",
            "v": "V faces: lon cell centers, lat -90.0..89.75; first/last rows zeroed",
        },
        "source_variables": source_vars,
    }
    path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def validate_output(path: Path) -> None:
    scripts_dir = Path(__file__).resolve().parents[2] / "Scripts"
    sys.path.insert(0, str(scripts_dir))
    from tc7_validate_initial_conditions import validate_npz

    summaries, warnings, unexpected = validate_npz(path)
    print(f"TC7 initial-condition file OK: {path}")
    for summary in summaries:
        print(summary)
    if unexpected:
        print("Unexpected arrays ignored by gendata_ref.py: %s" % ", ".join(unexpected))
    for warning in warnings:
        print("WARNING: %s" % warning)


def convert(args: argparse.Namespace) -> Path:
    xr = import_xarray()
    netcdf = args.netcdf.expanduser().resolve()
    output = args.output.expanduser().resolve()

    ds = xr.open_dataset(netcdf)
    try:
        data_names = list(ds.data_vars)
        z_name = find_dataset_name(data_names, VAR_CANDIDATES["geopotential"], "geopotential variable")
        u_name = find_dataset_name(data_names, VAR_CANDIDATES["u_wind"], "u-wind variable")
        v_name = find_dataset_name(data_names, VAR_CANDIDATES["v_wind"], "v-wind variable")

        z, z_lat, z_lon = extract_lat_lon_field(ds, z_name, args.pressure_level, args.time)
        u, u_lat, u_lon = extract_lat_lon_field(ds, u_name, args.pressure_level, args.time)
        v, v_lat, v_lon = extract_lat_lon_field(ds, v_name, args.pressure_level, args.time)
    finally:
        ds.close()

    grids = target_grids()
    height_m = interp_periodic_latlon(z, z_lat, z_lon, grids["lat_center"], grids["lon_center"]) / args.gravity
    eta_m = height_m - args.reference_depth
    u_m_s = interp_periodic_latlon(u, u_lat, u_lon, grids["lat_center"], grids["lon_u"])
    v_m_s = interp_periodic_latlon(v, v_lat, v_lon, grids["lat_v"], grids["lon_center"])
    v_m_s[0, :] = 0.0
    v_m_s[-1, :] = 0.0
    bathymetry_m = np.full((NY, NX), -args.reference_depth, dtype=np.float64)

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output,
        eta_m=eta_m,
        u_m_s=u_m_s,
        v_m_s=v_m_s,
        bathymetry_m=bathymetry_m,
    )
    write_metadata(output, args, {"geopotential": z_name, "u_wind": u_name, "v_wind": v_name})
    return output


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a TC7 500 hPa ERA5 NetCDF snapshot into tc7_initial_conditions.npz."
    )
    parser.add_argument("netcdf", type=Path, help="Downloaded ERA5 pressure-level NetCDF file.")
    parser.add_argument("--output", type=Path, default=default_output_path())
    parser.add_argument("--time", default="1978-12-21T00:00:00", help="UTC analysis time to select.")
    parser.add_argument("--pressure-level", type=float, default=500.0, help="Pressure level in hPa.")
    parser.add_argument("--gravity", type=float, default=DEFAULT_GRAVITY)
    parser.add_argument("--reference-depth", type=float, default=DEFAULT_REFERENCE_DEPTH)
    parser.add_argument("--no-validate", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    try:
        output = convert(args)
        print(f"Wrote {output}")
        print(f"Wrote {output.with_suffix('.metadata.json')}")
        if not args.no_validate:
            validate_output(output)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
