#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

NX = 1440
NY = 720
DEG_PER_CELL = 0.25
REFERENCE_DEPTH = 8000.0
TIME_ISO = "1978-12-21T00:00:00Z"
TIME_LABEL = "19781221_0000"
PRESSURE_LEVEL_HPA = 500.0
NOAA_NCSS_BASE = "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure"

SOURCES = {
    "hgt": {
        "dataset": "hgt.1978.nc",
        "variable": "hgt",
        "filename": f"ncep_hgt_{TIME_LABEL}_500hpa.nc",
        "units": "m",
    },
    "uwnd": {
        "dataset": "uwnd.1978.nc",
        "variable": "uwnd",
        "filename": f"ncep_uwnd_{TIME_LABEL}_500hpa.nc",
        "units": "m/s",
    },
    "vwnd": {
        "dataset": "vwnd.1978.nc",
        "variable": "vwnd",
        "filename": f"ncep_vwnd_{TIME_LABEL}_500hpa.nc",
        "units": "m/s",
    },
}


def default_raw_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "input" / "raw"


def default_output_path() -> Path:
    return default_raw_dir() / "tc7_initial_conditions.npz"


def ncss_url(dataset: str, variable: str) -> str:
    query = urllib.parse.urlencode(
        {
            "var": variable,
            "disableProjSubset": "on",
            "horizStride": "1",
            "time": TIME_ISO,
            "vertCoord": f"{PRESSURE_LEVEL_HPA:g}",
            "accept": "netcdf3",
        }
    )
    return f"{NOAA_NCSS_BASE}/{dataset}?{query}"


def download(url: str, output: Path, force: bool) -> None:
    if output.exists() and not force:
        print(f"Using existing {output}")
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": "MITgcm-TC7-prep/1.0"})
    print(f"Downloading {url}")
    with urllib.request.urlopen(request, timeout=120) as response, tmp.open("wb") as handle:
        handle.write(response.read())
    tmp.replace(output)
    print(f"Wrote {output}")


def read_slice(path: Path, variable: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with netcdf_file(path, "r", mmap=False) as nc:
        lat = np.asarray(nc.variables["lat"].data, dtype=np.float64).copy()
        lon = np.asarray(nc.variables["lon"].data, dtype=np.float64).copy()
        values = np.asarray(nc.variables[variable].data, dtype=np.float64).squeeze().copy()
        missing = getattr(nc.variables[variable], "missing_value", None)

    if values.ndim != 2:
        raise ValueError(f"{path}:{variable} did not reduce to a 2-D latitude/longitude field")
    if missing is not None:
        values[np.isclose(values, float(missing))] = np.nan
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{path}:{variable} contains non-finite or missing values")
    return values, lat, lon


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


def write_metadata(output: Path, raw_files: dict[str, Path]) -> None:
    metadata = {
        "source_dataset": "NOAA PSL NCEP/NCAR Reanalysis pressure-level analysis",
        "source_is_original_williamson_reference": False,
        "source_note": (
            "Downloaded from NOAA PSL THREDDS NetCDF Subset Service. This is a practical "
            "TC7 analysis source, not a confirmed copy of the original Williamson reference archive."
        ),
        "selected_time_utc": TIME_ISO,
        "pressure_level_hpa": PRESSURE_LEVEL_HPA,
        "height_variable": "hgt, geopotential height in meters",
        "eta_definition": f"eta_m = hgt_m - {REFERENCE_DEPTH:g}",
        "raw_files": {name: str(path) for name, path in raw_files.items()},
        "download_urls": {
            name: ncss_url(str(spec["dataset"]), str(spec["variable"]))
            for name, spec in SOURCES.items()
        },
        "target_grid": {
            "shape": [NY, NX],
            "spacing_degrees": DEG_PER_CELL,
            "eta": "cell centers: lon 0.125..359.875, lat -89.875..89.875",
            "u": "U faces: lon 0.0..359.75, lat cell centers",
            "v": "V faces: lon cell centers, lat -90.0..89.75; first/last rows zeroed",
        },
    }
    output.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


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


def prepare(args: argparse.Namespace) -> Path:
    raw_dir = args.raw_dir.expanduser().resolve()
    output = args.output.expanduser().resolve()

    raw_files = {}
    for name, spec in SOURCES.items():
        raw_file = raw_dir / str(spec["filename"])
        raw_files[name] = raw_file
        download(ncss_url(str(spec["dataset"]), str(spec["variable"])), raw_file, args.force_download)

    hgt, hgt_lat, hgt_lon = read_slice(raw_files["hgt"], "hgt")
    uwnd, uwnd_lat, uwnd_lon = read_slice(raw_files["uwnd"], "uwnd")
    vwnd, vwnd_lat, vwnd_lon = read_slice(raw_files["vwnd"], "vwnd")

    grids = target_grids()
    height_m = interp_periodic_latlon(hgt, hgt_lat, hgt_lon, grids["lat_center"], grids["lon_center"])
    eta_m = height_m - REFERENCE_DEPTH
    u_m_s = interp_periodic_latlon(uwnd, uwnd_lat, uwnd_lon, grids["lat_center"], grids["lon_u"])
    v_m_s = interp_periodic_latlon(vwnd, vwnd_lat, vwnd_lon, grids["lat_v"], grids["lon_center"])
    v_m_s[0, :] = 0.0
    v_m_s[-1, :] = 0.0
    bathymetry_m = np.full((NY, NX), -REFERENCE_DEPTH, dtype=np.float64)

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output, eta_m=eta_m, u_m_s=u_m_s, v_m_s=v_m_s, bathymetry_m=bathymetry_m)
    write_metadata(output, raw_files)
    print(f"Wrote {output}")
    print(f"Wrote {output.with_suffix('.metadata.json')}")
    return output


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download NOAA/NCEP 500 hPa analysis slices and prepare TC7 initial conditions."
    )
    parser.add_argument("--raw-dir", type=Path, default=default_raw_dir())
    parser.add_argument("--output", type=Path, default=default_output_path())
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--no-validate", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    try:
        output = prepare(args)
        if not args.no_validate:
            validate_output(output)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
