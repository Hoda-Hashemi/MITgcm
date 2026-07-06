#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
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
PRESSURE_LEVEL_HPA = 500.0
NOAA_NCSS_BASE = "https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis/pressure"

CASES = {
    "c1": {
        "label": "19781221_0000",
        "date_slug": "1978_12_21",
        "time_iso": "1978-12-21T00:00:00Z",
        "title": "0000 GMT 21 December 1978",
    },
    "c2": {
        "label": "19790116_0000",
        "date_slug": "1979_01_16",
        "time_iso": "1979-01-16T00:00:00Z",
        "title": "0000 GMT 16 January 1979",
    },
    "c3": {
        "label": "19790109_0000",
        "date_slug": "1979_01_09",
        "time_iso": "1979-01-09T00:00:00Z",
        "title": "0000 GMT 9 January 1979",
    },
}

VARIABLES = {
    "hgt": {"variable": "hgt", "units": "m"},
    "uwnd": {"variable": "uwnd", "units": "m/s"},
    "vwnd": {"variable": "vwnd", "units": "m/s"},
}


def default_raw_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "input" / "raw"


def default_output_path(case: dict[str, str]) -> Path:
    return default_raw_dir() / f"tc7_{case['label']}_initial_conditions.npz"


def parse_iso_time(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def ncss_url(dataset: str, variable: str, time_iso: str) -> str:
    query = urllib.parse.urlencode(
        {
            "var": variable,
            "disableProjSubset": "on",
            "horizStride": "1",
            "time": time_iso,
            "vertCoord": f"{PRESSURE_LEVEL_HPA:g}",
            "accept": "netcdf3",
        }
    )
    return f"{NOAA_NCSS_BASE}/{dataset}?{query}"


def dataset_name(variable: str, case: dict[str, str]) -> str:
    year = case["time_iso"][:4]
    return f"{variable}.{year}.nc"


def canonical_raw_filename(variable: str, case: dict[str, str]) -> str:
    return f"ncep_{variable}_{case['label']}_500hpa.nc"


def local_candidates(raw_dir: Path, variable: str, case: dict[str, str]) -> list[Path]:
    year = case["time_iso"][:4]
    date_slug = case["date_slug"]
    label = case["label"]
    names = [
        canonical_raw_filename(variable, case),
        f"{variable}.{date_slug}.nc",
        f"{variable}.{label}.nc",
        f"{variable}_{date_slug}.nc",
        f"{variable}_{label}.nc",
    ]
    if variable == "hgt":
        names.append(f"hgt.{year}.nc")
    else:
        names.append(f"{variable}.{year}.nc")
    return [raw_dir / name for name in names]


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


def decode_ncep_times(hours: np.ndarray, units: bytes | str) -> list[datetime]:
    text = units.decode("ascii", errors="ignore") if isinstance(units, bytes) else str(units)
    if "hours since 1800-01-01" not in text:
        raise ValueError(f"unsupported time units: {text}")
    base = datetime(1800, 1, 1, tzinfo=timezone.utc)
    return [base + timedelta(hours=float(value)) for value in np.asarray(hours).reshape(-1)]


def select_time_index(path: Path, nc: netcdf_file, expected_time: datetime) -> int | None:
    if "time" not in nc.variables:
        return None
    time_var = nc.variables["time"]
    decoded = decode_ncep_times(np.asarray(time_var.data).copy(), getattr(time_var, "units", ""))
    diffs = [abs((item - expected_time).total_seconds()) for item in decoded]
    index = int(np.argmin(diffs))
    if diffs[index] > 60.0:
        found = decoded[index].strftime("%Y-%m-%dT%H:%M:%SZ") if decoded else "none"
        raise ValueError(f"{path} time mismatch: expected {expected_time:%Y-%m-%dT%H:%M:%SZ}, found {found}")
    return index


def select_level_index(path: Path, nc: netcdf_file) -> int | None:
    if "level" not in nc.variables:
        return None
    levels = np.asarray(nc.variables["level"].data, dtype=np.float64).reshape(-1)
    index = int(np.argmin(np.abs(levels - PRESSURE_LEVEL_HPA)))
    if abs(float(levels[index]) - PRESSURE_LEVEL_HPA) > 1.0e-6:
        raise ValueError(f"{path} level mismatch: expected 500 hPa, found {levels[index]:g}")
    return index


def read_slice(path: Path, variable: str, expected_time_iso: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    expected_time = parse_iso_time(expected_time_iso)
    with netcdf_file(path, "r", mmap=False) as nc:
        if variable not in nc.variables:
            raise ValueError(f"{path} does not contain variable {variable}")

        lat = np.asarray(nc.variables["lat"].data, dtype=np.float64).copy()
        lon = np.asarray(nc.variables["lon"].data, dtype=np.float64).copy()
        var = nc.variables[variable]
        values = np.asarray(var.data, dtype=np.float64).copy()
        dims = list(getattr(var, "dimensions", ()))

        time_index = select_time_index(path, nc, expected_time)
        if time_index is not None and "time" in dims:
            axis = dims.index("time")
            values = np.take(values, time_index, axis=axis)
            dims.pop(axis)

        level_index = select_level_index(path, nc)
        if level_index is not None and "level" in dims:
            axis = dims.index("level")
            values = np.take(values, level_index, axis=axis)
            dims.pop(axis)

        missing = getattr(var, "missing_value", None)
        fill_value = getattr(var, "_FillValue", None)

    values = np.squeeze(values)
    if values.ndim != 2:
        raise ValueError(f"{path}:{variable} did not reduce to a 2-D latitude/longitude field")
    for marker in (missing, fill_value):
        if marker is not None:
            values[np.isclose(values, float(marker))] = np.nan
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


def write_metadata(output: Path, case_id: str, case: dict[str, str], raw_files: dict[str, Path]) -> None:
    metadata = {
        "case_id": case_id,
        "case_title": case["title"],
        "source_dataset": "NOAA PSL NCEP/NCAR Reanalysis pressure-level analysis",
        "source_is_original_williamson_reference": False,
        "source_note": (
            "Downloaded from NOAA PSL THREDDS NetCDF Subset Service. This is a practical "
            "TC7 analysis source, not a confirmed copy of the original Williamson reference archive."
        ),
        "selected_time_utc": case["time_iso"],
        "pressure_level_hpa": PRESSURE_LEVEL_HPA,
        "height_variable": "hgt, geopotential height in meters",
        "eta_definition": f"eta_m = hgt_m - {REFERENCE_DEPTH:g}",
        "raw_files": {name: str(path) for name, path in raw_files.items()},
        "download_urls": {
            name: ncss_url(dataset_name(name, case), str(spec["variable"]), case["time_iso"])
            for name, spec in VARIABLES.items()
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


def load_or_download_slice(
    raw_dir: Path,
    variable: str,
    case: dict[str, str],
    *,
    force_download: bool,
    no_download: bool,
) -> tuple[Path, np.ndarray, np.ndarray, np.ndarray]:
    last_error = None
    for candidate in local_candidates(raw_dir, variable, case):
        if not candidate.exists():
            continue
        try:
            values, lat, lon = read_slice(candidate, variable, case["time_iso"])
        except Exception as exc:
            last_error = exc
            continue
        print(f"Using {variable} from {candidate}")
        return candidate, values, lat, lon

    if no_download:
        detail = f" Last checked error: {last_error}" if last_error else ""
        raise FileNotFoundError(f"no usable local {variable} file for {case['label']}.{detail}")

    output = raw_dir / canonical_raw_filename(variable, case)
    url = ncss_url(dataset_name(variable, case), variable, case["time_iso"])
    download(url, output, force_download or output.exists())
    values, lat, lon = read_slice(output, variable, case["time_iso"])
    return output, values, lat, lon


def prepare_case(args: argparse.Namespace, case_id: str) -> Path:
    case = CASES[case_id]
    raw_dir = args.raw_dir.expanduser().resolve()
    if args.output and args.case == "all":
        raise ValueError("--output can only be used with one --case")
    output = args.output.expanduser().resolve() if args.output else default_output_path(case).resolve()

    raw_files = {}
    hgt, hgt_lat, hgt_lon = None, None, None
    uwnd, uwnd_lat, uwnd_lon = None, None, None
    vwnd, vwnd_lat, vwnd_lon = None, None, None
    for name in ("hgt", "uwnd", "vwnd"):
        raw_file, values, lat, lon = load_or_download_slice(
            raw_dir,
            name,
            case,
            force_download=args.force_download,
            no_download=args.no_download,
        )
        raw_files[name] = raw_file
        if name == "hgt":
            hgt, hgt_lat, hgt_lon = values, lat, lon
        elif name == "uwnd":
            uwnd, uwnd_lat, uwnd_lon = values, lat, lon
        else:
            vwnd, vwnd_lat, vwnd_lon = values, lat, lon

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
    write_metadata(output, case_id, case, raw_files)
    print(f"Wrote {output}")
    print(f"Wrote {output.with_suffix('.metadata.json')}")
    return output


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare the three Williamson TC7 500 hPa NCEP analysis cases."
    )
    parser.add_argument("--case", choices=["all", *CASES], default="all")
    parser.add_argument("--raw-dir", type=Path, default=default_raw_dir())
    parser.add_argument("--output", type=Path)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--no-download", action="store_true")
    parser.add_argument("--no-validate", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    case_ids = list(CASES) if args.case == "all" else [args.case]
    try:
        for case_id in case_ids:
            print(f"\nPreparing TC7 {case_id}: {CASES[case_id]['title']}")
            output = prepare_case(args, case_id)
            if not args.no_validate:
                validate_output(output)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
