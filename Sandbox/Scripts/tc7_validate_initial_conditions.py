#!/usr/bin/env python3
"""Validate the Williamson TC7 prepared initial-condition NPZ files."""

import argparse
import sys
from pathlib import Path

import numpy as np

NX = 1440
NY = 720
REQUIRED_ARRAYS = ("eta_m", "u_m_s", "v_m_s")
OPTIONAL_ARRAYS = ("bathymetry_m",)
DEFAULT_RELATIVE_PATH = (
    "vortexSphere_Williamson_TC7/input/raw"
)
DEFAULT_CASE_FILES = (
    "tc7_19781221_0000_initial_conditions.npz",
    "tc7_19790116_0000_initial_conditions.npz",
    "tc7_19790109_0000_initial_conditions.npz",
)


def default_input_path():
    sandbox_dir = Path(__file__).resolve().parent.parent
    return sandbox_dir / DEFAULT_RELATIVE_PATH


def describe_array(name, values):
    finite = np.isfinite(values)
    if not np.all(finite):
        bad_count = int(values.size - np.count_nonzero(finite))
        raise ValueError("%s contains %d non-finite value(s)" % (name, bad_count))

    return (
        "  %-12s shape=%s dtype=%s min=% .6e max=% .6e mean=% .6e"
        % (
            name,
            values.shape,
            values.dtype,
            float(np.min(values)),
            float(np.max(values)),
            float(np.mean(values)),
        )
    )


def validate_array(name, values):
    if values.shape != (NY, NX):
        raise ValueError("%s must have shape (%d, %d), got %s" % (name, NY, NX, values.shape))
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError("%s must be numeric, got dtype %s" % (name, values.dtype))
    return describe_array(name, values)


def collect_warnings(arrays):
    warnings = []
    v = arrays["v_m_s"]
    pole_v = max(float(np.max(np.abs(v[0, :]))), float(np.max(np.abs(v[-1, :]))))
    if pole_v > 0.0:
        warnings.append(
            "v_m_s has nonzero polar V-face rows; gendata_ref.py will zero first/last rows "
            "(max abs %.6e m/s)" % pole_v
        )

    eta_abs = float(np.max(np.abs(arrays["eta_m"])))
    if eta_abs > 10000.0:
        warnings.append(
            "eta_m has |value| > 10000 m; confirm this is free-surface anomaly, "
            "not raw geopotential or geopotential height"
        )

    for name in ("u_m_s", "v_m_s"):
        speed_abs = float(np.max(np.abs(arrays[name])))
        if speed_abs > 200.0:
            warnings.append(
                "%s has |value| > 200 m/s; confirm units are m/s" % name
            )

    bathy = arrays.get("bathymetry_m")
    if bathy is not None and float(np.max(bathy)) > 0.0:
        warnings.append(
            "bathymetry_m contains positive values; MITgcm bathymetry should be negative below sea level"
        )
    return warnings


def validate_npz(path):
    if not path.exists():
        raise FileNotFoundError(
            "%s does not exist. TC7 needs eta_m, u_m_s, and v_m_s arrays there."
            % path
        )
    metadata_path = path.with_suffix(".metadata.json")
    if not metadata_path.exists():
        raise FileNotFoundError("%s is missing matching metadata JSON" % path)

    summaries = []
    arrays = {}
    with np.load(path, allow_pickle=False) as data:
        names = set(data.files)
        missing = sorted(set(REQUIRED_ARRAYS) - names)
        if missing:
            raise ValueError("%s is missing required arrays: %s" % (path, ", ".join(missing)))

        unexpected = sorted(names - set(REQUIRED_ARRAYS) - set(OPTIONAL_ARRAYS))
        for name in REQUIRED_ARRAYS + OPTIONAL_ARRAYS:
            if name not in names:
                continue
            values = np.asarray(data[name])
            summaries.append(validate_array(name, values))
            arrays[name] = values

    warnings = collect_warnings(arrays)
    return summaries, warnings, unexpected


def expand_input_paths(paths):
    if not paths:
        paths = [default_input_path()]

    expanded = []
    for path in paths:
        path = path.expanduser()
        if path.is_dir():
            expected = [path / name for name in DEFAULT_CASE_FILES]
            if all(item.exists() for item in expected):
                expanded.extend(expected)
            else:
                expanded.extend(sorted(path.glob("tc7_*_initial_conditions.npz")))
        else:
            expanded.append(path)

    seen = set()
    unique = []
    for path in expanded:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    if not unique:
        raise FileNotFoundError("no TC7 initial-condition NPZ files found")
    return unique


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Validate TC7 prepared initial-condition NPZ files."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="NPZ file(s) or raw directory to validate; defaults to the TC7 raw input directory.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv or sys.argv[1:])

    try:
        paths = expand_input_paths(args.paths)
    except Exception as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1

    for path in paths:
        try:
            summaries, warnings, unexpected = validate_npz(path)
        except Exception as exc:
            print("ERROR: %s" % exc, file=sys.stderr)
            return 1

        print("TC7 initial-condition file OK: %s" % path)
        for summary in summaries:
            print(summary)
        if unexpected:
            print("Unexpected arrays ignored by gendata_ref.py: %s" % ", ".join(unexpected))
        for warning in warnings:
            print("WARNING: %s" % warning)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
