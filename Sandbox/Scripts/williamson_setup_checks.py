from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

import williamson_fields as wf

SHAPE = (wf.NY, wf.NX)
SIZE = wf.NY * wf.NX
F4_TOL = 1.0e-3
CORI_TOL = 1.0e-9


def read_f4(path: Path) -> np.ndarray:
    if not path.exists():
        raise ValueError(f"missing required file: {path}")
    field = np.fromfile(path, dtype=">f4")
    if field.size != SIZE:
        raise ValueError(f"{path} has {field.size} values, expected {SIZE}")
    field = field.reshape(SHAPE)
    if not np.isfinite(field).all():
        raise ValueError(f"{path} contains non-finite values")
    return field.astype(np.float64)


def read_data_value(data_path: Path, key: str) -> str | None:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*([^,!#]+)", re.I)
    for line in data_path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            return match.group(1).strip()
    return None


def parse_alpha(gendata_path: Path, explicit: float | None) -> float:
    if explicit is not None:
        return explicit
    text = gendata_path.read_text()
    match = re.search(r"^ALPHA_RAD\s*=\s*([^\n#]+)", text, re.M)
    if not match:
        return 0.0
    return float(match.group(1).strip())


def max_abs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.nanmax(np.abs(a - b)))


def assert_close(name: str, actual: np.ndarray, expected: np.ndarray, tol: float, failures: list[str]) -> None:
    err = max_abs(actual, expected)
    if err > tol:
        failures.append(f"{name} max abs error {err:.6g} exceeds tolerance {tol:g}")


def grid_corners() -> tuple[np.ndarray, np.ndarray]:
    lon = np.arange(wf.NX, dtype=np.float64) * wf.DLON_RAD
    lat = -0.5 * np.pi + np.arange(wf.NY, dtype=np.float64) * wf.DLAT_RAD
    return np.meshgrid(lon, lat, indexing="xy")


def rotated_mu(lon: np.ndarray, lat: np.ndarray, alpha: float) -> np.ndarray:
    return -np.cos(lon) * np.cos(lat) * np.sin(alpha) + np.sin(lat) * np.cos(alpha)


def check_rotated_coriolis(run_dir: Path, alpha: float, failures: list[str]) -> None:
    lon_c, lat_c = wf.grid_cell_centers()
    lon_g, lat_g = grid_corners()
    mu_c = rotated_mu(lon_c, lat_c, alpha)
    mu_g = rotated_mu(lon_g, lat_g, alpha)

    expected_c = 2.0 * wf.OMEGA * mu_c
    expected_g = 2.0 * wf.OMEGA * mu_g
    expected_cos = 2.0 * wf.OMEGA * np.sqrt(np.maximum(0.0, 1.0 - mu_c * mu_c))

    assert_close("fCoriC.bin", read_f4(run_dir / "fCoriC.bin"), expected_c, CORI_TOL, failures)
    assert_close("fCoriG.bin", read_f4(run_dir / "fCoriG.bin"), expected_g, CORI_TOL, failures)
    assert_close("fCorCs.bin", read_f4(run_dir / "fCorCs.bin"), expected_cos, CORI_TOL, failures)


def check_tc2(run_dir: Path, alpha: float) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    notes: list[str] = []
    select_cori = read_data_value(run_dir / "data", "selectCoriMap")
    if select_cori != "3":
        failures.append("TC2 must set selectCoriMap=3 so rotated-alpha Coriolis matches the rotated fields")

    eta = read_f4(run_dir / "eta_init.bin")
    u = read_f4(run_dir / "u_init.bin")
    v = read_f4(run_dir / "v_init.bin")
    bathy = read_f4(run_dir / "bathymetry_flat4000.bin")
    check_rotated_coriolis(run_dir, alpha, failures)

    lon_c, lat_c = wf.grid_cell_centers()
    lon_u, lat_u = wf.grid_u_faces()
    lon_v, lat_v = wf.grid_v_faces()
    assert_close("eta_init.bin", eta, wf.tc2_eta(lon_c, lat_c, alpha), F4_TOL, failures)
    assert_close("u_init.bin", u, wf.tc2_u(lon_u, lat_u, alpha), F4_TOL, failures)
    expected_v = wf.tc2_v(lon_v, lat_v, alpha)
    expected_v[0, :] = 0.0
    expected_v[-1, :] = 0.0
    assert_close("v_init.bin", v, expected_v, F4_TOL, failures)
    assert_close("bathymetry_flat4000.bin", bathy, wf.flat_bathymetry(wf.H0_SOLID), F4_TOL, failures)
    notes.append(f"alpha={alpha:.12g}; eta range {eta.min():.6g} to {eta.max():.6g}")
    return failures, notes


def check_tc3(run_dir: Path, alpha: float) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    notes: list[str] = []
    select_cori = read_data_value(run_dir / "data", "selectCoriMap")
    if select_cori != "3":
        failures.append("TC3 must set selectCoriMap=3 so rotated-alpha Coriolis matches the compact jet")

    eta = read_f4(run_dir / "eta_init.bin")
    u = read_f4(run_dir / "u_init.bin")
    v = read_f4(run_dir / "v_init.bin")
    bathy = read_f4(run_dir / "bathymetry_flat_tc3.bin")
    check_rotated_coriolis(run_dir, alpha, failures)

    lon_c, lat_c = wf.grid_cell_centers()
    expected_u, expected_v = wf.tc3_velocity_fields(alpha_rad=alpha)
    assert_close("eta_init.bin", eta, wf.tc3_eta(lon_c, lat_c, alpha), F4_TOL, failures)
    assert_close("u_init.bin", u, expected_u, F4_TOL, failures)
    assert_close("v_init.bin", v, expected_v, F4_TOL, failures)
    assert_close("bathymetry_flat_tc3.bin", bathy, wf.flat_bathymetry(wf.H0_SOLID), F4_TOL, failures)
    notes.append(f"alpha={alpha:.12g}; max wind {np.hypot(u, v).max():.6g} m/s")
    return failures, notes


def check_tc5(run_dir: Path, alpha: float) -> tuple[list[str], list[str]]:
    del alpha
    failures: list[str] = []
    notes: list[str] = []
    eta = read_f4(run_dir / "eta_init.bin")
    u = read_f4(run_dir / "u_init.bin")
    v = read_f4(run_dir / "v_init.bin")
    bathy = read_f4(run_dir / "bathymetry_mountain_tc5.bin")

    lon_c, lat_c = wf.grid_cell_centers()
    lon_u, lat_u = wf.grid_u_faces()
    lon_v, lat_v = wf.grid_v_faces()
    mountain = wf.tc5_mountain(lon_c, lat_c)
    geodrop = (wf.R_EARTH * wf.OMEGA * wf.TC5_U0 + 0.5 * wf.TC5_U0**2) * np.sin(lat_c) ** 2 / wf.G
    old_wrong_eta = mountain - geodrop

    assert_close("eta_init.bin", eta, -geodrop, F4_TOL, failures)
    assert_close("bathymetry_mountain_tc5.bin", bathy, wf.tc5_bathymetry(lon_c, lat_c), F4_TOL, failures)
    assert_close("u_init.bin", u, wf.tc5_u(lon_u, lat_u), F4_TOL, failures)
    assert_close("v_init.bin", v, wf.tc5_v(lon_v, lat_v), F4_TOL, failures)
    if max_abs(eta, old_wrong_eta) <= F4_TOL:
        failures.append("eta_init.bin still contains mountain height; mountain must live only in bathymetry")
    depth = eta - bathy
    if float(depth.min()) <= 100.0:
        failures.append(f"minimum wet depth is too small: {depth.min():.6g} m")
    notes.append(
        "eta excludes mountain; depth range "
        f"{depth.min():.6g} to {depth.max():.6g} m; mountain max {mountain.max():.6g} m"
    )
    return failures, notes


def check_tc6(run_dir: Path, alpha: float) -> tuple[list[str], list[str]]:
    del alpha
    failures: list[str] = []
    notes: list[str] = []
    eta = read_f4(run_dir / "eta_init.bin")
    u = read_f4(run_dir / "u_init.bin")
    v = read_f4(run_dir / "v_init.bin")
    bathy = read_f4(run_dir / "bathymetry_flat_tc6.bin")

    lon_c, lat_c = wf.grid_cell_centers()
    lon_u, lat_u = wf.grid_u_faces()
    lon_v, lat_v = wf.grid_v_faces()
    assert_close("eta_init.bin", eta, wf.tc6_eta(lon_c, lat_c), F4_TOL, failures)
    assert_close("u_init.bin", u, wf.tc6_u(lon_u, lat_u), F4_TOL, failures)
    expected_v = wf.tc6_v(lon_v, lat_v)
    expected_v[0, :] = 0.0
    expected_v[-1, :] = 0.0
    assert_close("v_init.bin", v, expected_v, F4_TOL, failures)
    assert_close("bathymetry_flat_tc6.bin", bathy, wf.flat_bathymetry(wf.TC6_H0), F4_TOL, failures)
    depth = eta - bathy
    if float(depth.min()) <= 100.0:
        failures.append(f"minimum wet depth is too small: {depth.min():.6g} m")
    notes.append(f"depth range {depth.min():.6g} to {depth.max():.6g} m")
    return failures, notes


CHECKS = {
    "TC2": check_tc2,
    "TC3": check_tc3,
    "TC5": check_tc5,
    "TC6": check_tc6,
}


def parse_args(case_code: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Validate generated {case_code} setup files before MITgcm submission.")
    parser.add_argument("--run-dir", type=Path, required=True, help="run directory after gendata_ref.py has been executed")
    parser.add_argument("--alpha", type=float, default=None, help="rotation angle in radians, when relevant")
    return parser.parse_args()


def main(case_code: str) -> int:
    args = parse_args(case_code)
    run_dir = args.run_dir.resolve()
    if case_code not in CHECKS:
        print(f"ERROR: unsupported case {case_code}", file=sys.stderr)
        return 2
    try:
        alpha = parse_alpha(run_dir / "gendata_ref.py", args.alpha)
        failures, notes = CHECKS[case_code](run_dir, alpha)
    except Exception as exc:
        print(f"{case_code} setup preflight: FAIL")
        print(f"  - {exc}")
        return 1

    if failures:
        print(f"{case_code} setup preflight: FAIL")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print(f"{case_code} setup preflight: PASS")
    for note in notes:
        print(f"  - {note}")
    return 0
