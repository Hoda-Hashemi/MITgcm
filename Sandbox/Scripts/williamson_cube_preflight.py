#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

import williamson_fields as wf
from williamson_cube import CS32_NC, GRID_FILES, CubeGrid, resolve_cs32_grid_dir


FIELD_NAMES = {
    "TC2": ("bathymetry_flat4000.bin", "eta_init.bin", "u_init.bin", "v_init.bin", "fCoriC.bin", "fCorCs.bin"),
    "TC5": ("bathymetry_mountain_tc5.bin", "eta_init.bin", "u_init.bin", "v_init.bin"),
    "TC7": ("bathymetry_tc7.bin", "eta_init.bin", "u_init.bin", "v_init.bin"),
}


def read_text_value(path: Path, key: str) -> str | None:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*([^,!#]+)", re.I)
    for line in path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            return match.group(1).strip()
    return None


def read_compact(run_dir: Path, name: str) -> np.ndarray:
    path = run_dir / name
    raw = np.fromfile(path, dtype=">f4")
    expected = 6 * CS32_NC * CS32_NC
    if raw.size != expected:
        raise ValueError(f"{name} has {raw.size} values, expected {expected}")
    compact = raw.reshape((6 * CS32_NC, CS32_NC)).astype(np.float64)
    field = np.stack(
        [compact[face * CS32_NC : (face + 1) * CS32_NC, :].T for face in range(6)]
    )
    if not np.isfinite(field).all():
        raise ValueError(f"{name} contains non-finite values")
    return field


def require_contains(path: Path, needle: str, description: str) -> None:
    if needle not in path.read_text():
        raise ValueError(f"{path} must contain {description}")


def parse_size(path: Path) -> dict[str, int]:
    params: dict[str, int] = {}
    pattern = re.compile(r"PARAMETER\s*\(\s*(sNx|sNy|nSx|nSy|nPx|nPy)\s*=\s*(\d+)\s*\)", re.I)
    for line in path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            params[match.group(1)] = int(match.group(2))
    return params


def check_case(run_dir: Path, case: str) -> list[str]:
    notes: list[str] = []
    for filename in ("data", "data.pkg", "eedata", "data.exch2"):
        if not (run_dir / filename).exists():
            raise ValueError(f"missing {filename}")
    for filename in GRID_FILES:
        if not (run_dir / filename).exists():
            raise ValueError(f"missing {filename}")
    for filename in FIELD_NAMES[case]:
        if not (run_dir / filename).exists():
            raise ValueError(f"missing {filename}")

    require_contains(run_dir / "data", "usingCurvilinearGrid = .TRUE.", "usingCurvilinearGrid = .TRUE.")
    require_contains(run_dir / "eedata", "useCubedSphereExchange=.TRUE.", "useCubedSphereExchange=.TRUE.")
    require_contains(run_dir / "data.exch2", "W2_mapIO", "EXCH2 map-IO setting")

    size_path = run_dir.parent / "code" / "SIZE.h"
    if not size_path.exists():
        size_path = run_dir.parent / "build" / "code" / "SIZE.h"
    if size_path.exists():
        params = parse_size(size_path)
        if params.get("nPx", 0) * params.get("nPy", 0) != 48:
            raise ValueError(f"SIZE.h must set nPx*nPy=48 for 48 MPI ranks, got {params}")
        if params.get("sNx") != 16 or params.get("sNy") != 8:
            raise ValueError(f"SIZE.h must use CS32 tile size 16x8, got {params}")

    grid = CubeGrid(resolve_cs32_grid_dir(run_dir))
    dt = float((read_text_value(run_dir / "data", "deltaT") or "0.0").replace("D", "E"))
    dx_min = float(min(np.min(grid.dxc), np.min(grid.dxf)))
    dy_min = float(min(np.min(grid.dyc), np.min(grid.dyf)))

    bathy = read_compact(run_dir, FIELD_NAMES[case][0])
    eta = read_compact(run_dir, "eta_init.bin")
    u = read_compact(run_dir, "u_init.bin")
    v = read_compact(run_dir, "v_init.bin")
    depth = eta - bathy
    if float(np.nanmin(depth)) <= 100.0:
        raise ValueError(f"minimum wet depth too small: {np.nanmin(depth):.6g} m")

    max_cfl = max(float(np.max(np.abs(u))) * dt / dx_min, float(np.max(np.abs(v))) * dt / dy_min)
    if max_cfl >= 0.5:
        raise ValueError(f"initial advective CFL too high for cube preflight: {max_cfl:.6g}")

    if case == "TC2":
        for face in range(1, 7):
            path = run_dir / f"fCoriG.face{face:03d}.bin"
            raw = np.fromfile(path, dtype=">f4")
            expected = (CS32_NC + 1) * (CS32_NC + 1)
            if raw.size != expected:
                raise ValueError(f"{path.name} has {raw.size} values, expected {expected}")
            if not np.isfinite(raw).all():
                raise ValueError(f"{path.name} contains non-finite values")

    if case == "TC5":
        expected_bathy = wf.tc5_bathymetry(grid.lon_c, grid.lat_c)
        bathy_error = float(np.nanmax(np.abs(bathy - expected_bathy)))
        if bathy_error > 1.0e-2:
            raise ValueError(f"TC5 bathymetry packing mismatch: max abs error {bathy_error:.6g} m")
        mountain = wf.TC5_H0 + bathy
        peak = np.unravel_index(np.nanargmax(mountain), mountain.shape)
        peak_lon = float(grid.xc[peak])
        peak_lat = float(grid.yc[peak])
        lon_error = abs(((peak_lon - 270.0 + 180.0) % 360.0) - 180.0)
        lat_error = abs(peak_lat - 30.0)
        if lon_error > 5.0 or lat_error > 5.0:
            raise ValueError(
                "TC5 mountain peak is misplaced: "
                f"lon={peak_lon:.3f}, lat={peak_lat:.3f}, "
                f"errors=({lon_error:.3f}, {lat_error:.3f}) deg"
            )

    notes.append(f"CS{CS32_NC} cube inputs finite; depth {np.nanmin(depth):.3f}..{np.nanmax(depth):.3f} m")
    notes.append(f"initial max speed {np.nanmax(np.hypot(u, v)):.6g} m/s; CFL {max_cfl:.6g}")
    if case == "TC5":
        notes.append(f"TC5 mountain peak checked near 270E, 30N")
    return notes


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Williamson cubed-sphere run inputs.")
    parser.add_argument("--case", choices=sorted(FIELD_NAMES), required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()

    try:
        notes = check_case(args.run_dir.resolve(), args.case)
    except Exception as exc:
        print(f"{args.case} cube preflight: FAIL")
        print(f"  - {exc}")
        return 1

    print(f"{args.case} cube preflight: PASS")
    for note in notes:
        print(f"  - {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
