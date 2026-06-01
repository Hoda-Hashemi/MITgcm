#!/usr/bin/env python3
"""Generate Williamson TC1 binaries for MITgcm.

This script is dependency-free on purpose so the case can be prepared on a
minimal Python installation.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from array import array
from pathlib import Path

EARTH_RADIUS = 6371e3
DAY = 86400.0
DEFAULT_NX = 1440
DEFAULT_NY = 720
DEFAULT_DEPTH = -4000.0
DEFAULT_ALPHA = 0.5 * math.pi
BELL_CENTER_LON_DEG = 270.0
BELL_CENTER_LAT_DEG = 0.0
BELL_HEIGHT_M = 1000.0
BELL_RADIUS_M = EARTH_RADIUS / 3.0

def read_alpha_from_mypackage(input_dir: Path) -> float | None:
    path = input_dir / "data.mypackage"
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8")
    match = re.search(r"myPa_param1\s*=\s*([+\-0-9.eEdD]+)", text, re.IGNORECASE)
    if match is None:
        return None
    return float(match.group(1).replace("D", "e").replace("d", "e"))

def parse_meta(meta_path: Path) -> dict[str, object]:
    text = meta_path.read_text(encoding="utf-8")
    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_vals = [int(v) for v in re.findall(r"[-+]?\d+", dim_block)]
    dims = []
    for index in range(n_dims):
        global_size, start, end = dim_vals[3 * index : 3 * index + 3]
        dims.append({"global": global_size, "start": start, "end": end, "n": end - start + 1})
    prec_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", text)
    nrecords_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", text)
    return {
        "dims": dims,
        "prec": prec_match.group(1).lower() if prec_match else "float64",
        "nrecords": int(nrecords_match.group(1)) if nrecords_match else 1,
    }

def read_mds_field(data_path: Path, meta_path: Path) -> tuple[array, int, int]:
    meta = parse_meta(meta_path)
    dims = meta["dims"]
    if len(dims) < 2:
        raise ValueError(f"{meta_path} does not describe a 2-D field")
    nx = int(dims[0]["global"])
    ny = int(dims[1]["global"])
    nvals = nx * ny

    prec = str(meta["prec"])
    typecode = "d" if "64" in prec or "real*8" in prec else "f"
    values = array(typecode)
    with open(data_path, "rb") as handle:
        values.fromfile(handle, nvals)
    if sys.byteorder != "big":
        values.byteswap()
    return values, nx, ny

def load_geometry(search_dirs: list[Path]) -> tuple[array | None, array | None, array | None, array | None, int, int]:
    for directory in search_dirs:
        xc_meta = directory / "XC.meta"
        yc_meta = directory / "YC.meta"
        if xc_meta.exists() and yc_meta.exists():
            xc, nx, ny = read_mds_field(directory / "XC.data", xc_meta)
            yc, nx2, ny2 = read_mds_field(directory / "YC.data", yc_meta)
            if nx != nx2 or ny != ny2:
                raise ValueError(f"Geometry shape mismatch in {directory}")
            xg = yg = None
            xg_meta = directory / "XG.meta"
            yg_meta = directory / "YG.meta"
            if xg_meta.exists() and yg_meta.exists():
                xg, nx3, ny3 = read_mds_field(directory / "XG.data", xg_meta)
                yg, nx4, ny4 = read_mds_field(directory / "YG.data", yg_meta)
                if nx != nx3 or ny != ny3 or nx != nx4 or ny != ny4:
                    raise ValueError(f"Velocity geometry shape mismatch in {directory}")
            return xc, yc, xg, yg, nx, ny

    return None, None, None, None, DEFAULT_NX, DEFAULT_NY

def write_binary(path: Path, values: array) -> None:
    out = array(values.typecode, values)
    if sys.byteorder != "big":
        out.byteswap()
    with open(path, "wb") as handle:
        out.tofile(handle)

def lonlat_to_unitvec(lon_deg: float, lat_deg: float) -> tuple[float, float, float]:
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    cos_lat = math.cos(lat)
    return (
        cos_lat * math.cos(lon),
        cos_lat * math.sin(lon),
        math.sin(lat),
    )

def rotate_vector(vec: tuple[float, float, float], axis: tuple[float, float, float], angle: float) -> tuple[float, float, float]:
    ax, ay, az = axis
    norm = math.sqrt(ax * ax + ay * ay + az * az)
    ax /= norm
    ay /= norm
    az /= norm

    vx, vy, vz = vec
    c = math.cos(angle)
    s = math.sin(angle)
    dot = ax * vx + ay * vy + az * vz
    cx = ay * vz - az * vy
    cy = az * vx - ax * vz
    cz = ax * vy - ay * vx
    return (
        vx * c + cx * s + ax * dot * (1.0 - c),
        vy * c + cy * s + ay * dot * (1.0 - c),
        vz * c + cz * s + az * dot * (1.0 - c),
    )

def point_lon_lat(index_i: int, index_j: int, nx: int, ny: int, kind: str) -> tuple[float, float]:
    dlon = 360.0 / nx
    dlat = 180.0 / ny
    if kind == "center":
        lon = (index_i + 0.5) * dlon
        lat = -90.0 + (index_j + 0.5) * dlat
    elif kind == "u":
        lon = index_i * dlon
        lat = -90.0 + (index_j + 0.5) * dlat
    elif kind == "v":
        lon = (index_i + 0.5) * dlon
        lat = -90.0 + index_j * dlat
    else:
        raise ValueError(kind)
    return lon, lat

def exact_eta_at(lon_deg: float, lat_deg: float, center_vec: tuple[float, float, float]) -> float:
    x, y, z = lonlat_to_unitvec(lon_deg, lat_deg)
    dot = max(-1.0, min(1.0, x * center_vec[0] + y * center_vec[1] + z * center_vec[2]))
    radius = EARTH_RADIUS * math.acos(dot)
    if radius >= BELL_RADIUS_M:
        return 0.0
    return 0.5 * BELL_HEIGHT_M * (1.0 + math.cos(math.pi * radius / BELL_RADIUS_M))

def solid_body_u(lon_deg: float, lat_deg: float, alpha: float) -> float:
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    u0 = 2.0 * math.pi * EARTH_RADIUS / (12.0 * DAY)
    return u0 * (math.cos(lat) * math.cos(alpha) + math.sin(lat) * math.cos(lon) * math.sin(alpha))

def solid_body_v(lon_deg: float, alpha: float) -> float:
    lon = math.radians(lon_deg)
    u0 = 2.0 * math.pi * EARTH_RADIUS / (12.0 * DAY)
    return -u0 * math.sin(lon) * math.sin(alpha)

def generate_regular_grid(nx: int, ny: int, alpha: float) -> tuple[array, array, array, array]:
    bathy = array("f", [DEFAULT_DEPTH]) * (nx * ny)
    eta = array("f")
    u = array("f")
    v = array("f")

    bell_center = lonlat_to_unitvec(BELL_CENTER_LON_DEG, BELL_CENTER_LAT_DEG)
    center_vec = bell_center

    for j in range(ny):
        for i in range(nx):
            lon_c, lat_c = point_lon_lat(i, j, nx, ny, "center")
            lon_u, lat_u = point_lon_lat(i, j, nx, ny, "u")
            lon_v, _ = point_lon_lat(i, j, nx, ny, "v")
            eta.append(exact_eta_at(lon_c, lat_c, center_vec))
            u.append(solid_body_u(lon_u, lat_u, alpha))
            v.append(solid_body_v(lon_v, alpha))

    return bathy, eta, u, v

def generate_from_geometry(
    xc: array,
    yc: array,
    xg: array | None,
    yg: array | None,
    nx: int,
    ny: int,
    alpha: float,
) -> tuple[array, array, array, array]:
    bathy = array("f", [DEFAULT_DEPTH]) * (nx * ny)
    eta = array("f")
    u = array("f")
    v = array("f")

    bell_center = lonlat_to_unitvec(BELL_CENTER_LON_DEG, BELL_CENTER_LAT_DEG)
    center_vec = bell_center

    for idx in range(nx * ny):
        lon_c = float(xc[idx])
        lat_c = float(yc[idx])
        lon_u = float(xg[idx]) if xg is not None else lon_c - 0.5 * (360.0 / nx)
        lat_u = lat_c
        lon_v = lon_c
        if yg is not None:
            lat_v = float(yg[idx])
        else:
            j = idx // nx
            lat_v = -90.0 + j * (180.0 / ny)
        eta.append(exact_eta_at(lon_c, lat_c, center_vec))
        u.append(solid_body_u(lon_u, lat_u, alpha))
        v.append(solid_body_v(lon_v, alpha))

    return bathy, eta, u, v

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, default=None, help="Rotation-axis angle in radians.")
    parser.add_argument(
        "--search-dir",
        action="append",
        default=[],
        help="Optional directory to search for XC/YC/XG/YG geometry files.",
    )
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    input_dir = Path(__file__).resolve().parent

    alpha = args.alpha
    if alpha is None:
        alpha = read_alpha_from_mypackage(input_dir)
    if alpha is None:
        alpha = DEFAULT_ALPHA

    search_dirs = [input_dir]
    search_dirs.extend(Path(p) for p in args.search_dir)
    search_dirs.extend([input_dir.parent / "run", input_dir.parent / "build"])
    xc, yc, xg, yg, nx, ny = load_geometry(search_dirs)

    if xc is None or yc is None:
        bathy, eta, u, v = generate_regular_grid(nx, ny, alpha)
        geometry_source = f"regular {nx}x{ny} lat-lon grid"
    else:
        bathy, eta, u, v = generate_from_geometry(xc, yc, xg, yg, nx, ny, alpha)
        geometry_source = f"geometry files for {nx}x{ny} grid"

    write_binary(input_dir / "bathymetry_flat.bin", bathy)
    write_binary(input_dir / "eta_init.bin", eta)
    write_binary(input_dir / "u_init.bin", u)
    write_binary(input_dir / "v_init.bin", v)

    print(f"TC1 generator: nx={nx}, ny={ny}, alpha={alpha:.15f} rad")
    print(f"Geometry source: {geometry_source}")
    print("Wrote bathymetry_flat.bin, eta_init.bin, u_init.bin, v_init.bin")

if __name__ == "__main__":
    main()

