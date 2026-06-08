#!/usr/bin/env python3
"""Generate atmospheric TC1 passive-tracer and velocity input binaries."""

import argparse
import math
import re
import sys
from array import array
from pathlib import Path

EARTH_RADIUS = 6371e3
DAY = 86400.0
NX = 1440
NY = 720
BELL_LON = 270.0
BELL_LAT = 0.0
BELL_HEIGHT = 1000.0
BELL_RADIUS = EARTH_RADIUS / 3.0
DEFAULT_ALPHA = 0.5 * math.pi
UEQ = 2.0 * math.pi * EARTH_RADIUS / (12.0 * DAY)


def read_alpha(input_dir: Path) -> float | None:
    path = input_dir / "data.mypackage"
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8")
    match = re.search(r"myPa_param1\s*=\s*([+\-0-9.eEdD]+)", text, re.IGNORECASE)
    if match is None:
        return None
    return float(match.group(1).replace("D", "e").replace("d", "e"))


def write_binary(path: Path, values: array) -> None:
    out = array(values.typecode, values)
    if sys.byteorder != "big":
        out.byteswap()
    with open(path, "wb") as handle:
        out.tofile(handle)


def unit_vector(lon_deg: float, lat_deg: float) -> tuple[float, float, float]:
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    cos_lat = math.cos(lat)
    return (
        cos_lat * math.cos(lon),
        cos_lat * math.sin(lon),
        math.sin(lat),
    )


BELL_CENTER = unit_vector(BELL_LON, BELL_LAT)


def psi(lon_deg: float, lat_deg: float, alpha: float) -> float:
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    return UEQ * EARTH_RADIUS * (
        math.cos(alpha) * math.sin(lat)
        - math.sin(alpha) * math.cos(lat) * math.cos(lon)
    )


def bell_eta(lon_deg: float, lat_deg: float) -> float:
    x, y, z = unit_vector(lon_deg, lat_deg)
    cx, cy, cz = BELL_CENTER
    dot = max(-1.0, min(1.0, x * cx + y * cy + z * cz))
    radius = EARTH_RADIUS * math.acos(dot)
    if radius >= BELL_RADIUS:
        return 0.0
    return 0.5 * BELL_HEIGHT * (1.0 + math.cos(math.pi * radius / BELL_RADIUS))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, default=None, help="Rotation-axis angle in radians.")
    args = parser.parse_args()

    input_dir = Path(__file__).resolve().parent
    alpha = args.alpha
    if alpha is None:
        alpha = read_alpha(input_dir)
    if alpha is None:
        alpha = DEFAULT_ALPHA

    dlon = 360.0 / NX
    dlat = 180.0 / NY
    dlon_rad = math.radians(dlon)
    dy = EARTH_RADIUS * math.radians(dlat)

    salt = array("d")
    u = array("d")
    v = array("d")

    for j in range(NY):
        lat_c = -90.0 + (j + 0.5) * dlat
        lat_s = -90.0 + j * dlat
        lat_n = lat_s + dlat
        lat_v = -90.0 + j * dlat
        dx = EARTH_RADIUS * math.cos(math.radians(lat_v)) * dlon_rad

        for i in range(NX):
            lon_c = (i + 0.5) * dlon
            lon_u = i * dlon
            lon_w = i * dlon
            lon_e = (i + 1) * dlon
            if lon_e >= 360.0:
                lon_e -= 360.0

            salt.append(bell_eta(lon_c, lat_c))
            u.append((psi(lon_u, lat_n, alpha) - psi(lon_u, lat_s, alpha)) / dy)
            if abs(dx) < 1.0e-14:
                v.append(0.0)
            else:
                v.append((psi(lon_w, lat_v, alpha) - psi(lon_e, lat_v, alpha)) / dx)

    write_binary(input_dir / "S_init.bin", salt)
    write_binary(input_dir / "u_init.bin", u)
    write_binary(input_dir / "v_init.bin", v)

    print(f"TC1 generator: alpha={alpha:.15f} rad")
    print("Wrote S_init.bin, u_init.bin, v_init.bin")


if __name__ == "__main__":
    main()
