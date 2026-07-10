#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent.parent.parent / "Scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import williamson_fields as wf
import williamson_cube as wc

NX = wf.NX
NY = wf.NY
R_EARTH = wf.R_EARTH
G = wf.G
H0 = wf.H0_SOLID
ALPHA_RAD = 1.0471975511965976

BATHY_FILE = "bathymetry_flat_tc3.bin"
ETA_FILE = "eta_init.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"
FCORI_C_FILE = "fCoriC.bin"
FCORI_G_FILE = "fCoriG.bin"
FCORI_COS_FILE = "fCorCs.bin"


def make_bathymetry() -> object:
    return wf.flat_bathymetry(H0)


def grid_corners() -> tuple[object, object]:
    lon = np.arange(NX, dtype=np.float64) * wf.DLON_RAD
    lat = -0.5 * np.pi + np.arange(NY, dtype=np.float64) * wf.DLAT_RAD
    return np.meshgrid(lon, lat, indexing="xy")


def rotated_mu(lon_rad: object, lat_rad: object, alpha_rad: float) -> object:
    return (
        -np.cos(lon_rad) * np.cos(lat_rad) * np.sin(alpha_rad)
        + np.sin(lat_rad) * np.cos(alpha_rad)
    )


def write_rotated_coriolis(base_dir: Path, alpha_rad: float = ALPHA_RAD) -> None:
    lon_c, lat_c = wf.grid_cell_centers()
    lon_g, lat_g = grid_corners()
    mu_c = rotated_mu(lon_c, lat_c, alpha_rad)
    mu_g = rotated_mu(lon_g, lat_g, alpha_rad)

    wf.write_field(base_dir, FCORI_C_FILE, 2.0 * wf.OMEGA * mu_c)
    wf.write_field(base_dir, FCORI_G_FILE, 2.0 * wf.OMEGA * mu_g)
    wf.write_field(
        base_dir,
        FCORI_COS_FILE,
        2.0 * wf.OMEGA * np.sqrt(np.maximum(0.0, 1.0 - mu_c * mu_c)),
    )


def make_eta(alpha_rad: float = ALPHA_RAD) -> object:
    lon, lat = wf.grid_cell_centers()
    return wf.tc3_eta(lon, lat, alpha_rad=alpha_rad)


def make_velocity_fields(alpha_rad: float = ALPHA_RAD) -> tuple[object, object]:
    return wf.tc3_velocity_fields(alpha_rad=alpha_rad)


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    if not (base_dir / ".use_latlon_grid").exists():
        wc.generate_tc3(base_dir, alpha_rad=ALPHA_RAD)
        print("")
        print("Done.")
        print(f"ALPHA_RAD = {ALPHA_RAD:.16g}")
        print(f"ALPHA_DEG = {np.rad2deg(ALPHA_RAD):.16g}")
        print(f"H0        = {H0:.12f} m")
        print(f"U0        = {wf.TC3_U0:.12f} m/s")
        print(f"XE        = {wf.TC3_XE:.12f}")
        return

    bathy = make_bathymetry()
    write_rotated_coriolis(base_dir, alpha_rad=ALPHA_RAD)
    eta = make_eta(alpha_rad=ALPHA_RAD)
    u, v = make_velocity_fields(alpha_rad=ALPHA_RAD)

    wf.write_field(base_dir, BATHY_FILE, bathy)
    wf.write_field(base_dir, ETA_FILE, eta)
    wf.write_field(base_dir, U_FILE, u)
    wf.write_field(base_dir, V_FILE, v)
    print("")
    print("Done.")
    print(f"ALPHA_RAD = {ALPHA_RAD:.16g}")
    print(f"ALPHA_DEG = {np.rad2deg(ALPHA_RAD):.16g}")
    print(f"H0        = {H0:.12f} m")
    print(f"U0        = {wf.TC3_U0:.12f} m/s")
    print(f"XE        = {wf.TC3_XE:.12f}")


if __name__ == "__main__":
    main()
