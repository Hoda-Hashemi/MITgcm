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


def make_bathymetry() -> object:
    return wf.flat_bathymetry(H0)


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
