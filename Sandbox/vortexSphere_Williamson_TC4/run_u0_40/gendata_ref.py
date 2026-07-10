#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import os
import sys

SCRIPTS_DIR = Path(__file__).resolve().parent.parent.parent / "Scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import williamson_fields as wf

NX = wf.NX
NY = wf.NY
R_EARTH = wf.R_EARTH
G = wf.G
H0 = wf.TC4_H0
U0 = float(os.environ.get("TC4_U0_VALUE", wf.TC4_U0))

BATHY_FILE = "bathymetry_flat_tc4.bin"
ETA_FILE = "eta_init.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"


def make_bathymetry() -> object:
    return wf.flat_bathymetry(H0)


def make_eta() -> object:
    lon, lat = wf.grid_cell_centers()
    return wf.tc4_eta(lon, lat, time_s=0.0, u0=U0)


def make_velocity_fields() -> tuple[object, object]:
    lon_u, lat_u = wf.grid_u_faces()
    lon_v, lat_v = wf.grid_v_faces()
    u = wf.tc4_u(lon_u, lat_u, time_s=0.0, u0=U0)
    v = wf.tc4_v(lon_v, lat_v, time_s=0.0, u0=U0)
    v[0, :] = 0.0
    v[-1, :] = 0.0
    return u, v


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    bathy = make_bathymetry()
    eta = make_eta()
    u, v = make_velocity_fields()

    wf.write_field(base_dir, BATHY_FILE, bathy)
    wf.write_field(base_dir, ETA_FILE, eta)
    wf.write_field(base_dir, U_FILE, u)
    wf.write_field(base_dir, V_FILE, v)
    print("")
    print("Done.")
    print("Williamson TC4 translating-low initial state generated.")
    print(f"U0        = {U0:.12f} m/s")
    print(f"H0        = {H0:.12f} m")


if __name__ == "__main__":
    main()
