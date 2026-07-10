#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys

SCRIPTS_DIR = Path(__file__).resolve().parent.parent.parent / "Scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import williamson_fields as wf
import williamson_cube as wc

NX = wf.NX
NY = wf.NY
R_EARTH = wf.R_EARTH
G = wf.G
H0 = wf.TC5_H0
ALPHA_RAD = 0.0

BATHY_FILE = "bathymetry_mountain_tc5.bin"
ETA_FILE = "eta_init.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"


def make_bathymetry() -> object:
    lon, lat = wf.grid_cell_centers()
    return wf.tc5_bathymetry(lon, lat)


def make_eta(alpha_rad: float = ALPHA_RAD) -> object:
    del alpha_rad
    lon, lat = wf.grid_cell_centers()
    # eta is the free-surface height. The mountain belongs in bathymetry,
    # so it must not be added here a second time.
    return wf.tc5_eta(lon, lat)


def make_velocity_fields(alpha_rad: float = ALPHA_RAD) -> tuple[object, object]:
    del alpha_rad
    lon_u, lat_u = wf.grid_u_faces()
    lon_v, lat_v = wf.grid_v_faces()
    u = wf.tc5_u(lon_u, lat_u)
    v = wf.tc5_v(lon_v, lat_v)
    v[0, :] = 0.0
    v[-1, :] = 0.0
    return u, v


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    if not (base_dir / ".use_latlon_grid").exists():
        wc.generate_tc5(base_dir)
        print("")
        print("Done.")
        print(f"H0              = {H0:.12f} m")
        print(f"mountain height = {wf.TC5_MOUNTAIN_HEIGHT:.12f} m")
        print(f"U0              = {wf.TC5_U0:.12f} m/s")
        return

    bathy = make_bathymetry()
    eta = make_eta()
    u, v = make_velocity_fields()

    wf.write_field(base_dir, BATHY_FILE, bathy)
    wf.write_field(base_dir, ETA_FILE, eta)
    wf.write_field(base_dir, U_FILE, u)
    wf.write_field(base_dir, V_FILE, v)
    print("")
    print("Done.")
    print(f"H0              = {H0:.12f} m")
    print(f"mountain height = {wf.TC5_MOUNTAIN_HEIGHT:.12f} m")
    print(f"U0              = {wf.TC5_U0:.12f} m/s")


if __name__ == "__main__":
    main()
