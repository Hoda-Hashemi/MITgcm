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
H0 = wf.TC7_REFERENCE_DEPTH
ALPHA_RAD = 0.0

BATHY_FILE = "bathymetry_tc7.bin"
ETA_FILE = "eta_init.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"
POLAR_TAPER_START_DEG = 85.0


def load_fields() -> tuple[object, object, object, object]:
    return wf.load_tc7_fields(Path(__file__).resolve().parent)


def polar_velocity_taper() -> object:
    lat = -0.5 * np.pi + (np.arange(NY, dtype=np.float64) + 0.5) * (np.pi / NY)
    abs_lat_deg = np.abs(np.rad2deg(lat))
    taper = np.ones(NY, dtype=np.float64)
    cap = abs_lat_deg > POLAR_TAPER_START_DEG
    x = (abs_lat_deg[cap] - POLAR_TAPER_START_DEG) / (90.0 - POLAR_TAPER_START_DEG)
    taper[cap] = np.cos(0.5 * np.pi * np.clip(x, 0.0, 1.0)) ** 2
    return taper[:, None]


def apply_polar_velocity_taper(u: object, v: object) -> tuple[object, object]:
    taper = polar_velocity_taper()
    u_tapered = np.asarray(u, dtype=np.float64) * taper
    v_tapered = np.asarray(v, dtype=np.float64) * taper
    v_tapered[0, :] = 0.0
    v_tapered[-1, :] = 0.0
    return u_tapered, v_tapered


def make_bathymetry() -> object:
    _eta, _u, _v, bathy = load_fields()
    return bathy


def make_eta(alpha_rad: float = ALPHA_RAD) -> object:
    del alpha_rad
    eta, u, v, _bathy = load_fields()
    eta, _u, _v = wf.preprocess_tc7_analysis_fields(eta, u, v)
    return eta


def make_velocity_fields(alpha_rad: float = ALPHA_RAD) -> tuple[object, object]:
    del alpha_rad
    eta, u, v, _bathy = load_fields()
    _eta, u, v = wf.preprocess_tc7_analysis_fields(eta, u, v)
    return apply_polar_velocity_taper(u, v)


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    raw_file = wf.resolve_tc7_raw_file(base_dir)
    if not (base_dir / ".use_latlon_grid").exists():
        wc.generate_tc7(base_dir)
        print("")
        print("Done.")
        print(f"TC7 input source: {raw_file}")
        print(
            "Applied TC7 large-scale filter: "
            f"zonal wavenumber <= {wf.TC7_LONGITUDE_WAVENUMBER_CUTOFF}, "
            f"{wf.TC7_SHAPIRO_SMOOTHING_PASSES} Shapiro passes"
        )
        return

    eta, u, v, bathy = load_fields()
    eta, u, v = wf.preprocess_tc7_analysis_fields(eta, u, v)
    u, v = apply_polar_velocity_taper(u, v)
    wf.write_field(base_dir, BATHY_FILE, bathy)
    wf.write_field(base_dir, ETA_FILE, eta)
    wf.write_field(base_dir, U_FILE, u)
    wf.write_field(base_dir, V_FILE, v)
    print("")
    print("Done.")
    print(f"TC7 input source: {raw_file}")
    print(
        "Applied TC7 large-scale filter: "
        f"zonal wavenumber <= {wf.TC7_LONGITUDE_WAVENUMBER_CUTOFF}, "
        f"{wf.TC7_SHAPIRO_SMOOTHING_PASSES} Shapiro passes"
    )
    print(f"Applied polar velocity taper poleward of {POLAR_TAPER_START_DEG:g} deg")


if __name__ == "__main__":
    main()
