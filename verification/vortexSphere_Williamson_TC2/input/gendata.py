#!/usr/bin/env python3
#%%
from pathlib import Path
import argparse
import numpy as np

NX = 1440
NY = 720

A = 6.37122e6
OMEGA = 7.292e-5
G = 9.80616

U0 = 2.0 * np.pi * A / (12.0 * 24.0 * 3600.0)
GH0 = 2.94e4
H0 = GH0 / G

BATHY_FILE = "bathymetry_tc2.bin"
ETA_FILE = "eta_init.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"

def cell_centers():
    lon = (np.arange(NX) + 0.5) * (2.0 * np.pi / NX)
    lat = -0.5 * np.pi + (np.arange(NY) + 0.5) * (np.pi / NY)
    return np.meshgrid(lon, lat, indexing="xy")

def u_faces():
    lon = np.arange(NX) * (2.0 * np.pi / NX)
    lat = -0.5 * np.pi + (np.arange(NY) + 0.5) * (np.pi / NY)
    return np.meshgrid(lon, lat, indexing="xy")

def v_faces():
    lon = (np.arange(NX) + 0.5) * (2.0 * np.pi / NX)
    lat = -0.5 * np.pi + np.arange(NY) * (np.pi / NY)
    return np.meshgrid(lon, lat, indexing="xy")

def write_field(path: Path, arr: np.ndarray):
    arr.astype(">f4").tofile(path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, required=True)
    args = ap.parse_args()

    alpha = args.alpha
    base = Path(__file__).resolve().parent

    #! flat bottom
    bathy = -H0 * np.ones((NY, NX), dtype=np.float64)
    write_field(base / BATHY_FILE, bathy)

    #! eta at tracer-cell centers
    lon_c, lat_c = cell_centers()
    q_c = -np.cos(lon_c) * np.cos(lat_c) * np.sin(alpha) + np.sin(lat_c) * np.cos(alpha)
    gh = GH0 - (A * OMEGA * U0 + 0.5 * U0**2) * q_c**2
    h = gh / G
    eta = h - H0

    #! u at west/east faces
    lon_u, lat_u = u_faces()
    u = U0 * (np.cos(lat_u) * np.cos(alpha) + np.cos(lon_u) * np.sin(lat_u) * np.sin(alpha))

    #! v at south/north faces
    lon_v, lat_v = v_faces()
    v = -U0 * np.sin(lon_v) * np.sin(alpha)

    write_field(base / ETA_FILE, eta)
    write_field(base / U_FILE, u)
    write_field(base / V_FILE, v)

    print(f"alpha = {alpha:.16f}")
    print(f"U0    = {U0:.8f} m/s")
    print(f"H0    = {H0:.8f} m")
    print(f"eta   : min={eta.min():.8e}, max={eta.max():.8e}")
    print(f"u     : min={u.min():.8e}, max={u.max():.8e}")
    print(f"v     : min={v.min():.8e}, max={v.max():.8e}")

if __name__ == "__main__":
    main()
#%%