#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np

NX = 1440
NY = 720

A = 6.37122e6
OMEGA = 7.292e-5
G = 9.80616

THETA_B = -np.pi / 6.0
THETA_E =  np.pi / 2.0
X_E = 0.3

U0 = 2.0 * np.pi * A / (12.0 * 24.0 * 3600.0)

H0 = 4000.0

BATHY_FILE = "bathymetry_tc3.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"
ETA_FILE = "eta_init.bin"

def write_field(path: Path, arr: np.ndarray):
    arr.astype(">f4").tofile(path)

def bfun(x):
    out = np.zeros_like(x)
    m = x > 0.0
    out[m] = np.exp(-1.0 / x[m])
    return out

def build_hprime_table():
    th = np.linspace(-0.5*np.pi, 0.5*np.pi, 20001)
    x = X_E * (th - THETA_B) / (THETA_E - THETA_B)
    up = U0 * bfun(x) * bfun(X_E - x) * np.exp(4.0 / X_E)
    f = 2.0 * OMEGA * np.sin(th)
    rhs = -(A / G) * (f * up + (up**2 * np.tan(th)) / A)

    hp = np.zeros_like(th)
    dth = th[1] - th[0]
    hp[1:] = np.cumsum(0.5 * (rhs[:-1] + rhs[1:]) * dth)
    hp -= hp.mean()
    return th, hp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, required=True)
    args = ap.parse_args()
    alpha = args.alpha

    base = Path(__file__).resolve().parent

    # tracer cell centers
    lon = (np.arange(NX) + 0.5) * (2.0 * np.pi / NX)
    lat = -0.5 * np.pi + (np.arange(NY) + 0.5) * (np.pi / NY)
    lon2d, lat2d = np.meshgrid(lon, lat, indexing="xy")

    # unprimed Cartesian position
    x = np.cos(lat2d) * np.cos(lon2d)
    y = np.cos(lat2d) * np.sin(lon2d)
    z = np.sin(lat2d)

    c = np.cos(alpha)
    s = np.sin(alpha)

    # rotated coordinates: primed Cartesian
    xp = c * x + s * z
    yp = y
    zp = -s * x + c * z

    theta_p = np.arcsin(np.clip(zp, -1.0, 1.0))
    lambda_p = np.arctan2(yp, xp)

    # compact-support zonal jet u'(theta'), v'=0
    xx = X_E * (theta_p - THETA_B) / (THETA_E - THETA_B)
    up = U0 * bfun(xx) * bfun(X_E - xx) * np.exp(4.0 / X_E)

    # primed east unit vector
    elon_p_x = -np.sin(lambda_p)
    elon_p_y =  np.cos(lambda_p)
    elon_p_z =  0.0 * lambda_p

    # velocity vector in primed Cartesian
    vp_x = up * elon_p_x
    vp_y = up * elon_p_y
    vp_z = up * elon_p_z

    # rotate velocity back to unprimed Cartesian
    vx = c * vp_x - s * vp_z
    vy = vp_y
    vz = s * vp_x + c * vp_z

    # unprimed local east/north unit vectors
    elon_x = -np.sin(lon2d)
    elon_y =  np.cos(lon2d)
    elon_z =  0.0 * lon2d

    elat_x = -np.sin(lat2d) * np.cos(lon2d)
    elat_y = -np.sin(lat2d) * np.sin(lon2d)
    elat_z =  np.cos(lat2d)

    u = vx * elon_x + vy * elon_y + vz * elon_z
    v = vx * elat_x + vy * elat_y + vz * elat_z

    # h'(theta') from steady balance ODE, then eta = h' anomaly
    th_tab, hp_tab = build_hprime_table()
    eta = np.interp(theta_p.ravel(), th_tab, hp_tab).reshape(NY, NX)
    eta -= eta.mean()

    bathy = -H0 * np.ones((NY, NX), dtype=np.float64)

    write_field(base / BATHY_FILE, bathy)
    write_field(base / U_FILE, u)
    write_field(base / V_FILE, v)
    write_field(base / ETA_FILE, eta)

    print(f"alpha = {alpha:.16f}")
    print(f"u:   min={u.min():.6e}, max={u.max():.6e}")
    print(f"v:   min={v.min():.6e}, max={v.max():.6e}")
    print(f"eta: min={eta.min():.6e}, max={eta.max():.6e}")

if __name__ == "__main__":
    main()
#%%