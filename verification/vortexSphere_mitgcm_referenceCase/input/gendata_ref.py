#!/usr/bin/env python3
"""Generate eta_init.bin for vortexSphere_mitgcm_referenceCase.

Creates a compact cosine free-surface anomaly (meters) on the 1440x720
global 0.25-degree lat-lon grid, applied only over ocean points.
"""
#%%
from pathlib import Path
import numpy as np

NX = 1440
NY = 720
AMP_M = 10.0
LAT0_DEG = 0.0
LON0_DEG = 180.0
COSINE_RADIUS = 0.04
BATHY_FILENAME = "bathymetry.bin"

def periodic_lon_distance_deg(lon_deg: np.ndarray, lon0_deg: float) -> np.ndarray:
    d = np.abs(lon_deg - lon0_deg)
    return np.minimum(d, 360.0 - d)

def load_ocean_mask(base_dir: Path) -> np.ndarray:
    """Return ocean mask from bathymetry.bin with automatic sign convention."""
    bathy_path = base_dir / BATHY_FILENAME
    bathy = np.fromfile(bathy_path, dtype=">f4").reshape(NY, NX)
    if np.nanmin(bathy) < 0.0:
        ocean = bathy < 0.0
        sign_note = "depth < 0"
    else:
        ocean = bathy > 0.0
        sign_note = "depth > 0"
    print(f"Loaded {bathy_path} ({sign_note} is ocean).")
    return ocean

def main() -> None:
    base_dir = Path(__file__).resolve().parent
    ocean_mask = load_ocean_mask(base_dir)

    lon = (np.arange(NX, dtype=np.float64) + 0.5) * (360.0 / NX)
    lat = -90.0 + (np.arange(NY, dtype=np.float64) + 0.5) * (180.0 / NY)

    lon2d, lat2d = np.meshgrid(lon, lat, indexing="xy")
    dlon = periodic_lon_distance_deg(lon2d, LON0_DEG)
    dlat = lat2d - LAT0_DEG

    # Same compact cosine shape used by adjustment.cs-32x32x1 gendata.m:
    # h = 0.5 + 0.5*cos(pi*min(R,R0)/R0), scaled by AMP_M.
    x = 0.25 * (dlon / 360.0)
    y = 0.25 * (dlat / 180.0)
    radius = np.sqrt(x**2 + y**2)
    eta = AMP_M * (
        0.5 + 0.5 * np.cos(np.pi * np.minimum(radius, COSINE_RADIUS) / COSINE_RADIUS)
    )

    # Depth-aware initialization: no anomaly over land.
    eta = np.where(ocean_mask, eta, 0.0)

    # Explicit double-check: anomaly should not exist over land.
    tol = max(1e-10, 1e-6 * np.nanmax(np.abs(eta)))
    signal = np.abs(eta) > tol
    land_mask = ~ocean_mask
    if np.any(signal & land_mask):
        print("FAIL: eta_init has non-zero values over land.")
    else:
        print("PASS: eta_init is ocean-only (no land leakage).")

    out = base_dir / "eta_init.bin"
    eta.astype(">f4").tofile(out)
    print(f"Wrote {out} with shape={eta.shape}, dtype='>f4'")
    print(f"eta min={eta.min():.6e} m, max={eta.max():.6e} m")
    print(f"ocean cells={int(ocean_mask.sum())}, land cells={int((~ocean_mask).sum())}")

if __name__ == "__main__":
    main()

# %%
