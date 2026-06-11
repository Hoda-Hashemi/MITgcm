#%%
from pathlib import Path
import numpy as np

# Grid / constants
NX = 1440
NY = 720
DLON_DEG = 360.0 / NX
DLAT_DEG = 180.0 / NY

R_EARTH = 6.371e6
OMEGA = 7.292e-5
G = 9.81

BATHY_FILENAME = "bathymetry.bin"

#!added for constant depth:
CONSTANT_DEPTH = -4000.0
BATHY_FILENAME = "bathymetry_flat4000.bin"
def write_flat_bathymetry(base_dir: Path):
    bathy = np.full((NY, NX), CONSTANT_DEPTH, dtype=">f4")
    bathy.tofile(base_dir / BATHY_FILENAME)
    print(f"Wrote {BATHY_FILENAME} with constant depth {CONSTANT_DEPTH} m.")

## Shared tools for both experiments:
def load_ocean_mask(base_dir: Path) -> np.ndarray:
    bathy = np.fromfile(base_dir / BATHY_FILENAME, dtype=">f4").reshape(NY, NX)

    if np.nanmin(bathy) < 0.0:
        ocean = bathy < 0.0
        print("Loaded bathymetry: ocean = depth < 0")
    else:
        ocean = bathy > 0.0
        print("Loaded bathymetry: ocean = depth > 0")

    return ocean

def grid_centers():
    lon = (np.arange(NX, dtype=np.float64) + 0.5) * DLON_DEG
    lat = -90.0 + (np.arange(NY, dtype=np.float64) + 0.5) * DLAT_DEG
    return np.meshgrid(lon, lat, indexing="xy")

def periodic_lon_distance_deg(lon_deg: np.ndarray, lon0_deg: float) -> np.ndarray:
    d = np.abs(lon_deg - lon0_deg)
    return np.minimum(d, 360.0 - d)

def spherical_distance_rad(lon_deg, lat_deg, lon0_deg, lat0_deg):
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    lon0 = np.deg2rad(lon0_deg)
    lat0 = np.deg2rad(lat0_deg)

    dlon = (lon - lon0 + np.pi) % (2.0 * np.pi) - np.pi

    cos_r = (
        np.sin(lat) * np.sin(lat0)
        + np.cos(lat) * np.cos(lat0) * np.cos(dlon)
    )

    return np.arccos(np.clip(cos_r, -1.0, 1.0))

def make_face_masks(ocean_mask: np.ndarray):
    u_mask = ocean_mask & np.roll(ocean_mask, 1, axis=1)

    v_mask = ocean_mask & np.roll(ocean_mask, 1, axis=0)
    v_mask[0, :] = False
    v_mask[-1, :] = False

    return u_mask, v_mask

def write_field(base_dir: Path, name: str, field: np.ndarray):
    field.astype(">f4").tofile(base_dir / name)
    print(f"Wrote {name}: min={field.min():.6e}, max={field.max():.6e}")

#! Experiment 1: Cosine bump -> gravity-wave adjustment
def make_cosine_bump(ocean_mask: np.ndarray):
    AMP_M = 10.0
    LAT0_DEG = 0.0
    LON0_DEG = 180.0
    COSINE_RADIUS = 0.04

    lon2d, lat2d = grid_centers()

    dlon = periodic_lon_distance_deg(lon2d, LON0_DEG)
    dlat = lat2d - LAT0_DEG

    x = 0.25 * (dlon / 360.0)
    y = 0.25 * (dlat / 180.0)
    radius = np.sqrt(x**2 + y**2)

    eta = AMP_M * (
        0.5
        + 0.5 * np.cos(np.pi * np.minimum(radius, COSINE_RADIUS) / COSINE_RADIUS)
    )

    eta = np.where(ocean_mask, eta, 0.0)

    u = np.zeros_like(eta)
    v = np.zeros_like(eta)

    return eta, u, v

#! Experiment 3: Gaussian Patch -> free-surface adjustment
def exp3_Gaussian_Patch(ocean_mask: np.ndarray):
    AMP_M = 0.1
    LAT0_DEG = 0.0
    LON0_DEG = 180.0
    SIGMA_DEG = 10.0
    SIGMA_M = R_EARTH * np.deg2rad(SIGMA_DEG)

    lon2d, lat2d = grid_centers()

    distance_m = R_EARTH * spherical_distance_rad(lon2d, lat2d, LON0_DEG, LAT0_DEG)
    eta = AMP_M * np.exp(-0.5 * (distance_m / SIGMA_M) ** 2)

    eta = np.where(ocean_mask, eta, 0.0)

    u = np.zeros_like(eta)
    v = np.zeros_like(eta)

    return eta, u, v

#! Experiment 2: Geostrophic Gaussian -> balanced eddy
def make_geostrophic_gaussian(ocean_mask: np.ndarray):
    AMP_M = 10.0
    LAT0_DEG = 30.0
    LON0_DEG = 180.0
    SIGMA_RAD = 0.04

    lon2d, lat2d = grid_centers()
    u_mask, v_mask = make_face_masks(ocean_mask)

    r = spherical_distance_rad(lon2d, lat2d, LON0_DEG, LAT0_DEG)

    eta = AMP_M * np.exp(-0.5 * (r / SIGMA_RAD) ** 2)
    eta = np.where(ocean_mask, eta, 0.0)

    f0 = 2.0 * OMEGA * np.sin(np.deg2rad(LAT0_DEG))
    if abs(f0) < 1.0e-10:
        raise ValueError("LAT0_DEG too close to equator for geostrophic balance.")

    psi = (G / f0) * eta

    dpsi_dlat = np.gradient(psi, np.deg2rad(DLAT_DEG), axis=0)
    dpsi_dlon = np.gradient(psi, np.deg2rad(DLON_DEG), axis=1)

    coslat = np.cos(np.deg2rad(lat2d))
    coslat = np.where(np.abs(coslat) < 1.0e-8, np.nan, coslat)

    u_c = -(1.0 / R_EARTH) * dpsi_dlat
    v_c = +(1.0 / (R_EARTH * coslat)) * dpsi_dlon

    u_c = np.nan_to_num(u_c)
    v_c = np.nan_to_num(v_c)

    u = 0.5 * (u_c + np.roll(u_c, 1, axis=1))
    v = 0.5 * (v_c + np.roll(v_c, 1, axis=0))

    u = np.where(u_mask, u, 0.0)
    v = np.where(v_mask, v, 0.0)

    return eta, u, v

# Main: choose experiment
def main():
    base_dir = Path(__file__).resolve().parent

    write_flat_bathymetry(base_dir)
    ocean_mask = load_ocean_mask(base_dir)

    #! Experiment 1: cosine bump, unbalanced, gravity waves
    # eta, u, v = make_cosine_bump(ocean_mask)

    #! Experiment 2: geostrophic Gaussian, balanced eddy
    # eta, u, v = make_geostrophic_gaussian(ocean_mask)

    #! Experiment 3: Gaussian Patch, amplitude=0.1 m, sigma=10 deg, distance in meters
    eta, u, v = exp3_Gaussian_Patch(ocean_mask)

    write_field(base_dir, "eta_init.bin", eta)
    write_field(base_dir, "u_init.bin", u)
    write_field(base_dir, "v_init.bin", v)

    print("Done.")

if __name__ == "__main__":
    main()

# %%
