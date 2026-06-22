#%%
from pathlib import Path
import numpy as np

# Grid / constants
NX = 1440
NY = 720
DLON_DEG = 360.0 / NX
DLAT_DEG = 180.0 / NY

R_EARTH = 6.371e6
DAY = 86400.0

# MITgcm ocean bathymetry is negative below sea level.
BATHYMETRY = -4000.0
BATHYMETRY_FILENAME = "bathymetry_flat4000.bin"

# TC1 parameters
U0 = 2.0 * np.pi * R_EARTH / (12.0 * DAY)
H0 = 1000.0               # passive-tracer amplitude
RADIUS_RAD = 1.0 / 3.0    # angular radius
LAT0_DEG = 0.0
LON0_DEG = 270.0          # 3*pi/2
ALPHA_RAD = 0.0
#4 cases: [0.0 , 0.05 ,1.5207963267948966 , 1.5707963267948966]

def write_flat_bathymetry(base_dir: Path):
    bathymetry = np.full((NY, NX), BATHYMETRY, dtype=">f4")
    bathymetry.tofile(base_dir / BATHYMETRY_FILENAME)
    print(
        f"Wrote {BATHYMETRY_FILENAME}: bathymetry={BATHYMETRY} m"
    )

def grid_centers():
    lon = (np.arange(NX, dtype=np.float64) + 0.5) * DLON_DEG
    lat = -90.0 + (np.arange(NY, dtype=np.float64) + 0.5) * DLAT_DEG
    return np.meshgrid(lon, lat, indexing="xy")

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

def write_field(base_dir: Path, name: str, field: np.ndarray):
    field.astype(">f4").tofile(base_dir / name)
    print(f"Wrote {name}: min={field.min():.6e}, max={field.max():.6e}")

def streamfunction(lon_rad, lat_rad, alpha_rad):
    return U0 * R_EARTH * (
        np.cos(alpha_rad) * np.sin(lat_rad)
        - np.sin(alpha_rad) * np.cos(lat_rad) * np.cos(lon_rad)
    )

def make_tc1_fields(alpha_rad: float = 0.0):
    lon2d, lat2d = grid_centers()

    # Cosine bell at tracer-cell centers.
    r = spherical_distance_rad(lon2d, lat2d, LON0_DEG, LAT0_DEG)

    tracer = np.zeros((NY, NX), dtype=np.float64)
    inside = r < RADIUS_RAD
    tracer[inside] = 0.5 * H0 * (
        1.0 + np.cos(np.pi * r[inside] / RADIUS_RAD)
    )

    # Derive C-grid face velocities from a streamfunction so their
    # finite-volume divergence vanishes to roundoff.
    lon_u = np.deg2rad(np.arange(NX, dtype=np.float64) * DLON_DEG)
    lat_s = np.deg2rad(-90.0 + np.arange(NY, dtype=np.float64) * DLAT_DEG)
    lat_n = lat_s + np.deg2rad(DLAT_DEG)
    dy = R_EARTH * np.deg2rad(DLAT_DEG)

    psi_n = streamfunction(lon_u[None, :], lat_n[:, None], alpha_rad)
    psi_s = streamfunction(lon_u[None, :], lat_s[:, None], alpha_rad)
    u = (psi_n - psi_s) / dy

    lon_w = np.deg2rad(np.arange(NX, dtype=np.float64) * DLON_DEG)
    lon_e = np.roll(lon_w, -1)
    lat_v = np.deg2rad(-90.0 + np.arange(NY, dtype=np.float64) * DLAT_DEG)
    dx = R_EARTH * np.cos(lat_v) * np.deg2rad(DLON_DEG)

    psi_w = streamfunction(lon_w[None, :], lat_v[:, None], alpha_rad)
    psi_e = streamfunction(lon_e[None, :], lat_v[:, None], alpha_rad)
    v = np.zeros((NY, NX), dtype=np.float64)
    non_polar = np.abs(dx) > 1.0e-14
    v[non_polar, :] = (
        psi_w[non_polar, :] - psi_e[non_polar, :]
    ) / dx[non_polar, None]

    return tracer, u, v

def main():
    base_dir = Path(__file__).resolve().parent

    write_flat_bathymetry(base_dir)

    tracer, u, v = make_tc1_fields(alpha_rad=ALPHA_RAD)

    write_field(base_dir, "S_init.bin", tracer)
    write_field(base_dir, "u_init.bin", u)
    write_field(base_dir, "v_init.bin", v)

    print("Done.")

if __name__ == "__main__":
    main()

# %%
