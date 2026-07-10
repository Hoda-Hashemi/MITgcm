#%%
from pathlib import Path
import sys
import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent.parent.parent / "Scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import williamson_cube as wc

# ============================================================
# Williamson TC2: Global steady-state nonlinear zonal geostrophic flow
# Clean MITgcm initializer for one-layer spherical shallow water
# ============================================================

# Grid
NX = 1440
NY = 720
DLON_RAD = 2.0 * np.pi / NX
DLAT_RAD = np.pi / NY

# Physical constants
R_EARTH = 6.371e6
OMEGA = 7.292e-5
G = 9.81

# Standard Williamson TC2 constants
DAY = 86400.0
U0 = 2.0 * np.pi * R_EARTH / (12.0 * DAY)   # one revolution in 12 days
GH0 = 2.94e4                                # m^2/s^2
H0 = GH0 / G                                # mean layer depth in meters

# Rotation angle. The job scripts patch ALPHA_RAD in the run-local copy
# before this generator writes the initial fields.
ALPHA_RAD = 0.05

# Output files
BATHY_FILE = "bathymetry_flat4000.bin"
ETA_FILE = "eta_init.bin"
U_FILE = "u_init.bin"
V_FILE = "v_init.bin"

def write_field(base_dir: Path, name: str, field: np.ndarray) -> None:
    field.astype(">f4").tofile(base_dir / name)
    print(f"Wrote {name}: shape={field.shape}, min={field.min():.6e}, max={field.max():.6e}")

def grid_cell_centers():
    """
    C-cell centers:
      lon_c = (i+1/2) dlon
      lat_c = -pi/2 + (j+1/2) dlat
    """
    lon_c = (np.arange(NX, dtype=np.float64) + 0.5) * DLON_RAD
    lat_c = -0.5 * np.pi + (np.arange(NY, dtype=np.float64) + 0.5) * DLAT_RAD
    return np.meshgrid(lon_c, lat_c, indexing="xy")

def grid_u_faces():
    """
    U faces:
      lon_u = i dlon
      lat_u = cell-center latitude
    """
    lon_u = np.arange(NX, dtype=np.float64) * DLON_RAD
    lat_u = -0.5 * np.pi + (np.arange(NY, dtype=np.float64) + 0.5) * DLAT_RAD
    return np.meshgrid(lon_u, lat_u, indexing="xy")

def grid_v_faces():
    """
    V faces:
      lon_v = cell-center longitude
      lat_v = -pi/2 + j dlat
    """
    lon_v = (np.arange(NX, dtype=np.float64) + 0.5) * DLON_RAD
    lat_v = -0.5 * np.pi + np.arange(NY, dtype=np.float64) * DLAT_RAD
    return np.meshgrid(lon_v, lat_v, indexing="xy")

def write_flat_bathymetry(base_dir: Path) -> None:
    bathy = np.full((NY, NX), -H0, dtype=">f4")
    bathy.tofile(base_dir / BATHY_FILE)
    print(f"Wrote {BATHY_FILE} with constant depth {-H0:.6f} m")

def make_tc2_eta(alpha_rad: float = ALPHA_RAD):
    """
    Cell-center free-surface anomaly eta = h - H0,
    with h from the Williamson TC2 analytic height field.
    """
    lon_c, lat_c = grid_cell_centers()

    mu = (
        -np.cos(lon_c) * np.cos(lat_c) * np.sin(alpha_rad)
        + np.sin(lat_c) * np.cos(alpha_rad)
    )

    gh = GH0 - (R_EARTH * OMEGA * U0 + 0.5 * U0**2) * mu**2
    h = gh / G
    eta = h - H0

    return eta

def make_tc2_u(alpha_rad: float = ALPHA_RAD):
    """
    Zonal-face velocity from the paper formula.
    """
    lon_u, lat_u = grid_u_faces()

    u = U0 * (
        np.cos(lat_u) * np.cos(alpha_rad)
        + np.cos(lon_u) * np.sin(lat_u) * np.sin(alpha_rad)
    )

    return u

def make_tc2_v(alpha_rad: float = ALPHA_RAD):
    """
    Meridional-face velocity from the paper formula.
    """
    lon_v, lat_v = grid_v_faces()

    v = -U0 * np.sin(lon_v) * np.sin(alpha_rad)

    # Pole rows are coordinate-singular V-faces on the lat-lon grid.
    # Set them to zero for a clean MITgcm initialization.
    v[0, :] = 0.0
    v[-1, :] = 0.0

    return v

def main():
    base_dir = Path(__file__).resolve().parent
    if not (base_dir / ".use_latlon_grid").exists():
        wc.generate_tc2(base_dir, alpha_rad=ALPHA_RAD)
        print("")
        print("Done.")
        print(f"ALPHA_RAD = {ALPHA_RAD:.16g}")
        print(f"ALPHA_DEG = {np.rad2deg(ALPHA_RAD):.16g}")
        print(f"U0        = {U0:.12f} m/s")
        print(f"H0        = {H0:.12f} m")
        print(f"GH0       = {GH0:.12f} m^2/s^2")
        return

    write_flat_bathymetry(base_dir)
    eta = make_tc2_eta(alpha_rad=ALPHA_RAD)
    u = make_tc2_u(alpha_rad=ALPHA_RAD)
    v = make_tc2_v(alpha_rad=ALPHA_RAD)

    write_field(base_dir, ETA_FILE, eta)
    write_field(base_dir, U_FILE, u)
    write_field(base_dir, V_FILE, v)

    print("")
    print("Done.")
    print(f"ALPHA_RAD = {ALPHA_RAD:.16g}")
    print(f"ALPHA_DEG = {np.rad2deg(ALPHA_RAD):.16g}")
    print(f"U0        = {U0:.12f} m/s")
    print(f"H0        = {H0:.12f} m")
    print(f"GH0       = {GH0:.12f} m^2/s^2")

if __name__ == "__main__":
    main()

# %%
