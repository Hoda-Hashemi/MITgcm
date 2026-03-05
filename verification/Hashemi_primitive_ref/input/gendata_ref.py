# #!/usr/bin/env python3
# """Generate eta_init.bin for Hashemi_primitive_ref.

# Creates a smooth, zero-mean Gaussian free-surface anomaly (meters)
# on the 1440x720 global 0.25-degree lat-lon grid.
# """

# from pathlib import Path
# import numpy as np

# NX = 1440
# NY = 720
# AMP_M = 0.20
# SIG_LAT_DEG = 10.0
# SIG_LON_DEG = 20.0
# LAT0_DEG = 0.0
# LON0_DEG = 180.0

# def periodic_lon_distance_deg(lon_deg: np.ndarray, lon0_deg: float) -> np.ndarray:
#     d = np.abs(lon_deg - lon0_deg)
#     return np.minimum(d, 360.0 - d)

# def main() -> None:
#     lon = (np.arange(NX, dtype=np.float64) + 0.5) * (360.0 / NX)
#     lat = -90.0 + (np.arange(NY, dtype=np.float64) + 0.5) * (180.0 / NY)

#     lon2d, lat2d = np.meshgrid(lon, lat, indexing="xy")
#     dlon = periodic_lon_distance_deg(lon2d, LON0_DEG)

#     eta = AMP_M * np.exp(
#         -(
#             ((lat2d - LAT0_DEG) ** 2) / (2.0 * SIG_LAT_DEG**2)
#             + (dlon**2) / (2.0 * SIG_LON_DEG**2)
#         )
#     )

#     # Remove mean to avoid artificial global mass bias.
#     eta -= np.mean(eta)

#     out = Path(__file__).resolve().parent / "eta_init.bin"
#     eta.astype(">f4").tofile(out)
#     print(f"Wrote {out} with shape={eta.shape}, dtype='>f4'")
#     print(f"eta min={eta.min():.6e} m, max={eta.max():.6e} m")

# if __name__ == "__main__":
#     main()

#!/usr/bin/env python3
"""Generate eta_init.bin for Hashemi_primitive_ref.

Creates a smooth, zero-mean Gaussian free-surface anomaly (meters)
on the 1440x720 global 0.25-degree lat-lon grid.
"""

from pathlib import Path
import numpy as np

NX = 1440
NY = 720
AMP_M = 0.20
SIG_LAT_DEG = 10.0
SIG_LON_DEG = 20.0
LAT0_DEG = 0.0
LON0_DEG = 180.0
COAST_BUFFER_CELLS = 1
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

def grow_land_mask(land: np.ndarray, ocean: np.ndarray, n_cells: int) -> np.ndarray:
    """
    Expand land mask into adjacent ocean cells by n_cells (4-neighbor stencil).
    Longitude is periodic; latitude is not.
    """
    if n_cells <= 0:
        return land.copy()

    blocked = land.copy()
    for _ in range(n_cells):
        north = np.zeros_like(blocked)
        north[1:, :] = blocked[:-1, :]
        south = np.zeros_like(blocked)
        south[:-1, :] = blocked[1:, :]
        west = np.roll(blocked, 1, axis=1)
        east = np.roll(blocked, -1, axis=1)
        adjacent = north | south | west | east
        blocked |= (ocean & adjacent)
    return blocked

def main() -> None:
    base_dir = Path(__file__).resolve().parent
    ocean_mask = load_ocean_mask(base_dir)
    land_mask = ~ocean_mask
    blocked_mask = grow_land_mask(land_mask, ocean_mask, COAST_BUFFER_CELLS)
    active_mask = ocean_mask & ~blocked_mask

    n_active = int(active_mask.sum())
    if n_active == 0:
        raise RuntimeError("No active ocean cells left after coastline buffering.")

    lon = (np.arange(NX, dtype=np.float64) + 0.5) * (360.0 / NX)
    lat = -90.0 + (np.arange(NY, dtype=np.float64) + 0.5) * (180.0 / NY)

    lon2d, lat2d = np.meshgrid(lon, lat, indexing="xy")
    dlon = periodic_lon_distance_deg(lon2d, LON0_DEG)

    eta = AMP_M * np.exp(
        -(
            ((lat2d - LAT0_DEG) ** 2) / (2.0 * SIG_LAT_DEG**2)
            + (dlon**2) / (2.0 * SIG_LON_DEG**2)
        )
    )

    # Keep anomaly only in active ocean cells (ocean minus coast buffer).
    eta = np.where(active_mask, eta, 0.0)

    # Remove mean over active ocean points only to avoid global mass bias.
    eta_mean_active = np.mean(eta[active_mask])
    eta[active_mask] -= eta_mean_active
    eta[~active_mask] = 0.0

    out = base_dir / "eta_init.bin"
    eta.astype(">f4").tofile(out)
    print(f"Wrote {out} with shape={eta.shape}, dtype='>f4'")
    print(f"eta min={eta.min():.6e} m, max={eta.max():.6e} m")
    print(
        f"active cells={n_active}, ocean cells={int(ocean_mask.sum())}, "
        f"coast buffer={COAST_BUFFER_CELLS} cell(s)"
    )

if __name__ == "__main__":
    main()
