#!/usr/bin/env python3
#%%
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert ETOPO5 CSV into MITgcm bathymetry.bin (big-endian float32)."
        )
    )
    parser.add_argument(
        "--input-csv",
        default="verification/vortexSphere_mitgcm_Poisson/resources/topography/etopo5.csv",
        help="Path to etopo5.csv (single column, header 'topo').",
    )
    parser.add_argument(
        "--output-bin",
        default="verification/vortexSphere_mitgcm_Poisson/input/bathymetry.bin",
        help="Output bathymetry.bin path.",
    )
    parser.add_argument(
        "--nx",
        type=int,
        default=1440,
        help="Target MITgcm Nx (default: 1440).",
    )
    parser.add_argument(
        "--ny",
        type=int,
        default=720,
        help="Target MITgcm Ny (default: 720).",
    )
    parser.add_argument(
        "--keep-land-positive",
        action="store_true",
        help=(
            "Keep land positive. Default clips land to 0 (ocean stays negative)."
        ),
    )
    # Ignore extra args (e.g., Jupyter's --f=...).
    args, _ = parser.parse_known_args(argv)
    return args

def load_flat_csv(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")
    values = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.float32)
    if values.ndim != 1:
        values = values.reshape(-1)
    return values

def _bilinear_interpolate_global(topo_lat_lon: np.ndarray, nx: int, ny: int) -> np.ndarray:
    """
    Bilinear interpolation from ETOPO5 [lat, lon] to MITgcm [ny, nx].
    Ocean should be negative (land already clipped to 0).
    """
    orig_nlat = topo_lat_lon.shape[0]
    topo_padded = np.vstack([topo_lat_lon, topo_lat_lon[-1:, :]])
    topo = topo_padded.T  # (lon, lat_padded)

    nlon_b, nlat_b = topo.shape
    dlat_b = np.pi / (orig_nlat - 1)
    dphi_b = 2.0 * np.pi / nlon_b

    dlat = np.pi / ny
    dphi = 2.0 * np.pi / nx

    depth = np.zeros((ny, nx), dtype=np.float32)

    for ilat in range(ny):
        lat0 = -0.5 * np.pi + 0.5 * dlat + dlat * ilat

        ilat_b = 1 + int((lat0 + 0.5 * np.pi - 0.5 * dlat_b) / dlat_b)
        ilat_t = ilat_b + 1
        if ilat_t > nlat_b:
            ilat_t = nlat_b

        lat_b = -0.5 * np.pi + 0.5 * dlat_b + dlat_b * (ilat_b - 1)
        lat_t = -0.5 * np.pi + 0.5 * dlat_b + dlat_b * (ilat_t - 1)
        alpha_lat = 0.0 if lat_t == lat_b else (lat0 - lat_b) / (lat_t - lat_b)

        for iphi in range(nx):
            phi0 = 0.5 * dphi + dphi * iphi

            iphi_l = 1 + int((phi0 - 0.5 * dphi_b) / dphi_b)
            iphi_r = iphi_l + 1

            if iphi_l > nlon_b:
                iphi_l -= nlon_b
            if iphi_r > nlon_b:
                iphi_r -= nlon_b
            if iphi_l < 1:
                iphi_l += nlon_b
            if iphi_r < 1:
                iphi_r += nlon_b

            ilat_b0, ilat_t0 = ilat_b - 1, ilat_t - 1
            iphi_l0, iphi_r0 = iphi_l - 1, iphi_r - 1

            phi_l = 0.5 * dphi_b + dphi_b * (iphi_l - 1)
            phi_r = 0.5 * dphi_b + dphi_b * (iphi_r - 1)
            alpha_lon = 0.0 if phi_r == phi_l else (phi0 - phi_l) / (phi_r - phi_l)

            z_b = topo[iphi_l0, ilat_b0]
            z_t = topo[iphi_l0, ilat_t0]
            z_phi_l = z_b + alpha_lat * (z_t - z_b)

            z_b = topo[iphi_r0, ilat_b0]
            z_t = topo[iphi_r0, ilat_t0]
            z_phi_r = z_b + alpha_lat * (z_t - z_b)

            depth[ilat, iphi] = z_phi_l + alpha_lon * (z_phi_r - z_phi_l)

    return depth

def main() -> int:
    args = parse_args()
    input_csv = Path(args.input_csv)
    output_bin = Path(args.output_bin)
    output_bin.parent.mkdir(parents=True, exist_ok=True)

    flat = load_flat_csv(input_csv)

    src_ny = 2161
    src_nx = 4320
    expected = src_ny * src_nx
    if flat.size != expected:
        raise ValueError(
            f"Unexpected ETOPO size: got {flat.size}, expected {expected} ({src_ny}x{src_nx})"
        )

    topo = flat.reshape(src_ny, src_nx)
    if not args.keep_land_positive:
        topo = np.minimum(topo, 0.0)

    bathy = _bilinear_interpolate_global(topo, args.nx, args.ny)

    if not args.keep_land_positive:
        bathy = np.minimum(bathy, 0.0)

    bathy.astype(">f4").tofile(output_bin)

    print(f"Wrote {output_bin}")
    print(f"Shape: {args.ny} x {args.nx}")
    print(f"Bytes: {output_bin.stat().st_size}")
    print(
        f"Range (m): min={float(np.nanmin(bathy)):.3f}, max={float(np.nanmax(bathy)):.3f}"
    )
    if not args.keep_land_positive:
        print("Land points clipped to 0.0 m (ocean remains negative).")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

# %%
