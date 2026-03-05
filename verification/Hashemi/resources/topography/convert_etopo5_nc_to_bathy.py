#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert ETOPO5 NetCDF into MITgcm bathymetry.bin "
            "(big-endian float32, R_low convention)."
        )
    )
    parser.add_argument(
        "--input-nc",
        default="verification/Hashemi/resources/topography/etopo5.nc",
        help="Path to ETOPO5 NetCDF file.",
    )
    parser.add_argument(
        "--output-bin",
        default="verification/Hashemi/input/bathymetry.bin",
        help="Output MITgcm bathymetry binary path.",
    )
    parser.add_argument("--nx", type=int, default=1440, help="Target MITgcm Nx.")
    parser.add_argument("--ny", type=int, default=720, help="Target MITgcm Ny.")
    parser.add_argument(
        "--lon-var",
        default="topo_lon",
        help="Longitude variable name in NetCDF.",
    )
    parser.add_argument(
        "--lat-var",
        default="topo_lat",
        help="Latitude variable name in NetCDF.",
    )
    parser.add_argument(
        "--topo-var",
        default="topo",
        help="Topography variable name in NetCDF.",
    )
    parser.add_argument(
        "--keep-land-positive",
        action="store_true",
        help=(
            "Keep positive topography over land. By default, land is clipped to 0 "
            "(ocean stays negative), which is usually safer for ocean-only runs."
        ),
    )
    return parser.parse_args()


def load_nc_fields(
    nc_path: Path, lon_var: str, lat_var: str, topo_var: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not nc_path.exists():
        raise FileNotFoundError(f"Input NetCDF not found: {nc_path}")

    with netcdf_file(str(nc_path), mode="r", mmap=False) as nc:
        if lon_var not in nc.variables:
            raise KeyError(f"Missing lon var '{lon_var}'")
        if lat_var not in nc.variables:
            raise KeyError(f"Missing lat var '{lat_var}'")
        if topo_var not in nc.variables:
            raise KeyError(f"Missing topo var '{topo_var}'")

        lon = np.array(nc.variables[lon_var][:], dtype=np.float64)
        lat = np.array(nc.variables[lat_var][:], dtype=np.float64)
        topo = np.array(nc.variables[topo_var][:], dtype=np.float32)

    if topo.ndim != 2:
        raise ValueError(f"Topography must be 2D, got shape {topo.shape}")
    if topo.shape != (lat.size, lon.size):
        raise ValueError(
            f"Topography shape {topo.shape} incompatible with lat/lon sizes "
            f"({lat.size}, {lon.size})"
        )

    return lon, lat, topo


def main() -> int:
    args = parse_args()
    input_nc = Path(args.input_nc)
    output_bin = Path(args.output_bin)
    output_bin.parent.mkdir(parents=True, exist_ok=True)

    src_lon, src_lat, topo = load_nc_fields(
        input_nc, args.lon_var, args.lat_var, args.topo_var
    )

    if not args.keep_land_positive:
        topo = np.minimum(topo, 0.0)

    dlon = 360.0 / float(args.nx)
    dlat = 180.0 / float(args.ny)
    tgt_lon = (np.arange(args.nx, dtype=np.float64) + 0.5) * dlon
    tgt_lat = -90.0 + (np.arange(args.ny, dtype=np.float64) + 0.5) * dlat

    lon_interp = np.empty((src_lat.size, args.nx), dtype=np.float32)
    for j in range(src_lat.size):
        lon_interp[j, :] = np.interp(tgt_lon, src_lon, topo[j, :], period=360.0)

    bathy = np.empty((args.ny, args.nx), dtype=np.float32)
    for i in range(args.nx):
        bathy[:, i] = np.interp(tgt_lat, src_lat, lon_interp[:, i])

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
