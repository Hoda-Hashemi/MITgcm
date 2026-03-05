#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert ETOPO5 flattened CSV (1 value/line) into MITgcm bathymetry.bin "
            "(big-endian float32, R_low convention)."
        )
    )
    parser.add_argument(
        "--input-csv",
        default="verification/Hashemi/resources/topography/etopo5.csv",
        help="Path to flattened etopo5.csv (single column, header 'topo').",
    )
    parser.add_argument(
        "--output-bin",
        default="verification/Hashemi/input/bathymetry.bin",
        help="Output MITgcm bathymetry binary path.",
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
            "Keep positive topography over land. By default, land is clipped to 0 "
            "(ocean stays negative), which is usually safer for ocean-only runs."
        ),
    )
    return parser.parse_args()


def load_flat_csv(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")
    values = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.float32)
    if values.ndim != 1:
        values = values.reshape(-1)
    return values


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

    src_lat = -90.0 + np.arange(src_ny, dtype=np.float64) / 12.0
    src_lon = np.arange(src_nx, dtype=np.float64) / 12.0

    dlon = 360.0 / float(args.nx)
    dlat = 180.0 / float(args.ny)
    tgt_lon = (np.arange(args.nx, dtype=np.float64) + 0.5) * dlon
    tgt_lat = -90.0 + (np.arange(args.ny, dtype=np.float64) + 0.5) * dlat

    lon_interp = np.empty((src_ny, args.nx), dtype=np.float32)
    for j in range(src_ny):
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

