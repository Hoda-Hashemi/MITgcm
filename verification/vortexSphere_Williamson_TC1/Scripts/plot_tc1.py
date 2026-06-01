#!/usr/bin/env python3
"""Plot the TC1 numerical and analytic fields."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

EARTH_RADIUS = 6.37122e6
DAY = 86400.0
OMEGA_12DAY = 2.0 * np.pi * EARTH_RADIUS / (12.0 * DAY)
BELT_CENTER_LON_DEG = 270.0
BELT_CENTER_LAT_DEG = 0.0
BELL_HEIGHT_M = 1000.0
BELL_RADIUS_M = EARTH_RADIUS / 3.0


def read_scalar_from_file(path: Path, name: str, default: float) -> float:
    if not path.exists():
        return default
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"{re.escape(name)}\s*=\s*([+\-0-9.eEdD]+)", text)
    if match is None:
        return default
    return float(match.group(1).replace("D", "e").replace("d", "e"))


def read_alpha(input_dir: Path, override: float | None) -> float:
    if override is not None:
        return override
    return read_scalar_from_file(input_dir / "data.mypackage", "myPa_param1", 0.5 * np.pi)


def read_delta_t(input_dir: Path) -> float:
    return read_scalar_from_file(input_dir / "data", "deltaT", 60.0)


def parse_meta(meta_path: Path) -> dict[str, object]:
    text = meta_path.read_text(encoding="utf-8")
    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_vals = [int(v) for v in re.findall(r"[-+]?\d+", dim_block)]
    dims = []
    for idx in range(n_dims):
        global_size, start, end = dim_vals[3 * idx : 3 * idx + 3]
        dims.append({"global": global_size, "start": start, "end": end, "n": end - start + 1})
    prec_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", text)
    nrecords_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", text)
    fld_block = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    field_names = [] if fld_block is None else [name.strip() for name in re.findall(r"'([^']+)'", fld_block.group(1))]
    return {
        "dims": dims,
        "prec": prec_match.group(1).lower() if prec_match else "float64",
        "nrecords": int(nrecords_match.group(1)) if nrecords_match else 1,
        "fields": field_names,
    }


def discover_latest_meta(run_dir: Path, var_name: str) -> Path:
    candidates = sorted(run_dir.glob(f"{var_name}.*.meta"))
    if candidates:
        def iteration(path: Path) -> int:
            match = re.search(r"\.(\d{10})\.meta$", path.name)
            if match is None:
                return -1
            return int(match.group(1))

        return max(candidates, key=iteration)

    static_meta = run_dir / f"{var_name}.meta"
    if static_meta.exists():
        return static_meta

    if not candidates:
        raise FileNotFoundError(f"No diagnostic files found for {var_name!r} in {run_dir}")


def read_mds_record(run_dir: Path, var_name: str, record: int = 0, iteration: int | None = None) -> tuple[np.ndarray, dict[str, object], Path, int]:
    if iteration is None:
        meta_path = discover_latest_meta(run_dir, var_name)
    else:
        meta_path = run_dir / f"{var_name}.{iteration:010d}.meta"
        if not meta_path.exists():
            raise FileNotFoundError(meta_path)

    meta = parse_meta(meta_path)
    data_path = meta_path.with_suffix(".data")
    prec = str(meta["prec"])
    dtype = ">f8" if "64" in prec or "real*8" in prec else ">f4"
    raw = np.fromfile(data_path, dtype=dtype)

    dims = meta["dims"]
    nx = dims[0]["global"]
    ny = dims[1]["global"]
    nrecords = int(meta["nrecords"])

    if nrecords == 1:
        field = raw.reshape((ny, nx))
    else:
        field = raw.reshape((nrecords, ny, nx))[record]

    match = re.search(r"\.(\d{10})\.meta$", meta_path.name)
    iteration_found = int(match.group(1)) if match else 0
    return field, meta, meta_path, iteration_found


def load_regular_geometry(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lon_c = (np.arange(nx, dtype=np.float64) + 0.5) * (360.0 / nx)
    lat_c = -90.0 + (np.arange(ny, dtype=np.float64) + 0.5) * (180.0 / ny)
    lon_u = np.arange(nx, dtype=np.float64) * (360.0 / nx)
    lat_u = lat_c.copy()
    lon_v = lon_c.copy()
    lat_v = -90.0 + np.arange(ny, dtype=np.float64) * (180.0 / ny)
    xc, yc = np.meshgrid(lon_c, lat_c, indexing="xy")
    xg, _ = np.meshgrid(lon_u, lat_u, indexing="xy")
    _, yg = np.meshgrid(lon_v, lat_v, indexing="xy")
    return xc, yc, xg, yg


def load_geometry(input_dir: Path, run_dir: Path, nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    for directory in (input_dir, run_dir, input_dir.parent / "build"):
        xc_meta = directory / "XC.meta"
        yc_meta = directory / "YC.meta"
        if xc_meta.exists() and yc_meta.exists():
            xc, _, _, _ = read_mds_record(directory, "XC")[:4]
            yc, _, _, _ = read_mds_record(directory, "YC")[:4]
            if (directory / "XG.meta").exists() and (directory / "YG.meta").exists():
                xg, _, _, _ = read_mds_record(directory, "XG")[:4]
                yg, _, _, _ = read_mds_record(directory, "YG")[:4]
                return xc, yc, xg, yg
            return xc, yc, *load_regular_geometry(nx, ny)[2:]
    return load_regular_geometry(nx, ny)


def lonlat_to_unitvec(lon_deg: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    return np.stack(
        (
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ),
        axis=-1,
    )


def rotate_vector(vec: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    return vec * c + np.cross(axis, vec) * s + axis * np.dot(axis, vec) * (1.0 - c)


def exact_eta(xc: np.ndarray, yc: np.ndarray, alpha: float, time_s: float) -> np.ndarray:
    grid_vec = lonlat_to_unitvec(xc, yc)
    bell_center = lonlat_to_unitvec(np.array(BELT_CENTER_LON_DEG), np.array(BELT_CENTER_LAT_DEG))
    axis = np.array([-math.sin(alpha), 0.0, math.cos(alpha)], dtype=np.float64)
    angle = OMEGA_12DAY * time_s / EARTH_RADIUS
    center_t = rotate_vector(bell_center, axis, angle)
    dot = np.clip(np.sum(grid_vec * center_t, axis=-1), -1.0, 1.0)
    radius = EARTH_RADIUS * np.arccos(dot)
    eta = np.zeros_like(radius)
    wet = radius < BELL_RADIUS_M
    eta[wet] = 0.5 * BELL_HEIGHT_M * (1.0 + np.cos(np.pi * radius[wet] / BELL_RADIUS_M))
    return eta


def exact_velocity(xg: np.ndarray, yc: np.ndarray, xc: np.ndarray, yg: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    u0 = OMEGA_12DAY
    lon_u = np.deg2rad(xg)
    lat_u = np.deg2rad(yc)
    lon_v = np.deg2rad(xc)
    u_exact = u0 * (np.cos(lat_u) * np.cos(alpha) + np.sin(lat_u) * np.cos(lon_u) * np.sin(alpha))
    v_exact = -u0 * np.sin(lon_v) * np.sin(alpha)
    return u_exact, v_exact


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iteration", type=int, default=None, help="Specific diagnostics iteration to plot.")
    parser.add_argument("--alpha", type=float, default=None, help="Override alpha in radians.")
    args = parser.parse_args()

    case_dir = Path(__file__).resolve().parent.parent
    input_dir = case_dir / "input"
    run_dir = case_dir / "run"
    out_dir = Path(__file__).resolve().parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    eta_num, eta_meta, _, iteration = read_mds_record(run_dir, "tc1Diag", record=0, iteration=args.iteration)
    if int(eta_meta["nrecords"]) < 3:
        raise ValueError("tc1Diag does not contain ETAN/UVEL/VVEL records")
    u_num, _, _, _ = read_mds_record(run_dir, "tc1Diag", record=1, iteration=iteration)
    v_num, _, _, _ = read_mds_record(run_dir, "tc1Diag", record=2, iteration=iteration)

    delta_t = read_delta_t(input_dir)
    alpha = read_alpha(input_dir, args.alpha)
    time_s = iteration * delta_t

    dims = eta_meta["dims"]
    nx = dims[0]["global"]
    ny = dims[1]["global"]
    xc, yc, xg, yg = load_geometry(input_dir, run_dir, nx, ny)

    eta_exact = exact_eta(xc, yc, alpha, time_s)
    u_exact, v_exact = exact_velocity(xg, yc, xc, yg, alpha)

    eta_err = eta_num - eta_exact
    vel_err = np.sqrt((u_num - u_exact) ** 2 + (v_num - v_exact) ** 2)

    eta_lim = max(float(np.max(np.abs(eta_exact))), float(np.max(np.abs(eta_num))))
    err_lim = max(float(np.max(np.abs(eta_err))), 1.0e-12)
    vel_err_lim = max(float(np.max(vel_err)), 1.0e-12)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    panels = [
        (axes[0, 0], eta_num, "Numerical ETAN", "m", 0.0, eta_lim, "viridis"),
        (axes[0, 1], eta_exact, "Exact ETAN", "m", 0.0, eta_lim, "viridis"),
        (axes[1, 0], eta_err, "ETAN error", "m", -err_lim, err_lim, "RdBu_r"),
        (axes[1, 1], vel_err, "Velocity error magnitude", "m/s", 0.0, vel_err_lim, "magma"),
    ]

    for ax, data, title, cbar_label, vmin, vmax, cmap in panels:
        im = ax.imshow(
            data,
            origin="lower",
            extent=(0.0, 360.0, -90.0, 90.0),
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(title)
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        fig.colorbar(im, ax=ax, shrink=0.9, label=cbar_label)

    fig.suptitle(
        f"Williamson TC1 | iteration {iteration} | t = {time_s / DAY:.6f} days | alpha = {alpha:.6f} rad"
    )

    png_path = out_dir / f"tc1_overview_it{iteration:010d}.png"
    pdf_path = out_dir / f"tc1_overview_it{iteration:010d}.pdf"
    fig.savefig(png_path, dpi=250, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=250, bbox_inches="tight")
    print(f"Saved {png_path}")
    print(f"Saved {pdf_path}")


if __name__ == "__main__":
    main()
