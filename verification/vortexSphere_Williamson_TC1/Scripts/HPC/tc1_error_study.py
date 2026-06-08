#!/usr/bin/env python3
"""TC1 error study for the Williamson cosine-bell advection test.

This script compares the MITgcm passive-salinity height field against the
analytic rigid-rotation
solution for Shallow Water Test Case 1 and writes:

- a final-time comparison figure with exact / model / signed error maps
- a time-series figure for normalized L1/L2/Linf error norms
- a CSV file with the full error history

The normalization follows Ullrich et al. (2010), JCP 229, Eqs. 81-86.
"""

# %%
from __future__ import annotations

import csv
import os
import math
import re
import tempfile
from pathlib import Path

_CACHE_ROOT = Path(tempfile.gettempdir()) / "mitgcm_tc1_error_study_cache"
_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg_cache"))

import matplotlib

def _running_in_ipykernel() -> bool:
    try:
        from IPython import get_ipython
    except Exception:
        return False
    shell = get_ipython()
    return shell is not None and shell.__class__.__name__ == "ZMQInteractiveShell"

NOTEBOOK_MODE = _running_in_ipykernel()
if not NOTEBOOK_MODE:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_DIR = SCRIPT_DIR.parents[1]
RUN_DIR = CASE_DIR / "run"
OUTPUT_DIR = SCRIPT_DIR.parent / "output" / RUN_DIR.name / "TC1Error"

EARTH_RADIUS = 6_371_000.0
DAY = 86_400.0
OMEGA = 2.0 * math.pi / (12.0 * DAY)
BELL_LON_DEG = 270.0
BELL_LAT_DEG = 0.0
BELL_HEIGHT = 1000.0
BELL_RADIUS = EARTH_RADIUS / 3.0
DEFAULT_ALPHA = 0.5 * math.pi
HEIGHT_CONTOUR_LEVELS = np.arange(160.0, 801.0, 160.0)
ERROR_CONTOUR_STEP = 10.0
MAX_ERROR_CONTOURS_PER_SIGN = 10
ZERO_CONTOUR_TOLERANCE = 0.1
STATIONARY_FIELD_TOLERANCE = 1.0

def _save_figure(fig: matplotlib.figure.Figure, path: Path) -> None:
    fig.savefig(path)
    if not NOTEBOOK_MODE:
        plt.close(fig)

def read_text_value(path: Path, pattern: str) -> float | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(pattern, text, re.IGNORECASE)
    if match is None:
        return None
    return float(match.group(1).replace("D", "e").replace("d", "e"))

def read_delta_t(run_dir: Path) -> float:
    for candidate in (run_dir / "data", run_dir.parent / "input" / "data"):
        value = read_text_value(candidate, r"deltaT\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return value
    return 60.0

def read_alpha(run_dir: Path) -> float:
    for candidate in (
        run_dir / "data.mypackage",
        run_dir.parent / "input" / "data.mypackage",
    ):
        value = read_text_value(
            candidate,
            r"myPa_param1\s*=\s*([+\-0-9.eEdD]+)",
        )
        if value is not None:
            return value
    return DEFAULT_ALPHA

def parse_mds_meta(meta_path: Path) -> dict[str, object]:
    text = meta_path.read_text(encoding="utf-8", errors="ignore")
    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_vals = [int(value) for value in re.findall(r"[-+]?\d+", dim_block)]
    dims = []
    for index in range(n_dims):
        global_size, start, end = dim_vals[3 * index : 3 * index + 3]
        dims.append(
            {
                "global": global_size,
                "start": start,
                "end": end,
                "n": end - start + 1,
            }
        )
    prec_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", text)
    nrecords_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", text)
    return {
        "n_dims": n_dims,
        "dims": dims,
        "prec": prec_match.group(1).lower() if prec_match else "float32",
        "nrecords": int(nrecords_match.group(1)) if nrecords_match else 1,
    }

def dtype_from_meta(meta: dict[str, object]) -> np.dtype:
    prec = str(meta["prec"])
    return np.dtype(">f8" if prec in {"float64", "real*8"} else ">f4")

def read_mds_field(run_dir: Path, name: str, iteration: int | None = None, record: int = 0) -> np.ndarray:
    run_dir = Path(run_dir)
    if iteration is None:
        meta_path = run_dir / f"{name}.meta"
        data_path = run_dir / f"{name}.data"
    else:
        meta_path = run_dir / f"{name}.{int(iteration):010d}.meta"
        data_path = meta_path.with_suffix(".data")
    if not meta_path.exists() or not data_path.exists():
        raise FileNotFoundError(f"Missing MDS field {name} in {run_dir}")

    meta = parse_mds_meta(meta_path)
    raw = np.fromfile(data_path, dtype=dtype_from_meta(meta))
    dims = meta["dims"]  # type: ignore[assignment]
    n_dims = int(meta["n_dims"])
    nrecords = int(meta["nrecords"])

    if n_dims == 2:
        nx = dims[0]["global"]  # type: ignore[index]
        ny = dims[1]["global"]  # type: ignore[index]
        shape = (nrecords, ny, nx) if nrecords > 1 else (ny, nx)
    elif n_dims == 3:
        nx = dims[0]["global"]  # type: ignore[index]
        ny = dims[1]["global"]  # type: ignore[index]
        nz = dims[2]["global"]  # type: ignore[index]
        shape = (nrecords, nz, ny, nx) if nrecords > 1 else (nz, ny, nx)
    else:
        raise ValueError(f"Unsupported nDims={n_dims} in {meta_path}")

    data = raw.reshape(shape)
    return data[record] if nrecords > 1 else data

def discover_iterations(run_dir: Path, field: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(field)}\.(\d{{10}})\.meta$")
    iterations = []
    for path in Path(run_dir).glob(f"{field}.*.meta"):
        match = pattern.match(path.name)
        if match is not None:
            iterations.append(int(match.group(1)))
    return sorted(set(iterations))

def unit_vector(lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    cos_lat = np.cos(lat)
    return cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)

def rotate_about_axis(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    alpha: float,
    angle: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ax = -math.sin(alpha)
    ay = 0.0
    az = math.cos(alpha)
    c = math.cos(angle)
    s = math.sin(angle)
    dot = ax * x + ay * y + az * z
    cx = ay * z - az * y
    cy = az * x - ax * z
    cz = ax * y - ay * x
    rx = x * c + cx * s + ax * dot * (1.0 - c)
    ry = y * c + cy * s + ay * dot * (1.0 - c)
    rz = z * c + cz * s + az * dot * (1.0 - c)
    return rx, ry, rz

def exact_height(xc: np.ndarray, yc: np.ndarray, t_sec: float, alpha: float) -> np.ndarray:
    x, y, z = unit_vector(xc, yc)
    x0, y0, z0 = rotate_about_axis(x, y, z, alpha, -OMEGA * t_sec)
    cx, cy, cz = unit_vector(BELL_LON_DEG, BELL_LAT_DEG)
    dot = np.clip(x0 * cx + y0 * cy + z0 * cz, -1.0, 1.0)
    radius = EARTH_RADIUS * np.arccos(dot)
    height = np.zeros_like(radius, dtype=np.float64)
    inside = radius < BELL_RADIUS
    height[inside] = 0.5 * BELL_HEIGHT * (
        1.0 + np.cos(math.pi * radius[inside] / BELL_RADIUS)
    )
    return height

def bell_center_at_time(alpha: float, t_sec: float) -> tuple[float, float]:
    x, y, z = unit_vector(BELL_LON_DEG, BELL_LAT_DEG)
    xr, yr, zr = rotate_about_axis(x, y, z, alpha, OMEGA * t_sec)
    lon = math.degrees(math.atan2(float(yr), float(xr))) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, float(zr)))))
    return lon, lat

def unwrap_longitudes(lons_deg: np.ndarray) -> np.ndarray:
    return np.degrees(np.unwrap(np.radians(np.asarray(lons_deg, dtype=np.float64))))

def wrap_longitude_delta(delta_deg: float) -> float:
    return ((delta_deg + 180.0) % 360.0) - 180.0

def great_circle_distance_km(
    lon1_deg: float,
    lat1_deg: float,
    lon2_deg: float,
    lat2_deg: float,
) -> float:
    lon1 = math.radians(lon1_deg)
    lat1 = math.radians(lat1_deg)
    lon2 = math.radians(lon2_deg)
    lat2 = math.radians(lat2_deg)
    sin_dlat = math.sin(0.5 * (lat2 - lat1))
    sin_dlon = math.sin(0.5 * (lon2 - lon1))
    a = sin_dlat * sin_dlat + math.cos(lat1) * math.cos(lat2) * sin_dlon * sin_dlon
    return EARTH_RADIUS * 2.0 * math.asin(min(1.0, math.sqrt(max(a, 0.0)))) / 1000.0

def local_window_grid(
    xc: np.ndarray,
    yc: np.ndarray,
    field: np.ndarray,
    center_lon_deg: float,
    center_lat_deg: float,
    half_window_deg: float = 30.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon = np.asarray(xc[0, :], dtype=np.float64)
    lat = np.asarray(yc[:, 0], dtype=np.float64)
    lon_rel = ((lon - center_lon_deg + 180.0) % 360.0) - 180.0
    lat_rel = lat - center_lat_deg

    lon_order = np.argsort(lon_rel)
    lat_order = np.argsort(lat_rel)
    lon_sorted = lon_rel[lon_order]
    lat_sorted = lat_rel[lat_order]

    lon_mask = np.abs(lon_sorted) <= half_window_deg
    lat_mask = np.abs(lat_sorted) <= half_window_deg

    lon_sel = lon_sorted[lon_mask]
    lat_sel = lat_sorted[lat_mask]
    field_sorted = np.asarray(field, dtype=np.float64)[np.ix_(lat_order, lon_order)]
    field_sel = field_sorted[np.ix_(lat_mask, lon_mask)]

    lon_grid, lat_grid = np.meshgrid(lon_sel, lat_sel)
    return lon_grid, lat_grid, field_sel

def get_weights(run_dir: Path, yc: np.ndarray) -> np.ndarray:
    for name in ("RAC", "Area"):
        meta_path = run_dir / f"{name}.meta"
        data_path = run_dir / f"{name}.data"
        if meta_path.exists() and data_path.exists():
            weights = np.asarray(read_mds_field(run_dir, name), dtype=np.float64)
            if weights.ndim == 3:
                weights = weights[0]
            return weights

    lat = np.deg2rad(np.asarray(yc, dtype=np.float64))
    return np.cos(lat)

def error_metrics(model: np.ndarray, exact: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    model = np.asarray(model, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)

    diff = model - exact
    abs_exact = np.abs(exact)
    abs_diff = np.abs(diff)
    max_exact = float(np.max(abs_exact))
    weighted_exact_l1 = float(np.sum(weights * abs_exact))
    weighted_exact_l2 = float(np.sum(weights * exact * exact))
    weighted_diff_l1 = float(np.sum(weights * abs_diff))
    weighted_diff_l2 = float(np.sum(weights * diff * diff))

    if max_exact == 0.0:
        max_exact = 1.0
    if weighted_exact_l1 == 0.0:
        weighted_exact_l1 = 1.0
    if weighted_exact_l2 == 0.0:
        weighted_exact_l2 = 1.0

    model_max = float(np.max(model))
    model_min = float(np.min(model))
    exact_max = float(np.max(exact))
    exact_min = float(np.min(exact))
    mass_exact = float(np.sum(weights * exact))
    mass_model = float(np.sum(weights * model))

    return {
        "l1": weighted_diff_l1 / weighted_exact_l1,
        "l2": math.sqrt(weighted_diff_l2 / weighted_exact_l2),
        "linf": float(np.max(abs_diff) / max_exact),
        "rel_max": (model_max - exact_max) / max_exact,
        "rel_min": (model_min - exact_min) / max_exact,
        "mass_rel": (mass_model - mass_exact) / abs(mass_exact) if mass_exact != 0.0 else 0.0,
        "model_max": model_max,
        "model_min": model_min,
        "exact_max": exact_max,
        "exact_min": exact_min,
    }

def field_extent(xc: np.ndarray, yc: np.ndarray) -> tuple[float, float, float, float]:
    lon = np.asarray(xc[0, :], dtype=np.float64)
    lat = np.asarray(yc[:, 0], dtype=np.float64)
    lon_step = float(np.median(np.diff(lon))) if lon.size > 1 else 1.0
    lat_step = float(np.median(np.diff(lat))) if lat.size > 1 else 1.0
    return (
        float(lon[0] - 0.5 * lon_step),
        float(lon[-1] + 0.5 * lon_step),
        float(lat[0] - 0.5 * lat_step),
        float(lat[-1] + 0.5 * lat_step),
    )

def plot_field(ax, field: np.ndarray, xc: np.ndarray, yc: np.ndarray, title: str, cmap: str, **kwargs):
    extent = field_extent(xc, yc)
    image = ax.imshow(
        field,
        origin="lower",
        extent=extent,
        interpolation="bilinear",
        cmap=cmap,
        aspect="auto",
        **kwargs,
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    return image

def save_final_comparison(
    output_dir: Path,
    xc: np.ndarray,
    yc: np.ndarray,
    model: np.ndarray,
    exact: np.ndarray,
    error: np.ndarray,
    it_used: int,
    days: float,
    alpha: float,
    filename: str,
    validation_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    center_lon, center_lat = bell_center_at_time(alpha, days * DAY)
    lon_grid, lat_grid, exact_local = local_window_grid(
        xc, yc, exact, center_lon, center_lat, half_window_deg=30.0
    )
    _, _, model_local = local_window_grid(
        xc, yc, model, center_lon, center_lat, half_window_deg=30.0
    )
    _, _, error_local = local_window_grid(
        xc, yc, error, center_lon, center_lat, half_window_deg=30.0
    )

    max_abs_error = float(np.max(np.abs(error_local)))
    available_error_levels = int(max_abs_error // ERROR_CONTOUR_STEP)
    max_error_level = min(
        available_error_levels,
        MAX_ERROR_CONTOURS_PER_SIGN,
    )
    positive_error_levels = ERROR_CONTOUR_STEP * np.arange(
        1, max_error_level + 1, dtype=np.float64
    )
    negative_error_levels = -positive_error_levels[::-1]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6), constrained_layout=True)

    ax = axes[0]
    ax.contour(
        lon_grid,
        lat_grid,
        exact_local,
        levels=HEIGHT_CONTOUR_LEVELS,
        colors="k",
        linestyles="--",
        linewidths=1.0,
    )
    ax.contour(
        lon_grid,
        lat_grid,
        model_local,
        levels=HEIGHT_CONTOUR_LEVELS,
        colors="k",
        linestyles="-",
        linewidths=0.8,
    )
    show_model_zero = float(np.min(model_local)) < -ZERO_CONTOUR_TOLERANCE
    if show_model_zero:
        model_zero = np.where(
            np.abs(model_local) < ZERO_CONTOUR_TOLERANCE,
            ZERO_CONTOUR_TOLERANCE,
            model_local,
        )
        ax.contour(
            lon_grid,
            lat_grid,
            model_zero,
            levels=[0.0],
            colors="k",
            linestyles=":",
            linewidths=1.0,
        )
    ax.set_title(f"TC1 height contours at t = {days:.2f} days (it = {it_used})")
    ax.set_xlabel("Longitude offset [deg]")
    ax.set_ylabel("Latitude offset [deg]")
    ax.set_xlim(-30.0, 30.0)
    ax.set_ylim(-30.0, 30.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(-30, 31, 15))
    ax.set_yticks(np.arange(-30, 31, 15))
    field_handles = [
        Line2D([0], [0], color="k", lw=1.0, linestyle="--", label="Exact"),
        Line2D([0], [0], color="k", lw=0.9, linestyle="-", label="Model"),
    ]
    if show_model_zero:
        field_handles.append(
            Line2D(
                [0],
                [0],
                color="k",
                lw=1.0,
                linestyle=":",
                label="Model zero contour",
            )
        )
    ax.legend(
        handles=field_handles,
        frameon=False,
        loc="upper right",
    )

    ax = axes[1]
    if negative_error_levels.size:
        ax.contour(
            lon_grid,
            lat_grid,
            error_local,
            levels=negative_error_levels,
            colors="k",
            linestyles="--",
            linewidths=0.8,
        )
    if positive_error_levels.size:
        ax.contour(
            lon_grid,
            lat_grid,
            error_local,
            levels=positive_error_levels,
            colors="k",
            linestyles="-",
            linewidths=0.8,
        )
    show_error_zero = (
        float(np.min(error_local)) < -ZERO_CONTOUR_TOLERANCE
        and float(np.max(error_local)) > ZERO_CONTOUR_TOLERANCE
    )
    if show_error_zero:
        error_zero = np.where(
            np.abs(error_local) < ZERO_CONTOUR_TOLERANCE,
            ZERO_CONTOUR_TOLERANCE,
            error_local,
        )
        ax.contour(
            lon_grid,
            lat_grid,
            error_zero,
            levels=[0.0],
            colors="k",
            linewidths=1.2,
        )
    if positive_error_levels.size == 0:
        ax.text(
            0.5,
            0.5,
            f"No +/-{ERROR_CONTOUR_STEP:g} m contours\n"
            f"max |error| = {max_abs_error:.3e} m",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    elif available_error_levels > MAX_ERROR_CONTOURS_PER_SIGN:
        displayed_limit = ERROR_CONTOUR_STEP * MAX_ERROR_CONTOURS_PER_SIGN
        ax.text(
            0.02,
            0.02,
            f"Contours capped at +/-{displayed_limit:g} m\n"
            f"max |error| = {max_abs_error:.3e} m",
            ha="left",
            va="bottom",
            transform=ax.transAxes,
        )
    ax.set_title("Signed error: model - exact")
    ax.set_xlabel("Longitude offset [deg]")
    ax.set_ylabel("Latitude offset [deg]")
    ax.set_xlim(-30.0, 30.0)
    ax.set_ylim(-30.0, 30.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(-30, 31, 15))
    ax.set_yticks(np.arange(-30, 31, 15))
    if positive_error_levels.size:
        error_handles = [
            Line2D([0], [0], color="k", lw=0.8, linestyle="-", label="Positive"),
            Line2D([0], [0], color="k", lw=0.8, linestyle="--", label="Negative"),
        ]
        if show_error_zero:
            error_handles.append(
                Line2D([0], [0], color="k", lw=1.2, linestyle="-", label="Zero")
            )
        ax.legend(
            handles=error_handles,
            frameon=False,
            loc="upper right",
        )

    fig.suptitle(
        f"Williamson TC1 comparison | {validation_label}",
        fontsize=13,
    )
    _save_figure(fig, output_dir / filename)

def save_metrics_plot(
    output_dir: Path,
    rows: list[dict[str, float]],
    validation_label: str,
) -> None:
    days = np.array([row["day"] for row in rows], dtype=np.float64)
    l1 = np.array([row["l1"] for row in rows], dtype=np.float64)
    l2 = np.array([row["l2"] for row in rows], dtype=np.float64)
    linf = np.array([row["linf"] for row in rows], dtype=np.float64)
    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)

    ax.plot(days, l1, color="k", lw=1.5, linestyle="-", label=r"$L_1$ Error")
    ax.plot(days, l2, color="k", lw=1.5, linestyle="--", label=r"$L_2$ Error")
    ax.plot(days, linf, color="k", lw=1.5, linestyle=":", label=r"$L_\infty$ Error")
    ax.set_title(f"TC1 normalized error growth | {validation_label}")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Normalized Errors")
    ax.set_xlim(float(days[0]), float(days[-1]))
    max_err = float(np.max(np.concatenate([l1, l2, linf])))
    ax.set_ylim(0.0, max_err * 1.08 if max_err > 0.0 else 1.0)
    ax.set_xticks(np.arange(0, max(13, int(days[-1]) + 1), 2))
    ax.grid(False)
    ax.legend(frameon=True, loc="upper left", fancybox=False, edgecolor="k")

    _save_figure(fig, output_dir / "tc1_error_metrics.pdf")

def write_validation_report(
    output_dir: Path,
    stationary_output: bool,
    max_field_change: float,
    model_peak_path_km: float,
    exact_peak_path_km: float,
    rows: list[dict[str, float]],
) -> None:
    status = "INVALID: height is stationary" if stationary_output else "MOVING: inspect errors"
    final = rows[-1]
    lines = [
        f"status: {status}",
        f"max_field_change_from_initial_m: {max_field_change:.12e}",
        f"model_peak_max_distance_from_start_km: {model_peak_path_km:.12e}",
        f"exact_peak_max_distance_from_start_km: {exact_peak_path_km:.12e}",
        f"final_l1: {final['l1']:.12e}",
        f"final_l2: {final['l2']:.12e}",
        f"final_linf: {final['linf']:.12e}",
    ]
    if stationary_output:
        lines.extend(
            [
                "",
                "The saved passive-salinity height field did not advect.",
                "A near-zero day-12 error is not a pass because the exact bell",
                "returns to its initial position after one full rotation.",
            ]
        )
    (output_dir / "tc1_validation.txt").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

def save_peak_plot(output_dir: Path, rows: list[dict[str, float]]) -> None:
    days = np.array([row["day"] for row in rows], dtype=np.float64)
    model_lon = unwrap_longitudes([row["model_lon"] for row in rows])
    exact_lon = unwrap_longitudes([row["exact_lon"] for row in rows])
    model_lat = np.array([row["model_lat"] for row in rows], dtype=np.float64)
    exact_lat = np.array([row["exact_lat"] for row in rows], dtype=np.float64)
    peak_dist = np.array([row["peak_distance_km"] for row in rows], dtype=np.float64)

    fig, axes = plt.subplots(3, 1, figsize=(8.8, 8.8), constrained_layout=True)

    ax = axes[0]
    ax.plot(days, exact_lon, color="k", lw=1.5, linestyle="--", label="Exact")
    ax.plot(days, model_lon, color="k", lw=1.3, linestyle="-", label="Model")
    ax.set_title("Bell peak longitude")
    ax.set_ylabel("Unwrapped lon [deg]")
    ax.legend(frameon=True, loc="upper left", fancybox=False, edgecolor="k")

    ax = axes[1]
    ax.plot(days, exact_lat, color="k", lw=1.5, linestyle="--", label="Exact")
    ax.plot(days, model_lat, color="k", lw=1.3, linestyle="-", label="Model")
    ax.set_title("Bell peak latitude")
    ax.set_ylabel("Latitude [deg]")

    ax = axes[2]
    ax.plot(days, peak_dist, color="k", lw=1.5)
    ax.set_title("Great-circle distance between model and exact peak")
    ax.set_ylabel("Distance [km]")
    ax.set_xlabel("Time [days]")

    _save_figure(fig, output_dir / "tc1_peak_trajectory.pdf")

def write_csv(output_dir: Path, rows: list[dict[str, float]]) -> None:
    csv_path = output_dir / "tc1_error_metrics.csv"
    fieldnames = [
        "iteration",
        "day",
        "l1",
        "l2",
        "linf",
        "rel_max",
        "rel_min",
        "mass_rel",
        "model_max",
        "model_min",
        "exact_max",
        "exact_min",
        "max_change_from_initial",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

def write_peak_csv(output_dir: Path, rows: list[dict[str, float]]) -> None:
    csv_path = output_dir / "tc1_peak_trajectory.csv"
    fieldnames = [
        "iteration",
        "day",
        "model_lon",
        "model_lat",
        "exact_lon",
        "exact_lat",
        "lon_offset_deg",
        "lat_offset_deg",
        "peak_distance_km",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

# %%
def run_error_study(
    run_dir: Path | None = None,
    alpha: float | None = None,
    output_dir: Path | None = None,
    field: str = "S",
) -> dict[str, object]:
    run_dir = RUN_DIR if run_dir is None else Path(run_dir)
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    alpha = read_alpha(run_dir) if alpha is None else alpha
    delta_t = read_delta_t(run_dir)
    output_dir = (
        output_dir.expanduser().resolve()
        if output_dir is not None
        else OUTPUT_DIR
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    iterations = discover_iterations(run_dir, field)
    if not iterations:
        raise FileNotFoundError(f"No {field} output found in {run_dir}")

    xc = np.asarray(read_mds_field(run_dir, "XC"), dtype=np.float64)
    yc = np.asarray(read_mds_field(run_dir, "YC"), dtype=np.float64)
    if xc.ndim == 3:
        xc = xc[0]
    if yc.ndim == 3:
        yc = yc[0]
    weights = np.asarray(get_weights(run_dir, yc), dtype=np.float64)
    if weights.ndim == 3:
        weights = weights[0]

    rows: list[dict[str, float]] = []
    peak_rows: list[dict[str, float]] = []
    initial_model: np.ndarray | None = None
    initial_exact: np.ndarray | None = None
    final_model: np.ndarray | None = None
    final_exact: np.ndarray | None = None
    diagnostic_model: np.ndarray | None = None
    diagnostic_exact: np.ndarray | None = None
    diagnostic_iteration = iterations[0]
    diagnostic_day = 0.0
    diagnostic_distance_km = -1.0
    max_field_change = 0.0

    for iteration in iterations:
        model = np.asarray(read_mds_field(run_dir, field, iteration), dtype=np.float64)
        if model.ndim == 3:
            if model.shape[0] != 1:
                raise ValueError(
                    f"TC1 expects one vertical level, found shape {model.shape}"
                )
            model = model[0]
        exact = exact_height(xc, yc, iteration * delta_t, alpha)
        metrics = error_metrics(model, exact, weights)
        day = iteration * delta_t / DAY

        model_peak_j, model_peak_i = np.unravel_index(int(np.argmax(model)), model.shape)
        model_peak_lon = float(xc[model_peak_j, model_peak_i])
        model_peak_lat = float(yc[model_peak_j, model_peak_i])
        exact_peak_lon, exact_peak_lat = bell_center_at_time(alpha, iteration * delta_t)

        peak_rows.append(
            {
                "iteration": float(iteration),
                "day": float(day),
                "model_lon": model_peak_lon,
                "model_lat": model_peak_lat,
                "exact_lon": exact_peak_lon,
                "exact_lat": exact_peak_lat,
                "lon_offset_deg": wrap_longitude_delta(model_peak_lon - exact_peak_lon),
                "lat_offset_deg": model_peak_lat - exact_peak_lat,
                "peak_distance_km": great_circle_distance_km(
                    model_peak_lon,
                    model_peak_lat,
                    exact_peak_lon,
                    exact_peak_lat,
                ),
            }
        )

        row = {
            "iteration": float(iteration),
            "day": float(day),
            **metrics,
        }
        if initial_model is None:
            initial_model = model
            initial_exact = exact
        field_change = float(np.max(np.abs(model - initial_model)))
        max_field_change = max(max_field_change, field_change)
        row["max_change_from_initial"] = field_change
        rows.append(row)
        exact_distance_from_start = great_circle_distance_km(
            BELL_LON_DEG,
            BELL_LAT_DEG,
            exact_peak_lon,
            exact_peak_lat,
        )
        if exact_distance_from_start > diagnostic_distance_km:
            diagnostic_distance_km = exact_distance_from_start
            diagnostic_model = model
            diagnostic_exact = exact
            diagnostic_iteration = iteration
            diagnostic_day = day
        final_model = model
        final_exact = exact

    assert final_model is not None and final_exact is not None
    assert initial_model is not None and initial_exact is not None
    assert diagnostic_model is not None and diagnostic_exact is not None

    model_start = peak_rows[0]
    model_peak_path_km = max(
        great_circle_distance_km(
            model_start["model_lon"],
            model_start["model_lat"],
            row["model_lon"],
            row["model_lat"],
        )
        for row in peak_rows
    )
    exact_peak_path_km = max(
        great_circle_distance_km(
            model_start["exact_lon"],
            model_start["exact_lat"],
            row["exact_lon"],
            row["exact_lat"],
        )
        for row in peak_rows
    )
    stationary_output = (
        max_field_change < STATIONARY_FIELD_TOLERANCE
        and model_peak_path_km < 100.0
        and exact_peak_path_km > 1000.0
    )
    validation_label = (
        "INVALID: stationary height"
        if stationary_output
        else "Height moves; check error norms"
    )

    write_csv(output_dir, rows)
    save_metrics_plot(output_dir, rows, validation_label)
    write_peak_csv(output_dir, peak_rows)
    save_peak_plot(output_dir, peak_rows)
    write_validation_report(
        output_dir,
        stationary_output,
        max_field_change,
        model_peak_path_km,
        exact_peak_path_km,
        rows,
    )
    final_it = iterations[-1]
    final_error = final_model - final_exact
    save_final_comparison(
        output_dir,
        xc,
        yc,
        final_model,
        final_exact,
        final_error,
        final_it,
        iterations[-1] * delta_t / DAY,
        alpha,
        "tc1_final_comparison.pdf",
        validation_label,
    )
    if diagnostic_iteration != final_it:
        save_final_comparison(
            output_dir,
            xc,
            yc,
            diagnostic_model,
            diagnostic_exact,
            diagnostic_model - diagnostic_exact,
            diagnostic_iteration,
            diagnostic_day,
            alpha,
            "tc1_worst_comparison.pdf",
            validation_label,
        )

    init_error = initial_model - initial_exact

    print(f"TC1 error study for {run_dir}")
    print(f"  alpha = {alpha:.6f} rad")
    print(f"  deltaT = {delta_t:.6f} s")
    print(f"  analyzed iterations = {iterations}")
    print(
        "  initial-time max abs error = "
        f"{np.max(np.abs(init_error)):.6e} m"
    )
    print(
        "  final-time max abs error = "
        f"{np.max(np.abs(final_error)):.6e} m"
    )
    print(
        "  final normalized errors: "
        f"L1={rows[-1]['l1']:.6e}, "
        f"L2={rows[-1]['l2']:.6e}, "
        f"Linf={rows[-1]['linf']:.6e}"
    )
    print(
        "  final relative extrema: "
        f"max={rows[-1]['rel_max']:.6e}, "
        f"min={rows[-1]['rel_min']:.6e}, "
        f"mass={rows[-1]['mass_rel']:.6e}"
    )
    model_lon_unwrapped = unwrap_longitudes([row["model_lon"] for row in peak_rows])
    exact_lon_unwrapped = unwrap_longitudes([row["exact_lon"] for row in peak_rows])
    print(
        "  peak-track span: "
        f"model_lon={np.ptp(model_lon_unwrapped):.6e} deg, "
        f"exact_lon={np.ptp(exact_lon_unwrapped):.6e} deg, "
        f"max_distance={max(row['peak_distance_km'] for row in peak_rows):.6e} km"
    )
    print(
        "  field validation: "
        f"{validation_label}; "
        f"max_change={max_field_change:.6e} m, "
        f"model_peak_path={model_peak_path_km:.6e} km"
    )
    print(f"  outputs written to {output_dir}")

    if np.max(np.abs(init_error)) > 1.0e-3:
        print(
            "  note: the initial model height does not exactly match the analytic "
            "initial condition; check input generation and file wiring."
        )

    if stationary_output:
        print(
            "  INVALID TC1 OUTPUT: the height field is stationary while the analytic bell "
            "travels around the sphere. A small day-12 error is misleading "
            "because the exact solution returns to its starting point."
        )

    return {
        "run_dir": run_dir,
        "output_dir": output_dir,
        "alpha": alpha,
        "delta_t": delta_t,
        "iterations": iterations,
        "rows": rows,
        "peak_rows": peak_rows,
        "stationary_output": stationary_output,
    }

RESULTS = run_error_study(RUN_DIR, output_dir=OUTPUT_DIR)

# %%
