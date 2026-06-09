#!/usr/bin/env python3
"""Williamson TC1 error analysis for MITgcm output.

This replaces:
    tc1_paper_plots.py
    tc1_paper_plots_fancy.py
    tc1_error_study.py

It writes only:
    1) paper-style cosine-bell overlay
    2) paper-style signed-error contours
    3) normalized-error line plot
    4) peak-trajectory line plot
    5) CSV, TXT, and LaTeX tables

Run:
    python tc1_error_analysis.py
    python tc1_error_analysis.py /path/to/run_alpha_0
"""
#%%
from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# ------------------------- easy settings -------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
CASE_DIR = SCRIPT_DIR.parents[1] if len(SCRIPT_DIR.parents) > 1 else SCRIPT_DIR

RUN_DIR = CASE_DIR / "run_alpha_0.05"    # edit this directly
FIELD = "S"                          # TC1 cosine-bell height stored as passive tracer

EARTH_RADIUS = 6_371_000.0
DAY = 86_400.0
OMEGA = 2.0 * math.pi / (12.0 * DAY)  # one revolution in 12 days
BELL_LON_DEG = 270.0
BELL_LAT_DEG = 0.0
BELL_HEIGHT = 1000.0
BELL_RADIUS = EARTH_RADIUS / 3.0
DEFAULT_ALPHA = 0.0

WINDOW_PADDING_DEG = 4.0
REFERENCE_HEIGHT_LEVELS = np.arange(160.0, 801.0, 160.0)
MODEL_HEIGHT_LEVELS = np.arange(80.0, 801.0, 80.0)
MODEL_TAIL_LEVELS = np.array([0.01, 0.02, 0.04])
THIN_LINE = 0.65
ZERO_LINE = 1.15
DPI = 400
# -----------------------------------------------------------------

def read_text_value(path: Path, pattern: str) -> float | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(pattern, text, re.IGNORECASE)
    if match is None:
        return None
    return float(match.group(1).replace("D", "e").replace("d", "e"))

def read_delta_t(run_dir: Path) -> float:
    for path in (run_dir / "data", run_dir.parent / "input" / "data"):
        value = read_text_value(path, r"deltaT\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return value
    return 60.0

def read_alpha(run_dir: Path) -> float:
    for path in (run_dir / "data.mypackage", run_dir.parent / "input" / "data.mypackage"):
        value = read_text_value(path, r"myPa_param1\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return value
    return DEFAULT_ALPHA

def parse_mds_meta(meta_path: Path) -> dict[str, object]:
    text = meta_path.read_text(encoding="utf-8", errors="ignore")
    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_values = [int(v) for v in re.findall(r"[-+]?\d+", dim_block)]
    dims = []
    for i in range(n_dims):
        n_global, start, end = dim_values[3 * i : 3 * i + 3]
        dims.append({"global": n_global, "start": start, "end": end})
    prec = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", text)
    nrecords = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", text)
    return {
        "n_dims": n_dims,
        "dims": dims,
        "prec": prec.group(1).lower() if prec else "float32",
        "nrecords": int(nrecords.group(1)) if nrecords else 1,
    }

def dtype_from_meta(meta: dict[str, object]) -> np.dtype:
    prec = str(meta["prec"])
    return np.dtype(">f8" if prec in {"float64", "real*8"} else ">f4")

def read_mds_field(run_dir: Path, name: str, iteration: int | None = None) -> np.ndarray:
    if iteration is None:
        meta_path = run_dir / f"{name}.meta"
    else:
        meta_path = run_dir / f"{name}.{iteration:010d}.meta"
    data_path = meta_path.with_suffix(".data")

    if not meta_path.exists() or not data_path.exists():
        raise FileNotFoundError(f"Missing {name} in {run_dir}")

    meta = parse_mds_meta(meta_path)
    raw = np.fromfile(data_path, dtype=dtype_from_meta(meta))
    dims = meta["dims"]
    n_dims = int(meta["n_dims"])
    nrecords = int(meta["nrecords"])

    if n_dims == 2:
        nx = dims[0]["global"]
        ny = dims[1]["global"]
        shape = (nrecords, ny, nx) if nrecords > 1 else (ny, nx)
    elif n_dims == 3:
        nx = dims[0]["global"]
        ny = dims[1]["global"]
        nz = dims[2]["global"]
        shape = (nrecords, nz, ny, nx) if nrecords > 1 else (nz, ny, nx)
    else:
        raise ValueError(f"Unsupported nDims={n_dims}: {meta_path}")

    data = raw.reshape(shape)
    return np.squeeze(data[0] if nrecords > 1 else data)

def to_2d(field: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    if field.ndim == 3 and field.shape[0] == 1:
        return field[0]
    if field.ndim != 2:
        raise ValueError(f"Expected 2D or 1-level 3D field, got {field.shape}")
    return field

def discover_iterations(run_dir: Path, field: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(field)}\.(\d{{10}})\.meta$")
    out = []
    for path in run_dir.glob(f"{field}.*.meta"):
        match = pattern.match(path.name)
        if match:
            out.append(int(match.group(1)))
    return sorted(set(out))

def unit_vector(lon_deg: np.ndarray | float, lat_deg: np.ndarray | float):
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    cos_lat = np.cos(lat)
    return cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)

def rotate_about_axis(x, y, z, alpha: float, angle: float):
    ax = -math.sin(alpha)
    ay = 0.0
    az = math.cos(alpha)
    c = math.cos(angle)
    s = math.sin(angle)
    dot = ax * x + ay * y + az * z
    cx = ay * z - az * y
    cy = az * x - ax * z
    cz = ax * y - ay * x
    return (
        x * c + cx * s + ax * dot * (1.0 - c),
        y * c + cy * s + ay * dot * (1.0 - c),
        z * c + cz * s + az * dot * (1.0 - c),
    )

def exact_height(xc: np.ndarray, yc: np.ndarray, t_sec: float, alpha: float) -> np.ndarray:
    x, y, z = unit_vector(xc, yc)
    x0, y0, z0 = rotate_about_axis(x, y, z, alpha, -OMEGA * t_sec)
    cx, cy, cz = unit_vector(BELL_LON_DEG, BELL_LAT_DEG)
    dot = np.clip(x0 * cx + y0 * cy + z0 * cz, -1.0, 1.0)
    radius = EARTH_RADIUS * np.arccos(dot)
    height = np.zeros_like(radius, dtype=np.float64)
    inside = radius < BELL_RADIUS
    height[inside] = 0.5 * BELL_HEIGHT * (1.0 + np.cos(math.pi * radius[inside] / BELL_RADIUS))
    return height

def bell_center_at_time(alpha: float, t_sec: float) -> tuple[float, float]:
    x, y, z = unit_vector(BELL_LON_DEG, BELL_LAT_DEG)
    xr, yr, zr = rotate_about_axis(x, y, z, alpha, OMEGA * t_sec)
    lon = math.degrees(math.atan2(float(yr), float(xr))) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, float(zr)))))
    return lon, lat

def unwrap_longitudes(lons_deg) -> np.ndarray:
    return np.degrees(np.unwrap(np.radians(np.asarray(lons_deg, dtype=np.float64))))

def wrap_longitude_delta(delta_deg: float) -> float:
    return ((delta_deg + 180.0) % 360.0) - 180.0

def great_circle_distance_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lon1 = math.radians(lon1)
    lat1 = math.radians(lat1)
    lon2 = math.radians(lon2)
    lat2 = math.radians(lat2)
    sin_dlat = math.sin(0.5 * (lat2 - lat1))
    sin_dlon = math.sin(0.5 * (lon2 - lon1))
    a = sin_dlat**2 + math.cos(lat1) * math.cos(lat2) * sin_dlon**2
    return 2.0 * EARTH_RADIUS * math.asin(min(1.0, math.sqrt(max(a, 0.0)))) / 1000.0

def get_weights(run_dir: Path, yc: np.ndarray) -> np.ndarray:
    for name in ("RAC", "Area"):
        if (run_dir / f"{name}.meta").exists():
            return to_2d(read_mds_field(run_dir, name))
    return np.cos(np.deg2rad(yc))

def error_metrics(model: np.ndarray, exact: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    diff = model - exact
    abs_exact = np.abs(exact)
    abs_diff = np.abs(diff)

    max_exact = max(float(np.max(abs_exact)), 1.0)
    exact_l1 = max(float(np.sum(weights * abs_exact)), 1.0)
    exact_l2 = max(float(np.sum(weights * exact * exact)), 1.0)

    mass_exact = float(np.sum(weights * exact))
    mass_model = float(np.sum(weights * model))

    return {
        "normalized_l1": float(np.sum(weights * abs_diff) / exact_l1),
        "normalized_l2": float(math.sqrt(np.sum(weights * diff * diff) / exact_l2)),
        "normalized_linf": float(np.max(abs_diff) / max_exact),
        "relative_mass_error": float((mass_model - mass_exact) / abs(mass_exact)) if mass_exact != 0.0 else 0.0,
        "model_max": float(np.max(model)),
        "model_min": float(np.min(model)),
        "exact_max": float(np.max(exact)),
        "exact_min": float(np.min(exact)),
    }

def local_window(xc: np.ndarray, yc: np.ndarray, field: np.ndarray, center_lon: float, center_lat: float, half_width: float):
    lon_axis = ((xc[0, :] - center_lon + 180.0) % 360.0) - 180.0
    lat_axis = yc[:, 0] - center_lat
    lon_order = np.argsort(lon_axis)
    lat_order = np.argsort(lat_axis)
    lon_sorted = lon_axis[lon_order]
    lat_sorted = lat_axis[lat_order]
    lon_mask = np.abs(lon_sorted) <= half_width
    lat_mask = np.abs(lat_sorted) <= half_width
    field_sorted = field[np.ix_(lat_order, lon_order)]
    lon, lat = np.meshgrid(lon_sorted[lon_mask], lat_sorted[lat_mask])
    return lon, lat, field_sorted[np.ix_(lat_mask, lon_mask)]

def support_window(xc: np.ndarray, yc: np.ndarray, exact: np.ndarray, center_lon: float, center_lat: float) -> float:
    support = exact > 0.0
    lon = ((xc[support] - center_lon + 180.0) % 360.0) - 180.0
    lat = yc[support] - center_lat
    return float(max(np.max(np.abs(lon)), np.max(np.abs(lat))) + WINDOW_PADDING_DEG)

def levels_inside(levels: np.ndarray, field: np.ndarray) -> np.ndarray:
    lo = float(np.nanmin(field))
    hi = float(np.nanmax(field))
    return np.array([v for v in levels if lo <= v <= hi], dtype=np.float64)

def nice_error_step(max_abs: float) -> float:
    if max_abs <= 0.0 or not math.isfinite(max_abs):
        return 1.0
    raw = max_abs / 6.0
    scale = 10.0 ** math.floor(math.log10(raw))
    frac = raw / scale
    if frac <= 1.0:
        return scale
    if frac <= 2.0:
        return 2.0 * scale
    if frac <= 5.0:
        return 5.0 * scale
    return 10.0 * scale

def plot_cosine_bell_overlay(output_dir: Path, lon, lat, model, exact, day: float) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 5.2))

    c_model = None
    model_levels = levels_inside(MODEL_HEIGHT_LEVELS, model)
    if model_levels.size:
        c_model = ax.contour(lon, lat, model, levels=model_levels, colors="black", linestyles="-", linewidths=THIN_LINE)
        ax.clabel(c_model, levels=c_model.levels[1::2], fmt="%g", fontsize=7, inline=True)

    tail_levels = levels_inside(MODEL_TAIL_LEVELS, model)
    if tail_levels.size:
        c_tail = ax.contour(lon, lat, model, levels=tail_levels, colors="black", linestyles=":", linewidths=THIN_LINE)
        ax.clabel(c_tail, fmt="%g", fontsize=6, inline=True)

    exact_levels = levels_inside(REFERENCE_HEIGHT_LEVELS, exact)
    if exact_levels.size:
        c_exact = ax.contour(lon, lat, exact, levels=exact_levels, colors="black", linestyles="--", linewidths=THIN_LINE)
        ax.clabel(c_exact, levels=c_exact.levels[::2], fmt="%g", fontsize=7, inline=True)

    ax.legend(
        handles=[
            Line2D([0], [0], color="black", linestyle="--", lw=THIN_LINE, label="Exact"),
            Line2D([0], [0], color="black", linestyle="-", lw=THIN_LINE, label="MITgcm"),
            Line2D([0], [0], color="black", linestyle=":", lw=THIN_LINE, label="MITgcm tails"),
        ],
        frameon=False,
        loc="upper right",
        fontsize=8,
    )
    ax.set_title(f"Cosine Bell Advection - MITgcm - {day:.0f} days", fontsize=10)
    ax.set_xlabel("Longitude offset [deg]")
    ax.set_ylabel("Latitude offset [deg]")
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(direction="out", top=True, right=True)
    fig.tight_layout()

    for suffix in ("pdf", "png"):
        fig.savefig(output_dir / f"tc1_cosine_bell_overlay.{suffix}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def plot_signed_error(output_dir: Path, lon, lat, error) -> None:
    max_abs = float(np.nanmax(np.abs(error)))
    step = nice_error_step(max_abs)
    max_level = int(math.floor(max_abs / step))
    pos = step * np.arange(1, max_level + 1)
    neg = -pos[::-1]

    fig, ax = plt.subplots(figsize=(5.2, 5.2))

    if neg.size:
        c_neg = ax.contour(lon, lat, error, levels=neg, colors="black", linestyles="--", linewidths=THIN_LINE)
        ax.clabel(c_neg, levels=c_neg.levels[::2], fmt="%g", fontsize=7, inline=True)
    if pos.size:
        c_pos = ax.contour(lon, lat, error, levels=pos, colors="black", linestyles="-", linewidths=THIN_LINE)
        ax.clabel(c_pos, fmt="%g", fontsize=7, inline=True)

    zero_tol = max(1.0e-10, 0.01 * step)
    if float(np.nanmin(error)) < -zero_tol and float(np.nanmax(error)) > zero_tol:
        error_zero = np.where(np.abs(error) < zero_tol, zero_tol, error)
        ax.contour(lon, lat, error_zero, levels=[0.0], colors="black", linestyles="-", linewidths=ZERO_LINE)

    ax.legend(
        handles=[
            Line2D([0], [0], color="black", linestyle="-", lw=THIN_LINE, label="Positive"),
            Line2D([0], [0], color="black", linestyle="--", lw=THIN_LINE, label="Negative"),
            Line2D([0], [0], color="black", linestyle="-", lw=ZERO_LINE, label="Zero"),
        ],
        frameon=False,
        loc="upper right",
        fontsize=8,
    )
    ax.set_title(f"Signed Error - MITgcm ({step:g} m contours)", fontsize=10)
    ax.set_xlabel("Longitude offset [deg]")
    ax.set_ylabel("Latitude offset [deg]")
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(direction="out", top=True, right=True)
    fig.tight_layout()

    for suffix in ("pdf", "png"):
        fig.savefig(output_dir / f"tc1_signed_error_contours.{suffix}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def plot_error_metrics(output_dir: Path, rows: list[dict[str, float]]) -> None:
    days = np.array([r["day"] for r in rows])
    l1 = np.array([r["normalized_l1"] for r in rows])
    l2 = np.array([r["normalized_l2"] for r in rows])
    linf = np.array([r["normalized_linf"] for r in rows])

    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    ax.plot(days, l1, color="black", lw=1.5, linestyle="-", label=r"$L_1$")
    ax.plot(days, l2, color="black", lw=1.5, linestyle="--", label=r"$L_2$")
    ax.plot(days, linf, color="black", lw=1.5, linestyle=":", label=r"$L_\infty$")
    ax.set_title("TC1 normalized error growth")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Normalized error")
    ax.set_xlim(float(days[0]), float(days[-1]))
    max_error = float(np.max([np.max(l1), np.max(l2), np.max(linf)]))
    ax.set_ylim(0.0, 1.08 * max_error if max_error > 0.0 else 1.0)
    ax.set_xticks(np.arange(0, max(13, int(days[-1]) + 1), 2))
    ax.legend(frameon=True, loc="upper left", fancybox=False, edgecolor="black")
    fig.savefig(output_dir / "tc1_error_metrics.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_peak_trajectory(output_dir: Path, rows: list[dict[str, float]]) -> None:
    days = np.array([r["day"] for r in rows])
    model_lon = unwrap_longitudes([r["model_lon"] for r in rows])
    exact_lon = unwrap_longitudes([r["exact_lon"] for r in rows])
    model_lat = np.array([r["model_lat"] for r in rows])
    exact_lat = np.array([r["exact_lat"] for r in rows])
    distance = np.array([r["peak_distance_km"] for r in rows])

    fig, axes = plt.subplots(3, 1, figsize=(8.6, 8.6), constrained_layout=True)

    axes[0].plot(days, exact_lon, color="black", lw=1.5, linestyle="--", label="Exact")
    axes[0].plot(days, model_lon, color="black", lw=1.2, linestyle="-", label="MITgcm")
    axes[0].set_title("Bell peak longitude")
    axes[0].set_ylabel("Unwrapped lon [deg]")
    axes[0].legend(frameon=True, loc="upper left", fancybox=False, edgecolor="black")

    axes[1].plot(days, exact_lat, color="black", lw=1.5, linestyle="--", label="Exact")
    axes[1].plot(days, model_lat, color="black", lw=1.2, linestyle="-")
    axes[1].set_title("Bell peak latitude")
    axes[1].set_ylabel("Latitude [deg]")

    axes[2].plot(days, distance, color="black", lw=1.5)
    axes[2].set_title("Great-circle distance between MITgcm and exact peak")
    axes[2].set_xlabel("Time [days]")
    axes[2].set_ylabel("Distance [km]")

    fig.savefig(output_dir / "tc1_peak_trajectory.pdf", bbox_inches="tight")
    plt.close(fig)

def write_tables(output_dir: Path, rows: list[dict[str, float]]) -> None:
    csv_path = output_dir / "tc1_error_table.csv"
    txt_path = output_dir / "tc1_error_table.txt"
    tex_path = output_dir / "tc1_error_table.tex"

    columns = [
        "iteration", "day", "normalized_l1", "normalized_l2", "normalized_linf",
        "mean_error", "mean_abs_error", "rmse", "max_abs_error",
        "min_error", "max_error", "relative_mass_error",
    ]

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in columns})

    with txt_path.open("w", encoding="utf-8") as handle:
        handle.write("Williamson TC1 error table\n")
        handle.write("error = MITgcm - exact\n\n")
        handle.write(f"{'day':>6s} {'L1':>12s} {'L2':>12s} {'Linf':>12s} {'RMSE[m]':>12s} {'mass':>12s}\n")
        for r in rows:
            handle.write(
                f"{r['day']:6.2f} {r['normalized_l1']:12.5e} {r['normalized_l2']:12.5e} "
                f"{r['normalized_linf']:12.5e} {r['rmse']:12.5e} {r['relative_mass_error']:12.5e}\n"
            )

    with tex_path.open("w", encoding="utf-8") as handle:
        handle.write("% Requires: \\usepackage[table]{xcolor}\n")
        handle.write("% Requires: \\usepackage{booktabs}\n")
        handle.write("\\begin{table}[htbp]\n")
        handle.write("\\centering\n")
        handle.write("\\rowcolors{2}{gray!10}{white}\n")
        handle.write("\\begin{tabular}{rrrrrr}\n")
        handle.write("\\toprule\n")
        handle.write("\\rowcolor{gray!25}\n")
        handle.write("day & $L_1$ & $L_2$ & $L_\\infty$ & RMSE [m] & mass err. \\\\ \n")
        handle.write("\\midrule\n")
        for r in rows:
            handle.write(
                f"{r['day']:.2f} & {r['normalized_l1']:.3e} & {r['normalized_l2']:.3e} & "
                f"{r['normalized_linf']:.3e} & {r['rmse']:.3e} & {r['relative_mass_error']:.3e} \\\\ \n"
            )
        handle.write("\\bottomrule\n")
        handle.write("\\end{tabular}\n")
        handle.write("\\caption{Williamson TC1 error metrics. Error is MITgcm minus exact cosine-bell height.}\n")
        handle.write("\\end{table}\n")

def analyze(run_dir: Path) -> None:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    alpha = read_alpha(run_dir)
    output_dir = SCRIPT_DIR.parent / "output" / run_dir.name / "TC1Error"
    output_dir.mkdir(parents=True, exist_ok=True)

    iterations = discover_iterations(run_dir, FIELD)
    if not iterations:
        raise FileNotFoundError(f"No {FIELD} output found in {run_dir}")

    xc = to_2d(read_mds_field(run_dir, "XC"))
    yc = to_2d(read_mds_field(run_dir, "YC"))
    weights = get_weights(run_dir, yc)

    rows: list[dict[str, float]] = []
    peak_rows: list[dict[str, float]] = []
    initial_model: np.ndarray | None = None
    max_change_from_initial = 0.0

    for iteration in iterations:
        model = to_2d(read_mds_field(run_dir, FIELD, iteration))
        exact = exact_height(xc, yc, iteration * delta_t, alpha)
        error = model - exact
        metrics = error_metrics(model, exact, weights)
        day = iteration * delta_t / DAY

        valid = np.isfinite(error) & np.isfinite(weights) & (weights > 0.0)
        area = weights[valid]
        diff = error[valid]
        area_sum = float(np.sum(area))

        if initial_model is None:
            initial_model = model.copy()
        max_change_from_initial = max(max_change_from_initial, float(np.max(np.abs(model - initial_model))))

        model_j, model_i = np.unravel_index(int(np.argmax(model)), model.shape)
        model_lon = float(xc[model_j, model_i])
        model_lat = float(yc[model_j, model_i])
        exact_lon, exact_lat = bell_center_at_time(alpha, iteration * delta_t)

        rows.append({
            "iteration": float(iteration),
            "day": float(day),
            **metrics,
            "mean_error": float(np.sum(area * diff) / area_sum),
            "mean_abs_error": float(np.sum(area * np.abs(diff)) / area_sum),
            "rmse": float(np.sqrt(np.sum(area * diff * diff) / area_sum)),
            "max_abs_error": float(np.max(np.abs(diff))),
            "min_error": float(np.min(diff)),
            "max_error": float(np.max(diff)),
            "max_change_from_initial": max_change_from_initial,
        })

        peak_rows.append({
            "iteration": float(iteration),
            "day": float(day),
            "model_lon": model_lon,
            "model_lat": model_lat,
            "exact_lon": exact_lon,
            "exact_lat": exact_lat,
            "lon_offset_deg": wrap_longitude_delta(model_lon - exact_lon),
            "lat_offset_deg": model_lat - exact_lat,
            "peak_distance_km": great_circle_distance_km(model_lon, model_lat, exact_lon, exact_lat),
        })

    final_iteration = iterations[-1]
    final_day = final_iteration * delta_t / DAY
    final_model = to_2d(read_mds_field(run_dir, FIELD, final_iteration))
    final_exact = exact_height(xc, yc, final_iteration * delta_t, alpha)
    final_error = final_model - final_exact
    center_lon, center_lat = bell_center_at_time(alpha, final_iteration * delta_t)
    half_width = support_window(xc, yc, final_exact, center_lon, center_lat)
    lon, lat, final_model_local = local_window(xc, yc, final_model, center_lon, center_lat, half_width)
    _, _, final_exact_local = local_window(xc, yc, final_exact, center_lon, center_lat, half_width)
    _, _, final_error_local = local_window(xc, yc, final_error, center_lon, center_lat, half_width)

    plot_cosine_bell_overlay(output_dir, lon, lat, final_model_local, final_exact_local, final_day)
    plot_signed_error(output_dir, lon, lat, final_error_local)
    plot_error_metrics(output_dir, rows)
    plot_peak_trajectory(output_dir, peak_rows)
    write_tables(output_dir, rows)

    print(f"TC1 error analysis for: {run_dir}")
    print(f"alpha = {alpha:.8g} rad, deltaT = {delta_t:.8g} s")
    print(f"iterations = {iterations[0]} ... {iterations[-1]}  ({len(iterations)} outputs)")
    print(
        "final errors: "
        f"L1={rows[-1]['normalized_l1']:.6e}, "
        f"L2={rows[-1]['normalized_l2']:.6e}, "
        f"Linf={rows[-1]['normalized_linf']:.6e}, "
        f"mass={rows[-1]['relative_mass_error']:.6e}"
    )
    print(f"outputs written to: {output_dir}")

if __name__ == "__main__":
    analyze(Path(RUN_DIR))

# %%

