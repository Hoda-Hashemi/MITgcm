#!/usr/bin/env python3
#%%
from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from github_pages import build_site
from postprocessing_ErrorAnalysis import (
    DAY,
    ErrorAnalysisResult,
    ErrorAnalysisSpec,
    ErrorFieldSpec,
    analyze_error,
)
from postprocessing_Quantities import PostprocessingSpec, analyze_run
from run_log import write_run_log
from shared import save_figure_variants
from snapshot_plots import SnapshotSpec, run_snapshots

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC1"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC1"

TC1_EARTH_RADIUS = 6_371_000.0
TC1_OMEGA = 2.0 * math.pi / (12.0 * DAY)
TC1_BELL_LON_DEG = 270.0
TC1_BELL_LAT_DEG = 0.0
TC1_BELL_HEIGHT = 1000.0
TC1_BELL_RADIUS = TC1_EARTH_RADIUS / 3.0
TC1_WINDOW_PADDING_DEG = 4.0
TC1_REFERENCE_LEVELS = np.arange(160.0, 801.0, 160.0)
TC1_MODEL_LEVELS = np.arange(80.0, 801.0, 80.0)
TC1_TAIL_LEVELS = np.array([0.01, 0.02, 0.04])
THIN_LINE = 0.65
ZERO_LINE = 1.15

alpha_case = "run_alpha_0"
RUN_DIRS = [
    Path(item).expanduser()
    for item in os.environ.get("TC1_RUN_DIRS", "").split(os.pathsep)
    if item
] or [CASE_DIR / alpha_case]

MAKE_SNAPSHOTS = True
MAKE_ERROR_ANALYSIS = True
MAKE_POSTPROCESSING = True

POSTPROCESSING_SPEC = PostprocessingSpec(
    case_code=CASE_CODE,
    eta_candidates=("Eta", "ETAN"),
    u_candidates=("U", "UVEL"),
    v_candidates=("V", "VVEL"),
    conserved_candidates=("S", "SALT"),
    conserved_label="S",
    conserved_units=r"psu m$^2$",
    compute_kinetic_energy_if_missing=True,
    compute_vorticity_if_missing=True,
    compute_potential_vorticity_if_missing=True,
)

SNAPSHOT_FIELDS = (
    SnapshotSpec("S", "tracer", "passive tracer", "psu", 0.0, 1000.0, center_zero=True),
    SnapshotSpec("ETAN", "etan", "ETAN", "m", center_zero=True),
    SnapshotSpec("PsiVEL", "psi", "PsiVEL", r"m$^3$ s$^{-1}$", center_zero=True),
    SnapshotSpec("PhiVEL", "phi", "PhiVEL", r"m$^2$ s$^{-1}$", center_zero=True),
)

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

def tc1_exact_height(xc: np.ndarray, yc: np.ndarray, t_sec: float, alpha: float) -> np.ndarray:
    x, y, z = unit_vector(xc, yc)
    x0, y0, z0 = rotate_about_axis(x, y, z, alpha, -TC1_OMEGA * t_sec)
    cx, cy, cz = unit_vector(TC1_BELL_LON_DEG, TC1_BELL_LAT_DEG)
    dot = np.clip(x0 * cx + y0 * cy + z0 * cz, -1.0, 1.0)
    radius = TC1_EARTH_RADIUS * np.arccos(dot)
    height = np.zeros_like(radius, dtype=np.float64)
    inside = radius < TC1_BELL_RADIUS
    height[inside] = 0.5 * TC1_BELL_HEIGHT * (
        1.0 + np.cos(math.pi * radius[inside] / TC1_BELL_RADIUS)
    )
    return height

def tc1_bell_center(alpha: float, t_sec: float) -> tuple[float, float]:
    x, y, z = unit_vector(TC1_BELL_LON_DEG, TC1_BELL_LAT_DEG)
    xr, yr, zr = rotate_about_axis(x, y, z, alpha, TC1_OMEGA * t_sec)
    lon = math.degrees(math.atan2(float(yr), float(xr))) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, float(zr)))))
    return lon, lat

def local_window(
    xc: np.ndarray,
    yc: np.ndarray,
    field: np.ndarray,
    center_lon: float,
    center_lat: float,
    half_width: float,
):
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

def tc1_support_window(xc, yc, exact, center_lon: float, center_lat: float) -> float:
    support = exact > 0.0
    if not np.any(support):
        return 30.0
    lon = ((xc[support] - center_lon + 180.0) % 360.0) - 180.0
    lat = yc[support] - center_lat
    return float(max(np.max(np.abs(lon)), np.max(np.abs(lat))) + TC1_WINDOW_PADDING_DEG)

def levels_inside(levels: np.ndarray, field: np.ndarray) -> np.ndarray:
    lo = float(np.nanmin(field))
    hi = float(np.nanmax(field))
    return np.array([value for value in levels if lo <= value <= hi], dtype=np.float64)

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

def plot_tc1_cosine_bell_overlay(output_dir: Path, lon, lat, model, exact, day: float, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    model_levels = levels_inside(TC1_MODEL_LEVELS, model)
    if model_levels.size:
        c_model = ax.contour(
            lon, lat, model, levels=model_levels, colors="black", linestyles="-", linewidths=THIN_LINE
        )
        ax.clabel(c_model, levels=c_model.levels[1::2], fmt="%g", fontsize=7, inline=True)
    tail_levels = levels_inside(TC1_TAIL_LEVELS, model)
    if tail_levels.size:
        c_tail = ax.contour(
            lon, lat, model, levels=tail_levels, colors="black", linestyles=":", linewidths=THIN_LINE
        )
        ax.clabel(c_tail, fmt="%g", fontsize=6, inline=True)
    exact_levels = levels_inside(TC1_REFERENCE_LEVELS, exact)
    if exact_levels.size:
        c_exact = ax.contour(
            lon, lat, exact, levels=exact_levels, colors="black", linestyles="--", linewidths=THIN_LINE
        )
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
    save_figure_variants(fig, output_dir / "tc1_cosine_bell_overlay.pdf", dpi=dpi)
    plt.close(fig)

def plot_tc1_signed_error(output_dir: Path, lon, lat, error, dpi: int) -> None:
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
    ax.set_title(f"Signed Error - MITgcm ({step:g} psu contours)", fontsize=10)
    ax.set_xlabel("Longitude offset [deg]")
    ax.set_ylabel("Latitude offset [deg]")
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(direction="out", top=True, right=True)
    fig.tight_layout()
    save_figure_variants(fig, output_dir / "tc1_signed_error_contours.pdf", dpi=dpi)
    plt.close(fig)

def tc1_error_spec(case_code: str, field_candidates: tuple[str, ...]) -> ErrorAnalysisSpec:
    return ErrorAnalysisSpec(
        case_code=case_code,
        fields=(
            ErrorFieldSpec(
                name="tracer",
                label="tracer",
                candidates=field_candidates,
                reference=tc1_exact_height,
                normalize=True,
                mass_conservation=True,
            ),
        ),
        title=f"{case_code} normalized tracer error growth",
        table_title=f"Williamson {case_code} error table",
        plot_stem="tc1_error_metrics",
        table_stem="tc1_error_table",
        legacy_single_field_columns=True,
        extra_saved=("tc1_cosine_bell_overlay", "tc1_signed_error_contours"),
    )

def write_tc1_final_plots(result: ErrorAnalysisResult, *, dpi: int) -> None:
    final_model = result.final_models["tracer"]
    final_exact = result.final_references["tracer"]
    final_error = result.final_errors["tracer"]
    center_lon, center_lat = tc1_bell_center(result.alpha, result.final_iteration * result.delta_t)
    half_width = tc1_support_window(result.xc, result.yc, final_exact, center_lon, center_lat)
    lon, lat, model_local = local_window(result.xc, result.yc, final_model, center_lon, center_lat, half_width)
    _, _, exact_local = local_window(result.xc, result.yc, final_exact, center_lon, center_lat, half_width)
    _, _, error_local = local_window(result.xc, result.yc, final_error, center_lon, center_lat, half_width)

    plot_tc1_cosine_bell_overlay(result.output_dir, lon, lat, model_local, exact_local, result.final_day, dpi)
    plot_tc1_signed_error(result.output_dir, lon, lat, error_local, dpi)

def analyze_tc1_error(
    run_dir: Path,
    *,
    case_code: str = CASE_CODE,
    field: str = "S",
    dpi: int = 400,
) -> Path:
    field_candidates = ("S", "SALT") if field == "S" else (field,)
    result = analyze_error(run_dir, tc1_error_spec(case_code, field_candidates), dpi=dpi)
    write_tc1_final_plots(result, dpi=dpi)
    return result.output_dir

def main() -> None:
    for run_dir in RUN_DIRS:
        if MAKE_SNAPSHOTS:
            run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        if MAKE_ERROR_ANALYSIS:
            analyze_tc1_error(run_dir)
        if MAKE_POSTPROCESSING:
            analyze_run(run_dir, spec=POSTPROCESSING_SPEC)
    build_site(["testcase1"])
    for run_dir in RUN_DIRS:
        write_run_log(run_dir, case_code=CASE_CODE)

if __name__ == "__main__":
    main()

# %%

