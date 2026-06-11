from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from mitgcm_io import (
    cell_weights,
    center_to_shape,
    discover_iterations,
    first_existing_field,
    read_delta_t,
    read_mds_field,
    read_package_alpha,
    to_2d,
)
from shared import diagnosis_output_dir, infer_alpha_label, save_figure_variants, write_manifest

DAY = 86_400.0

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

TC2_R_EARTH = 6.371e6
TC2_OMEGA = 7.292e-5
TC2_G = 9.81
TC2_U0 = 2.0 * math.pi * TC2_R_EARTH / (12.0 * DAY)
TC2_GH0 = 2.94e4
TC2_H0 = TC2_GH0 / TC2_G


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


def weighted_error_metrics(model: np.ndarray, exact: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    diff = model - exact
    abs_exact = np.abs(exact)
    abs_diff = np.abs(diff)
    exact_l1 = max(float(np.sum(weights * abs_exact)), 1.0)
    exact_l2 = max(float(np.sum(weights * exact * exact)), 1.0)
    max_exact = max(float(np.max(abs_exact)), 1.0)
    mass_exact = float(np.sum(weights * exact))
    mass_model = float(np.sum(weights * model))
    return {
        "normalized_l1": float(np.sum(weights * abs_diff) / exact_l1),
        "normalized_l2": float(math.sqrt(np.sum(weights * diff * diff) / exact_l2)),
        "normalized_linf": float(np.max(abs_diff) / max_exact),
        "relative_mass_error": (
            float((mass_model - mass_exact) / abs(mass_exact)) if mass_exact != 0.0 else 0.0
        ),
        "model_max": float(np.max(model)),
        "model_min": float(np.min(model)),
        "exact_max": float(np.max(exact)),
        "exact_min": float(np.min(exact)),
    }


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
        c_model = ax.contour(lon, lat, model, levels=model_levels, colors="black", linestyles="-", linewidths=THIN_LINE)
        ax.clabel(c_model, levels=c_model.levels[1::2], fmt="%g", fontsize=7, inline=True)
    tail_levels = levels_inside(TC1_TAIL_LEVELS, model)
    if tail_levels.size:
        c_tail = ax.contour(lon, lat, model, levels=tail_levels, colors="black", linestyles=":", linewidths=THIN_LINE)
        ax.clabel(c_tail, fmt="%g", fontsize=6, inline=True)
    exact_levels = levels_inside(TC1_REFERENCE_LEVELS, exact)
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
    ax.set_title(f"Signed Error - MITgcm ({step:g} m contours)", fontsize=10)
    ax.set_xlabel("Longitude offset [deg]")
    ax.set_ylabel("Latitude offset [deg]")
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(direction="out", top=True, right=True)
    fig.tight_layout()
    save_figure_variants(fig, output_dir / "tc1_signed_error_contours.pdf", dpi=dpi)
    plt.close(fig)


def plot_tc1_error_metrics(output_dir: Path, rows: list[dict[str, float]]) -> None:
    days = np.array([row["day"] for row in rows])
    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    ax.plot(days, [row["normalized_l1"] for row in rows], color="black", lw=1.5, linestyle="-", label=r"$L_1$")
    ax.plot(days, [row["normalized_l2"] for row in rows], color="black", lw=1.5, linestyle="--", label=r"$L_2$")
    ax.plot(days, [row["normalized_linf"] for row in rows], color="black", lw=1.5, linestyle=":", label=r"$L_\infty$")
    ax.set_title("TC1 normalized error growth")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Normalized error")
    ax.set_xlim(float(days[0]), float(days[-1]))
    max_error = max(float(row["normalized_l1"]) for row in rows)
    max_error = max(max_error, max(float(row["normalized_l2"]) for row in rows))
    max_error = max(max_error, max(float(row["normalized_linf"]) for row in rows))
    ax.set_ylim(0.0, 1.08 * max_error if max_error > 0.0 else 1.0)
    ax.set_xticks(np.arange(0, max(13, int(days[-1]) + 1), 2))
    ax.legend(frameon=True, loc="upper left", fancybox=False, edgecolor="black")
    save_figure_variants(fig, output_dir / "tc1_error_metrics.pdf")
    plt.close(fig)


def write_tc1_error_tables(output_dir: Path, rows: list[dict[str, float]]) -> None:
    columns = [
        "iteration", "day", "normalized_l1", "normalized_l2", "normalized_linf",
        "mean_error", "mean_abs_error", "rmse", "max_abs_error",
        "min_error", "max_error", "relative_mass_error",
    ]
    with (output_dir / "tc1_error_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in columns})
    with (output_dir / "tc1_error_table.txt").open("w", encoding="utf-8") as handle:
        handle.write("Williamson TC1 error table\nerror = MITgcm - exact\n\n")
        handle.write(f"{'day':>6s} {'L1':>12s} {'L2':>12s} {'Linf':>12s} {'RMSE[m]':>12s} {'mass':>12s}\n")
        for row in rows:
            handle.write(
                f"{row['day']:6.2f} {row['normalized_l1']:12.5e} {row['normalized_l2']:12.5e} "
                f"{row['normalized_linf']:12.5e} {row['rmse']:12.5e} {row['relative_mass_error']:12.5e}\n"
            )
    with (output_dir / "tc1_error_table.tex").open("w", encoding="utf-8") as handle:
        handle.write("% Requires: \\usepackage[table]{xcolor}\n% Requires: \\usepackage{booktabs}\n")
        handle.write("\\begin{table}[htbp]\n\\centering\n\\rowcolors{2}{gray!10}{white}\n")
        handle.write("\\begin{tabular}{rrrrrr}\n\\toprule\n\\rowcolor{gray!25}\n")
        handle.write("day & $L_1$ & $L_2$ & $L_\\infty$ & RMSE [m] & mass err. \\\\ \n\\midrule\n")
        for row in rows:
            handle.write(
                f"{row['day']:.2f} & {row['normalized_l1']:.3e} & {row['normalized_l2']:.3e} & "
                f"{row['normalized_linf']:.3e} & {row['rmse']:.3e} & {row['relative_mass_error']:.3e} \\\\ \n"
            )
        handle.write("\\bottomrule\n\\end{tabular}\n")
        handle.write("\\caption{Williamson TC1 error metrics. Error is MITgcm minus exact cosine-bell height.}\n\\end{table}\n")


def analyze_tc1_error(run_dir: Path, *, case_code: str = "TC1", field: str = "S", dpi: int = 400) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    alpha = read_package_alpha(run_dir)
    alpha_label = infer_alpha_label(run_dir, case_code)
    output_dir = diagnosis_output_dir(case_code, "error", alpha_label)
    output_dir.mkdir(parents=True, exist_ok=True)

    iterations = discover_iterations(run_dir, field)
    if not iterations:
        raise FileNotFoundError(f"No {field} output found in {run_dir}")

    xc = to_2d(read_mds_field(run_dir, "XC"))
    yc = to_2d(read_mds_field(run_dir, "YC"))
    weights = cell_weights(run_dir, yc)
    rows: list[dict[str, float]] = []
    initial_model: np.ndarray | None = None
    max_change = 0.0

    for iteration in iterations:
        t_sec = iteration * delta_t
        model = to_2d(read_mds_field(run_dir, field, iteration))
        exact = tc1_exact_height(xc, yc, t_sec, alpha)
        error = model - exact
        valid = np.isfinite(error) & np.isfinite(weights) & (weights > 0.0)
        area = weights[valid]
        diff = error[valid]
        area_sum = float(np.sum(area))
        if initial_model is None:
            initial_model = model.copy()
        max_change = max(max_change, float(np.max(np.abs(model - initial_model))))

        rows.append({
            "iteration": float(iteration),
            "day": float(t_sec / DAY),
            **weighted_error_metrics(model, exact, weights),
            "mean_error": float(np.sum(area * diff) / area_sum),
            "mean_abs_error": float(np.sum(area * np.abs(diff)) / area_sum),
            "rmse": float(np.sqrt(np.sum(area * diff * diff) / area_sum)),
            "max_abs_error": float(np.max(np.abs(diff))),
            "min_error": float(np.min(diff)),
            "max_error": float(np.max(diff)),
            "max_change_from_initial": max_change,
        })

    final_iteration = iterations[-1]
    final_t = final_iteration * delta_t
    final_day = final_t / DAY
    final_model = to_2d(read_mds_field(run_dir, field, final_iteration))
    final_exact = tc1_exact_height(xc, yc, final_t, alpha)
    center_lon, center_lat = tc1_bell_center(alpha, final_t)
    half_width = tc1_support_window(xc, yc, final_exact, center_lon, center_lat)
    lon, lat, model_local = local_window(xc, yc, final_model, center_lon, center_lat, half_width)
    _, _, exact_local = local_window(xc, yc, final_exact, center_lon, center_lat, half_width)
    _, _, error_local = local_window(xc, yc, final_model - final_exact, center_lon, center_lat, half_width)

    plot_tc1_cosine_bell_overlay(output_dir, lon, lat, model_local, exact_local, final_day, dpi)
    plot_tc1_signed_error(output_dir, lon, lat, error_local, dpi)
    plot_tc1_error_metrics(output_dir, rows)
    write_tc1_error_tables(output_dir, rows)
    write_manifest(
        output_dir,
        {
            "alpha": alpha_label,
            "case": case_code,
            "product": "error_analysis",
            "saved": [
                "tc1_cosine_bell_overlay",
                "tc1_signed_error_contours",
                "tc1_error_metrics",
                "tc1_error_table",
            ],
            "source_run": str(run_dir),
        },
    )
    print(f"TC1 error analysis for: {run_dir}")
    print(f"alpha = {alpha:.8g} rad, deltaT = {delta_t:.8g} s")
    print(f"iterations = {iterations[0]} ... {iterations[-1]} ({len(iterations)} outputs)")
    print(
        "final errors: "
        f"L1={rows[-1]['normalized_l1']:.6e}, "
        f"L2={rows[-1]['normalized_l2']:.6e}, "
        f"Linf={rows[-1]['normalized_linf']:.6e}, "
        f"mass={rows[-1]['relative_mass_error']:.6e}"
    )
    print(f"outputs written to: {output_dir}")
    return output_dir


def exact_tc2_eta(xc: np.ndarray, yc: np.ndarray, alpha: float) -> np.ndarray:
    lon = np.deg2rad(xc)
    lat = np.deg2rad(yc)
    mu = -np.cos(lon) * np.cos(lat) * np.sin(alpha) + np.sin(lat) * np.cos(alpha)
    gh = TC2_GH0 - (TC2_R_EARTH * TC2_OMEGA * TC2_U0 + 0.5 * TC2_U0**2) * mu**2
    return gh / TC2_G - TC2_H0


def exact_tc2_u(xc: np.ndarray, yc: np.ndarray, alpha: float) -> np.ndarray:
    lon = np.deg2rad(xc)
    lat = np.deg2rad(yc)
    return TC2_U0 * (np.cos(lat) * np.cos(alpha) + np.cos(lon) * np.sin(lat) * np.sin(alpha))


def exact_tc2_v(xc: np.ndarray, yc: np.ndarray, alpha: float) -> np.ndarray:
    return -TC2_U0 * np.sin(np.deg2rad(xc)) * np.sin(alpha)


def relative_l2(num: np.ndarray, exact: np.ndarray, weights: np.ndarray) -> float:
    denom = float(np.sum(weights * exact * exact))
    return 0.0 if denom <= 0.0 else float(np.sqrt(np.sum(weights * (num - exact) ** 2) / denom))


def relative_linf(num: np.ndarray, exact: np.ndarray) -> float:
    denom = float(np.max(np.abs(exact)))
    return 0.0 if denom <= 0.0 else float(np.max(np.abs(num - exact)) / denom)


def absolute_l2(num: np.ndarray, exact: np.ndarray, weights: np.ndarray) -> float:
    denom = float(np.sum(weights))
    return 0.0 if denom <= 0.0 else float(np.sqrt(np.sum(weights * (num - exact) ** 2) / denom))


def absolute_linf(num: np.ndarray, exact: np.ndarray) -> float:
    return float(np.max(np.abs(num - exact)))


def plot_tc2_error_metrics(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> None:
    days = np.array([row["day"] for row in rows], dtype=np.float64)
    fig, ax = plt.subplots(figsize=(8.6, 5.0), constrained_layout=True)
    ax.semilogy(days, [row["l2_eta"] for row in rows], "-", color="black", lw=1.5, label=r"$L_2(\eta)$")
    ax.semilogy(days, [row["linf_eta"] for row in rows], "--", color="black", lw=1.5, label=r"$L_\infty(\eta)$")
    ax.semilogy(days, [row["l2_u"] for row in rows], "-.", color="black", lw=1.5, label=r"$L_2(u)$")
    ax.semilogy(days, [row["linf_u"] for row in rows], ":", color="black", lw=1.5, label=r"$L_\infty(u)$")
    ax.semilogy(days, [row["l2_v"] for row in rows], "-", color="gray", lw=1.2, label=r"$L_2(v)$ abs")
    ax.semilogy(days, [row["linf_v"] for row in rows], "--", color="gray", lw=1.2, label=r"$L_\infty(v)$ abs")
    ax.set_title("Williamson TC2 steady-state error norms")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Error")
    ax.set_xlim(float(days[0]), float(days[-1]))
    ax.legend(frameon=True, loc="upper left", fancybox=False, edgecolor="black")
    save_figure_variants(fig, output_dir / "TC2_error_norms.pdf", dpi=dpi)
    plt.close(fig)


def write_tc2_error_tables(output_dir: Path, rows: list[dict[str, float]]) -> None:
    columns = [
        "iteration", "day", "l2_eta", "linf_eta", "l2_u", "linf_u",
        "l2_v", "linf_v", "eta_min", "eta_max", "u_min", "u_max", "v_min", "v_max",
    ]
    with (output_dir / "TC2_error_norms.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in columns})
    with (output_dir / "TC2_error_norms.txt").open("w", encoding="utf-8") as handle:
        handle.write("Williamson TC2 steady-state error norms\n\n")
        handle.write(
            f"{'day':>6s} {'L2(eta)':>14s} {'Linf(eta)':>14s} "
            f"{'L2(u)':>14s} {'Linf(u)':>14s} {'L2(v)abs':>14s} {'Linf(v)abs':>14s}\n"
        )
        for row in rows:
            handle.write(
                f"{row['day']:6.2f} {row['l2_eta']:14.6e} {row['linf_eta']:14.6e} "
                f"{row['l2_u']:14.6e} {row['linf_u']:14.6e} "
                f"{row['l2_v']:14.6e} {row['linf_v']:14.6e}\n"
            )
    with (output_dir / "TC2_error_norms.tex").open("w", encoding="utf-8") as handle:
        handle.write("% Requires: \\usepackage[table]{xcolor}\n% Requires: \\usepackage{booktabs}\n")
        handle.write("\\begin{table}[htbp]\n\\centering\n\\rowcolors{2}{gray!10}{white}\n")
        handle.write("\\begin{tabular}{rrrrrr}\n\\toprule\n\\rowcolor{gray!25}\n")
        handle.write("day & $L_2(\\eta)$ & $L_\\infty(\\eta)$ & $L_2(u)$ & $L_\\infty(u)$ & $L_2(v)$ abs \\\\\n\\midrule\n")
        for row in rows:
            handle.write(
                f"{row['day']:.2f} & {row['l2_eta']:.3e} & {row['linf_eta']:.3e} & "
                f"{row['l2_u']:.3e} & {row['linf_u']:.3e} & {row['l2_v']:.3e} \\\\\n"
            )
        handle.write("\\bottomrule\n\\end{tabular}\n")
        handle.write("\\caption{Williamson TC2 steady-state error norms for free surface and velocity components.}\n\\end{table}\n")


def analyze_tc2_error(
    run_dir: Path,
    *,
    case_code: str = "TC2",
    eta_field: str = "Eta",
    u_candidates: tuple[str, ...] = ("U", "UVEL", "UVELMASS"),
    v_candidates: tuple[str, ...] = ("V", "VVEL", "VVELMASS"),
    dpi: int = 400,
) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    alpha = read_package_alpha(run_dir)
    alpha_label = infer_alpha_label(run_dir, case_code)
    output_dir = diagnosis_output_dir(case_code, "error", alpha_label)
    output_dir.mkdir(parents=True, exist_ok=True)

    eta_iterations = discover_iterations(run_dir, eta_field)
    u_name = first_existing_field(run_dir, u_candidates)
    v_name = first_existing_field(run_dir, v_candidates)
    if not eta_iterations or u_name is None or v_name is None:
        raise FileNotFoundError("Missing Eta/U/V output for TC2 error analysis")
    iterations = sorted(
        set(eta_iterations)
        & set(discover_iterations(run_dir, u_name))
        & set(discover_iterations(run_dir, v_name))
    )
    if not iterations:
        raise FileNotFoundError("Eta/U/V iterations do not match")

    xc = to_2d(read_mds_field(run_dir, "XC"))
    yc = to_2d(read_mds_field(run_dir, "YC"))
    weights = cell_weights(run_dir, yc)
    eta_exact = exact_tc2_eta(xc, yc, alpha)
    u_exact = exact_tc2_u(xc, yc, alpha)
    v_exact = exact_tc2_v(xc, yc, alpha)
    rows: list[dict[str, float]] = []

    for iteration in iterations:
        eta_num = to_2d(read_mds_field(run_dir, eta_field, iteration))
        u_num = center_to_shape(read_mds_field(run_dir, u_name, iteration), xc.shape)
        v_num = center_to_shape(read_mds_field(run_dir, v_name, iteration), xc.shape)
        rows.append({
            "iteration": float(iteration),
            "day": float(iteration * delta_t / DAY),
            "l2_eta": relative_l2(eta_num, eta_exact, weights),
            "linf_eta": relative_linf(eta_num, eta_exact),
            "l2_u": relative_l2(u_num, u_exact, weights),
            "linf_u": relative_linf(u_num, u_exact),
            "l2_v": absolute_l2(v_num, v_exact, weights),
            "linf_v": absolute_linf(v_num, v_exact),
            "eta_min": float(np.min(eta_num)),
            "eta_max": float(np.max(eta_num)),
            "u_min": float(np.min(u_num)),
            "u_max": float(np.max(u_num)),
            "v_min": float(np.min(v_num)),
            "v_max": float(np.max(v_num)),
        })

    plot_tc2_error_metrics(output_dir, rows, dpi)
    write_tc2_error_tables(output_dir, rows)
    write_manifest(
        output_dir,
        {
            "alpha": alpha_label,
            "case": case_code,
            "product": "error_analysis",
            "saved": ["TC2_error_norms", "TC2_error_norms.csv", "TC2_error_norms.txt", "TC2_error_norms.tex"],
            "source_run": str(run_dir),
        },
    )
    print(f"TC2 steady-state analysis for: {run_dir}")
    print(f"alpha = {alpha:.8g} rad, deltaT = {delta_t:.8g} s")
    print(f"iterations = {iterations[0]} ... {iterations[-1]} ({len(iterations)} outputs)")
    print(
        "final errors: "
        f"L2(eta)={rows[-1]['l2_eta']:.6e}, "
        f"Linf(eta)={rows[-1]['linf_eta']:.6e}, "
        f"L2(u)={rows[-1]['l2_u']:.6e}, "
        f"Linf(u)={rows[-1]['linf_u']:.6e}, "
        f"L2(v)abs={rows[-1]['l2_v']:.6e}, "
        f"Linf(v)abs={rows[-1]['linf_v']:.6e}"
    )
    print(f"outputs written to: {output_dir}")
    return output_dir
