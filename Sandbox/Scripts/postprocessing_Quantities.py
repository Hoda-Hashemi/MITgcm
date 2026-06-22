#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

try:
    import numpy as np
except ImportError:
    raise SystemExit("numpy is required; run with ./.venv/bin/python")

from mitgcm_io import (
    center_to_shape,
    discover_iterations,
    first_existing_field,
    lon_lat,
    read_data_value,
    read_delta_t,
    read_mds_field,
    to_2d,
)
from shared import detect_case_code, diagnosis_output_dir, infer_alpha_label, save_figure_variants, write_manifest

DAY = 86_400.0
SCRIPT_DIR = Path(__file__).resolve().parent
SANDBOX_DIR = SCRIPT_DIR.parent

GRID_COLOR = "#D8D8D8"
GRID_ALPHA = 0.45
GRID_LINEWIDTH = 0.6
COLOR_PRIMARY = "#090708"
COLOR_MEAN = "#FF0B0B"
COLOR_RMS = "#3F5F9D"
COLOR_BAND = "#3F5F9D"
COLOR_ZERO = "#8E8E8E"

ETA_NAMES = ("Eta", "ETAN")
U_NAMES = ("U", "UVEL", "UVELMASS")
V_NAMES = ("V", "VVEL", "VVELMASS")

CSV_COLUMNS = (
    "iteration",
    "day",
    "volume_m3",
    "mass_kg",
    "conserved_scalar_integral",
    "conserved_scalar_rel_change",
    "kinetic_energy_j",
    "potential_energy_j",
    "mechanical_energy_j",
    "ke_area_mean_j_m-2",
    "zeta_mean_s-1",
    "zeta_rms_s-1",
    "zeta_min_s-1",
    "zeta_max_s-1",
    "pv_mean_s-1_m-1",
    "pv_rms_s-1_m-1",
    "pv_min_s-1_m-1",
    "pv_max_s-1_m-1",
)


@dataclass(frozen=True)
class PostprocessingSpec:
    case_code: str | None = None
    eta_candidates: tuple[str, ...] = ETA_NAMES
    u_candidates: tuple[str, ...] = U_NAMES
    v_candidates: tuple[str, ...] = V_NAMES
    kinetic_energy_candidates: tuple[str, ...] = ("momKE", "K.E.", "KE")
    vorticity_candidates: tuple[str, ...] = ("momVort3", "Vort3")
    conserved_candidates: tuple[str, ...] = ()
    conserved_label: str = "conserved scalar"
    conserved_units: str = ""
    compute_kinetic_energy_if_missing: bool = True
    compute_vorticity_if_missing: bool = True
    compute_potential_vorticity_if_missing: bool = True


def default_postprocessing_spec(case_code: str) -> PostprocessingSpec:
    if case_code == "TC1":
        return PostprocessingSpec(
            case_code=case_code,
            conserved_candidates=("S", "SALT"),
            conserved_label="S",
            conserved_units=r"psu m$^2$",
        )
    if case_code == "TC1_prime":
        return PostprocessingSpec(
            case_code=case_code,
            eta_candidates=("ETAN", "Eta"),
            u_candidates=("UVEL", "U", "UVELMASS"),
            v_candidates=("VVEL", "V", "VVELMASS"),
        )
    return PostprocessingSpec(case_code=case_code)


def finite_sum(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return float(np.sum(finite)) if finite.size else float("nan")


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.sum(values[mask] * weights[mask]) / np.sum(weights[mask]))


def weighted_rms(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.sqrt(np.sum(values[mask] * values[mask] * weights[mask]) / np.sum(weights[mask])))


def finite_min(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.min(finite)) if finite.size else float("nan")


def finite_max(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.max(finite)) if finite.size else float("nan")


def load_pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plots; use ./.venv/bin/python or add --no-plots") from exc
    return plt


def grid_area(run_dir: Path, shape: tuple[int, int], xc: np.ndarray, yc: np.ndarray) -> np.ndarray:
    if (run_dir / "RAC.meta").exists():
        return to_2d(read_mds_field(run_dir, "RAC"))

    radius = read_data_value(run_dir, "rSphere", 6_371_000.0)
    lon = np.deg2rad(xc[0, :])
    lat = np.deg2rad(yc[:, 0])
    dlon = float(np.median(np.diff(lon)))
    dlat = float(np.median(np.diff(lat)))
    lat_s = np.clip(lat - 0.5 * dlat, -0.5 * np.pi, 0.5 * np.pi)
    lat_n = np.clip(lat + 0.5 * dlat, -0.5 * np.pi, 0.5 * np.pi)
    area_1d = radius * radius * dlon * (np.sin(lat_n) - np.sin(lat_s))
    return np.broadcast_to(area_1d[:, None], shape).copy()


def resting_thickness(run_dir: Path, shape: tuple[int, int]) -> np.ndarray:
    if (run_dir / "Depth.meta").exists():
        depth = np.abs(to_2d(read_mds_field(run_dir, "Depth")))
        if depth.shape == shape:
            return depth

    hfac = np.ones(shape, dtype=np.float64)
    if (run_dir / "hFacC.meta").exists():
        hfac = center_to_shape(read_mds_field(run_dir, "hFacC"), shape)
    return read_data_value(run_dir, "delR", 1.0) * hfac


def lat_derivative(values: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    edge_order = 2 if lat_rad.size > 2 else 1
    return np.gradient(values, lat_rad, axis=0, edge_order=edge_order)


def relative_vorticity(
    u: np.ndarray,
    v: np.ndarray,
    xc: np.ndarray,
    yc: np.ndarray,
    radius: float,
) -> np.ndarray:
    u = center_to_shape(u, xc.shape)
    v = center_to_shape(v, xc.shape)
    lon_rad = np.deg2rad(xc[0, :])
    lat_rad = np.deg2rad(yc[:, 0])
    dlon = float(np.median(np.diff(lon_rad)))

    theta = np.deg2rad(yc)
    cos_theta = np.cos(theta)
    dv_dlambda = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2.0 * dlon)
    ducos_dtheta = lat_derivative(u * cos_theta, lat_rad)
    denom = radius * np.where(np.abs(cos_theta) > 1.0e-12, cos_theta, np.nan)
    return (dv_dlambda - ducos_dtheta) / denom


def potential_vorticity(zeta: np.ndarray, h: np.ndarray, yc: np.ndarray, omega: float) -> np.ndarray:
    coriolis = 2.0 * omega * np.sin(np.deg2rad(yc))
    return np.where(h > 0.0, (zeta + coriolis) / h, np.nan)


def common_iterations(
    run_dir: Path,
    spec: PostprocessingSpec,
) -> tuple[str, str | None, str | None, str | None, str | None, str | None, list[int]]:
    eta_name = first_existing_field(run_dir, spec.eta_candidates)
    if eta_name is None:
        raise FileNotFoundError("need eta output field")

    kinetic_name = first_existing_field(run_dir, spec.kinetic_energy_candidates)
    vorticity_name = first_existing_field(run_dir, spec.vorticity_candidates)
    need_velocity = (
        (kinetic_name is None and spec.compute_kinetic_energy_if_missing)
        or (vorticity_name is None and spec.compute_vorticity_if_missing)
    )
    u_name = first_existing_field(run_dir, spec.u_candidates) if need_velocity else None
    v_name = first_existing_field(run_dir, spec.v_candidates) if need_velocity else None
    if need_velocity and (u_name is None or v_name is None):
        raise FileNotFoundError("need U and V output fields to compute missing KE/vorticity")

    common = set(discover_iterations(run_dir, eta_name))
    if u_name is not None:
        common &= set(discover_iterations(run_dir, u_name))
    if v_name is not None:
        common &= set(discover_iterations(run_dir, v_name))
    conserved_name = first_existing_field(run_dir, spec.conserved_candidates) if spec.conserved_candidates else None
    if conserved_name is not None:
        common &= set(discover_iterations(run_dir, conserved_name))
    if kinetic_name is not None:
        common &= set(discover_iterations(run_dir, kinetic_name))
    if vorticity_name is not None:
        common &= set(discover_iterations(run_dir, vorticity_name))

    iterations = sorted(common)
    if not iterations:
        raise FileNotFoundError("postprocessing fields have no common output iterations")
    return eta_name, u_name, v_name, conserved_name, kinetic_name, vorticity_name, iterations


def row_for_iteration(
    run_dir: Path,
    iteration: int,
    fields: tuple[str, str | None, str | None, str | None, str | None, str | None],
    xc: np.ndarray,
    yc: np.ndarray,
    area: np.ndarray,
    h0: np.ndarray,
    rho0: float,
    gravity: float,
    omega: float,
    radius: float,
    delta_t: float,
    compute_kinetic_energy_if_missing: bool,
    compute_vorticity_if_missing: bool,
    compute_potential_vorticity_if_missing: bool,
) -> dict[str, float]:
    eta_name, u_name, v_name, conserved_name, kinetic_name, vorticity_name = fields
    eta = center_to_shape(read_mds_field(run_dir, eta_name, iteration), area.shape)

    h = h0 + eta
    potential_energy = 0.5 * rho0 * gravity * finite_sum(eta * eta * area)

    need_velocity = (
        (kinetic_name is None and compute_kinetic_energy_if_missing)
        or (vorticity_name is None and compute_vorticity_if_missing)
    )
    u = v = None
    if need_velocity:
        if u_name is None or v_name is None:
            raise FileNotFoundError("need U and V output fields to compute missing KE/vorticity")
        u = center_to_shape(read_mds_field(run_dir, u_name, iteration), area.shape)
        v = center_to_shape(read_mds_field(run_dir, v_name, iteration), area.shape)

    if kinetic_name is not None:
        ke_specific = center_to_shape(read_mds_field(run_dir, kinetic_name, iteration), area.shape)
        ke_column = rho0 * h * ke_specific
        kinetic = finite_sum(ke_column * area)
        ke_area_mean = weighted_mean(ke_column, area)
    elif compute_kinetic_energy_if_missing and u is not None and v is not None:
        ke_column = 0.5 * rho0 * h * (u * u + v * v)
        kinetic = finite_sum(ke_column * area)
        ke_area_mean = weighted_mean(ke_column, area)
    else:
        kinetic = float("nan")
        ke_area_mean = float("nan")

    if vorticity_name is not None:
        zeta = center_to_shape(read_mds_field(run_dir, vorticity_name, iteration), area.shape)
    elif compute_vorticity_if_missing and u is not None and v is not None:
        zeta = relative_vorticity(u, v, xc, yc, radius)
    else:
        zeta = np.full(area.shape, np.nan, dtype=np.float64)

    if compute_potential_vorticity_if_missing:
        pv = potential_vorticity(zeta, h, yc, omega)
    else:
        pv = np.full(area.shape, np.nan, dtype=np.float64)

    conserved_integral = float("nan")
    if conserved_name is not None:
        conserved = center_to_shape(read_mds_field(run_dir, conserved_name, iteration), area.shape)
        conserved_integral = finite_sum(conserved * area)

    volume = finite_sum(h * area)
    return {
        "iteration": int(iteration),
        "day": float(iteration * delta_t / DAY),
        "volume_m3": volume,
        "mass_kg": rho0 * volume,
        "conserved_scalar_integral": conserved_integral,
        "conserved_scalar_rel_change": float("nan"),
        "kinetic_energy_j": kinetic,
        "potential_energy_j": potential_energy,
        "mechanical_energy_j": kinetic + potential_energy,
        "ke_area_mean_j_m-2": ke_area_mean,
        "zeta_mean_s-1": weighted_mean(zeta, area),
        "zeta_rms_s-1": weighted_rms(zeta, area),
        "zeta_min_s-1": finite_min(zeta),
        "zeta_max_s-1": finite_max(zeta),
        "pv_mean_s-1_m-1": weighted_mean(pv, area),
        "pv_rms_s-1_m-1": weighted_rms(pv, area),
        "pv_min_s-1_m-1": finite_min(pv),
        "pv_max_s-1_m-1": finite_max(pv),
    }


def fill_relative_change(rows: list[dict[str, float]], key: str, out_key: str) -> None:
    reference = next(
        (row[key] for row in rows if np.isfinite(row[key]) and abs(row[key]) > 0.0),
        float("nan"),
    )
    for row in rows:
        value = row[key]
        row[out_key] = (
            float((value - reference) / abs(reference))
            if np.isfinite(value) and np.isfinite(reference)
            else float("nan")
        )


def write_csv_table(output_dir: Path, rows: list[dict[str, float]]) -> Path:
    path = output_dir / "postprocessing_values.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in CSV_COLUMNS})
    return path


def style_timeseries_axis(ax) -> None:
    ax.grid(True, color=GRID_COLOR, linewidth=GRID_LINEWIDTH, alpha=GRID_ALPHA)
    ax.set_axisbelow(True)
    ax.tick_params(direction="out", top=True, right=True)


def legend_inside(ax) -> None:
    ax.legend(frameon=False, fontsize=8, loc="best")


def plot_mass(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)
    mass = np.array([row["mass_kg"] for row in rows], dtype=np.float64)
    reference = next((value for value in mass if np.isfinite(value) and abs(value) > 0.0), float("nan"))
    mass_delta = mass - reference
    mass_relative = np.where(np.isfinite(reference), mass_delta / abs(reference), np.nan)

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.8), sharex=True, constrained_layout=True)
    axes[0].plot(day, mass_relative, color=COLOR_PRIMARY, lw=1.4, label="mass drift")
    axes[0].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[0].set_title("Mass drift")
    axes[0].set_ylabel(r"$\Delta M/M_0$")
    axes[1].plot(day, mass_delta, color=COLOR_RMS, lw=1.25, label="mass anomaly")
    axes[1].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[1].set_title("Mass anomaly")
    axes[1].set_ylabel(r"$\Delta M$ [kg]")
    axes[1].set_xlabel("Time [days]")
    for ax in axes:
        style_timeseries_axis(ax)
        legend_inside(ax)
    paths = save_figure_variants(fig, output_dir / "postprocessing_mass.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_conserved_integral(
    output_dir: Path,
    rows: list[dict[str, float]],
    label: str,
    units: str,
    dpi: int,
) -> list[Path]:
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)
    integral = np.array([row["conserved_scalar_integral"] for row in rows], dtype=np.float64)
    relative = np.array([row["conserved_scalar_rel_change"] for row in rows], dtype=np.float64)
    if not np.any(np.isfinite(integral)):
        return []

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 6.0), sharex=True, constrained_layout=True)
    unit_text = f" [{units}]" if units else ""
    axes[0].plot(day, integral, color=COLOR_PRIMARY, lw=1.4, label=f"{label} integral")
    axes[0].set_title(f"{label} integral")
    axes[0].set_ylabel(r"$I$" + unit_text)
    axes[1].plot(day, relative, color=COLOR_MEAN, lw=1.25, label="relative change")
    axes[1].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[1].set_ylabel("Relative change")
    axes[1].set_xlabel("Time [days]")
    for ax in axes:
        style_timeseries_axis(ax)
        legend_inside(ax)
    paths = save_figure_variants(fig, output_dir / "postprocessing_tracer_integral.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_energies(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)
    panels = (
        ("Kinetic energy", "kinetic_energy_j", r"$E_K$ [J]", COLOR_MEAN),
        ("Potential energy", "potential_energy_j", r"$E_P$ [J]", COLOR_RMS),
        ("Mechanical energy", "mechanical_energy_j", r"$E_M$ [J]", COLOR_PRIMARY),
    )

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 7.2), sharex=True, constrained_layout=True)
    for ax, (title, key, ylabel, color) in zip(axes, panels):
        ax.plot(day, [row[key] for row in rows], color=color, lw=1.35, label=ylabel.split()[0])
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        style_timeseries_axis(ax)
        legend_inside(ax)
    axes[-1].set_xlabel("Time [days]")
    paths = save_figure_variants(fig, output_dir / "postprocessing_energies.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_min_mean_rms(
    ax,
    day: np.ndarray,
    mean: np.ndarray,
    rms: np.ndarray,
    minimum: np.ndarray,
    maximum: np.ndarray,
    title: str,
    ylabel: str,
) -> None:
    markevery = max(1, len(day) // 12)
    ax.fill_between(day, minimum, maximum, color=COLOR_BAND, alpha=0.14, label="min-max")
    ax.plot(
        day,
        mean,
        color=COLOR_MEAN,
        marker="o",
        markersize=3.0,
        markerfacecolor="white",
        markeredgewidth=0.8,
        markevery=markevery,
        lw=1.3,
        label="mean",
    )
    ax.plot(
        day,
        rms,
        color=COLOR_RMS,
        marker="s",
        markersize=3.0,
        markerfacecolor="white",
        markeredgewidth=0.8,
        markevery=markevery,
        lw=1.2,
        ls="--",
        label="rms",
    )
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    legend_inside(ax)
    style_timeseries_axis(ax)


def plot_vorticity_pv(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)

    z_mean = np.array([row["zeta_mean_s-1"] for row in rows], dtype=np.float64)
    z_min = np.array([row["zeta_min_s-1"] for row in rows], dtype=np.float64)
    z_max = np.array([row["zeta_max_s-1"] for row in rows], dtype=np.float64)
    z_rms = np.array([row["zeta_rms_s-1"] for row in rows], dtype=np.float64)
    q_mean = np.array([row["pv_mean_s-1_m-1"] for row in rows], dtype=np.float64)
    q_min = np.array([row["pv_min_s-1_m-1"] for row in rows], dtype=np.float64)
    q_max = np.array([row["pv_max_s-1_m-1"] for row in rows], dtype=np.float64)
    q_rms = np.array([row["pv_rms_s-1_m-1"] for row in rows], dtype=np.float64)

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 6.4), sharex=True, constrained_layout=True)
    plot_min_mean_rms(
        axes[0],
        day,
        z_mean,
        z_rms,
        z_min,
        z_max,
        "Relative vorticity",
        r"$\zeta$ [s$^{-1}$]",
    )
    plot_min_mean_rms(
        axes[1],
        day,
        q_mean,
        q_rms,
        q_min,
        q_max,
        "Potential vorticity",
        r"$Q$ [s$^{-1}$ m$^{-1}$]",
    )
    axes[-1].set_xlabel("Time [days]")
    paths = save_figure_variants(fig, output_dir / "postprocessing_vorticity_pv.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_timeseries(
    output_dir: Path,
    rows: list[dict[str, float]],
    spec: PostprocessingSpec,
    dpi: int,
) -> list[Path]:
    saved: list[Path] = []
    saved.extend(plot_mass(output_dir, rows, dpi))
    saved.extend(plot_conserved_integral(output_dir, rows, spec.conserved_label, spec.conserved_units, dpi))
    saved.extend(plot_vorticity_pv(output_dir, rows, dpi))
    saved.extend(plot_energies(output_dir, rows, dpi))
    return saved


def analyze_run(
    run_dir: Path,
    *,
    spec: PostprocessingSpec | None = None,
    make_plots: bool = True,
    dpi: int = 220,
) -> Path:
    run_dir = run_dir.expanduser().resolve()
    detected_case = detect_case_code(run_dir)
    spec = spec or default_postprocessing_spec(detected_case)
    case_code = spec.case_code or detected_case
    alpha = infer_alpha_label(run_dir, case_code)
    output_dir = diagnosis_output_dir(case_code, "postprocessing", alpha)
    output_dir.mkdir(parents=True, exist_ok=True)

    eta_name, u_name, v_name, conserved_name, kinetic_name, vorticity_name, iterations = common_iterations(run_dir, spec)
    xc, yc = lon_lat(run_dir)
    area = grid_area(run_dir, xc.shape, xc, yc)
    h0 = resting_thickness(run_dir, area.shape)
    rho0 = read_data_value(run_dir, "rhoConst", 1000.0)
    gravity = read_data_value(run_dir, "gravity", 9.81)
    omega = read_data_value(run_dir, "omega", 7.292e-5)
    radius = read_data_value(run_dir, "rSphere", 6_371_000.0)
    delta_t = read_delta_t(run_dir)

    rows = [
        row_for_iteration(
            run_dir,
            iteration,
            (eta_name, u_name, v_name, conserved_name, kinetic_name, vorticity_name),
            xc,
            yc,
            area,
            h0,
            rho0,
            gravity,
            omega,
            radius,
            delta_t,
            spec.compute_kinetic_energy_if_missing,
            spec.compute_vorticity_if_missing,
            spec.compute_potential_vorticity_if_missing,
        )
        for iteration in iterations
    ]
    fill_relative_change(rows, "conserved_scalar_integral", "conserved_scalar_rel_change")

    table = write_csv_table(output_dir, rows)
    saved = [table.name]
    if make_plots:
        plots = plot_timeseries(output_dir, rows, spec, dpi)
        saved.extend(path.name for path in plots)

    write_manifest(
        output_dir,
        {
            "alpha": alpha,
            "case": case_code,
            "product": "postprocessing",
            "source_run": str(run_dir),
            "fields": {
                "eta": eta_name,
                "u": u_name,
                "v": v_name,
                "conserved_scalar": conserved_name,
                "kinetic_energy": kinetic_name or "computed_from_velocity",
                "vorticity": vorticity_name or "computed_from_velocity",
            },
            "iterations": [int(value) for value in iterations],
            "saved": saved,
        },
    )

    print(
        f"{case_code} alpha={alpha}: {len(iterations)} iterations, "
        f"mass={rows[-1]['mass_kg']:.6e} kg, "
        f"ME={rows[-1]['mechanical_energy_j']:.6e} J"
    )
    if conserved_name is not None:
        print(
            f"  {conserved_name} integral={rows[-1]['conserved_scalar_integral']:.6e}, "
            f"relative_change={rows[-1]['conserved_scalar_rel_change']:.6e}"
        )
    print(f"  wrote {output_dir}")
    return output_dir


def default_run_dirs() -> list[Path]:
    return sorted(SANDBOX_DIR.glob("vortexSphere_Williamson_TC*/run_alpha_*"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Williamson MITgcm mass, KE, vorticity, and PV diagnostics.")
    parser.add_argument("run_dirs", nargs="*", type=Path, help="MITgcm run directories. Defaults to all Williamson runs.")
    parser.add_argument("--no-plots", action="store_true", help="Write CSV only.")
    parser.add_argument("--dpi", type=int, default=220, help="Figure DPI.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dirs = args.run_dirs or default_run_dirs()
    if not run_dirs:
        raise SystemExit("no run directories found")

    failures = 0
    for run_dir in run_dirs:
        try:
            analyze_run(run_dir, make_plots=not args.no_plots, dpi=args.dpi)
        except Exception as exc:
            failures += 1
            print(f"failed {run_dir}: {exc}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
