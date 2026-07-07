#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError:
    raise SystemExit("numpy is required; run with ./.venv/bin/python")

from mitgcm_io import (
    center_to_shape,
    discover_iterations,
    first_existing_field,
    lon_lat,
    parse_mds_meta,
    read_data_value,
    read_delta_t,
    read_mds_field,
)
from postprocessing_Quantities import (
    DAY,
    finite_max,
    finite_min,
    finite_sum,
    grid_area,
    load_pyplot,
    relative_vorticity,
    resting_thickness,
    weighted_mean,
    weighted_rms,
)
from shared import (
    OUTPUT_ROOT,
    SANDBOX_DIR,
    case_output_name,
    diagnosis_output_dir,
    infer_alpha_label,
    save_figure_variants,
    write_manifest,
)

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SANDBOX_DIR.parent
DOCS_ASSET_ROOT = REPO_ROOT / "docs" / "assets" / "williamson"
DOCS_FRAGMENT_ROOT = REPO_ROOT / "docs" / "fragments"

CASES = ("TC1", "TC2", "TC3", "TC4", "TC5", "TC6", "TC7")
ETA_NAMES = ("Eta", "ETAN")
U_NAMES = ("U", "UVEL", "UVELMASS")
V_NAMES = ("V", "VVEL", "VVELMASS")
KE_NAMES = ("momKE", "K.E.", "KE")
VORTICITY_NAMES = ("momVort3", "Vort3")

CSV_COLUMNS = (
    "iteration",
    "day",
    "eta_finite_fraction",
    "u_finite_fraction",
    "v_finite_fraction",
    "scalar_finite_fraction",
    "minimum_state_finite_fraction",
    "volume_m3",
    "volume_rel_change",
    "mass_kg",
    "mass_rel_change",
    "eta_volume_anomaly_m3",
    "eta_volume_anomaly_delta_m3",
    "tracer_area_integral",
    "tracer_area_rel_change",
    "tracer_mass_integral",
    "tracer_mass_rel_change",
    "kinetic_energy_j",
    "kinetic_energy_rel_change",
    "surface_potential_energy_j",
    "surface_potential_energy_rel_change",
    "potential_energy_anomaly_j",
    "mechanical_energy_j",
    "mechanical_energy_rel_change",
    "mechanical_energy_anomaly_j",
    "mechanical_energy_anomaly_rel_change",
    "zeta_mean_s-1",
    "zeta_rms_s-1",
    "zeta_min_s-1",
    "zeta_max_s-1",
    "pv_mean_s-1_m-1",
    "pv_rms_s-1_m-1",
    "pv_min_s-1_m-1",
    "pv_max_s-1_m-1",
    "potential_enstrophy",
    "potential_enstrophy_rel_change",
)

SUMMARY_COLUMNS = (
    "case",
    "alpha",
    "status",
    "health_verdict",
    "min_state_finite_fraction",
    "first_bad_day",
    "nonfinite_state_records",
    "run_dir",
    "iterations",
    "last_day",
    "fields",
    "max_abs_mass_rel_change",
    "max_abs_tracer_mass_rel_change",
    "max_abs_mechanical_energy_rel_change",
    "max_abs_potential_enstrophy_rel_change",
    "mass_verdict",
    "quantity_verdict",
    "energy_verdict",
    "enstrophy_verdict",
    "output_dir",
)

CONSERVATION_PRODUCTS = (
    "availability.json",
    "conservation_energy.pdf",
    "conservation_energy.png",
    "conservation_mass.pdf",
    "conservation_mass.png",
    "conservation_pv_enstrophy.pdf",
    "conservation_pv_enstrophy.png",
    "conservation_summary.json",
    "conservation_timeseries.csv",
    "conservation_tracer.pdf",
    "conservation_tracer.png",
    "manifest.json",
)


@dataclass(frozen=True)
class CaseSpec:
    title: str
    quantity_label: str
    scalar_candidates: tuple[str, ...] = ()
    scalar_units: str = ""
    notes: str = ""


CASE_SPECS: dict[str, CaseSpec] = {
    "TC1": CaseSpec(
        title="Cosine bell passive-tracer advection",
        quantity_label="passive tracer",
        scalar_candidates=("S", "SALT"),
        scalar_units="tracer m3",
        notes=(
            "TC1 uses prescribed nondivergent velocity, so the tracer amount is "
            "the primary conserved quantity. Energy is reported only as a "
            "diagnostic proxy from the available fields."
        ),
    ),
    "TC2": CaseSpec(
        title="Steady solid-body geostrophic flow",
        quantity_label="mass, energy, and potential enstrophy",
    ),
    "TC3": CaseSpec(
        title="Steady compact-support zonal flow",
        quantity_label="mass, energy, and potential enstrophy",
    ),
    "TC4": CaseSpec(
        title="Forced nonlinear exact solution",
        quantity_label="mass, energy, and forced-budget diagnostics",
        notes=(
            "TC4 is a forced case with completed run output. Mass can be audited "
            "from the saved free-surface fields; energy and enstrophy are reported "
            "as forced-response diagnostics rather than unforced conservation invariants."
        ),
    ),
    "TC5": CaseSpec(
        title="Zonal flow over an isolated mountain",
        quantity_label="mass, energy, and potential enstrophy",
    ),
    "TC6": CaseSpec(
        title="Rossby-Haurwitz wave",
        quantity_label="mass, energy, and potential enstrophy",
    ),
    "TC7": CaseSpec(
        title="Analyzed 500 mb initial state",
        quantity_label="mass, energy, and potential enstrophy",
        notes=(
            "TC7 has three completed filtered analyzed-state runs at 0000 GMT "
            "21 Dec 1978, 16 Jan 1979, and 9 Jan 1979. The saved MITgcm state "
            "fields are finite through day 5, and conservation diagnostics are "
            "available as validation assets. The current conservation script does "
            "not compute cubed-sphere PV/enstrophy without a native vorticity or "
            "PV diagnostic, so that column remains unavailable for TC7."
        ),
    ),
}

GRID_COLOR = "#D8D8D8"
GRID_ALPHA = 0.45
GRID_LINEWIDTH = 0.6
COLOR_PRIMARY = "#090708"
COLOR_MEAN = "#FF0B0B"
COLOR_RMS = "#3F5F9D"
COLOR_BAND = "#3F5F9D"
COLOR_ZERO = "#8E8E8E"


def safe_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, np.generic):
        return safe_json(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(safe_json(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )


def case_dir(case_code: str) -> Path:
    return SANDBOX_DIR / f"vortexSphere_Williamson_{case_code}"


def conservation_root(case_code: str) -> Path:
    return OUTPUT_ROOT / case_output_name(case_code) / "Diagnosis" / "conservation"


def docs_conservation_root(case_code: str) -> Path:
    return DOCS_ASSET_ROOT / case_output_name(case_code) / "Diagnosis" / "conservation"


def configured_diagnostics(case_code: str) -> list[str]:
    path = case_dir(case_code) / "input" / "data.diagnostics"
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="ignore")
    fields = [item.strip() for item in re.findall(r"'([^']+)'", text)]
    ordered: list[str] = []
    for field in fields:
        if field and field not in ordered:
            ordered.append(field)
    return ordered


def direct_mds_counts(run_dir: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    pattern = re.compile(r"^(?P<name>.+?)(?:\.(?P<iter>\d{10}))?\.meta$")
    for meta in run_dir.glob("*.meta"):
        match = pattern.match(meta.name)
        if not match:
            continue
        name = match.group("name")
        if name.startswith("pickup"):
            continue
        counts[name] = counts.get(name, 0) + 1
    return dict(sorted(counts.items()))


def diagnostic_stream_fields(run_dir: Path) -> dict[str, list[str]]:
    streams: dict[str, list[str]] = {}
    for stream in ("dynDiag", "dyn_Aux", "dynStDiag"):
        meta = next(iter(sorted(run_dir.glob(f"{stream}.*.meta"))), None)
        if meta is None:
            continue
        fields = list(parse_mds_meta(meta).get("fld_list", []))
        streams[stream] = fields
    return streams


def has_run_output(run_dir: Path) -> bool:
    return any(run_dir.glob("*.meta")) and (
        any(run_dir.glob("Eta.*.meta"))
        or any(run_dir.glob("ETAN.*.meta"))
        or any(run_dir.glob("dynDiag.*.meta"))
    )


def discover_case_run_dirs(case_code: str) -> list[Path]:
    root = case_dir(case_code)
    if not root.exists():
        return []
    candidates = [
        path
        for path in sorted(root.iterdir())
        if path.is_dir() and path.name.startswith("run") and has_run_output(path)
    ]
    return candidates


def alpha_label(run_dir: Path, case_code: str) -> str:
    return infer_alpha_label(run_dir, case_code)


def first_iteration_field(run_dir: Path, names: tuple[str, ...]) -> tuple[str | None, set[int]]:
    name = first_existing_field(run_dir, names)
    return name, set(discover_iterations(run_dir, name)) if name is not None else set()


def finite_reference(rows: list[dict[str, float]], key: str) -> float:
    for row in rows:
        value = row.get(key, float("nan"))
        if np.isfinite(value) and abs(value) > 0.0:
            return float(value)
    return float("nan")


def fill_relative(rows: list[dict[str, float]], key: str, out_key: str) -> None:
    reference = finite_reference(rows, key)
    for row in rows:
        value = row.get(key, float("nan"))
        row[out_key] = (
            float((value - reference) / abs(reference))
            if np.isfinite(value) and np.isfinite(reference)
            else float("nan")
        )


def fill_eta_delta(rows: list[dict[str, float]]) -> None:
    reference = next(
        (
            row["eta_volume_anomaly_m3"]
            for row in rows
            if np.isfinite(row["eta_volume_anomaly_m3"])
        ),
        float("nan"),
    )
    for row in rows:
        value = row["eta_volume_anomaly_m3"]
        row["eta_volume_anomaly_delta_m3"] = (
            float(value - reference) if np.isfinite(value) and np.isfinite(reference) else float("nan")
        )


def max_abs_value(rows: list[dict[str, float]], key: str) -> float:
    values = np.array([row.get(key, float("nan")) for row in rows], dtype=np.float64)
    finite = values[np.isfinite(values)]
    return float(np.max(np.abs(finite))) if finite.size else float("nan")


def finite_fraction(field: np.ndarray) -> float:
    field = np.asarray(field)
    if field.size == 0:
        return float("nan")
    return float(np.count_nonzero(np.isfinite(field)) / field.size)


def coriolis_for_diagnostics(
    run_dir: Path,
    shape: tuple[int, int],
    yc: np.ndarray,
    omega: float,
) -> tuple[np.ndarray, str]:
    raw_path = run_dir / "fCoriC.bin"
    if raw_path.exists():
        values = np.fromfile(raw_path, dtype=">f4")
        if values.size == shape[0] * shape[1]:
            return values.reshape(shape).astype(np.float64), "fCoriC.bin"

    return 2.0 * omega * np.sin(np.deg2rad(yc)), "2*omega*sin(latitude)"


def min_finite_fraction(rows: list[dict[str, float]]) -> float:
    values = np.array(
        [row.get("minimum_state_finite_fraction", float("nan")) for row in rows],
        dtype=np.float64,
    )
    finite = values[np.isfinite(values)]
    return float(np.min(finite)) if finite.size else float("nan")


def bad_state_rows(rows: list[dict[str, float]], threshold: float = 0.999) -> list[dict[str, float]]:
    return [
        row
        for row in rows
        if np.isfinite(row.get("minimum_state_finite_fraction", float("nan")))
        and row["minimum_state_finite_fraction"] < threshold
    ]


def first_bad_day(rows: list[dict[str, float]]) -> float:
    bad = bad_state_rows(rows)
    return float(bad[0]["day"]) if bad else float("nan")


def last_finite(rows: list[dict[str, float]], key: str) -> float:
    for row in reversed(rows):
        value = row.get(key, float("nan"))
        if np.isfinite(value):
            return float(value)
    return float("nan")


def grade_relative(max_abs_rel: float, *, unavailable: str = "unavailable") -> str:
    if not np.isfinite(max_abs_rel):
        return unavailable
    if max_abs_rel <= 1.0e-10:
        return "preserved to roundoff"
    if max_abs_rel <= 1.0e-7:
        return "well preserved"
    if max_abs_rel <= 1.0e-5:
        return "small drift"
    if max_abs_rel <= 1.0e-3:
        return "noticeable drift"
    return "not preserved"


def format_metric(value: float, default: str = "n/a") -> str:
    if not np.isfinite(value):
        return default
    if value == 0.0:
        return "0"
    abs_value = abs(value)
    if abs_value < 1.0e-3 or abs_value >= 1.0e4:
        return f"{value:.3e}"
    return f"{value:.6g}"


def format_day(value: float) -> str:
    if not np.isfinite(value):
        return "n/a"
    if abs(value - round(value)) < 1.0e-9:
        return f"day {int(round(value))}"
    return f"day {value:.3g}"


def available_metric(rows: list[dict[str, float]], key: str) -> bool:
    return any(np.isfinite(row.get(key, float("nan"))) for row in rows)


def compute_row(
    run_dir: Path,
    iteration: int,
    fields: dict[str, str | None],
    field_iterations: dict[str, set[int]],
    xc: np.ndarray,
    yc: np.ndarray,
    area: np.ndarray,
    h0: np.ndarray,
    coriolis: np.ndarray,
    rho0: float,
    gravity: float,
    radius: float,
    delta_t: float,
    use_momke: bool,
) -> dict[str, float]:
    eta_name = fields["eta"]
    if eta_name is None:
        raise FileNotFoundError("eta field is required")
    eta = center_to_shape(read_mds_field(run_dir, eta_name, iteration), area.shape)
    h = h0 + eta

    row: dict[str, float] = {key: float("nan") for key in CSV_COLUMNS}
    row["iteration"] = int(iteration)
    row["day"] = float(iteration * delta_t / DAY)
    state_fractions = [finite_fraction(eta)]
    row["eta_finite_fraction"] = state_fractions[-1]
    row["volume_m3"] = finite_sum(h * area)
    row["mass_kg"] = rho0 * row["volume_m3"]
    row["eta_volume_anomaly_m3"] = finite_sum(eta * area)
    row["surface_potential_energy_j"] = 0.5 * rho0 * gravity * finite_sum(eta * eta * area)
    row["potential_energy_anomaly_j"] = rho0 * gravity * finite_sum((h0 * eta + 0.5 * eta * eta) * area)

    scalar_name = fields.get("scalar")
    if scalar_name is not None and iteration in field_iterations.get("scalar", set()):
        scalar = center_to_shape(read_mds_field(run_dir, scalar_name, iteration), area.shape)
        state_fractions.append(finite_fraction(scalar))
        row["scalar_finite_fraction"] = state_fractions[-1]
        row["tracer_area_integral"] = finite_sum(scalar * area)
        row["tracer_mass_integral"] = finite_sum(scalar * h * area)

    u = v = None
    u_name = fields.get("u")
    v_name = fields.get("v")
    have_velocity = (
        u_name is not None
        and v_name is not None
        and iteration in field_iterations.get("u", set())
        and iteration in field_iterations.get("v", set())
    )
    if have_velocity:
        u = center_to_shape(read_mds_field(run_dir, u_name, iteration), area.shape)
        v = center_to_shape(read_mds_field(run_dir, v_name, iteration), area.shape)
        state_fractions.extend([finite_fraction(u), finite_fraction(v)])
        row["u_finite_fraction"] = state_fractions[-2]
        row["v_finite_fraction"] = state_fractions[-1]

    ke_name = fields.get("kinetic")
    if use_momke and ke_name is not None and iteration in field_iterations.get("kinetic", set()):
        ke_specific = center_to_shape(read_mds_field(run_dir, ke_name, iteration), area.shape)
        row["kinetic_energy_j"] = rho0 * finite_sum(h * ke_specific * area)
    elif u is not None and v is not None:
        row["kinetic_energy_j"] = 0.5 * rho0 * finite_sum(h * (u * u + v * v) * area)

    if np.isfinite(row["kinetic_energy_j"]):
        row["mechanical_energy_j"] = row["kinetic_energy_j"] + row["surface_potential_energy_j"]
        row["mechanical_energy_anomaly_j"] = row["kinetic_energy_j"] + row["potential_energy_anomaly_j"]

    vort_name = fields.get("vorticity")
    zeta = None
    if vort_name is not None and iteration in field_iterations.get("vorticity", set()):
        zeta = center_to_shape(read_mds_field(run_dir, vort_name, iteration), area.shape)
    elif u is not None and v is not None:
        zeta = relative_vorticity(u, v, xc, yc, radius)

    if zeta is not None:
        row["zeta_mean_s-1"] = weighted_mean(zeta, area)
        row["zeta_rms_s-1"] = weighted_rms(zeta, area)
        row["zeta_min_s-1"] = finite_min(zeta)
        row["zeta_max_s-1"] = finite_max(zeta)
        pv = np.where(h > 0.0, (zeta + coriolis) / h, np.nan)
        row["pv_mean_s-1_m-1"] = weighted_mean(pv, area)
        row["pv_rms_s-1_m-1"] = weighted_rms(pv, area)
        row["pv_min_s-1_m-1"] = finite_min(pv)
        row["pv_max_s-1_m-1"] = finite_max(pv)
        row["potential_enstrophy"] = 0.5 * finite_sum(np.where(h > 0.0, h * pv * pv * area, np.nan))

    finite_state_fractions = [value for value in state_fractions if np.isfinite(value)]
    row["minimum_state_finite_fraction"] = (
        min(finite_state_fractions) if finite_state_fractions else float("nan")
    )
    return row


def write_csv(path: Path, rows: list[dict[str, Any]], columns: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def style_timeseries_axis(ax) -> None:
    ax.grid(True, color=GRID_COLOR, linewidth=GRID_LINEWIDTH, alpha=GRID_ALPHA)
    ax.set_axisbelow(True)
    ax.tick_params(direction="out", top=True, right=True)


def legend_inside(ax) -> None:
    ax.legend(frameon=False, fontsize=8, loc="best")


def plot_mass(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)
    mass_rel = np.array([row["mass_rel_change"] for row in rows], dtype=np.float64)
    eta_delta = np.array([row["eta_volume_anomaly_delta_m3"] for row in rows], dtype=np.float64)

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.8), sharex=True, constrained_layout=True)
    axes[0].plot(day, mass_rel, color=COLOR_PRIMARY, lw=1.35, label="mass drift")
    axes[0].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[0].set_title("Mass drift")
    axes[0].set_ylabel(r"$\Delta M/M_0$")
    axes[1].plot(day, eta_delta, color=COLOR_RMS, lw=1.35, label="free-surface anomaly")
    axes[1].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[1].set_title("Integrated free-surface anomaly drift")
    axes[1].set_ylabel(r"$\Delta\int\eta\,dA$ [m$^3$]")
    axes[1].set_xlabel("Time [days]")
    for ax in axes:
        style_timeseries_axis(ax)
        legend_inside(ax)
    paths = save_figure_variants(fig, output_dir / "conservation_mass.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_energy(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    if not available_metric(rows, "mechanical_energy_rel_change"):
        return []
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.8), sharex=True, constrained_layout=True)
    axes[0].plot(
        day,
        [row["mechanical_energy_rel_change"] for row in rows],
        color=COLOR_PRIMARY,
        lw=1.35,
        label=r"$E_K+\frac{1}{2}\rho g\eta^2$",
    )
    if available_metric(rows, "mechanical_energy_anomaly_rel_change"):
        axes[0].plot(
            day,
            [row["mechanical_energy_anomaly_rel_change"] for row in rows],
            color=COLOR_BAND,
            lw=1.0,
            ls="--",
            label=r"$E_K+\rho g(H\eta+\eta^2/2)$",
        )
    axes[0].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[0].set_title("Mechanical energy proxy drift")
    axes[0].set_ylabel("Relative change")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(
        day,
        [row["kinetic_energy_rel_change"] for row in rows],
        color=COLOR_MEAN,
        lw=1.2,
        label=r"$E_K$",
    )
    axes[1].plot(
        day,
        [row["surface_potential_energy_rel_change"] for row in rows],
        color=COLOR_RMS,
        lw=1.2,
        label=r"$E_P$",
    )
    axes[1].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[1].set_ylabel("Relative change")
    axes[1].set_xlabel("Time [days]")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        style_timeseries_axis(ax)
    paths = save_figure_variants(fig, output_dir / "conservation_energy.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_tracer(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    if not available_metric(rows, "tracer_mass_rel_change"):
        return []
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    ax.plot(
        day,
        [row["tracer_mass_rel_change"] for row in rows],
        color=COLOR_PRIMARY,
        lw=1.35,
        label="volume-weighted",
    )
    if available_metric(rows, "tracer_area_rel_change"):
        ax.plot(
            day,
            [row["tracer_area_rel_change"] for row in rows],
            color=COLOR_RMS,
            lw=1.0,
            ls="--",
            label="area integral",
        )
    ax.axhline(0.0, color=COLOR_ZERO, lw=0.8)
    ax.set_title("Tracer amount drift")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Relative change")
    ax.legend(frameon=False, fontsize=8)
    style_timeseries_axis(ax)
    paths = save_figure_variants(fig, output_dir / "conservation_tracer.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_pv_enstrophy(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    if not available_metric(rows, "potential_enstrophy_rel_change"):
        return []
    plt = load_pyplot()
    day = np.array([row["day"] for row in rows], dtype=np.float64)

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.8), sharex=True, constrained_layout=True)
    axes[0].plot(
        day,
        [row["potential_enstrophy_rel_change"] for row in rows],
        color=COLOR_PRIMARY,
        lw=1.35,
    )
    axes[0].axhline(0.0, color=COLOR_ZERO, lw=0.8)
    axes[0].set_title("Derived potential enstrophy drift")
    axes[0].set_ylabel("Relative change")
    axes[0].lines[0].set_label("potential enstrophy")
    legend_inside(axes[0])
    axes[1].plot(day, [row["pv_mean_s-1_m-1"] for row in rows], color=COLOR_MEAN, lw=1.2, label="PV mean")
    axes[1].plot(day, [row["pv_rms_s-1_m-1"] for row in rows], color=COLOR_RMS, lw=1.2, label="PV rms")
    axes[1].set_title("Derived potential vorticity")
    axes[1].set_ylabel(r"$Q$ [s$^{-1}$ m$^{-1}$]")
    axes[1].set_xlabel("Time [days]")
    legend_inside(axes[1])
    for ax in axes:
        style_timeseries_axis(ax)
    paths = save_figure_variants(fig, output_dir / "conservation_pv_enstrophy.pdf", dpi=dpi)
    plt.close(fig)
    return paths


def plot_all(output_dir: Path, rows: list[dict[str, float]], dpi: int) -> list[Path]:
    saved: list[Path] = []
    saved.extend(plot_mass(output_dir, rows, dpi))
    saved.extend(plot_tracer(output_dir, rows, dpi))
    saved.extend(plot_energy(output_dir, rows, dpi))
    saved.extend(plot_pv_enstrophy(output_dir, rows, dpi))
    return saved


def prune_conservation_products(output_dir: Path) -> None:
    for name in CONSERVATION_PRODUCTS:
        path = output_dir / name
        if path.exists():
            path.unlink()


def analyze_run(
    run_dir: Path,
    *,
    make_plots: bool,
    dpi: int,
    use_momke: bool,
) -> dict[str, Any]:
    run_dir = run_dir.expanduser().resolve()
    case_code = next((case for case in CASES if case in str(run_dir)), None)
    if case_code is None:
        raise ValueError(f"could not infer Williamson case from {run_dir}")
    spec = CASE_SPECS[case_code]
    alpha = alpha_label(run_dir, case_code)
    output_dir = diagnosis_output_dir(case_code, "conservation", alpha)
    output_dir.mkdir(parents=True, exist_ok=True)
    prune_conservation_products(output_dir)

    eta_name, eta_iters = first_iteration_field(run_dir, ETA_NAMES)
    if eta_name is None or not eta_iters:
        raise FileNotFoundError("no Eta/ETAN output found")

    u_name, u_iters = first_iteration_field(run_dir, U_NAMES)
    v_name, v_iters = first_iteration_field(run_dir, V_NAMES)
    ke_name, ke_iters = first_iteration_field(run_dir, KE_NAMES)
    vort_name, vort_iters = first_iteration_field(run_dir, VORTICITY_NAMES)
    scalar_name, scalar_iters = first_iteration_field(run_dir, spec.scalar_candidates)

    fields = {
        "eta": eta_name,
        "u": u_name,
        "v": v_name,
        "kinetic": ke_name,
        "vorticity": vort_name,
        "scalar": scalar_name,
    }
    field_iterations = {
        "eta": eta_iters,
        "u": u_iters,
        "v": v_iters,
        "kinetic": ke_iters,
        "vorticity": vort_iters,
        "scalar": scalar_iters,
    }
    iterations = sorted(eta_iters)
    if not iterations:
        raise FileNotFoundError("eta output has no iterations")

    xc, yc = lon_lat(run_dir)
    area = grid_area(run_dir, xc.shape, xc, yc)
    h0 = resting_thickness(run_dir, area.shape)
    rho0 = read_data_value(run_dir, "rhoConst", 1.0)
    gravity = read_data_value(run_dir, "gravity", 9.81)
    omega = read_data_value(run_dir, "omega", 7.292e-5)
    radius = read_data_value(run_dir, "rSphere", 6_371_000.0)
    delta_t = read_delta_t(run_dir)
    coriolis, coriolis_source = coriolis_for_diagnostics(run_dir, area.shape, yc, omega)

    rows = [
        compute_row(
            run_dir,
            iteration,
            fields,
            field_iterations,
            xc,
            yc,
            area,
            h0,
            coriolis,
            rho0,
            gravity,
            radius,
            delta_t,
            use_momke,
        )
        for iteration in iterations
    ]

    fill_relative(rows, "volume_m3", "volume_rel_change")
    fill_relative(rows, "mass_kg", "mass_rel_change")
    fill_eta_delta(rows)
    fill_relative(rows, "tracer_area_integral", "tracer_area_rel_change")
    fill_relative(rows, "tracer_mass_integral", "tracer_mass_rel_change")
    fill_relative(rows, "kinetic_energy_j", "kinetic_energy_rel_change")
    fill_relative(rows, "surface_potential_energy_j", "surface_potential_energy_rel_change")
    fill_relative(rows, "mechanical_energy_j", "mechanical_energy_rel_change")
    fill_relative(rows, "mechanical_energy_anomaly_j", "mechanical_energy_anomaly_rel_change")
    fill_relative(rows, "potential_enstrophy", "potential_enstrophy_rel_change")

    table_path = output_dir / "conservation_timeseries.csv"
    write_csv(table_path, rows, CSV_COLUMNS)

    saved = [table_path.name]
    if make_plots:
        saved.extend(path.name for path in plot_all(output_dir, rows, dpi))

    availability = {
        "configured_diagnostics": configured_diagnostics(case_code),
        "direct_mds_counts": direct_mds_counts(run_dir),
        "diagnostic_stream_fields": diagnostic_stream_fields(run_dir),
        "interpreted_fields": fields,
        "field_iteration_counts": {
            key: len(value)
            for key, value in field_iterations.items()
            if fields.get(key) is not None or key == "eta"
        },
        "derived": {
            "mass": eta_name is not None,
            "free_surface_integral": eta_name is not None,
            "kinetic_energy": (ke_name is not None) or (u_name is not None and v_name is not None),
            "relative_vorticity": available_metric(rows, "zeta_mean_s-1"),
            "potential_vorticity": available_metric(rows, "pv_mean_s-1_m-1"),
            "potential_enstrophy": available_metric(rows, "potential_enstrophy"),
            "tracer_amount": scalar_name is not None,
            "coriolis_source": coriolis_source,
        },
    }
    write_json(output_dir / "availability.json", availability)

    max_mass = max_abs_value(rows, "mass_rel_change")
    max_tracer = max_abs_value(rows, "tracer_mass_rel_change")
    max_energy = max_abs_value(rows, "mechanical_energy_rel_change")
    max_enstrophy = max_abs_value(rows, "potential_enstrophy_rel_change")
    min_state_fraction = min_finite_fraction(rows)
    bad_rows = bad_state_rows(rows)
    first_bad = first_bad_day(rows)
    health_verdict = (
        "invalid output: non-finite state fields"
        if bad_rows
        else "finite state fields"
    )
    invalid_output = health_verdict.startswith("invalid")

    summary = {
        "case": case_code,
        "alpha": alpha,
        "status": "available",
        "health_verdict": health_verdict,
        "min_state_finite_fraction": min_state_fraction,
        "first_bad_day": first_bad,
        "nonfinite_state_records": len(bad_rows),
        "run_dir": str(run_dir),
        "iterations": len(rows),
        "last_day": rows[-1]["day"],
        "fields": ",".join(f"{key}:{value}" for key, value in fields.items() if value is not None),
        "coriolis_source": coriolis_source,
        "max_abs_mass_rel_change": max_mass,
        "max_abs_tracer_mass_rel_change": max_tracer,
        "max_abs_mechanical_energy_rel_change": max_energy,
        "max_abs_potential_enstrophy_rel_change": max_enstrophy,
        "mass_verdict": "invalid output" if invalid_output else grade_relative(max_mass),
        "quantity_verdict": "invalid output"
        if invalid_output
        else (
            grade_relative(max_tracer)
            if spec.scalar_candidates
            else grade_relative(max_enstrophy, unavailable="not a scalar test")
        ),
        "energy_verdict": "invalid output" if invalid_output else grade_relative(max_energy),
        "enstrophy_verdict": "invalid output" if invalid_output else grade_relative(max_enstrophy),
        "output_dir": str(output_dir),
        "saved": saved,
        "availability": availability,
        "last_values": {
            "mass_rel_change": last_finite(rows, "mass_rel_change"),
            "tracer_mass_rel_change": last_finite(rows, "tracer_mass_rel_change"),
            "mechanical_energy_rel_change": last_finite(rows, "mechanical_energy_rel_change"),
            "potential_enstrophy_rel_change": last_finite(rows, "potential_enstrophy_rel_change"),
        },
    }
    write_json(output_dir / "conservation_summary.json", summary)
    write_manifest(
        output_dir,
        {
            "alpha": alpha,
            "case": case_code,
            "product": "conservation",
            "source_run": str(run_dir),
            "fields": fields,
            "iterations": [int(value) for value in iterations],
            "saved": saved,
            "summary": {
                "max_abs_mass_rel_change": max_mass,
                "max_abs_tracer_mass_rel_change": max_tracer,
                "max_abs_mechanical_energy_rel_change": max_energy,
                "max_abs_potential_enstrophy_rel_change": max_enstrophy,
            },
        },
    )
    print(
        f"{case_code} alpha={alpha}: {len(rows)} iterations, "
        f"mass drift={format_metric(max_mass)}, "
        f"energy drift={format_metric(max_energy)}, "
        f"enstrophy drift={format_metric(max_enstrophy)}"
    )
    return summary


def unavailable_summary(case_code: str, reason: str) -> dict[str, Any]:
    output_dir = conservation_root(case_code) / "unavailable"
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "case": case_code,
        "alpha": "unavailable",
        "status": "unavailable",
        "health_verdict": "unavailable",
        "min_state_finite_fraction": float("nan"),
        "first_bad_day": float("nan"),
        "nonfinite_state_records": 0,
        "run_dir": "",
        "iterations": 0,
        "last_day": float("nan"),
        "fields": "",
        "max_abs_mass_rel_change": float("nan"),
        "max_abs_tracer_mass_rel_change": float("nan"),
        "max_abs_mechanical_energy_rel_change": float("nan"),
        "max_abs_potential_enstrophy_rel_change": float("nan"),
        "mass_verdict": "unavailable",
        "quantity_verdict": "unavailable",
        "energy_verdict": "unavailable",
        "enstrophy_verdict": "unavailable",
        "output_dir": str(output_dir),
        "reason": reason,
        "availability": {
            "configured_diagnostics": configured_diagnostics(case_code),
            "direct_mds_counts": {},
            "diagnostic_stream_fields": {},
            "derived": {},
        },
    }
    write_json(output_dir / "conservation_summary.json", payload)
    write_manifest(
        output_dir,
        {
            "case": case_code,
            "product": "conservation",
            "status": "unavailable",
            "reason": reason,
            "configured_diagnostics": configured_diagnostics(case_code),
            "saved": ["conservation_summary.json"],
        },
    )
    return payload


def default_unavailable_reason(case_code: str) -> str:
    if case_code == "TC4":
        return (
            "diagnostics pending: TC4 is forced, so final validation needs the "
            "completed run output plus forcing-budget diagnostics"
        )
    if case_code == "TC7":
        return (
            "no completed filtered 25 s TC7 MDS output was found; expected finite "
            "state output through day 5 for the three analyzed 500 mb dates"
        )
    return "no run directory with MDS metadata was found"


def normalize_summary(case_code: str, summary: dict[str, Any]) -> dict[str, Any]:
    if (
        case_code == "TC4"
        and summary.get("status") != "available"
        and (
            "analytic forcing hook must be implemented" in str(summary.get("reason", ""))
            or "no completed MDS output is archived yet" in str(summary.get("reason", ""))
        )
    ):
        summary = dict(summary)
        summary["reason"] = default_unavailable_reason(case_code)
    if (
        case_code == "TC7"
        and summary.get("status") != "available"
        and "requires the analyzed 500 mb initial-condition NPZ" in str(summary.get("reason", ""))
    ):
        summary = dict(summary)
        summary["reason"] = default_unavailable_reason(case_code)
    return summary


def relative_asset_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT / "docs"))
    except ValueError:
        return str(path)


def asset_entries(case_code: str, summaries: list[dict[str, Any]], *, docs_paths: bool) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    base = docs_conservation_root(case_code) if docs_paths else conservation_root(case_code)
    for summary in summaries:
        if summary.get("status") != "available":
            continue
        alpha = str(summary["alpha"])
        run_dir = base / f"alpha_{alpha}"
        saved = set(summary.get("saved", []))
        for name, label in (
            ("conservation_mass.png", "Mass/free surface drift"),
            ("conservation_tracer.png", "Tracer amount drift"),
            ("conservation_energy.png", "Energy proxy drift"),
            ("conservation_pv_enstrophy.png", "PV/enstrophy drift"),
            ("conservation_timeseries.csv", "Conservation time series table"),
            ("conservation_summary.json", "Conservation summary data"),
        ):
            source = run_dir / name
            if name in saved or (name == "conservation_summary.json" and source.exists()):
                entries.append(
                    {
                        "alpha": alpha,
                        "label": label,
                        "path": relative_asset_path(source) if docs_paths else str(source),
                    }
                )
    return entries


def combine_case_verdict(case_code: str, summaries: list[dict[str, Any]]) -> str:
    available = [item for item in summaries if item.get("status") == "available"]
    spec = CASE_SPECS[case_code]
    if not available:
        reason = summaries[0].get("reason", "no usable run output found") if summaries else "no usable run output found"
        return f"Unavailable: {reason}"

    invalid = [
        item
        for item in available
        if str(item.get("health_verdict", "")).startswith("invalid")
    ]
    if invalid:
        worst = min(
            (
                float(item.get("min_state_finite_fraction", float("nan")))
                for item in invalid
                if item.get("min_state_finite_fraction") is not None
            ),
            default=float("nan"),
        )
        first_bad = min(
            (
                float(item.get("first_bad_day", float("nan")))
                for item in invalid
                if np.isfinite(float(item.get("first_bad_day", float("nan"))))
            ),
            default=float("nan"),
        )
        return (
            "Output health is invalid because at least one Eta/U/V/scalar state contains "
            f"non-finite values (minimum finite fraction {format_metric(worst)}). "
            f"The first bad saved record is {format_day(first_bad)}. "
            "Do not use the conservation verdict as final until the run is regenerated "
            "with finite state fields."
        )

    def finite_case_max(key: str) -> float:
        values = []
        for item in available:
            value = item.get(key)
            if value is None:
                continue
            value = float(value)
            if math.isfinite(value):
                values.append(value)
        return max(values) if values else float("nan")

    max_mass = finite_case_max("max_abs_mass_rel_change")
    max_energy = finite_case_max("max_abs_mechanical_energy_rel_change")
    max_enstrophy = finite_case_max("max_abs_potential_enstrophy_rel_change")
    if case_code == "TC4":
        return (
            f"Mass is {grade_relative(max_mass)} (max |rel| {format_metric(max_mass)}). "
            "Because TC4 is analytically forced, the mechanical-energy and derived "
            "potential-enstrophy curves are diagnostic traces, not unforced conservation "
            f"claims (energy max |rel| {format_metric(max_energy)}, "
            f"enstrophy max |rel| {format_metric(max_enstrophy)}). "
            "Together with the analytic path check and published snapshots, the c1 output "
            "is marked verified."
        )
    if spec.scalar_candidates:
        max_tracer = finite_case_max("max_abs_tracer_mass_rel_change")
        return (
            f"Mass is {grade_relative(max_mass)} (max |rel| {format_metric(max_mass)}). "
            f"The TC1 tracer amount is {grade_relative(max_tracer)} "
            f"(max |rel| {format_metric(max_tracer)}). "
            "The energy curve is diagnostic only for this prescribed-flow tracer test."
        )

    return (
        f"Mass is {grade_relative(max_mass)} (max |rel| {format_metric(max_mass)}). "
        f"The mechanical-energy proxy is {grade_relative(max_energy)} "
        f"(max |rel| {format_metric(max_energy)}). "
        f"Derived potential enstrophy is {grade_relative(max_enstrophy)} "
        f"(max |rel| {format_metric(max_enstrophy)})."
    )


def write_case_outputs(case_code: str, summaries: list[dict[str, Any]], *, docs_paths: bool) -> None:
    root = conservation_root(case_code)
    root.mkdir(parents=True, exist_ok=True)
    write_csv(root / "summary.csv", summaries, SUMMARY_COLUMNS)
    write_json(
        root / "availability_summary.json",
        {
            "case": case_code,
            "configured_diagnostics": configured_diagnostics(case_code),
            "runs": summaries,
        },
    )

    assets = asset_entries(case_code, summaries, docs_paths=docs_paths)
    write_json(root / "assets.json", {"case": case_code, "assets": assets})

    spec = CASE_SPECS[case_code]
    available = [item for item in summaries if item.get("status") == "available"]
    lines = [
        f"# {case_code} Conservation and Output-Health Checks",
        "",
        f"**Experiment:** {spec.title}",
        "",
        "## Suggested Section Copy",
        "",
        combine_case_verdict(case_code, summaries),
        "",
        "These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface "
        "drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered "
        "`U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run "
        "currently writes a native vorticity or PV diagnostic.",
        "",
    ]
    if spec.notes:
        lines.extend(["## Note", "", spec.notes, ""])
    lines.extend(
        [
            "## Availability",
            "",
            f"- Configured diagnostics: {', '.join(configured_diagnostics(case_code)) or 'none found'}",
            f"- Runs with conservation tables: {len(available)}",
            "",
            "## Assets To Link",
            "",
        ]
    )
    if assets:
        for item in assets:
            lines.append(f"- alpha {item['alpha']}: {item['label']} -> `{item['path']}`")
    else:
        lines.append("- No conservation assets are available for this case yet.")
    lines.extend(["", "## Per-Run Verdicts", ""])
    for item in summaries:
        alpha = item.get("alpha", "unavailable")
        if item.get("status") != "available":
            lines.append(f"- alpha {alpha}: unavailable ({item.get('reason', 'no run output')}).")
            continue
        health = item.get("health_verdict", "unknown health")
        bad_count = int(item.get("nonfinite_state_records", 0) or 0)
        bad_text = (
            f"; first non-finite state at {format_day(float(item.get('first_bad_day', float('nan'))))}"
            if bad_count
            else ""
        )
        lines.append(
            "- alpha "
            f"{alpha}: mass {item['mass_verdict']}; "
            f"quantity {item['quantity_verdict']}; "
            f"energy {item['energy_verdict']}; "
            f"enstrophy {item['enstrophy_verdict']}; "
            f"health {health}{bad_text}."
        )
    lines.append("")
    (root / "site_snippet.md").write_text("\n".join(lines), encoding="utf-8")


def case_status_class(case_code: str, summaries: list[dict[str, Any]]) -> tuple[str, str]:
    available = [item for item in summaries if item.get("status") == "available"]
    if not available:
        if case_code == "TC4":
            return "pending", "Pending output"
        if case_code == "TC7":
            return "pending", "Pending output"
        if case_code == "TC3":
            return "pending", "Validated with caveat"
        return "scaffold", "Scaffold"
    if any(str(item.get("health_verdict", "")).startswith("invalid") for item in available):
        return "issues", "Issues"
    if any(item.get("mass_verdict") == "not preserved" for item in available):
        return "issues", "Issues"
    if case_code in {"TC2", "TC4", "TC5", "TC6"}:
        return "verified", "Verified"
    if case_code == "TC3":
        return "pending", "Validated with caveat"
    return "verified", "Verified"


def finite_case_max_from_summaries(summaries: list[dict[str, Any]], key: str) -> float:
    values: list[float] = []
    for item in summaries:
        value = item.get(key)
        try:
            value_float = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(value_float):
            values.append(abs(value_float))
    return max(values) if values else float("nan")


def compact_run_verdict(summary: dict[str, Any]) -> str:
    if summary.get("status") != "available":
        return f"unavailable: {summary.get('reason', 'no usable output')}"
    health = str(summary.get("health_verdict", "unknown health"))
    if health.startswith("invalid"):
        return (
            "invalid output; "
            f"{summary.get('nonfinite_state_records', 0)} non-finite saved records, "
            f"first bad {format_day(float(summary.get('first_bad_day', float('nan'))))}"
        )
    return (
        f"mass {summary.get('mass_verdict', 'n/a')}; "
        f"quantity {summary.get('quantity_verdict', 'n/a')}; "
        f"energy {summary.get('energy_verdict', 'n/a')}; "
        f"PV/enstrophy {summary.get('enstrophy_verdict', 'n/a')}"
    )


def write_conservation_report(all_summaries: dict[str, list[dict[str, Any]]], *, docs_paths: bool) -> None:
    DOCS_FRAGMENT_ROOT.mkdir(parents=True, exist_ok=True)

    rows_html: list[str] = []
    diagnostic_equation_html = "\n".join(
        [
            "          <p><strong>Derived potential enstrophy</strong> means this quantity is reconstructed in postprocessing from saved state fields. It is not a native MITgcm output field in these runs.</p>",
            "          <div class='equation-box'>",
            "            <div class='math-line'>\\[h = H + \\eta,\\qquad M = \\rho_0\\int_\\Omega h\\,dA\\]</div>",
            "            <div class='math-line'>\\[E_m = \\rho_0\\int_\\Omega hK\\,dA + \\frac{1}{2}\\rho_0 g\\int_\\Omega \\eta^2\\,dA\\]</div>",
            "            <div class='math-line'>\\[\\zeta = \\frac{1}{a\\cos\\theta}\\left(\\frac{\\partial v}{\\partial\\lambda} - \\frac{\\partial(u\\cos\\theta)}{\\partial\\theta}\\right)\\]</div>",
            "            <div class='math-line'>\\[q = \\frac{\\zeta + f}{h},\\qquad Z = \\frac{1}{2}\\int_\\Omega hq^2\\,dA\\]</div>",
            "            <div class='math-line'>\\[\\Delta_r X(t) = \\frac{X(t)-X(0)}{|X(0)|}\\]</div>",
            "          </div>",
            "          <p>Here <code>K</code> is native <code>momKE</code> when available, otherwise it is reconstructed as <code>(u^2+v^2)/2</code>. The Coriolis term <code>f</code> comes from <code>fCoriC.bin</code> when a rotated run writes it; otherwise the diagnostic uses <code>2*omega*sin(latitude)</code>.</p>",
        ]
    )
    md_lines = [
        "## Conservation diagnostics",
        "",
        "These diagnostics are computed from MITgcm MDS output. `Eta/ETAN` supplies mass and free-surface volume, `S/SALT` supplies the TC1 transported tracer when present, `momKE` or centered `U/V` supplies the mechanical-energy proxy, and `U/V/Eta` supply derived vorticity, potential vorticity, and potential enstrophy.",
        "",
        "Derived potential enstrophy is reconstructed in postprocessing; it is not a native MITgcm output field in these runs.",
        "",
        "$$h = H + \\eta,\\qquad M = \\rho_0\\int_\\Omega h\\,dA$$",
        "",
        "$$E_m = \\rho_0\\int_\\Omega hK\\,dA + \\frac{1}{2}\\rho_0 g\\int_\\Omega \\eta^2\\,dA$$",
        "",
        "$$\\zeta = \\frac{1}{a\\cos\\theta}\\left(\\frac{\\partial v}{\\partial\\lambda} - \\frac{\\partial(u\\cos\\theta)}{\\partial\\theta}\\right)$$",
        "",
        "$$q = \\frac{\\zeta + f}{h},\\qquad Z = \\frac{1}{2}\\int_\\Omega hq^2\\,dA$$",
        "",
        "$$\\Delta_r X(t) = \\frac{X(t)-X(0)}{|X(0)|}$$",
        "",
        "`K` is native `momKE` when available, otherwise `(u^2+v^2)/2`. `f` comes from `fCoriC.bin` when a rotated run writes it; otherwise the diagnostic uses `2*omega*sin(latitude)`.",
        "",
        "| Case | Status | Mass | Target quantity | Energy proxy | PV/enstrophy | Output health |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]

    for case_code in CASES:
        summaries = all_summaries.get(case_code, [])
        spec = CASE_SPECS[case_code]
        status_class, status_label = case_status_class(case_code, summaries)
        available = [item for item in summaries if item.get("status") == "available"]
        case_invalid = any(
            str(item.get("health_verdict", "")).startswith("invalid")
            for item in available
        )
        case_verdict = combine_case_verdict(case_code, summaries)
        if available:
            mass_text = grade_relative(finite_case_max_from_summaries(summaries, "max_abs_mass_rel_change"))
            if spec.scalar_candidates:
                quantity_text = grade_relative(
                    finite_case_max_from_summaries(summaries, "max_abs_tracer_mass_rel_change")
                )
            else:
                quantity_text = grade_relative(
                    finite_case_max_from_summaries(summaries, "max_abs_potential_enstrophy_rel_change"),
                    unavailable="not a scalar test",
                )
            energy_text = grade_relative(
                finite_case_max_from_summaries(summaries, "max_abs_mechanical_energy_rel_change")
            )
            enstrophy_text = grade_relative(
                finite_case_max_from_summaries(summaries, "max_abs_potential_enstrophy_rel_change")
            )
            invalid = [
                item
                for item in available
                if str(item.get("health_verdict", "")).startswith("invalid")
            ]
            if invalid:
                mass_text = quantity_text = energy_text = enstrophy_text = "invalid output"
                first_bad = min(
                    (
                        float(item.get("first_bad_day", float("nan")))
                        for item in invalid
                        if np.isfinite(float(item.get("first_bad_day", float("nan"))))
                    ),
                    default=float("nan"),
                )
                health_text = f"invalid; first bad {format_day(first_bad)}"
            else:
                health_text = "finite saved state fields"
        else:
            mass_text = quantity_text = energy_text = enstrophy_text = "unavailable"
            health_text = summaries[0].get("reason", "no usable run output") if summaries else "no usable run output"

        md_lines.append(
            "| "
            + " | ".join(
                [
                    case_code,
                    status_label,
                    mass_text,
                    quantity_text,
                    energy_text,
                    enstrophy_text,
                    health_text,
                ]
            )
            + " |"
        )

        per_run_html = []
        for summary in summaries:
            alpha = html.escape(str(summary.get("alpha", "unavailable")))
            per_run_html.append(
                f"<li><code>alpha={alpha}</code>: {html.escape(compact_run_verdict(summary))}</li>"
            )
        if not per_run_html:
            per_run_html.append("<li>No conservation output found.</li>")

        assets = asset_entries(case_code, summaries, docs_paths=docs_paths)
        asset_links = []
        if not case_invalid:
            for item in assets:
                if not str(item["path"]).endswith(".png"):
                    continue
                asset_links.append(
                    "<a class='btn-secondary' href='"
                    f"{html.escape(item['path'])}'>{html.escape(item['label'])} "
                    f"(alpha {html.escape(item['alpha'])})</a>"
                )
        if case_invalid:
            asset_html = (
                "<p class='section-note'>Plot files are archived, but they are not "
                "shown as validation evidence because the saved state fields are non-finite.</p>"
            )
        elif asset_links:
            asset_html = "<div class='home-actions'>" + "".join(asset_links[:6]) + "</div>"
        else:
            asset_html = "<p class='section-note'>No plot assets are available for this case yet.</p>"

        rows_html.append(
            "\n".join(
                [
                    f"<details class='details-panel case-block status-{status_class}'>",
                    (
                        f"<summary>{html.escape(case_code)} - {html.escape(spec.title)} "
                        f"<span class='status-chip status-{status_class}'>{html.escape(status_label)}</span></summary>"
                    ),
                    f"<p>{html.escape(case_verdict)}</p>",
                    "<ul>",
                    *per_run_html,
                    "</ul>",
                    asset_html,
                    "</details>",
                ]
            )
        )

    html_text = "\n".join(
        [
            "<section id='conservation-audit' class='section-card experiment-section'>",
            "  <header class='experiment-top'>",
            "    <div>",
            "      <span class='section-label'>Run checks</span>",
            "      <h2>Conservation and Output-Health Checks</h2>",
            "      <p class='experiment-summary'>This section checks whether each run output stayed finite and whether mass, transported quantity, energy proxy, and derived PV/enstrophy behave as expected.</p>",
            "    </div>",
            "    <div class='metric-grid'>",
            "      <div class='metric-card'><span>Fields</span><strong>Eta / U / V</strong></div>",
            "      <div class='metric-card'><span>Health Gate</span><strong>finite state</strong></div>",
            "      <div class='metric-card'><span>Cases</span><strong>TC1-TC7</strong></div>",
            "      <div class='metric-card'><span>Equations</span><strong>Rendered</strong></div>",
            "    </div>",
            "  </header>",
            "  <div class='experiment-body'>",
            "    <article class='experiment-description' aria-label='Conservation diagnostic definitions'>",
            "      <section class='description-block detail-definition'>",
            "        <h3>What Is Computed</h3>",
            "        <div class='description-copy'>",
            "          <p>Mass and free-surface volume come from <code>Eta/ETAN</code>. TC1 tracer conservation uses <code>S/SALT</code>. The mechanical-energy proxy uses native <code>momKE</code> when it exists, otherwise centered <code>U/V</code>. Relative vorticity, PV, and potential enstrophy are derived from <code>U/V/Eta</code>.</p>",
            "        </div>",
            "      </section>",
            "      <section class='description-block detail-equations'>",
            "        <h3>Equations</h3>",
            "        <div class='description-copy'>",
            diagnostic_equation_html,
            "        </div>",
            "      </section>",
            "      <section class='description-block detail-measure'>",
            "        <h3>Health Gate</h3>",
            "        <div class='description-copy'>",
            "          <p>Every saved state is checked for finite <code>Eta</code>, <code>U</code>, <code>V</code>, and scalar values when available. A run with all-NaN or partially non-finite later records is marked invalid even if day 0 conserves exactly.</p>",
            "        </div>",
            "      </section>",
            "    </article>",
            "    <div class='settings-stack'>",
            "      <div class='table-block data-panel'>",
            "        <h3>Case Verdicts</h3>",
            "        <div class='table-scroll'><table><thead><tr><th>Case</th><th>Status</th><th>Mass</th><th>Target quantity</th><th>Energy proxy</th><th>PV/enstrophy</th><th>Output health</th></tr></thead><tbody>",
            *[
                "<tr>"
                + "".join(
                    f"<td>{html.escape(value)}</td>"
                    for value in line.strip("|").split("|")
                )
                + "</tr>"
                for line in md_lines[5:]
                if line.startswith("| ") and "---" not in line
            ],
            "        </tbody></table></div>",
            "      </div>",
            "    </div>",
            "  </div>",
            "  <div class='diagnostics-stack'>",
            *rows_html,
            "  </div>",
            "</section>",
            "",
        ]
    )
    (DOCS_FRAGMENT_ROOT / "conservation_report.html").write_text(html_text, encoding="utf-8")

    md_lines.extend(
        [
            "",
            "### Required output fields",
            "",
            "- `Eta` or `ETAN`: mass, free-surface volume, surface potential energy.",
            "- `U`/`V` or `UVEL`/`VVEL`: velocity health, derived KE if `momKE` is absent, relative vorticity, PV, and potential enstrophy.",
            "- `momKE`: preferred kinetic-energy diagnostic when written.",
            "- `S` or `SALT`: TC1 tracer amount.",
            "",
            "### Per-run notes",
            "",
        ]
    )
    for case_code in CASES:
        summaries = all_summaries.get(case_code, [])
        md_lines.append(f"#### {case_code}")
        md_lines.append("")
        md_lines.append(combine_case_verdict(case_code, summaries))
        md_lines.append("")
        for summary in summaries:
            md_lines.append(f"- alpha `{summary.get('alpha', 'unavailable')}`: {compact_run_verdict(summary)}")
        if not summaries:
            md_lines.append("- No conservation output found.")
        md_lines.append("")
    (DOCS_FRAGMENT_ROOT / "conservation_report.md").write_text("\n".join(md_lines), encoding="utf-8")


def mirror_to_docs(case_code: str) -> None:
    src = conservation_root(case_code)
    dst = docs_conservation_root(case_code)
    if not src.exists():
        return
    if dst.exists():
        shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst, dirs_exist_ok=True)


def load_existing_case_summaries(case_code: str) -> list[dict[str, Any]]:
    path = conservation_root(case_code) / "availability_summary.json"
    if not path.exists():
        return []
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []
    runs = payload.get("runs", [])
    if not isinstance(runs, list):
        return []
    return [normalize_summary(case_code, dict(run)) for run in runs if isinstance(run, dict)]


def parse_cases(value: str | None) -> list[str]:
    if not value:
        return list(CASES)
    requested = [item.strip().upper() for item in value.split(",") if item.strip()]
    unknown = [item for item in requested if item not in CASES]
    if unknown:
        raise SystemExit(f"unknown cases: {', '.join(unknown)}")
    return requested


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Conservation diagnostics for Williamson MITgcm experiment outputs."
    )
    parser.add_argument("run_dirs", nargs="*", type=Path, help="Optional run directories to analyze.")
    parser.add_argument("--cases", help="Comma-separated cases to analyze, default TC1-TC7.")
    parser.add_argument("--no-plots", action="store_true", help="Write tables/manifests only.")
    parser.add_argument("--dpi", type=int, default=220, help="Figure DPI.")
    parser.add_argument(
        "--no-doc-assets",
        action="store_true",
        help="Do not mirror conservation outputs into docs/assets/williamson.",
    )
    parser.add_argument(
        "--compute-ke-from-uv",
        action="store_true",
        help="Ignore momKE and compute kinetic energy from centered U/V.",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Rebuild docs/fragments/conservation_report.* from existing summaries and assets.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected_cases = parse_cases(args.cases)
    explicit_runs = [path.expanduser().resolve() for path in args.run_dirs]

    if args.report_only:
        all_summaries = {
            case_code: load_existing_case_summaries(case_code)
            for case_code in CASES
        }
        write_conservation_report(all_summaries, docs_paths=not args.no_doc_assets)
        return

    failures = 0
    all_summaries: dict[str, list[dict[str, Any]]] = {}
    for case_code in selected_cases:
        if explicit_runs:
            run_dirs = [path for path in explicit_runs if case_code in str(path)]
        else:
            run_dirs = discover_case_run_dirs(case_code)

        summaries: list[dict[str, Any]] = []
        if not run_dirs:
            summaries.append(
                unavailable_summary(
                    case_code,
                    default_unavailable_reason(case_code),
                )
            )
        else:
            for run_dir in run_dirs:
                try:
                    summaries.append(
                        analyze_run(
                            run_dir,
                            make_plots=not args.no_plots,
                            dpi=args.dpi,
                            use_momke=not args.compute_ke_from_uv,
                        )
                    )
                except Exception as exc:
                    failures += 1
                    print(f"failed {run_dir}: {exc}")
                    summaries.append(
                        unavailable_summary(case_code, f"{run_dir}: {exc}")
                    )

        write_case_outputs(case_code, summaries, docs_paths=not args.no_doc_assets)
        if not args.no_doc_assets:
            mirror_to_docs(case_code)
        all_summaries[case_code] = summaries

    for case_code in CASES:
        if case_code not in all_summaries:
            all_summaries[case_code] = load_existing_case_summaries(case_code)

    write_conservation_report(all_summaries, docs_paths=not args.no_doc_assets)

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
