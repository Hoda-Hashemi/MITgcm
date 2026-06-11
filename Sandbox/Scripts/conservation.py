from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mitgcm_io import (
    center_to_shape,
    discover_iterations,
    first_existing_field,
    read_data_value,
    read_delta_t,
    read_mds_field,
    to_2d,
)
from shared import diagnosis_output_dir, infer_alpha_label, save_figure_variants, write_manifest

DAY = 86_400.0


@dataclass(frozen=True)
class ConservationSpec:
    case_code: str
    output_kind: str
    report_title: str
    mass_stem: str
    mass_title: str
    mass_ylabel: str
    energy_stem: str
    energy_title: str
    energy_ylabel: str
    table_stem: str
    table_title: str
    mass_column_label: str
    kinetic_column_label: str
    potential_column_label: str
    energy_column_label: str
    eta_candidates: tuple[str, ...]
    u_candidates: tuple[str, ...]
    v_candidates: tuple[str, ...]
    cleanup_stems: tuple[str, ...] = ()


def scalar_grid_value(run_dir: Path, field: str, default: float) -> float:
    if not (run_dir / f"{field}.meta").exists():
        return default
    values = np.asarray(read_mds_field(run_dir, field), dtype=np.float64)
    if values.size == 0:
        return default
    return float(np.squeeze(values).ravel()[0])


def grid_field(run_dir: Path, field: str, shape: tuple[int, int], default: float = 1.0) -> np.ndarray:
    if not (run_dir / f"{field}.meta").exists():
        return np.full(shape, default, dtype=np.float64)
    values = np.asarray(read_mds_field(run_dir, field), dtype=np.float64)
    if values.shape == shape:
        return values
    return center_to_shape(values, shape)


def finite_sum(values: np.ndarray) -> float:
    return float(np.sum(np.where(np.isfinite(values), values, 0.0)))


def resting_thickness(run_dir: Path, shape: tuple[int, int]) -> np.ndarray:
    if (run_dir / "Depth.meta").exists():
        depth = to_2d(read_mds_field(run_dir, "Depth"))
        if depth.shape == shape:
            return depth
    hfac_c = grid_field(run_dir, "hFacC", shape)
    thickness = scalar_grid_value(run_dir, "DRF", read_data_value(run_dir, "delR", 1.0))
    return thickness * hfac_c


def mass_integral(h: np.ndarray, area: np.ndarray) -> float:
    return finite_sum(h * area)


def kinetic_energy(
    u: np.ndarray,
    v: np.ndarray,
    h: np.ndarray,
    rho0: float,
    area: np.ndarray,
) -> float:
    u_center = center_to_shape(u, area.shape)
    v_center = center_to_shape(v, area.shape)
    return 0.5 * rho0 * finite_sum(h * (u_center * u_center + v_center * v_center) * area)


def potential_energy(eta: np.ndarray, rho0: float, gravity: float, area: np.ndarray) -> float:
    return 0.5 * rho0 * gravity * finite_sum(eta * eta * area)


def remove_output_stems(output_dir: Path, stems: tuple[str, ...]) -> None:
    for stem in stems:
        for suffix in ("pdf", "png"):
            path = output_dir / f"{stem}.{suffix}"
            if path.exists():
                path.unlink()


def plot_mass(output_dir: Path, rows: list[dict[str, float]], dpi: int, spec: ConservationSpec) -> None:
    days = np.array([row["day"] for row in rows], dtype=np.float64)
    mass = np.array([row["mass"] for row in rows], dtype=np.float64)

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.plot(days, mass, color="black", lw=1.5)
    ax.set_title(spec.mass_title)
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(spec.mass_ylabel)
    ax.tick_params(direction="out", top=True, right=True)

    save_figure_variants(fig, output_dir / f"{spec.mass_stem}.pdf", dpi=dpi)
    plt.close(fig)


def plot_mechanical_energy(output_dir: Path, rows: list[dict[str, float]], dpi: int, spec: ConservationSpec) -> None:
    days = np.array([row["day"] for row in rows], dtype=np.float64)
    mechanical = np.array([row["mechanical_energy"] for row in rows], dtype=np.float64)

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.plot(days, mechanical, color="black", lw=1.5)
    ax.set_title(spec.energy_title)
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(spec.energy_ylabel)
    ax.tick_params(direction="out", top=True, right=True)

    save_figure_variants(fig, output_dir / f"{spec.energy_stem}.pdf", dpi=dpi)
    plt.close(fig)


def write_tables(output_dir: Path, rows: list[dict[str, float]], spec: ConservationSpec) -> None:
    columns = [
        "iteration",
        "day",
        "mass",
        "kinetic_energy",
        "potential_energy",
        "mechanical_energy",
    ]
    with (output_dir / f"{spec.table_stem}.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in columns})

    with (output_dir / f"{spec.table_stem}.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"{spec.table_title}\n\n")
        handle.write(
            f"{'day':>6s} {spec.mass_column_label:>14s} {spec.kinetic_column_label:>14s} "
            f"{spec.potential_column_label:>14s} {spec.energy_column_label:>14s}\n"
        )
        for row in rows:
            handle.write(
                f"{row['day']:6.2f} {row['mass']:14.6e} {row['kinetic_energy']:14.6e} "
                f"{row['potential_energy']:14.6e} {row['mechanical_energy']:14.6e}\n"
            )


def analyze_conservation(
    run_dir: Path,
    spec: ConservationSpec,
    *,
    dpi: int = 300,
) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    alpha_label = infer_alpha_label(run_dir, spec.case_code)
    output_dir = diagnosis_output_dir(spec.case_code, spec.output_kind, alpha_label)
    output_dir.mkdir(parents=True, exist_ok=True)

    eta_name = first_existing_field(run_dir, spec.eta_candidates)
    u_name = first_existing_field(run_dir, spec.u_candidates)
    v_name = first_existing_field(run_dir, spec.v_candidates)
    if eta_name is None:
        raise FileNotFoundError(f"Missing Eta output from {spec.eta_candidates}")
    if u_name is None or v_name is None:
        raise FileNotFoundError("Missing U/V output required for conservation diagnostics")

    iterations = set(discover_iterations(run_dir, eta_name))
    iterations &= set(discover_iterations(run_dir, u_name))
    iterations &= set(discover_iterations(run_dir, v_name))
    iterations = sorted(iterations)
    if not iterations:
        raise FileNotFoundError("Eta/U/V iterations do not overlap")

    area = to_2d(read_mds_field(run_dir, "RAC"))
    h0 = resting_thickness(run_dir, area.shape)
    rho0 = read_data_value(run_dir, "rhoConst", 1000.0)
    gravity = read_data_value(run_dir, "gravity", 9.81)

    rows: list[dict[str, float]] = []

    for iteration in iterations:
        eta = to_2d(read_mds_field(run_dir, eta_name, iteration))
        u = np.asarray(read_mds_field(run_dir, u_name, iteration), dtype=np.float64)
        v = np.asarray(read_mds_field(run_dir, v_name, iteration), dtype=np.float64)
        h = h0 + eta

        mass = mass_integral(h, area)
        ke = kinetic_energy(u, v, h, rho0, area)
        pe = potential_energy(eta, rho0, gravity, area)
        me = ke + pe

        rows.append({
            "iteration": float(iteration),
            "day": float(iteration * delta_t / DAY),
            "mass": mass,
            "kinetic_energy": ke,
            "potential_energy": pe,
            "mechanical_energy": me,
        })

    remove_output_stems(output_dir, spec.cleanup_stems)
    plot_mass(output_dir, rows, dpi, spec)
    plot_mechanical_energy(output_dir, rows, dpi, spec)
    write_tables(output_dir, rows, spec)
    write_manifest(
        output_dir,
        {
            "alpha": alpha_label,
            "case": spec.case_code,
            "product": spec.output_kind,
            "saved": [
                spec.mass_stem,
                spec.energy_stem,
                spec.table_stem,
            ],
            "source_run": str(run_dir),
            "fields": {
                "eta": eta_name,
                "u": u_name,
                "v": v_name,
            },
        },
    )

    print(f"{spec.report_title} for: {run_dir}")
    print(f"iterations = {iterations[0]} ... {iterations[-1]} ({len(iterations)} outputs)")
    print(
        "final values: "
        f"mass={rows[-1]['mass']:.6e}, "
        f"mechanical_energy={rows[-1]['mechanical_energy']:.6e}"
    )
    print(f"outputs written to: {output_dir}")
    return output_dir
