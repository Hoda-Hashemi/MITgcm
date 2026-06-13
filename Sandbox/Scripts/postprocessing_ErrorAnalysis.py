from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

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

ReferenceFunction = Callable[[np.ndarray, np.ndarray, float, float], np.ndarray]


@dataclass(frozen=True)
class ErrorFieldSpec:
    name: str
    candidates: tuple[str, ...]
    reference: ReferenceFunction | None
    label: str | None = None
    normalize: bool = True
    mass_conservation: bool = False
    reference_kind: str = "exact"


@dataclass(frozen=True)
class ErrorAnalysisSpec:
    case_code: str
    fields: tuple[ErrorFieldSpec, ...]
    title: str
    table_title: str
    plot_stem: str = "error_norms"
    table_stem: str = "error_table"
    alpha_label: str | None = None
    log_y: bool = False
    legacy_single_field_columns: bool = False
    extra_saved: tuple[str, ...] = ()


@dataclass(frozen=True)
class ErrorAnalysisResult:
    output_dir: Path
    rows: list[dict[str, float]]
    final_iteration: int
    final_day: float
    alpha: float
    delta_t: float
    xc: np.ndarray
    yc: np.ndarray
    weights: np.ndarray
    field_names: dict[str, str]
    final_models: dict[str, np.ndarray]
    final_references: dict[str, np.ndarray]
    final_errors: dict[str, np.ndarray]


def weighted_error_norms(
    model: np.ndarray,
    reference: np.ndarray,
    weights: np.ndarray,
    *,
    normalize: bool = True,
) -> dict[str, float]:
    model = np.asarray(model, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    mask = np.isfinite(model) & np.isfinite(reference) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return {"l1": float("nan"), "l2": float("nan"), "linf": float("nan")}

    diff = model[mask] - reference[mask]
    ref = reference[mask]
    weight = weights[mask]
    abs_diff = np.abs(diff)

    if normalize:
        l1_denom = max(float(np.sum(weight * np.abs(ref))), 1.0)
        l2_denom = max(float(np.sum(weight * ref * ref)), 1.0)
        linf_denom = max(float(np.max(np.abs(ref))), 1.0)
        return {
            "l1": float(np.sum(weight * abs_diff) / l1_denom),
            "l2": float(math.sqrt(np.sum(weight * diff * diff) / l2_denom)),
            "linf": float(np.max(abs_diff) / linf_denom),
        }

    area = float(np.sum(weight))
    if area <= 0.0:
        return {"l1": float("nan"), "l2": float("nan"), "linf": float("nan")}
    return {
        "l1": float(np.sum(weight * abs_diff) / area),
        "l2": float(math.sqrt(np.sum(weight * diff * diff) / area)),
        "linf": float(np.max(abs_diff)),
    }


def weighted_error_statistics(
    model: np.ndarray,
    reference: np.ndarray,
    weights: np.ndarray,
) -> dict[str, float]:
    mask = np.isfinite(model) & np.isfinite(reference) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return {
            "mean_error": float("nan"),
            "mean_abs_error": float("nan"),
            "rmse": float("nan"),
            "max_abs_error": float("nan"),
            "min_error": float("nan"),
            "max_error": float("nan"),
        }

    diff = model[mask] - reference[mask]
    weight = weights[mask]
    area = float(np.sum(weight))
    return {
        "mean_error": float(np.sum(weight * diff) / area),
        "mean_abs_error": float(np.sum(weight * np.abs(diff)) / area),
        "rmse": float(math.sqrt(np.sum(weight * diff * diff) / area)),
        "max_abs_error": float(np.max(np.abs(diff))),
        "min_error": float(np.min(diff)),
        "max_error": float(np.max(diff)),
    }


def weighted_integral(values: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.sum(values[mask] * weights[mask]))


def finite_min(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.min(finite)) if finite.size else float("nan")


def finite_max(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.max(finite)) if finite.size else float("nan")


def read_centered_field(run_dir: Path, field: str, iteration: int, shape: tuple[int, int]) -> np.ndarray:
    return center_to_shape(read_mds_field(run_dir, field, iteration), shape)


def common_error_iterations(
    run_dir: Path,
    field_specs: tuple[ErrorFieldSpec, ...],
) -> tuple[dict[str, str], list[int]]:
    field_names: dict[str, str] = {}
    common: set[int] | None = None
    for spec in field_specs:
        field_name = first_existing_field(run_dir, spec.candidates)
        if field_name is None:
            raise FileNotFoundError(f"Missing output field for {spec.name}: {spec.candidates}")
        field_names[spec.name] = field_name
        iterations = set(discover_iterations(run_dir, field_name))
        common = iterations if common is None else common & iterations

    iterations = sorted(common or [])
    if not iterations:
        raise FileNotFoundError("Error-analysis fields have no common output iterations")
    return field_names, iterations


def reference_field(
    spec: ErrorFieldSpec,
    model: np.ndarray,
    initial_fields: dict[str, np.ndarray],
    xc: np.ndarray,
    yc: np.ndarray,
    t_sec: float,
    alpha: float,
) -> np.ndarray:
    if spec.reference is not None:
        return center_to_shape(spec.reference(xc, yc, t_sec, alpha), xc.shape)
    if spec.reference_kind == "initial":
        return initial_fields[spec.name]
    raise ValueError(f"No reference field defined for {spec.name}")


def add_legacy_single_field_columns(
    row: dict[str, float],
    field_name: str,
    *,
    mass_conservation: bool,
) -> None:
    row["normalized_l1"] = row[f"l1_{field_name}"]
    row["normalized_l2"] = row[f"l2_{field_name}"]
    row["normalized_linf"] = row[f"linf_{field_name}"]
    for key in ("mean_error", "mean_abs_error", "rmse", "max_abs_error", "min_error", "max_error"):
        row[key] = row[f"{key}_{field_name}"]
    if mass_conservation:
        row["relative_mass_error"] = row[f"relative_mass_error_{field_name}"]


def csv_columns(spec: ErrorAnalysisSpec) -> list[str]:
    columns = ["iteration", "day"]
    for field in spec.fields:
        name = field.name
        columns.extend(
            [
                f"l1_{name}",
                f"l2_{name}",
                f"linf_{name}",
                f"mean_error_{name}",
                f"mean_abs_error_{name}",
                f"rmse_{name}",
                f"max_abs_error_{name}",
                f"min_error_{name}",
                f"max_error_{name}",
                f"{name}_min",
                f"{name}_max",
                f"reference_min_{name}",
                f"reference_max_{name}",
            ]
        )
        if field.mass_conservation:
            columns.append(f"relative_mass_error_{name}")

    if spec.legacy_single_field_columns and len(spec.fields) == 1:
        columns.extend(
            [
                "normalized_l1",
                "normalized_l2",
                "normalized_linf",
                "mean_error",
                "mean_abs_error",
                "rmse",
                "max_abs_error",
                "min_error",
                "max_error",
            ]
        )
        if spec.fields[0].mass_conservation:
            columns.append("relative_mass_error")
    return columns


def write_error_tables(output_dir: Path, spec: ErrorAnalysisSpec, rows: list[dict[str, float]]) -> None:
    columns = csv_columns(spec)
    with (output_dir / f"{spec.table_stem}.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, float("nan")) for key in columns})

    metric_columns = [
        key
        for key in columns
        if key == "day" or key.startswith(("l1_", "l2_", "linf_", "relative_mass_error"))
    ]
    with (output_dir / f"{spec.table_stem}.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"{spec.table_title}\nerror = model - reference\n\n")
        handle.write(" ".join(f"{key:>18s}" for key in metric_columns) + "\n")
        for row in rows:
            handle.write(" ".join(f"{row.get(key, float('nan')):18.8e}" for key in metric_columns) + "\n")

    with (output_dir / f"{spec.table_stem}.tex").open("w", encoding="utf-8") as handle:
        handle.write("% Requires: \\usepackage[table]{xcolor}\n% Requires: \\usepackage{booktabs}\n")
        handle.write("\\begin{table}[htbp]\n\\centering\n\\rowcolors{2}{gray!10}{white}\n")
        handle.write("\\begin{tabular}{" + "r" * len(metric_columns) + "}\n\\toprule\n")
        handle.write(" & ".join(metric_columns).replace("_", "\\_") + " \\\\\n\\midrule\n")
        for row in rows:
            handle.write(
                " & ".join(f"{row.get(key, float('nan')):.3e}" for key in metric_columns)
                + " \\\\\n"
            )
        handle.write("\\bottomrule\n\\end{tabular}\n")
        handle.write(f"\\caption{{{spec.table_title}.}}\n\\end{{table}}\n")


def plot_error_norms(
    output_dir: Path,
    spec: ErrorAnalysisSpec,
    rows: list[dict[str, float]],
    dpi: int,
) -> None:
    days = np.array([row["day"] for row in rows], dtype=np.float64)
    metric_specs = (
        ("l1", r"$L_1$", "tab:blue"),
        ("l2", r"$L_2$", "tab:orange"),
        ("linf", r"$L_\infty$", "tab:green"),
    )
    fig, axes = plt.subplots(3, 1, figsize=(8.6, 7.2), sharex=True, constrained_layout=True)

    for ax, (metric, metric_label, color) in zip(axes, metric_specs):
        for field in spec.fields:
            label = field.label or field.name
            key = f"{metric}_{field.name}"
            if key not in rows[0]:
                continue
            values = np.array([row[key] for row in rows], dtype=np.float64)
            line_label = f"tracked {label}"
            if spec.log_y:
                values = np.maximum(values, 1.0e-30)
                ax.semilogy(days, values, color=color, lw=1.35, label=line_label)
            else:
                ax.plot(days, values, color=color, lw=1.45, label=line_label)
        ax.set_title(metric_label)
        ax.set_ylabel("error norm")
        if days.size > 1:
            ax.set_xlim(float(days[0]), float(days[-1]))
        ax.grid(True, color="0.72", linewidth=0.6, alpha=0.28)
        ax.legend(frameon=True, loc="best", fancybox=False, edgecolor="black", fontsize=8)
        ax.tick_params(direction="out", top=True, right=True)

    fig.suptitle(spec.title)
    axes[-1].set_xlabel("Time [days]")
    save_figure_variants(fig, output_dir / f"{spec.plot_stem}.pdf", dpi=dpi)
    plt.close(fig)


def analyze_error(
    run_dir: Path,
    spec: ErrorAnalysisSpec,
    *,
    dpi: int = 400,
) -> ErrorAnalysisResult:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    alpha = read_package_alpha(run_dir)
    alpha_label = infer_alpha_label(run_dir, spec.case_code, spec.alpha_label)
    output_dir = diagnosis_output_dir(spec.case_code, "error", alpha_label)
    output_dir.mkdir(parents=True, exist_ok=True)

    field_names, iterations = common_error_iterations(run_dir, spec.fields)
    xc = to_2d(read_mds_field(run_dir, "XC"))
    yc = to_2d(read_mds_field(run_dir, "YC"))
    weights = cell_weights(run_dir, yc)

    initial_fields = {
        field.name: read_centered_field(run_dir, field_names[field.name], iterations[0], xc.shape)
        for field in spec.fields
        if field.reference_kind == "initial"
    }

    rows: list[dict[str, float]] = []
    final_models: dict[str, np.ndarray] = {}
    final_references: dict[str, np.ndarray] = {}
    final_errors: dict[str, np.ndarray] = {}

    for iteration in iterations:
        t_sec = iteration * delta_t
        row: dict[str, float] = {"iteration": float(iteration), "day": float(t_sec / DAY)}

        for field in spec.fields:
            model = read_centered_field(run_dir, field_names[field.name], iteration, xc.shape)
            reference = reference_field(field, model, initial_fields, xc, yc, t_sec, alpha)
            error = model - reference
            norms = weighted_error_norms(model, reference, weights, normalize=field.normalize)
            stats = weighted_error_statistics(model, reference, weights)

            row[f"l1_{field.name}"] = norms["l1"]
            row[f"l2_{field.name}"] = norms["l2"]
            row[f"linf_{field.name}"] = norms["linf"]
            for key, value in stats.items():
                row[f"{key}_{field.name}"] = value
            row[f"{field.name}_min"] = finite_min(model)
            row[f"{field.name}_max"] = finite_max(model)
            row[f"reference_min_{field.name}"] = finite_min(reference)
            row[f"reference_max_{field.name}"] = finite_max(reference)

            if field.mass_conservation:
                reference_mass = weighted_integral(reference, weights)
                model_mass = weighted_integral(model, weights)
                row[f"relative_mass_error_{field.name}"] = (
                    float((model_mass - reference_mass) / abs(reference_mass))
                    if reference_mass != 0.0
                    else 0.0
                )

            if iteration == iterations[-1]:
                final_models[field.name] = model
                final_references[field.name] = reference
                final_errors[field.name] = error

        if spec.legacy_single_field_columns and len(spec.fields) == 1:
            add_legacy_single_field_columns(
                row,
                spec.fields[0].name,
                mass_conservation=spec.fields[0].mass_conservation,
            )
        rows.append(row)

    plot_error_norms(output_dir, spec, rows, dpi)
    write_error_tables(output_dir, spec, rows)
    write_manifest(
        output_dir,
        {
            "alpha": alpha_label,
            "case": spec.case_code,
            "product": "error_analysis",
            "saved": [
                spec.plot_stem,
                f"{spec.table_stem}.csv",
                f"{spec.table_stem}.txt",
                f"{spec.table_stem}.tex",
                *spec.extra_saved,
            ],
            "source_run": str(run_dir),
        },
    )

    final_parts: list[str] = []
    for field in spec.fields:
        final_parts.extend(
            [
                f"L1({field.name})={rows[-1][f'l1_{field.name}']:.6e}",
                f"L2({field.name})={rows[-1][f'l2_{field.name}']:.6e}",
                f"Linf({field.name})={rows[-1][f'linf_{field.name}']:.6e}",
            ]
        )
        if field.mass_conservation:
            final_parts.append(
                f"mass({field.name})={rows[-1][f'relative_mass_error_{field.name}']:.6e}"
            )

    print(f"{spec.case_code} error analysis for: {run_dir}")
    print(f"alpha = {alpha:.8g} rad, deltaT = {delta_t:.8g} s")
    print(f"iterations = {iterations[0]} ... {iterations[-1]} ({len(iterations)} outputs)")
    print("final errors: " + ", ".join(final_parts))
    print(f"outputs written to: {output_dir}")

    return ErrorAnalysisResult(
        output_dir=output_dir,
        rows=rows,
        final_iteration=iterations[-1],
        final_day=float(iterations[-1] * delta_t / DAY),
        alpha=alpha,
        delta_t=delta_t,
        xc=xc,
        yc=yc,
        weights=weights,
        field_names=field_names,
        final_models=final_models,
        final_references=final_references,
        final_errors=final_errors,
    )
