from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from mitgcm_io import (
    center_to_shape,
    day_number,
    discover_iterations,
    finite_min_max,
    first_existing_field,
    grid_extent,
    lon_lat,
    read_delta_t,
    read_mds_field,
    time_text,
    to_2d,
)
from shared import infer_alpha_label, resolve_color_limits, save_figure_variants, snapshot_output_dir, write_manifest

CS_FACE_SIZE = 32
CS_NET_POSITIONS = {
    0: (1, 0),
    1: (1, 1),
    2: (1, 2),
    3: (1, 3),
    4: (0, 1),
    5: (2, 1),
}

plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
})


@dataclass(frozen=True)
class SnapshotSpec:
    field: str
    folder: str
    display: str
    units: str
    vmin: float | None = None
    vmax: float | None = None
    center_zero: bool = False


@dataclass
class SnapshotResult:
    folder: str
    field: str
    display: str
    iterations: list[int]
    saved_iterations: list[int]
    skipped_nonfinite_iterations: list[int]

    @property
    def wrote_any(self) -> bool:
        return bool(self.saved_iterations)

    @property
    def invalid(self) -> bool:
        return bool(self.skipped_nonfinite_iterations)


def plot_snapshot(
    field: np.ndarray,
    xc: np.ndarray,
    yc: np.ndarray,
    title: str,
    units: str,
    out_path: Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    center_zero: bool = False,
    dpi: int = 220,
) -> None:
    field = to_2d(field)
    vmin, vmax = resolve_color_limits(field, vmin, vmax, symmetric=center_zero)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if center_zero else None

    fig, ax = plt.subplots(figsize=(9.0, 5.8), constrained_layout=True)
    image = ax.imshow(
        field,
        origin="lower",
        extent=grid_extent(xc, yc),
        interpolation="bilinear",
        aspect="auto",
        cmap="seismic",
        norm=norm,
        vmin=None if norm else vmin,
        vmax=None if norm else vmax,
    )
    ax.set_title(title, pad=18)
    ax.set_xlabel(r"longitude $\lambda$ [deg]")
    ax.set_ylabel(r"latitude $\theta$ [deg]")
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    cbar = fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.10)
    cbar.set_label(units)
    save_figure_variants(fig, out_path, dpi=dpi, formats=("png",))
    plt.close(fig)


def is_cs_compact(field: np.ndarray) -> bool:
    field = np.squeeze(np.asarray(field))
    return field.shape in {
        (CS_FACE_SIZE * 6, CS_FACE_SIZE),
        (CS_FACE_SIZE, CS_FACE_SIZE * 6),
        (CS_FACE_SIZE, 6, CS_FACE_SIZE),
        (6, CS_FACE_SIZE, CS_FACE_SIZE),
    }


def compact_to_faces(field: np.ndarray) -> np.ndarray:
    arr = np.squeeze(np.asarray(field, dtype=np.float64))
    if arr.shape == (6, CS_FACE_SIZE, CS_FACE_SIZE):
        return arr
    if arr.shape == (CS_FACE_SIZE, 6, CS_FACE_SIZE):
        return np.moveaxis(arr, 1, 0)
    if arr.shape == (CS_FACE_SIZE, CS_FACE_SIZE * 6):
        return np.asarray(
            [arr[:, index * CS_FACE_SIZE : (index + 1) * CS_FACE_SIZE] for index in range(6)]
        )
    if arr.shape == (CS_FACE_SIZE * 6, CS_FACE_SIZE):
        return np.asarray(
            [arr[index * CS_FACE_SIZE : (index + 1) * CS_FACE_SIZE, :] for index in range(6)]
        )
    raise ValueError(f"cannot convert cubed-sphere field with shape {arr.shape} to six faces")


def angular_distance_degrees(left: float, right: float) -> float:
    return abs(((left - right + 180.0) % 360.0) - 180.0)


def ordered_cube_faces(field: np.ndarray, xc: np.ndarray, yc: np.ndarray) -> np.ndarray:
    faces = compact_to_faces(field)
    lon_faces = compact_to_faces(xc)
    lat_faces = compact_to_faces(yc)

    center = CS_FACE_SIZE // 2
    center_lon = [float(lon_faces[index, center, center]) for index in range(6)]
    center_lat = [float(lat_faces[index, center, center]) for index in range(6)]

    north = int(np.argmax(center_lat))
    south = int(np.argmin(center_lat))
    remaining = [index for index in range(6) if index not in {north, south}]

    middle: list[int] = []
    for target_lon in (0.0, 90.0, 180.0, -90.0):
        chosen = min(remaining, key=lambda index: angular_distance_degrees(center_lon[index], target_lon))
        middle.append(chosen)
        remaining.remove(chosen)

    return faces[middle + [north, south]]


def faces_to_net(faces: np.ndarray) -> np.ndarray:
    net = np.full((CS_FACE_SIZE * 3, CS_FACE_SIZE * 4), np.nan, dtype=np.float64)
    for face_index, (row, col) in CS_NET_POSITIONS.items():
        net[
            row * CS_FACE_SIZE : (row + 1) * CS_FACE_SIZE,
            col * CS_FACE_SIZE : (col + 1) * CS_FACE_SIZE,
        ] = faces[face_index]
    return net


def plot_cube_net_snapshot(
    field: np.ndarray,
    xc: np.ndarray,
    yc: np.ndarray,
    title: str,
    units: str,
    out_path: Path,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    center_zero: bool = False,
    dpi: int = 220,
) -> None:
    ordered = ordered_cube_faces(field, xc, yc)
    net = faces_to_net(ordered)
    vmin, vmax = resolve_color_limits(net, vmin, vmax, symmetric=center_zero)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if center_zero else None

    fig, ax = plt.subplots(figsize=(8.8, 6.4), constrained_layout=True)
    cmap = plt.get_cmap("seismic").copy()
    cmap.set_bad("white")
    image = ax.imshow(
        np.ma.masked_invalid(net),
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
        vmin=None if norm else vmin,
        vmax=None if norm else vmax,
    )
    for xline in range(0, CS_FACE_SIZE * 4 + 1, CS_FACE_SIZE):
        ax.axvline(xline - 0.5, color="0.25", linewidth=0.8)
    for yline in range(0, CS_FACE_SIZE * 3 + 1, CS_FACE_SIZE):
        ax.axhline(yline - 0.5, color="0.25", linewidth=0.8)
    labels = {
        "face 5": (CS_FACE_SIZE * 1 + 2, CS_FACE_SIZE * 2 + 25),
        "face 1": (2, CS_FACE_SIZE * 1 + 25),
        "face 2": (CS_FACE_SIZE * 1 + 2, CS_FACE_SIZE * 1 + 25),
        "face 3": (CS_FACE_SIZE * 2 + 2, CS_FACE_SIZE * 1 + 25),
        "face 4": (CS_FACE_SIZE * 3 + 2, CS_FACE_SIZE * 1 + 25),
        "face 6": (CS_FACE_SIZE * 1 + 2, 25),
    }
    for label, (xpos, ypos) in labels.items():
        ax.text(xpos, ypos, label, color="0.2", fontsize=9, alpha=0.75)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, pad=14)
    cbar = fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.10, fraction=0.08)
    cbar.set_label(units)
    save_figure_variants(fig, out_path, dpi=dpi, formats=("png",))
    plt.close(fig)


def clear_generated_snapshot_images(output_root: Path, folder: str) -> None:
    field_dir = output_root / folder
    if not field_dir.exists():
        return
    for pattern in ("*.pdf", "*.png"):
        for path in field_dir.glob(pattern):
            path.unlink()


def save_scalar_series(
    case_code: str,
    run_dir: Path,
    output_root: Path,
    spec: SnapshotSpec,
    delta_t: float,
    *,
    dpi: int = 220,
) -> SnapshotResult:
    iterations = discover_iterations(run_dir, spec.field)
    if not iterations:
        print(f"skip {spec.display}: no {spec.field} files")
        return SnapshotResult(spec.folder, spec.field, spec.display, [], [], [])

    xc, yc = lon_lat(run_dir)
    cube_plot = is_cs_compact(xc) and is_cs_compact(yc)
    case = case_code.lower()
    saved_iterations: list[int] = []
    skipped_nonfinite_iterations: list[int] = []
    for iteration in iterations:
        field = to_2d(read_mds_field(run_dir, spec.field, iteration))
        value_range = finite_min_max(field)
        if value_range is None:
            print(f"skip {spec.display} iteration {iteration}: no finite values")
            skipped_nonfinite_iterations.append(iteration)
            continue
        missing = int(field.size - np.count_nonzero(np.isfinite(field)))
        if missing:
            print(f"warning {spec.display} iteration {iteration}: {missing} non-finite values")
        field_min, field_max = value_range
        title = (
            f"{case_code} {spec.display} | iteration {iteration} | "
            f"time {time_text(iteration, delta_t)} | "
            f"min {field_min:.3e} | max {field_max:.3e}"
        )
        out = output_root / spec.folder / f"{case}_{spec.folder}_day_{day_number(iteration, delta_t):05.2f}_iter_{iteration:010d}.pdf"
        if cube_plot and is_cs_compact(field):
            plot_cube_net_snapshot(
                field,
                xc,
                yc,
                title,
                spec.units,
                out,
                vmin=spec.vmin,
                vmax=spec.vmax,
                center_zero=spec.center_zero,
                dpi=dpi,
            )
        else:
            plot_snapshot(
                field,
                xc,
                yc,
                title,
                spec.units,
                out,
                vmin=spec.vmin,
                vmax=spec.vmax,
                center_zero=spec.center_zero,
                dpi=dpi,
            )
        saved_iterations.append(iteration)
    if saved_iterations:
        print(f"wrote {spec.display}: {output_root / spec.folder}")
    else:
        print(f"skip {spec.display}: no finite snapshots")
    return SnapshotResult(
        spec.folder,
        spec.field,
        spec.display,
        [int(value) for value in iterations],
        saved_iterations,
        skipped_nonfinite_iterations,
    )


def save_velocity_magnitude(
    case_code: str,
    run_dir: Path,
    output_root: Path,
    delta_t: float,
    *,
    u_candidates: tuple[str, ...] = ("U", "UVEL", "UVELMASS"),
    v_candidates: tuple[str, ...] = ("V", "VVEL", "VVELMASS"),
    dpi: int = 220,
) -> SnapshotResult:
    u_name = first_existing_field(run_dir, u_candidates)
    v_name = first_existing_field(run_dir, v_candidates)
    if u_name is None or v_name is None:
        print("skip velocity magnitude: no U/V files")
        return SnapshotResult("velocity_magnitude", "velocity_magnitude", "velocity magnitude", [], [], [])

    iterations = sorted(set(discover_iterations(run_dir, u_name)) & set(discover_iterations(run_dir, v_name)))
    if not iterations:
        print("skip velocity magnitude: U/V iterations do not match")
        return SnapshotResult("velocity_magnitude", "velocity_magnitude", "velocity magnitude", [], [], [])

    xc, yc = lon_lat(run_dir)
    cube_plot = is_cs_compact(xc) and is_cs_compact(yc)
    case = case_code.lower()
    saved_iterations: list[int] = []
    skipped_nonfinite_iterations: list[int] = []
    for iteration in iterations:
        u = center_to_shape(read_mds_field(run_dir, u_name, iteration), xc.shape)
        v = center_to_shape(read_mds_field(run_dir, v_name, iteration), xc.shape)
        speed = np.sqrt(u * u + v * v)
        value_range = finite_min_max(speed)
        if value_range is None:
            print(f"skip velocity magnitude iteration {iteration}: no finite values")
            skipped_nonfinite_iterations.append(iteration)
            continue
        missing = int(speed.size - np.count_nonzero(np.isfinite(speed)))
        if missing:
            print(f"warning velocity magnitude iteration {iteration}: {missing} non-finite values")
        speed_min, speed_max = value_range
        title = (
            f"{case_code} velocity magnitude | iteration {iteration} | "
            f"time {time_text(iteration, delta_t)} | "
            f"min {speed_min:.3e} | max {speed_max:.3e}"
        )
        out = output_root / "velocity_magnitude" / f"{case}_velocity_magnitude_day_{day_number(iteration, delta_t):05.2f}_iter_{iteration:010d}.pdf"
        if cube_plot and is_cs_compact(speed):
            plot_cube_net_snapshot(speed, xc, yc, title, r"m s$^{-1}$", out, vmin=0.0, dpi=dpi)
        else:
            plot_snapshot(speed, xc, yc, title, r"m s$^{-1}$", out, vmin=0.0, dpi=dpi)
        saved_iterations.append(iteration)
    if saved_iterations:
        print(f"wrote velocity magnitude: {output_root / 'velocity_magnitude'}")
    else:
        print("skip velocity magnitude: no finite snapshots")
    return SnapshotResult(
        "velocity_magnitude",
        "velocity_magnitude",
        "velocity magnitude",
        [int(value) for value in iterations],
        saved_iterations,
        skipped_nonfinite_iterations,
    )


def run_snapshots(
    case_code: str,
    run_dir: Path,
    specs: tuple[SnapshotSpec, ...],
    *,
    save_velocity: bool = True,
    dpi: int = 220,
    fail_on_invalid: bool = True,
) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    alpha = infer_alpha_label(run_dir, case_code)
    output_root = snapshot_output_dir(case_code, alpha)
    output_root.mkdir(parents=True, exist_ok=True)

    saved = []
    results: list[SnapshotResult] = []
    for spec in specs:
        clear_generated_snapshot_images(output_root, spec.folder)
        result = save_scalar_series(case_code, run_dir, output_root, spec, delta_t, dpi=dpi)
        results.append(result)
        if result.wrote_any:
            saved.append(spec.folder)
    if save_velocity:
        clear_generated_snapshot_images(output_root, "velocity_magnitude")
        result = save_velocity_magnitude(case_code, run_dir, output_root, delta_t, dpi=dpi)
        results.append(result)
        if result.wrote_any:
            saved.append("velocity_magnitude")

    invalid = [result for result in results if result.invalid]
    field_results = {
        result.folder: {
            "field": result.field,
            "display": result.display,
            "iterations": [int(value) for value in result.iterations],
            "saved_iterations": [int(value) for value in result.saved_iterations],
            "skipped_nonfinite_iterations": [
                int(value) for value in result.skipped_nonfinite_iterations
            ],
        }
        for result in results
    }
    write_manifest(
        output_root,
        {
            "alpha": alpha,
            "case": case_code,
            "product": "snapshots",
            "saved": saved,
            "source_run": str(run_dir),
            "valid": not invalid,
            "field_results": field_results,
        },
    )
    print(f"outputs written to: {output_root}")
    if fail_on_invalid and invalid:
        details = "; ".join(
            f"{result.display}: {len(result.skipped_nonfinite_iterations)} non-finite iteration(s)"
            for result in invalid
        )
        raise RuntimeError(f"{case_code} snapshots are invalid: {details}")
    return output_root
