#!/usr/bin/env python3
"""Generate, verify, and clean advect_cs Williamson TC1 assets."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

REPO_DIR = Path(__file__).resolve().parents[3]
SETUP_ROOT = REPO_DIR / "Sandbox" / "MITgcm_Williamson_TC1" / "existingTutorials" / "advect_cs"
OUTPUT_ROOT = (
    REPO_DIR
    / "Sandbox"
    / "output"
    / "existingTutorials"
    / "test1"
    / "MITGCM_Williamson_TC1"
    / "advect_cs"
)

DAY = 86400.0
DELTA_T = 2700.0
# Match ini_vel.F: one full solid-body revolution in 12 days.
ROTATION_RATE = 2.0 * math.pi / (12.0 * DAY)
REQUESTED_SNAPSHOT_DAYS = (0, 3, 6, 9, 10, 12)
SNAPSHOT_FIELDS = ("tracer", "theta", "eta", "velocity_magnitude", "tracer_error")
ALPHA_ORDER = ("0", "0.05", "1.52", "1.57")
FACE_SIZE = 32

RUN_PRODUCT_PATTERNS = (
    "*.data",
    "*.meta",
    "STDOUT*",
    "STDERR*",
    "pickup*",
    "available_diagnostics.log",
    "mitgcmuv",
)


def alpha_sort_key(path: Path) -> Tuple[int, str]:
    label = path.name.removeprefix("alpha_").removeprefix("run_alpha_")
    if label in ALPHA_ORDER:
        return (ALPHA_ORDER.index(label), label)
    return (len(ALPHA_ORDER), label)


def alpha_from_path(path: Path) -> str:
    for part in path.parts:
        if part.startswith("alpha_"):
            return part.removeprefix("alpha_")
        if part.startswith("run_alpha_"):
            return part.removeprefix("run_alpha_")
    return "unknown"


def alpha_value(run_dir: Path) -> float:
    alpha_file = run_dir / "tc1_alpha.txt"
    if alpha_file.exists():
        return float(alpha_file.read_text().strip())
    label = alpha_from_path(run_dir)
    if label == "1.52":
        return math.pi / 2.0 - 0.05
    if label == "1.57":
        return math.pi / 2.0
    return float(label)


def safe_remove_tree(path: Path) -> None:
    path = path.resolve()
    root = OUTPUT_ROOT.resolve()
    if root not in path.parents and path != root:
        raise ValueError(f"refusing unsafe removal target: {path}")
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def save_figure(fig: plt.Figure, path_without_suffix: Path, dpi: int = 180) -> List[str]:
    path_without_suffix.parent.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    for suffix in (".pdf", ".png"):
        out = Path(f"{path_without_suffix}{suffix}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        written.append(str(out))
    plt.close(fig)
    return written


def parse_data_value(data_path: Path, key: str, default: float) -> float:
    if not data_path.exists():
        return default
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*([^,]+)", re.IGNORECASE)
    for line in data_path.read_text(errors="ignore").splitlines():
        clean = line.split("#", 1)[0]
        match = pattern.search(clean)
        if match:
            return float(match.group(1).replace("d", "e").replace("D", "e"))
    return default


def parse_meta(meta_path: Path) -> Tuple[List[int], List[int], List[int], str, int, List[str]]:
    text = meta_path.read_text(errors="ignore")
    dim_match = re.search(r"dimList\s*=\s*\[([^\]]+)\]", text, re.S)
    if not dim_match:
        raise ValueError(f"missing dimList in {meta_path}")
    dim_nums = [int(value) for value in re.findall(r"[-+]?\d+", dim_match.group(1))]
    triples = [dim_nums[index : index + 3] for index in range(0, len(dim_nums), 3)]
    dims = [triple[0] for triple in triples]
    starts = [triple[1] for triple in triples]
    ends = [triple[2] for triple in triples]
    prec_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'", text)
    data_prec = prec_match.group(1).lower() if prec_match else "float64"
    rec_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)", text)
    nrecords = int(rec_match.group(1)) if rec_match else 1
    fields = [field.strip() for field in re.findall(r"'([^']+)'", text.split("fldList", 1)[-1])]
    return dims, starts, ends, data_prec, nrecords, fields


def read_mds(prefix: Path, record: Optional[int] = None) -> np.ndarray:
    meta_path = Path(f"{prefix}.meta")
    data_path = Path(f"{prefix}.data")
    _dims, starts, ends, data_prec, nrecords, _fields = parse_meta(meta_path)
    dtype = ">f8" if "64" in data_prec else ">f4"
    raw = np.fromfile(data_path, dtype=dtype)
    dims = [end - start + 1 for start, end in zip(starts, ends)]
    shape = tuple(reversed(dims))
    expected = int(np.prod(shape)) * nrecords
    if raw.size != expected:
        raise ValueError(f"{data_path} has {raw.size} values, expected {expected}")
    if nrecords > 1:
        arr = raw.reshape((nrecords, *shape))
        if record is None:
            return arr
        return np.squeeze(arr[record])
    return np.squeeze(raw.reshape(shape))


def tiled_meta_paths(run_dir: Path, name: str, iteration: int) -> List[Path]:
    return sorted(run_dir.glob(f"{name}.{iteration:010d}.[0-9][0-9][0-9].[0-9][0-9][0-9].meta"))


def read_mds_tiled(paths: Sequence[Path], record: Optional[int] = None) -> np.ndarray:
    if not paths:
        raise FileNotFoundError("no tiled MDS files")
    dims, _starts, _ends, _prec, nrecords, _fields = parse_meta(paths[0])
    global_shape = tuple(reversed(dims))
    if nrecords > 1 and record is None:
        assembled = np.full((nrecords, *global_shape), np.nan, dtype=np.float64)
    else:
        assembled = np.full(global_shape, np.nan, dtype=np.float64)
    for meta_path in paths:
        _dims, starts, ends, _prec, _nrecords, _fields = parse_meta(meta_path)
        prefix = meta_path.with_suffix("")
        tile = read_mds(prefix, record=record)
        if len(starts) != 2:
            raise ValueError(f"expected 2D tiled MDS field: {meta_path}")
        x_slice = slice(starts[0] - 1, ends[0])
        y_slice = slice(starts[1] - 1, ends[1])
        if nrecords > 1 and record is None:
            assembled[:, y_slice, x_slice] = tile
        else:
            assembled[y_slice, x_slice] = tile
    return np.squeeze(assembled)


def read_mds_named(run_dir: Path, name: str, iteration: int, record: Optional[int] = None) -> Optional[np.ndarray]:
    prefix = mds_prefix(run_dir, name, iteration)
    if prefix is not None:
        return read_mds(prefix, record=record)
    tile_paths = tiled_meta_paths(run_dir, name, iteration)
    if tile_paths:
        return read_mds_tiled(tile_paths, record=record)
    return None


def read_init_faces(path: Path) -> np.ndarray:
    raw = np.fromfile(path, dtype=">f8")
    if raw.size != FACE_SIZE * 6 * FACE_SIZE:
        raise ValueError(f"unexpected init size for {path}: {raw.size}")
    return np.moveaxis(raw.reshape((FACE_SIZE, 6, FACE_SIZE), order="F"), 1, 0)


def read_grid_faces(run_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    lon_faces = []
    lat_faces = []
    for face in range(1, 7):
        path = run_dir / f"grid_cs32.face{face:03d}.bin"
        raw = np.fromfile(path, dtype=">f8")
        if raw.size < 2 * 33 * 33:
            raise ValueError(f"unexpected grid size for {path}: {raw.size}")
        lon = raw[: 33 * 33].reshape((33, 33), order="F")[:FACE_SIZE, :FACE_SIZE]
        lat = raw[33 * 33 : 2 * 33 * 33].reshape((33, 33), order="F")[:FACE_SIZE, :FACE_SIZE]
        lon_faces.append(lon)
        lat_faces.append(lat)
    return np.asarray(lon_faces), np.asarray(lat_faces)


def compact_to_faces(field: np.ndarray) -> np.ndarray:
    arr = np.squeeze(np.asarray(field, dtype=np.float64))
    if arr.shape == (6, FACE_SIZE, FACE_SIZE):
        return arr
    if arr.shape == (FACE_SIZE, 6, FACE_SIZE):
        return np.moveaxis(arr, 1, 0)
    if arr.shape == (FACE_SIZE, FACE_SIZE * 6):
        return np.asarray(
            [arr[:, index * FACE_SIZE : (index + 1) * FACE_SIZE] for index in range(6)]
        )
    if arr.shape == (FACE_SIZE * 3, FACE_SIZE * 2):
        faces = []
        for row in range(3):
            for col in range(2):
                faces.append(arr[row * FACE_SIZE : (row + 1) * FACE_SIZE, col * FACE_SIZE : (col + 1) * FACE_SIZE])
        return np.asarray(faces)
    if arr.shape == (FACE_SIZE * 2, FACE_SIZE * 3):
        faces = []
        for col in range(3):
            for row in range(2):
                faces.append(arr[row * FACE_SIZE : (row + 1) * FACE_SIZE, col * FACE_SIZE : (col + 1) * FACE_SIZE])
        return np.asarray(faces[:6])
    raise ValueError(f"cannot convert field with shape {arr.shape} to six faces")


def faces_to_net(faces: np.ndarray) -> np.ndarray:
    net = np.full((FACE_SIZE * 3, FACE_SIZE * 4), np.nan, dtype=np.float64)
    positions = {
        0: (1, 0),
        1: (1, 1),
        2: (1, 2),
        3: (1, 3),
        4: (0, 1),
        5: (2, 1),
    }
    for face_index, (row, col) in positions.items():
        net[row * FACE_SIZE : (row + 1) * FACE_SIZE, col * FACE_SIZE : (col + 1) * FACE_SIZE] = faces[face_index]
    return net


def unit_advect(lon_deg: np.ndarray, lat_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    return np.cos(lat) * np.sin(lon), -np.cos(lat) * np.cos(lon), np.sin(lat)


def rotate_about_axis(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, alpha: float, angle: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ax, ay, az = 0.0, math.sin(alpha), math.cos(alpha)
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


def advect_bell_exact(
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    alpha: float,
    seconds: float,
    center_lon: float,
    center_lat: float,
) -> np.ndarray:
    x, y, z = unit_advect(lon_deg, lat_deg)
    xb, yb, zb = rotate_about_axis(x, y, z, alpha, -ROTATION_RATE * seconds)
    cx, cy, cz = unit_advect(center_lon, center_lat)
    chord = np.sqrt((xb - cx) ** 2 + (yb - cy) ** 2 + (zb - cz) ** 2)
    bell = 1.0 + 0.5 * (1.0 + np.cos(math.pi * np.minimum(chord / 0.3, 1.0)))
    return bell


def analytic_velocity_faces(lon_deg: np.ndarray, lat_deg: np.ndarray, alpha: float) -> np.ndarray:
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    u0 = ROTATION_RATE * 6_371_000.0
    u = u0 * (np.cos(lat) * math.cos(alpha) + np.cos(lon) * np.sin(lat) * math.sin(alpha))
    v = -u0 * np.sin(lon) * math.sin(alpha)
    return np.sqrt(u * u + v * v)


def mds_prefix(run_dir: Path, name: str, iteration: int) -> Optional[Path]:
    prefix = run_dir / f"{name}.{iteration:010d}"
    if Path(f"{prefix}.meta").exists() and Path(f"{prefix}.data").exists():
        return prefix
    return None


def read_diag_field(run_dir: Path, iteration: int, wanted: str) -> Optional[np.ndarray]:
    prefix = mds_prefix(run_dir, "ts_Diag", iteration)
    if prefix is not None:
        meta_path = Path(f"{prefix}.meta")
    else:
        tile_paths = tiled_meta_paths(run_dir, "ts_Diag", iteration)
        if not tile_paths:
            return None
        meta_path = tile_paths[0]
    _dims, _starts, _ends, _prec, _nrecords, fields = parse_meta(meta_path)
    wanted = wanted.upper()
    for index, field in enumerate(fields):
        if field.upper() == wanted:
            return read_mds_named(run_dir, "ts_Diag", iteration, record=index)
    fallback = {"THETA": 0, "SALT": 1}.get(wanted)
    if fallback is not None:
        return read_mds_named(run_dir, "ts_Diag", iteration, record=fallback)
    return None


def read_model_field(run_dir: Path, iteration: int, field: str, lon: np.ndarray, lat: np.ndarray, alpha: float) -> np.ndarray:
    if iteration == 0:
        if field == "tracer":
            return read_init_faces(run_dir / "S.init")
        if field == "theta":
            return read_init_faces(run_dir / "T.init")
        if field == "eta":
            return np.zeros((6, FACE_SIZE, FACE_SIZE), dtype=np.float64)
        if field == "velocity_magnitude":
            return analytic_velocity_faces(lon, lat, alpha)

    candidate_names = {
        "tracer": ("S", "SALT"),
        "theta": ("T", "THETA"),
        "eta": ("Eta", "ETAN"),
        "u": ("U", "UVEL"),
        "v": ("V", "VVEL"),
    }
    if field in ("tracer", "theta", "eta"):
        for name in candidate_names[field]:
            data = read_mds_named(run_dir, name, iteration)
            if data is not None:
                return compact_to_faces(data)
        diag = read_diag_field(run_dir, iteration, "SALT" if field == "tracer" else "THETA")
        if diag is not None:
            return compact_to_faces(diag)
        if field == "eta":
            return np.zeros((6, FACE_SIZE, FACE_SIZE), dtype=np.float64)
    if field == "velocity_magnitude":
        u_field = None
        v_field = None
        for name in candidate_names["u"]:
            data = read_mds_named(run_dir, name, iteration)
            if data is not None:
                u_field = compact_to_faces(data)
                break
        for name in candidate_names["v"]:
            data = read_mds_named(run_dir, name, iteration)
            if data is not None:
                v_field = compact_to_faces(data)
                break
        if u_field is not None and v_field is not None:
            return np.sqrt(u_field * u_field + v_field * v_field)
        return analytic_velocity_faces(lon, lat, alpha)
    raise FileNotFoundError(f"missing {field} at iteration {iteration} in {run_dir}")


def available_iterations(run_dir: Path) -> List[int]:
    iterations = {0}
    for path in run_dir.glob("*.meta"):
        match = re.search(r"\.(\d{10})(?:\.\d{3}\.\d{3})?\.meta$", path.name)
        if match:
            iterations.add(int(match.group(1)))
    return sorted(iterations)


def requested_iterations(delta_t: float) -> Dict[int, int]:
    return {day: int(round(day * DAY / delta_t)) for day in REQUESTED_SNAPSHOT_DAYS}


def plot_net(
    faces: np.ndarray,
    output_base: Path,
    title: str,
    units: str,
    *,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    center_zero: bool = False,
) -> List[str]:
    net = faces_to_net(faces)
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    cmap = plt.get_cmap("seismic").copy()
    cmap.set_bad("white")
    norm = None
    if center_zero:
        limit = float(np.nanmax(np.abs(net)))
        if limit == 0.0:
            limit = 1.0
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    image = ax.imshow(np.ma.masked_invalid(net), origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    for xline in range(0, FACE_SIZE * 4 + 1, FACE_SIZE):
        ax.axvline(xline - 0.5, color="0.25", linewidth=0.8)
    for yline in range(0, FACE_SIZE * 3 + 1, FACE_SIZE):
        ax.axhline(yline - 0.5, color="0.25", linewidth=0.8)
    labels = {
        "face 5": (FACE_SIZE * 1 + 2, FACE_SIZE * 2 + 25),
        "face 1": (2, FACE_SIZE * 1 + 25),
        "face 2": (FACE_SIZE * 1 + 2, FACE_SIZE * 1 + 25),
        "face 3": (FACE_SIZE * 2 + 2, FACE_SIZE * 1 + 25),
        "face 4": (FACE_SIZE * 3 + 2, FACE_SIZE * 1 + 25),
        "face 6": (FACE_SIZE * 1 + 2, 25),
    }
    for label, (xpos, ypos) in labels.items():
        ax.text(xpos, ypos, label, color="0.2", fontsize=9, alpha=0.75)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    cbar = fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.12, fraction=0.08)
    cbar.set_label(units)
    return save_figure(fig, output_base)


def write_experiment_file(path: Path, alpha_label: str, alpha: float, snapshot_iterations: Sequence[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "experiment: MITGCM_Williamson_TC1/existingTutorials/advect_cs",
                "source_verification: verification/advect_cs",
                "test_case: Williamson TC1 passive tracer analogue",
                f"alpha: {alpha}",
                f"alpha_label: {alpha_label}",
                "tutorial_name: advect_cs",
                "logic: passive tracer S is advected by the TC1 solid-body velocity; momentum is not stepped.",
                f"requested_snapshot_days: {', '.join(str(day) for day in REQUESTED_SNAPSHOT_DAYS)}",
                f"snapshot_iterations: {', '.join(str(item) for item in snapshot_iterations)}",
                "plot_style: unfolded cubed-sphere net with seismic colormap",
                "",
            ]
        )
    )


def clean_alpha_assets(alpha_label: str) -> None:
    for target in (
        OUTPUT_ROOT / "Snapshots" / f"alpha_{alpha_label}",
        OUTPUT_ROOT / "Diagnosis" / "error" / f"alpha_{alpha_label}",
        OUTPUT_ROOT / "Diagnosis" / "postprocessing" / f"alpha_{alpha_label}",
    ):
        safe_remove_tree(target)


def process_run(run_dir: Path) -> Dict[str, object]:
    alpha_label = alpha_from_path(run_dir)
    alpha = alpha_value(run_dir)
    delta_t = parse_data_value(run_dir / "data", "deltaT", DELTA_T)
    day_to_iter = requested_iterations(delta_t)
    lon_faces, lat_faces = read_grid_faces(run_dir)
    iterations = available_iterations(run_dir)
    outputs: List[str] = []
    clean_alpha_assets(alpha_label)

    snap_root = OUTPUT_ROOT / "Snapshots" / f"alpha_{alpha_label}"
    error_root = OUTPUT_ROOT / "Diagnosis" / "error" / f"alpha_{alpha_label}"
    post_root = OUTPUT_ROOT / "Diagnosis" / "postprocessing" / f"alpha_{alpha_label}"

    metric_rows: List[Dict[str, float]] = []
    integral_rows: List[Dict[str, float]] = []
    snapshot_iterations: List[int] = []

    for day, iteration in day_to_iter.items():
        seconds = day * DAY
        exact = advect_bell_exact(lon_faces, lat_faces, alpha, seconds, 180.0, 35.0)
        tracer = read_model_field(run_dir, iteration, "tracer", lon_faces, lat_faces, alpha)
        theta = read_model_field(run_dir, iteration, "theta", lon_faces, lat_faces, alpha)
        eta = read_model_field(run_dir, iteration, "eta", lon_faces, lat_faces, alpha)
        velocity = read_model_field(run_dir, iteration, "velocity_magnitude", lon_faces, lat_faces, alpha)
        error = tracer - exact
        snapshot_iterations.append(iteration)
        fields = (
            ("tracer", "S", tracer, "psu", 1.0, 2.0, False),
            ("theta", "T", theta, "degC", 1.0, 2.0, False),
            ("eta", "Eta", eta, "m", None, None, True),
            ("velocity_magnitude", "|u|", velocity, "m s^-1", None, None, False),
            ("tracer_error", "S error", error, "psu", None, None, True),
        )
        for folder, label, data, units, vmin, vmax, center_zero in fields:
            outputs.extend(
                plot_net(
                    data,
                    snap_root / folder / f"advect_cs_alpha_{alpha_label}_day_{day:02d}",
                    f"advect_cs {label} | alpha={alpha_label} | day={day:.2f}",
                    units,
                    vmin=vmin,
                    vmax=vmax,
                    center_zero=center_zero,
                )
            )

    for iteration in iterations:
        day = iteration * delta_t / DAY
        seconds = iteration * delta_t
        exact = advect_bell_exact(lon_faces, lat_faces, alpha, seconds, 180.0, 35.0)
        tracer = read_model_field(run_dir, iteration, "tracer", lon_faces, lat_faces, alpha)
        error = tracer - exact
        metric_rows.append(
            {
                "day": day,
                "iteration": float(iteration),
                "l1": float(np.nanmean(np.abs(error))),
                "l2": float(np.sqrt(np.nanmean(error * error))),
                "linf": float(np.nanmax(np.abs(error))),
                "tracer_min": float(np.nanmin(tracer)),
                "tracer_max": float(np.nanmax(tracer)),
            }
        )
        integral_rows.append(
            {
                "day": day,
                "iteration": float(iteration),
                "tracer_mean": float(np.nanmean(tracer)),
                "tracer_sum": float(np.nansum(tracer)),
            }
        )

    error_root.mkdir(parents=True, exist_ok=True)
    with (error_root / "advect_cs_error_table.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metric_rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metric_rows)
    with (post_root / "postprocessing_values.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(integral_rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(integral_rows)

    days = np.array([row["day"] for row in metric_rows])
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogy(days, [row["l1"] for row in metric_rows], label="L1")
    ax.semilogy(days, [row["l2"] for row in metric_rows], label="L2")
    ax.semilogy(days, [row["linf"] for row in metric_rows], label="Linf")
    ax.set_xlabel("day")
    ax.set_ylabel("error norm [psu]")
    ax.set_title(f"advect_cs error norms | alpha={alpha_label}")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    outputs.extend(save_figure(fig, error_root / "advect_cs_error_metrics"))

    final_day = REQUESTED_SNAPSHOT_DAYS[-1]
    final_iter = day_to_iter[final_day]
    final_seconds = final_day * DAY
    final_exact = advect_bell_exact(lon_faces, lat_faces, alpha, final_seconds, 180.0, 35.0)
    final_model = read_model_field(run_dir, final_iter, "tracer", lon_faces, lat_faces, alpha)
    outputs.extend(plot_net(final_model, error_root / "advect_cs_final_tracer_model", f"advect_cs final model tracer | alpha={alpha_label}", "psu", vmin=1.0, vmax=2.0))
    outputs.extend(plot_net(final_exact, error_root / "advect_cs_final_tracer_exact", f"advect_cs final exact tracer | alpha={alpha_label}", "psu", vmin=1.0, vmax=2.0))
    outputs.extend(plot_net(final_model - final_exact, error_root / "advect_cs_final_tracer_error", f"advect_cs final tracer error | alpha={alpha_label}", "psu", center_zero=True))

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot([row["day"] for row in integral_rows], [row["tracer_sum"] for row in integral_rows], marker="o", ms=3)
    ax.set_xlabel("day")
    ax.set_ylabel("sum(S)")
    ax.set_title(f"advect_cs tracer integral | alpha={alpha_label}")
    ax.grid(True, alpha=0.25)
    outputs.extend(save_figure(fig, post_root / "postprocessing_tracer_integral"))

    velocity = read_model_field(run_dir, final_iter, "velocity_magnitude", lon_faces, lat_faces, alpha)
    outputs.extend(plot_net(velocity, post_root / "postprocessing_velocity", f"advect_cs velocity magnitude | alpha={alpha_label}", "m s^-1"))

    write_experiment_file(snap_root / "experiment.txt", alpha_label, alpha, snapshot_iterations)
    write_experiment_file(error_root / "experiment.txt", alpha_label, alpha, snapshot_iterations)
    write_experiment_file(post_root / "experiment.txt", alpha_label, alpha, snapshot_iterations)

    manifest = {
        "alpha": alpha_label,
        "alpha_value": alpha,
        "delta_t": delta_t,
        "experiment": "MITGCM_Williamson_TC1/existingTutorials/advect_cs",
        "logic": "Reads advect_cs MDS output and compares S against the advect_cs initial bell transported by the same TC1 solid-body rotation.",
        "iterations": iterations,
        "requested_snapshot_days": REQUESTED_SNAPSHOT_DAYS,
        "snapshot_iterations": snapshot_iterations,
        "outputs": outputs,
    }
    for folder in (snap_root, error_root, post_root):
        (folder / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def run_dirs() -> List[Path]:
    return sorted(
        [path for path in OUTPUT_ROOT.glob("alpha_*/run_alpha_*") if path.is_dir()],
        key=alpha_sort_key,
    )


def count_outputs() -> Dict[str, int]:
    return {
        "png": len(list(OUTPUT_ROOT.rglob("*.png"))),
        "pdf": len(list(OUTPUT_ROOT.rglob("*.pdf"))),
        "csv": len(list(OUTPUT_ROOT.rglob("*.csv"))),
        "json": len(list(OUTPUT_ROOT.rglob("*.json"))),
    }


def existing_snapshot_days(alpha: str, field: str) -> List[int]:
    field_dir = OUTPUT_ROOT / "Snapshots" / f"alpha_{alpha}" / field
    days: set[int] = set()
    for path in field_dir.glob("*.png"):
        match = re.search(r"_day_([0-9]+)", path.name)
        if match:
            days.add(int(match.group(1)))
    return sorted(days)


def missing_snapshot_days() -> Dict[str, Dict[str, List[int]]]:
    missing: Dict[str, Dict[str, List[int]]] = {}
    alpha_dirs = sorted((OUTPUT_ROOT / "Snapshots").glob("alpha_*"), key=alpha_sort_key)
    for alpha_dir in alpha_dirs:
        if not alpha_dir.is_dir():
            continue
        alpha = alpha_dir.name.removeprefix("alpha_")
        for field in SNAPSHOT_FIELDS:
            existing = set(existing_snapshot_days(alpha, field))
            field_missing = [day for day in REQUESTED_SNAPSHOT_DAYS if day not in existing]
            if field_missing:
                missing.setdefault(alpha, {})[field] = field_missing
    return missing


def clean_run_products() -> int:
    removed = 0
    for run_dir in sorted(OUTPUT_ROOT.glob("alpha_*/run_alpha_*"), key=alpha_sort_key):
        if not run_dir.is_dir():
            continue
        resolved = run_dir.resolve()
        root = OUTPUT_ROOT.resolve()
        if root not in resolved.parents:
            raise ValueError(f"refusing unsafe cleanup target: {run_dir}")
        for path in list(run_dir.iterdir()):
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()
            removed += 1
    return removed


def existing_manifests() -> List[Dict[str, object]]:
    manifests: List[Dict[str, object]] = []
    paths = sorted((OUTPUT_ROOT / "Snapshots").glob("alpha_*/manifest.json"), key=lambda item: alpha_sort_key(item.parent))
    for path in paths:
        manifests.append(json.loads(path.read_text()))
    return manifests


def write_summary(manifests: Sequence[Dict[str, object]], removed: int) -> Dict[str, object]:
    summary_path = OUTPUT_ROOT / "postprocess_summary.json"
    if not manifests:
        manifests = existing_manifests()
    if not manifests and summary_path.exists():
        previous = json.loads(summary_path.read_text())
        manifests = previous.get("runs", [])
    summary = {
        "output_root": str(OUTPUT_ROOT),
        "counts": count_outputs(),
        "requested_snapshot_days": REQUESTED_SNAPSHOT_DAYS,
        "missing_snapshot_days": missing_snapshot_days(),
        "removed_run_products": removed,
        "runs": list(manifests),
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--clean-run-products", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    manifests: List[Dict[str, object]] = []
    if args.generate:
        dirs = run_dirs()
        if not dirs:
            raise SystemExit(f"no run directories found under {OUTPUT_ROOT}")
        for run_dir in dirs:
            manifests.append(process_run(run_dir))

    removed = clean_run_products() if args.clean_run_products else 0
    summary = write_summary(manifests, removed)
    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(f"output_root: {summary['output_root']}")
        print(f"counts: {summary['counts']}")
        print(f"requested_snapshot_days: {', '.join(str(day) for day in REQUESTED_SNAPSHOT_DAYS)}")
        print(f"missing_snapshot_days: {summary['missing_snapshot_days'] or 'none'}")
        print(f"removed_run_products: {removed}")


if __name__ == "__main__":
    main()
