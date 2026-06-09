#!/usr/bin/env python3
"""Williamson TC1 snapshots from MITgcm MDS output.

Outputs are written to:
    ../output/<run_name>/Eta
    ../output/<run_name>/VelocityMagnitude
    ../output/<run_name>/Psi
    ../output/<run_name>/Phi

Run:
    set RUN_DIR below, then run this file
"""
#%%
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
})

# ------------------------- easy settings -------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
CASE_DIR = SCRIPT_DIR.parents[1] if len(SCRIPT_DIR.parents) > 1 else SCRIPT_DIR

RUN_DIR = CASE_DIR / "run_alpha_1.57"   # edit this directly
ETA_FIELD = "S"                   # use "S" or "ETAN"
U_CANDIDATES = ("U", "UVEL", "UVELMASS")
V_CANDIDATES = ("V", "VVEL", "VVELMASS")
PSI_FIELD = "PsiVEL"
PHI_FIELD = "PhiVEL"

ETA_LIMITS = (0.0, 1000.0)            # cosine bell: 0 ... 1000
CMAP_POSITIVE = "seismic"
CMAP_SIGNED = "seismic"
DPI = 220
# -----------------------------------------------------------------

DIAGNOSTIC_STREAM_FOR_FIELD = {
    "ETAN": "dynDiag",
    "PhiVEL": "dyn_Aux",
    "PsiVEL": "dyn_Aux",
}

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
    fld_list_match = re.search(r"fldList\s*=\s*\{(.*?)\}", text, re.S)
    fld_list = []
    if fld_list_match is not None:
        fld_list = [field.strip() for field in re.findall(r"'([^']+)'", fld_list_match.group(1))]
    return {
        "n_dims": n_dims,
        "dims": dims,
        "prec": prec.group(1).lower() if prec else "float32",
        "nrecords": int(nrecords.group(1)) if nrecords else 1,
        "fld_list": fld_list,
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

    record = 0
    if not meta_path.exists() or not data_path.exists():
        diag_stream = DIAGNOSTIC_STREAM_FOR_FIELD.get(name.strip())
        if diag_stream is None:
            raise FileNotFoundError(f"Missing {name} in {run_dir}")
        if iteration is None:
            meta_path = run_dir / f"{diag_stream}.meta"
        else:
            meta_path = run_dir / f"{diag_stream}.{iteration:010d}.meta"
        data_path = meta_path.with_suffix(".data")
        if not meta_path.exists() or not data_path.exists():
            raise FileNotFoundError(f"Missing {name} / {diag_stream} in {run_dir}")
        meta = parse_mds_meta(meta_path)
        fld_list = [field.strip() for field in meta.get("fld_list", [])]
        if name.strip() not in fld_list:
            raise FileNotFoundError(f"{name} not found in {diag_stream} in {run_dir}")
        record = fld_list.index(name.strip())
    else:
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
    return np.squeeze(data[record] if nrecords > 1 else data)

def discover_iterations(run_dir: Path, field: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(field)}\.(\d{{10}})\.meta$")
    out = []
    for path in run_dir.glob(f"{field}.*.meta"):
        match = pattern.match(path.name)
        if match:
            out.append(int(match.group(1)))
    if out:
        return sorted(set(out))

    diag_stream = DIAGNOSTIC_STREAM_FOR_FIELD.get(field.strip())
    if diag_stream is None:
        return []
    diag_pattern = re.compile(rf"^{re.escape(diag_stream)}\.(\d{{10}})\.meta$")
    for path in run_dir.glob(f"{diag_stream}.*.meta"):
        match = diag_pattern.match(path.name)
        if not match:
            continue
        meta = parse_mds_meta(path)
        fld_list = [name.strip() for name in meta.get("fld_list", [])]
        if field.strip() in fld_list:
            out.append(int(match.group(1)))
    return sorted(set(out))

def first_existing_field(run_dir: Path, names: tuple[str, ...]) -> str | None:
    for name in names:
        if discover_iterations(run_dir, name):
            return name
    return None

def to_2d(field: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    if field.ndim == 3 and field.shape[0] == 1:
        return field[0]
    if field.ndim != 2:
        raise ValueError(f"Expected 2D or 1-level 3D field, got {field.shape}")
    return field

def center_to_shape(field: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Average a C-grid face field to tracer-cell shape if needed."""
    field = to_2d(field)
    ny, nx = shape
    if field.shape == shape:
        return field
    if field.shape == (ny, nx + 1):
        return 0.5 * (field[:, :-1] + field[:, 1:])
    if field.shape == (ny + 1, nx):
        return 0.5 * (field[:-1, :] + field[1:, :])
    return field[:ny, :nx]

def lon_lat(run_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    xc = to_2d(read_mds_field(run_dir, "XC"))
    yc = to_2d(read_mds_field(run_dir, "YC"))
    return xc, yc

def extent_from_grid(xc: np.ndarray, yc: np.ndarray) -> tuple[float, float, float, float]:
    lon = xc[0, :]
    lat = yc[:, 0]
    dlon = float(np.median(np.diff(lon))) if lon.size > 1 else 1.0
    dlat = float(np.median(np.diff(lat))) if lat.size > 1 else 1.0
    return (lon[0] - 0.5 * dlon, lon[-1] + 0.5 * dlon,
            lat[0] - 0.5 * dlat, lat[-1] + 0.5 * dlat)

def time_text(iteration: int, delta_t: float) -> str:
    seconds = int(round(iteration * delta_t))
    day = seconds // 86400
    seconds -= day * 86400
    hour = seconds // 3600
    seconds -= hour * 3600
    minute = seconds // 60
    return f"{day:d} d {hour:02d} h {minute:02d} min"

def day_number(iteration: int, delta_t: float) -> float:
    return iteration * delta_t / 86400.0

def finite_min_max(field: np.ndarray) -> tuple[float, float] | None:
    finite = np.isfinite(field)
    if not np.any(finite):
        return None
    values = field[finite]
    return float(np.min(values)), float(np.max(values))

def plot_snapshot(
    field: np.ndarray,
    xc: np.ndarray,
    yc: np.ndarray,
    title: str,
    units: str,
    out_path: Path,
    vmin: float | None = None,
    vmax: float | None = None,
    symmetric: bool = False,
) -> None:
    field = to_2d(field)
    finite = np.isfinite(field)
    if symmetric:
        limit = float(np.nanmax(np.abs(field[finite]))) if np.any(finite) else 1.0
        vmin, vmax = -limit, limit
    elif vmin is None or vmax is None:
        vmin = float(np.nanmin(field[finite])) if np.any(finite) else 0.0
        vmax = float(np.nanmax(field[finite])) if np.any(finite) else 1.0
        if vmin == vmax:
            vmax = vmin + 1.0

    fig, ax = plt.subplots(figsize=(9.0, 5.8), constrained_layout=True)
    image = ax.imshow(
        field,
        origin="lower",
        extent=extent_from_grid(xc, yc),
        interpolation="bilinear",
        aspect="auto",
        cmap=CMAP_SIGNED if symmetric else CMAP_POSITIVE,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title, pad=18)
    ax.set_xlabel(r"longitude $\lambda$ [deg]")
    ax.set_ylabel(r"latitude $\theta$ [deg]")
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    cbar = fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.10)
    cbar.set_label(units)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def save_scalar_series(
    run_dir: Path,
    output_root: Path,
    field_name: str,
    folder: str,
    display: str,
    units: str,
    delta_t: float,
    vmin: float | None = None,
    vmax: float | None = None,
    symmetric: bool = False,
) -> None:
    iterations = discover_iterations(run_dir, field_name)
    if not iterations:
        print(f"skip {display}: no {field_name} files")
        return

    xc, yc = lon_lat(run_dir)
    for iteration in iterations:
        field = read_mds_field(run_dir, field_name, iteration)
        if field.ndim == 3 and field.shape[0] == 1:
            field = field[0]
        value_range = finite_min_max(field)
        if value_range is None:
            print(
                f"skip {display} iteration {iteration}: "
                "field contains no finite values"
            )
            continue
        non_finite_count = int(field.size - np.count_nonzero(np.isfinite(field)))
        if non_finite_count:
            print(
                f"warning {display} iteration {iteration}: "
                f"{non_finite_count} non-finite values"
            )
        field_min, field_max = value_range
        title = (
            f"TC1 {display} | iteration {iteration} | "
            f"time {time_text(iteration, delta_t)} | "
            f"min {field_min:.3e} | max {field_max:.3e}"
        )
        out = output_root / folder / f"tc1_{folder.lower()}_day_{day_number(iteration, delta_t):05.2f}_iter_{iteration:010d}.pdf"
        plot_snapshot(field, xc, yc, title, units, out, vmin, vmax, symmetric)
    print(f"wrote {display}: {output_root / folder}")

def save_velocity_magnitude(run_dir: Path, output_root: Path, delta_t: float) -> None:
    u_name = first_existing_field(run_dir, U_CANDIDATES)
    v_name = first_existing_field(run_dir, V_CANDIDATES)
    if u_name is None or v_name is None:
        print("skip velocity magnitude: no U/V files")
        return

    iterations = sorted(set(discover_iterations(run_dir, u_name)) & set(discover_iterations(run_dir, v_name)))
    if not iterations:
        print("skip velocity magnitude: U/V iterations do not match")
        return

    xc, yc = lon_lat(run_dir)
    shape = xc.shape
    for iteration in iterations:
        u = center_to_shape(read_mds_field(run_dir, u_name, iteration), shape)
        v = center_to_shape(read_mds_field(run_dir, v_name, iteration), shape)
        speed = np.sqrt(u * u + v * v)
        title = (
            f"TC1 velocity magnitude | iteration {iteration} | "
            f"time {time_text(iteration, delta_t)} | "
            f"min {np.nanmin(speed):.3e} | max {np.nanmax(speed):.3e}"
        )
        out = output_root / "VelocityMagnitude" / f"tc1_velocity_magnitude_day_{day_number(iteration, delta_t):05.2f}_iter_{iteration:010d}.pdf"
        plot_snapshot(speed, xc, yc, title, r"m s$^{-1}$", out, vmin=0.0, vmax=None, symmetric=True)
    print(f"wrote velocity magnitude: {output_root / 'VelocityMagnitude'}")

def run_snapshots(run_dir: Path) -> None:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    delta_t = read_delta_t(run_dir)
    output_root = SCRIPT_DIR.parent / "output" / run_dir.name
    output_root.mkdir(parents=True, exist_ok=True)

    save_scalar_series(
        run_dir, output_root, ETA_FIELD, "Eta",
        "eta / passive tracer height", "m", delta_t,
        vmin=ETA_LIMITS[0], vmax=ETA_LIMITS[1], symmetric=False,
    )
    save_scalar_series(
        run_dir, output_root, "ETAN", "ETAN",
        "ETAN", "m", delta_t,
    )
    save_velocity_magnitude(run_dir, output_root, delta_t)
    save_scalar_series(run_dir, output_root, PSI_FIELD, "Psi", "PsiVEL", r"m$^3$ s$^{-1}$", delta_t, symmetric=True)
    save_scalar_series(run_dir, output_root, PHI_FIELD, "Phi", "PhiVEL", r"m$^2$ s$^{-1}$", delta_t, symmetric=True)

    print(f"outputs written to: {output_root}")

if __name__ == "__main__":
    run_snapshots(Path(RUN_DIR))

# %%
