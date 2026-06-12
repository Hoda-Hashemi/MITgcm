from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from shared import alpha_from_data_mypackage, alpha_from_gendata, read_text_value

DIAGNOSTIC_STREAMS = {
    "ETAN": "dynDiag",
    "PhiVEL": "dyn_Aux",
    "PsiVEL": "dyn_Aux",
}


def read_data_value(run_dir: Path, name: str, default: float | None = None) -> float:
    for path in (run_dir / "data", run_dir.parent / "input" / "data"):
        value = read_text_value(path, rf"\b{re.escape(name)}\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return value
    if default is None:
        raise ValueError(f"{name} not found for {run_dir}")
    return default


def read_delta_t(run_dir: Path, default: float = 60.0) -> float:
    return read_data_value(run_dir, "deltaT", default)


def read_package_alpha(run_dir: Path, default: float = 0.0) -> float:
    value = alpha_from_data_mypackage(run_dir)
    if value is not None:
        return value
    value = alpha_from_gendata(run_dir)
    if value is not None:
        return value
    return default


def parse_mds_meta(meta_path: Path) -> dict[str, object]:
    text = meta_path.read_text(encoding="utf-8", errors="ignore")
    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_values = [int(value) for value in re.findall(r"[-+]?\d+", dim_block)]
    dims = []
    for index in range(n_dims):
        start = 3 * index
        n_global, tile_start, tile_end = dim_values[start : start + 3]
        dims.append({"global": n_global, "start": tile_start, "end": tile_end})

    precision = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", text)
    nrecords = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", text)
    fld_list = re.search(r"fldList\s*=\s*\{(.*?)\}", text, re.S)
    fields = re.findall(r"'([^']+)'", fld_list.group(1)) if fld_list else []

    return {
        "n_dims": n_dims,
        "dims": dims,
        "prec": precision.group(1).lower() if precision else "float32",
        "nrecords": int(nrecords.group(1)) if nrecords else 1,
        "fld_list": [field.strip() for field in fields],
    }


def dtype_from_meta(meta: dict[str, object]) -> np.dtype:
    precision = str(meta["prec"])
    return np.dtype(">f8" if precision in {"float64", "real*8"} else ">f4")


def read_mds_field(run_dir: Path, name: str, iteration: int | None = None) -> np.ndarray:
    stem = run_dir / name if iteration is None else run_dir / f"{name}.{iteration:010d}"
    meta_path = Path(f"{stem}.meta")
    data_path = Path(f"{stem}.data")
    record = 0

    if not meta_path.exists() or not data_path.exists():
        stream = DIAGNOSTIC_STREAMS.get(name.strip())
        if stream is None:
            raise FileNotFoundError(f"Missing {name} in {run_dir}")
        stem = run_dir / stream if iteration is None else run_dir / f"{stream}.{iteration:010d}"
        meta_path = Path(f"{stem}.meta")
        data_path = Path(f"{stem}.data")
        if not meta_path.exists() or not data_path.exists():
            raise FileNotFoundError(f"Missing {name} / {stream} in {run_dir}")
        meta = parse_mds_meta(meta_path)
        fields = list(meta.get("fld_list", []))
        if name.strip() not in fields:
            raise FileNotFoundError(f"{name} not found in {stream} for {run_dir}")
        record = fields.index(name.strip())
    else:
        meta = parse_mds_meta(meta_path)

    raw = np.fromfile(data_path, dtype=dtype_from_meta(meta))
    dims = list(meta["dims"])
    n_dims = int(meta["n_dims"])
    nrecords = int(meta["nrecords"])

    if n_dims == 2:
        nx = int(dims[0]["global"])
        ny = int(dims[1]["global"])
        shape = (nrecords, ny, nx) if nrecords > 1 else (ny, nx)
    elif n_dims == 3:
        nx = int(dims[0]["global"])
        ny = int(dims[1]["global"])
        nz = int(dims[2]["global"])
        shape = (nrecords, nz, ny, nx) if nrecords > 1 else (nz, ny, nx)
    else:
        raise ValueError(f"Unsupported nDims={n_dims}: {meta_path}")

    data = raw.reshape(shape)
    return np.squeeze(data[record] if nrecords > 1 else data)


def discover_iterations(run_dir: Path, field: str) -> list[int]:
    pattern = re.compile(rf"^{re.escape(field)}\.(\d{{10}})\.meta$")
    iterations = [
        int(match.group(1))
        for path in run_dir.glob(f"{field}.*.meta")
        if (match := pattern.match(path.name))
    ]
    if iterations:
        return sorted(set(iterations))

    stream = DIAGNOSTIC_STREAMS.get(field.strip())
    if stream is None:
        return []

    pattern = re.compile(rf"^{re.escape(stream)}\.(\d{{10}})\.meta$")
    for path in run_dir.glob(f"{stream}.*.meta"):
        match = pattern.match(path.name)
        if not match:
            continue
        fields = list(parse_mds_meta(path).get("fld_list", []))
        if field.strip() in fields:
            iterations.append(int(match.group(1)))
    return sorted(set(iterations))


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
        raise ValueError(f"Expected 2D or one-level 3D field, got {field.shape}")
    return field


def center_to_shape(field: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    field = to_2d(field)
    ny, nx = shape
    if field.shape == shape:
        return field
    if field.shape == (ny, nx + 1):
        return 0.5 * (field[:, :-1] + field[:, 1:])
    if field.shape == (ny + 1, nx):
        return 0.5 * (field[:-1, :] + field[1:, :])
    raise ValueError(f"Cannot center field of shape {field.shape} to {shape}")


def lon_lat(run_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    return to_2d(read_mds_field(run_dir, "XC")), to_2d(read_mds_field(run_dir, "YC"))


def grid_extent(xc: np.ndarray, yc: np.ndarray) -> tuple[float, float, float, float]:
    lon = xc[0, :]
    lat = yc[:, 0]
    dlon = float(np.median(np.diff(lon))) if lon.size > 1 else 1.0
    dlat = float(np.median(np.diff(lat))) if lat.size > 1 else 1.0
    return (
        float(lon[0] - 0.5 * dlon),
        float(lon[-1] + 0.5 * dlon),
        float(lat[0] - 0.5 * dlat),
        float(lat[-1] + 0.5 * dlat),
    )


def time_text(iteration: int, delta_t: float) -> str:
    seconds = int(round(iteration * delta_t))
    day, seconds = divmod(seconds, 86400)
    hour, seconds = divmod(seconds, 3600)
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


def cell_weights(run_dir: Path, yc: np.ndarray) -> np.ndarray:
    for name in ("RAC", "Area"):
        if (run_dir / f"{name}.meta").exists():
            return to_2d(read_mds_field(run_dir, name))
    return np.cos(np.deg2rad(yc))
