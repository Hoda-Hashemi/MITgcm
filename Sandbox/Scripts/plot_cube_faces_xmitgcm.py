#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from mitgcm_io import discover_iterations, read_delta_t, read_mds_field
from snapshot_plots import CS_FACE_SIZE, compact_to_faces, faces_to_net
from shared import resolve_color_limits


DEFAULT_FIELD = "Eta"
DEFAULT_CMAP = "seismic"


def find_repo_root(path: Path) -> Path:
    for candidate in (path.resolve(), *path.resolve().parents):
        if (candidate / "verification").exists() and (candidate / "Sandbox").exists():
            return candidate
    raise FileNotFoundError(f"could not infer MITgcm root from {path}")


def clean_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def default_output_dir(run_dir: Path) -> Path:
    repo_root = find_repo_root(run_dir)
    return repo_root / "Sandbox" / "output" / "xmitgcm_cube_faces" / run_dir.parent.name / run_dir.name


def latest_iteration(run_dir: Path, field: str) -> int:
    iterations = discover_iterations(run_dir, field)
    if not iterations:
        raise FileNotFoundError(f"no iterations found for {field} in {run_dir}")
    return iterations[-1]


def simple_extra_dims(field: str) -> list[str]:
    field_upper = field.upper()
    if field_upper in {"U", "UVEL", "UVELMASS"}:
        return ["j", "i_g"]
    if field_upper in {"V", "VVEL", "VVELMASS"}:
        return ["j_g", "i"]
    return ["j", "i"]


def xarray_dataarray_to_faces(data_array) -> np.ndarray:
    data_array = data_array.squeeze(drop=True)
    if "face" in data_array.dims:
        y_dim = next((dim for dim in ("j", "j_g") if dim in data_array.dims), None)
        x_dim = next((dim for dim in ("i", "i_g") if dim in data_array.dims), None)
        if y_dim is None or x_dim is None:
            raise ValueError(f"cannot identify x/y dims for {data_array.name}: {data_array.dims}")
        return np.asarray(data_array.transpose("face", y_dim, x_dim).values, dtype=np.float64)
    return compact_to_faces(np.asarray(data_array.values, dtype=np.float64))


def read_faces_xmitgcm(run_dir: Path, field: str, iteration: int, face_size: int) -> np.ndarray:
    faces, _lon_faces, _lat_faces = read_faces_grid_xmitgcm(run_dir, field, iteration, face_size)
    return faces


def read_faces_grid_xmitgcm(
    run_dir: Path,
    field: str,
    iteration: int,
    face_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        from xmitgcm import open_mdsdataset
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "xmitgcm is not installed in this Python environment.\n"
            "Install it first, for example:\n"
            "  python -m pip install xmitgcm xarray dask xgcm"
        ) from exc

    extra_variables = {
        field: {
            "dims": simple_extra_dims(field),
            "attrs": {"long_name": field, "units": ""},
        }
    }
    dataset = open_mdsdataset(
        str(run_dir),
        grid_dir=str(run_dir),
        iters=[iteration],
        prefix=[field],
        read_grid=True,
        geometry="cs",
        nx=face_size,
        delta_t=read_delta_t(run_dir, default=1.0),
        default_dtype=np.dtype("float64"),
        ignore_unknown_vars=False,
        extra_variables=extra_variables,
    )
    if field in dataset:
        data_array = dataset[field]
    elif field.upper() in dataset:
        data_array = dataset[field.upper()]
    else:
        names = ", ".join(sorted(dataset.data_vars))
        raise KeyError(f"{field} not present in xmitgcm dataset. Available variables: {names}")
    return (
        xarray_dataarray_to_faces(data_array),
        xarray_dataarray_to_faces(dataset["XC"]),
        xarray_dataarray_to_faces(dataset["YC"]),
    )


def read_faces_raw(run_dir: Path, field: str, iteration: int) -> np.ndarray:
    return compact_to_faces(read_mds_field(run_dir, field, iteration))


def interpolate_faces_to_latlon(
    run_dir: Path,
    faces: np.ndarray,
    *,
    lon_faces: np.ndarray | None = None,
    lat_faces: np.ndarray | None = None,
    lon_res: float,
    lat_res: float,
    method: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        from scipy.interpolate import griddata
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "scipy is required for --layout latlon.\n"
            "Install it first, for example:\n"
            "  python -m pip install scipy"
        ) from exc

    if lon_faces is None or lat_faces is None:
        lon_faces = compact_to_faces(read_mds_field(run_dir, "XC"))
        lat_faces = compact_to_faces(read_mds_field(run_dir, "YC"))
    lon = np.mod(lon_faces.ravel(), 360.0)
    lat = lat_faces.ravel()
    values = np.asarray(faces, dtype=np.float64).ravel()
    valid = np.isfinite(values) & np.isfinite(lon) & np.isfinite(lat)

    lon_valid = lon[valid]
    lat_valid = lat[valid]
    values_valid = values[valid]
    points = np.column_stack(
        [
            np.concatenate([lon_valid - 360.0, lon_valid, lon_valid + 360.0]),
            np.tile(lat_valid, 3),
        ]
    )
    tiled_values = np.tile(values_valid, 3)

    lon_target = np.arange(0.0, 360.0 + 0.5 * lon_res, lon_res)
    lat_target = np.arange(-90.0, 90.0 + 0.5 * lat_res, lat_res)
    lon_grid, lat_grid = np.meshgrid(lon_target, lat_target)
    interpolated = griddata(points, tiled_values, (lon_grid, lat_grid), method=method)
    if np.isnan(interpolated).any() and method != "nearest":
        fill = griddata(points, tiled_values, (lon_grid, lat_grid), method="nearest")
        interpolated = np.where(np.isfinite(interpolated), interpolated, fill)
    return lon_grid, lat_grid, interpolated


def plot_latlon(
    run_dir: Path,
    faces: np.ndarray,
    out_path: Path,
    *,
    lon_faces: np.ndarray | None = None,
    lat_faces: np.ndarray | None = None,
    title: str,
    units: str,
    cmap_name: str,
    center_zero: bool,
    vmin: float | None,
    vmax: float | None,
    lon_res: float,
    lat_res: float,
    interp_method: str,
    dpi: int,
) -> None:
    lon_grid, lat_grid, field = interpolate_faces_to_latlon(
        run_dir,
        faces,
        lon_faces=lon_faces,
        lat_faces=lat_faces,
        lon_res=lon_res,
        lat_res=lat_res,
        method=interp_method,
    )
    vmin, vmax = resolve_color_limits(field, vmin, vmax, symmetric=center_zero)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if center_zero else None
    fig, ax = plt.subplots(figsize=(10.5, 5.6), constrained_layout=True)
    image = ax.pcolormesh(
        lon_grid,
        lat_grid,
        field,
        shading="auto",
        cmap=cmap_name,
        norm=norm,
        vmin=None if norm else vmin,
        vmax=None if norm else vmax,
    )
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xlabel(r"longitude $\lambda$ [deg]")
    ax.set_ylabel(r"latitude $\theta$ [deg]")
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.set_title(title)
    fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.10).set_label(units)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_panels(
    faces: np.ndarray,
    out_path: Path,
    *,
    title: str,
    units: str,
    layout: str,
    cmap_name: str,
    center_zero: bool,
    vmin: float | None,
    vmax: float | None,
    dpi: int,
) -> None:
    faces = np.asarray(faces, dtype=np.float64)
    if faces.shape[0] != 6:
        raise ValueError(f"expected six cube faces; got {faces.shape}")

    vmin, vmax = resolve_color_limits(faces, vmin, vmax, symmetric=center_zero)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if center_zero else None
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("white")

    if layout == "net":
        fig, ax = plt.subplots(figsize=(8.8, 6.4), constrained_layout=True)
        image = ax.imshow(
            np.ma.masked_invalid(faces_to_net(faces)),
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
            vmin=None if norm else vmin,
            vmax=None if norm else vmax,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title)
        fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.08).set_label(units)
    else:
        if layout == "strip":
            nrows, ncols, figsize = 1, 6, (15.0, 3.4)
        else:
            nrows, ncols, figsize = 2, 3, (10.0, 6.8)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
        axes_flat = np.ravel(axes)
        image = None
        for index, ax in enumerate(axes_flat):
            image = ax.imshow(
                np.ma.masked_invalid(faces[index]),
                origin="lower",
                interpolation="nearest",
                cmap=cmap,
                norm=norm,
                vmin=None if norm else vmin,
                vmax=None if norm else vmax,
            )
            ax.set_title(f"face {index + 1}")
            ax.set_xticks([])
            ax.set_yticks([])
        fig.suptitle(title)
        fig.colorbar(image, ax=axes_flat.tolist(), orientation="horizontal", pad=0.05).set_label(units)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot MITgcm CS32 compact cubed-sphere output as six separate faces."
    )
    parser.add_argument("run_dir", type=Path, help="MITgcm run directory containing MDS .meta/.data files")
    parser.add_argument("--field", default=DEFAULT_FIELD, help=f"MDS field/prefix to plot, default: {DEFAULT_FIELD}")
    parser.add_argument("--iter", type=int, default=None, help="MITgcm iteration; default: latest for field")
    parser.add_argument("--reader", choices=("xmitgcm", "raw"), default="xmitgcm")
    parser.add_argument("--layout", choices=("panels", "strip", "net", "latlon"), default="panels")
    parser.add_argument("--units", default="", help="Colorbar units label")
    parser.add_argument("--cmap", default=DEFAULT_CMAP)
    parser.add_argument("--center-zero", action="store_true", help="Use a zero-centered diverging norm")
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--lon-res", type=float, default=1.0, help="latlon output longitude spacing in degrees")
    parser.add_argument("--lat-res", type=float, default=1.0, help="latlon output latitude spacing in degrees")
    parser.add_argument("--interp-method", choices=("linear", "nearest", "cubic"), default="linear")
    parser.add_argument("--face-size", type=int, default=CS_FACE_SIZE)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    iteration = args.iter if args.iter is not None else latest_iteration(run_dir, args.field)

    lon_faces = None
    lat_faces = None
    if args.reader == "xmitgcm":
        faces, lon_faces, lat_faces = read_faces_grid_xmitgcm(run_dir, args.field, iteration, args.face_size)
    else:
        faces = read_faces_raw(run_dir, args.field, iteration)

    output_dir = args.output_dir.resolve() if args.output_dir else default_output_dir(run_dir)
    out_name = f"{clean_name(args.field)}_iter_{iteration:010d}_{args.reader}_{args.layout}.png"
    out_path = output_dir / out_name
    title = f"{args.field} iter {iteration:010d} ({args.reader}, {args.layout})"
    if args.layout == "latlon":
        plot_latlon(
            run_dir,
            faces,
            out_path,
            lon_faces=lon_faces,
            lat_faces=lat_faces,
            title=title,
            units=args.units,
            cmap_name=args.cmap,
            center_zero=args.center_zero,
            vmin=args.vmin,
            vmax=args.vmax,
            lon_res=args.lon_res,
            lat_res=args.lat_res,
            interp_method=args.interp_method,
            dpi=args.dpi,
        )
    else:
        plot_panels(
            faces,
            out_path,
            title=title,
            units=args.units,
            layout=args.layout,
            cmap_name=args.cmap,
            center_zero=args.center_zero,
            vmin=args.vmin,
            vmax=args.vmax,
            dpi=args.dpi,
        )
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
