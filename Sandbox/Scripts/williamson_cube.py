from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

import williamson_fields as wf

CS32_NC = 32
GRID_FILES = tuple(f"tile{face:03d}.mitgrid" for face in range(1, 7))


def _find_mitgcm_root(anchor: Path) -> Path:
    for path in (anchor.resolve(), *anchor.resolve().parents):
        if (path / "verification" / "aim.5l_cs" / "input" / "tile001.mitgrid").exists():
            return path
    raise FileNotFoundError("could not locate MITgcm root with CS32 mitgrid files")


def resolve_cs32_grid_dir(base_dir: Path) -> Path:
    env = os.environ.get("WILLIAMSON_CS_GRID_DIR")
    if env:
        path = Path(env).expanduser().resolve()
        if (path / "tile001.mitgrid").exists():
            return path
        raise FileNotFoundError(f"WILLIAMSON_CS_GRID_DIR is missing tile001.mitgrid: {path}")

    local = base_dir.resolve()
    if (local / "tile001.mitgrid").exists():
        return local
    return _find_mitgcm_root(base_dir) / "verification" / "aim.5l_cs" / "input"


def _deg_to_rad(lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon = np.deg2rad(np.mod(lon_deg, 360.0))
    lat = np.deg2rad(lat_deg)
    return lon, lat


def _sph_to_cart(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    cos_lat = np.cos(lat)
    return np.stack(
        (cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)),
        axis=-1,
    )


def _cart_to_lonlat(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xyz = _normalise(xyz)
    lon = np.mod(np.arctan2(xyz[..., 1], xyz[..., 0]), 2.0 * np.pi)
    lat = np.arcsin(np.clip(xyz[..., 2], -1.0, 1.0))
    return lon, lat


def _normalise(xyz: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(xyz, axis=-1, keepdims=True)
    return xyz / np.maximum(norm, 1.0e-300)


def _axis_tangent(points: np.ndarray, axis: int) -> np.ndarray:
    tangent = np.empty_like(points)
    if axis == 0:
        tangent[1:-1, ...] = points[2:, ...] - points[:-2, ...]
        tangent[0, ...] = points[1, ...] - points[0, ...]
        tangent[-1, ...] = points[-1, ...] - points[-2, ...]
    else:
        tangent[:, 1:-1, ...] = points[:, 2:, ...] - points[:, :-2, ...]
        tangent[:, 0, ...] = points[:, 1, ...] - points[:, 0, ...]
        tangent[:, -1, ...] = points[:, -1, ...] - points[:, -2, ...]
    tangent = tangent - np.sum(tangent * points, axis=-1, keepdims=True) * points
    return _normalise(tangent)


def _east_north_unit(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    east = np.stack((-np.sin(lon), np.cos(lon), np.zeros_like(lon)), axis=-1)
    north = np.stack(
        (-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)),
        axis=-1,
    )
    return east, north


class CubeGrid:
    def __init__(self, grid_dir: Path, nc: int = CS32_NC) -> None:
        self.grid_dir = grid_dir
        self.nc = nc
        self.xc = np.empty((6, nc, nc), dtype=np.float64)
        self.yc = np.empty((6, nc, nc), dtype=np.float64)
        self.xg = np.empty((6, nc + 1, nc + 1), dtype=np.float64)
        self.yg = np.empty((6, nc + 1, nc + 1), dtype=np.float64)
        self.dxc = np.empty((6, nc, nc), dtype=np.float64)
        self.dyc = np.empty((6, nc, nc), dtype=np.float64)
        self.dxf = np.empty((6, nc, nc), dtype=np.float64)
        self.dyf = np.empty((6, nc, nc), dtype=np.float64)

        for face, name in enumerate(GRID_FILES):
            path = grid_dir / name
            raw = np.fromfile(path, dtype=">f8")
            expected = (nc + 1) * (nc + 1) * 16
            if raw.size != expected:
                raise ValueError(f"{path} has {raw.size} values, expected {expected}")
            data = raw.reshape((nc + 1, nc + 1, 16), order="F")
            self.xc[face] = data[:nc, :nc, 0]
            self.yc[face] = data[:nc, :nc, 1]
            self.dxf[face] = data[:nc, :nc, 2]
            self.dyf[face] = data[:nc, :nc, 3]
            self.xg[face] = data[:, :, 5]
            self.yg[face] = data[:, :, 6]
            self.dxc[face] = data[:nc, :nc, 10]
            self.dyc[face] = data[:nc, :nc, 11]

        self.lon_c, self.lat_c = _deg_to_rad(self.xc, self.yc)
        self.lon_g, self.lat_g = _deg_to_rad(self.xg, self.yg)
        self.rg = _sph_to_cart(self.lon_g, self.lat_g)

        u_points_all = _normalise(self.rg[:, :-1, :-1, :] + self.rg[:, :-1, 1:, :])
        v_points_all = _normalise(self.rg[:, :-1, :-1, :] + self.rg[:, 1:, :-1, :])
        self.r_u = u_points_all
        self.r_v = v_points_all
        self.lon_u, self.lat_u = _cart_to_lonlat(self.r_u)
        self.lon_v, self.lat_v = _cart_to_lonlat(self.r_v)
        self.e_u = _axis_tangent(u_points_all, axis=1)
        self.e_v = _axis_tangent(v_points_all, axis=2)


def copy_cs32_grid_files(run_dir: Path) -> None:
    grid_dir = resolve_cs32_grid_dir(run_dir)
    for name in GRID_FILES:
        src = grid_dir / name
        dst = run_dir / name
        if not dst.exists() or dst.stat().st_size != src.stat().st_size:
            dst.write_bytes(src.read_bytes())


def write_compact(base_dir: Path, name: str, faces: np.ndarray, dtype: str = ">f4") -> None:
    faces = np.asarray(faces, dtype=np.float64)
    if faces.shape != (6, CS32_NC, CS32_NC):
        raise ValueError(f"{name}: expected shape {(6, CS32_NC, CS32_NC)}, got {faces.shape}")
    compact = np.transpose(faces, (1, 0, 2))
    np.asfortranarray(compact).astype(dtype).tofile(base_dir / name)
    print(
        f"Wrote {name}: cube=CS{CS32_NC}, shape={faces.shape}, "
        f"min={np.nanmin(faces):.6e}, max={np.nanmax(faces):.6e}"
    )


def write_corner_faces(base_dir: Path, stem: str, faces: np.ndarray) -> None:
    faces = np.asarray(faces, dtype=np.float64)
    if faces.shape != (6, CS32_NC + 1, CS32_NC + 1):
        raise ValueError(f"{stem}: expected corner shape {(6, CS32_NC + 1, CS32_NC + 1)}, got {faces.shape}")
    for face in range(6):
        name = f"{stem}.face{face + 1:03d}.bin"
        np.asfortranarray(faces[face]).astype(">f4").tofile(base_dir / name)
        print(f"Wrote {name}: shape={faces[face].shape}")


def interp_latlon(field: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (wf.NY, wf.NX):
        raise ValueError(f"expected lat-lon field {(wf.NY, wf.NX)}, got {field.shape}")

    x = np.mod(lon, 2.0 * np.pi) / wf.DLON_RAD - 0.5
    y = (lat + 0.5 * np.pi) / wf.DLAT_RAD - 0.5
    y = np.clip(y, 0.0, wf.NY - 1.0)

    i0 = np.floor(x).astype(np.int64)
    j0 = np.floor(y).astype(np.int64)
    fx = x - i0
    fy = y - j0
    i0 = np.mod(i0, wf.NX)
    i1 = np.mod(i0 + 1, wf.NX)
    j1 = np.minimum(j0 + 1, wf.NY - 1)

    return (
        (1.0 - fx) * (1.0 - fy) * field[j0, i0]
        + fx * (1.0 - fy) * field[j0, i1]
        + (1.0 - fx) * fy * field[j1, i0]
        + fx * fy * field[j1, i1]
    )


def project_east_north_to_cube(
    lon: np.ndarray,
    lat: np.ndarray,
    eastward: np.ndarray,
    northward: np.ndarray,
    target_axis: np.ndarray,
) -> np.ndarray:
    east, north = _east_north_unit(lon, lat)
    vector = eastward[..., None] * east + northward[..., None] * north
    return np.sum(vector * target_axis, axis=-1)


def velocity_faces_from_functions(
    grid: CubeGrid,
    u_func,
    v_func,
    *,
    alpha_rad: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    u_e = u_func(grid.lon_u, grid.lat_u, alpha_rad)
    v_n = v_func(grid.lon_u, grid.lat_u, alpha_rad)
    u_cube = project_east_north_to_cube(grid.lon_u, grid.lat_u, u_e, v_n, grid.e_u)

    u_e = u_func(grid.lon_v, grid.lat_v, alpha_rad)
    v_n = v_func(grid.lon_v, grid.lat_v, alpha_rad)
    v_cube = project_east_north_to_cube(grid.lon_v, grid.lat_v, u_e, v_n, grid.e_v)
    return u_cube, v_cube


def velocity_faces_from_latlon(
    grid: CubeGrid,
    u_east: np.ndarray,
    v_north: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    u_e = interp_latlon(u_east, grid.lon_u, grid.lat_u)
    v_n = interp_latlon(v_north, grid.lon_u, grid.lat_u)
    u_cube = project_east_north_to_cube(grid.lon_u, grid.lat_u, u_e, v_n, grid.e_u)

    u_e = interp_latlon(u_east, grid.lon_v, grid.lat_v)
    v_n = interp_latlon(v_north, grid.lon_v, grid.lat_v)
    v_cube = project_east_north_to_cube(grid.lon_v, grid.lat_v, u_e, v_n, grid.e_v)
    return u_cube, v_cube


def generate_tc2(base_dir: Path, alpha_rad: float) -> None:
    copy_cs32_grid_files(base_dir)
    grid = CubeGrid(resolve_cs32_grid_dir(base_dir))

    eta = wf.tc2_eta(grid.lon_c, grid.lat_c, alpha_rad)
    bathy = np.full_like(eta, -wf.H0_SOLID)
    u, v = velocity_faces_from_functions(grid, wf.tc2_u, wf.tc2_v, alpha_rad=alpha_rad)

    mu_c = -np.cos(grid.lon_c) * np.cos(grid.lat_c) * np.sin(alpha_rad) + np.sin(grid.lat_c) * np.cos(alpha_rad)
    mu_g = -np.cos(grid.lon_g) * np.cos(grid.lat_g) * np.sin(alpha_rad) + np.sin(grid.lat_g) * np.cos(alpha_rad)
    write_compact(base_dir, "fCoriC.bin", 2.0 * wf.OMEGA * mu_c)
    write_compact(base_dir, "fCorCs.bin", 2.0 * wf.OMEGA * np.sqrt(np.maximum(0.0, 1.0 - mu_c * mu_c)))
    write_corner_faces(base_dir, "fCoriG", 2.0 * wf.OMEGA * mu_g)

    write_compact(base_dir, "bathymetry_flat4000.bin", bathy)
    write_compact(base_dir, "eta_init.bin", eta)
    write_compact(base_dir, "u_init.bin", u)
    write_compact(base_dir, "v_init.bin", v)
    print(f"Cube grid: CS{CS32_NC}; alpha={alpha_rad:.16g}; max speed={np.nanmax(np.hypot(u, v)):.6e} m/s")


def generate_tc5(base_dir: Path) -> None:
    copy_cs32_grid_files(base_dir)
    grid = CubeGrid(resolve_cs32_grid_dir(base_dir))
    bathy = wf.tc5_bathymetry(grid.lon_c, grid.lat_c)
    eta = wf.tc5_eta(grid.lon_c, grid.lat_c)
    u, v = velocity_faces_from_functions(
        grid,
        lambda lon, lat, _alpha: wf.tc5_u(lon, lat),
        lambda lon, lat, _alpha: wf.tc5_v(lon, lat),
    )
    write_compact(base_dir, "bathymetry_mountain_tc5.bin", bathy)
    write_compact(base_dir, "eta_init.bin", eta)
    write_compact(base_dir, "u_init.bin", u)
    write_compact(base_dir, "v_init.bin", v)
    print(f"Cube grid: CS{CS32_NC}; depth min={np.nanmin(eta - bathy):.6e} m")


def generate_tc7(base_dir: Path) -> None:
    copy_cs32_grid_files(base_dir)
    grid = CubeGrid(resolve_cs32_grid_dir(base_dir))
    eta_ll, u_ll, v_ll, bathy_ll = wf.load_tc7_fields(base_dir)
    eta_ll, u_ll, v_ll = wf.preprocess_tc7_analysis_fields(eta_ll, u_ll, v_ll)
    eta = interp_latlon(eta_ll, grid.lon_c, grid.lat_c)
    bathy = interp_latlon(bathy_ll, grid.lon_c, grid.lat_c)
    u, v = velocity_faces_from_latlon(grid, u_ll, v_ll)
    write_compact(base_dir, "bathymetry_tc7.bin", bathy)
    write_compact(base_dir, "eta_init.bin", eta)
    write_compact(base_dir, "u_init.bin", u)
    write_compact(base_dir, "v_init.bin", v)
    print(f"Cube grid: CS{CS32_NC}; max speed={np.nanmax(np.hypot(u, v)):.6e} m/s")


def add_scripts_path(anchor: Path) -> None:
    scripts_dir = _find_mitgcm_root(anchor) / "Sandbox" / "Scripts"
    scripts_text = str(scripts_dir)
    if scripts_text not in sys.path:
        sys.path.insert(0, scripts_text)
