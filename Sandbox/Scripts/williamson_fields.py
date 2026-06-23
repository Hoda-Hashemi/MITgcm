from __future__ import annotations

from pathlib import Path

import numpy as np

NX = 1440
NY = 720
R_EARTH = 6.371e6
OMEGA = 7.292e-5
G = 9.81
DAY = 86400.0
DLON_RAD = 2.0 * np.pi / NX
DLAT_RAD = np.pi / NY

U0_SOLID = 2.0 * np.pi * R_EARTH / (12.0 * DAY)
GH0_SOLID = 2.94e4
H0_SOLID = GH0_SOLID / G

TC3_U0 = U0_SOLID
TC3_LAT0 = -np.pi / 6.0
TC3_LAT1 = np.pi / 2.0
TC3_LAT_GRID = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 20001)

TC5_U0 = 20.0
TC5_H0 = 5400.0
TC5_MOUNTAIN_HEIGHT = 2000.0
TC5_MOUNTAIN_RADIUS = np.pi / 9.0
TC5_MOUNTAIN_LON = 1.5 * np.pi
TC5_MOUNTAIN_LAT = np.pi / 6.0

TC6_WAVE_NUMBER = 4
TC6_K = 7.848e-6
TC6_WAVE_OMEGA = 7.848e-6
TC6_H0 = 8000.0

TC7_REFERENCE_DEPTH = 8000.0
TC7_LONGITUDE_WAVENUMBER_CUTOFF = 42
TC7_SHAPIRO_SMOOTHING_PASSES = 8


def add_scripts_path(anchor: Path) -> None:
    """Let run-local gendata copies import this module from Sandbox/Scripts."""
    import sys

    case_dir = anchor.resolve().parent
    if case_dir.name in {"input", "run_standard", "run_analysis"} or case_dir.name.startswith("run_"):
        sandbox_dir = case_dir.parent.parent
    else:
        sandbox_dir = case_dir.parent
    scripts_dir = sandbox_dir / "Scripts"
    scripts_text = str(scripts_dir)
    if scripts_text not in sys.path:
        sys.path.insert(0, scripts_text)


def write_field(base_dir: Path, name: str, field: np.ndarray) -> None:
    field.astype(">f4").tofile(base_dir / name)
    print(f"Wrote {name}: shape={field.shape}, min={np.nanmin(field):.6e}, max={np.nanmax(field):.6e}")


def grid_cell_centers(nx: int = NX, ny: int = NY) -> tuple[np.ndarray, np.ndarray]:
    lon = (np.arange(nx, dtype=np.float64) + 0.5) * (2.0 * np.pi / nx)
    lat = -0.5 * np.pi + (np.arange(ny, dtype=np.float64) + 0.5) * (np.pi / ny)
    return np.meshgrid(lon, lat, indexing="xy")


def grid_u_faces(nx: int = NX, ny: int = NY) -> tuple[np.ndarray, np.ndarray]:
    lon = np.arange(nx, dtype=np.float64) * (2.0 * np.pi / nx)
    lat = -0.5 * np.pi + (np.arange(ny, dtype=np.float64) + 0.5) * (np.pi / ny)
    return np.meshgrid(lon, lat, indexing="xy")


def grid_v_faces(nx: int = NX, ny: int = NY) -> tuple[np.ndarray, np.ndarray]:
    lon = (np.arange(nx, dtype=np.float64) + 0.5) * (2.0 * np.pi / nx)
    lat = -0.5 * np.pi + np.arange(ny, dtype=np.float64) * (np.pi / ny)
    return np.meshgrid(lon, lat, indexing="xy")


def rotated_latitude(lon_rad: np.ndarray, lat_rad: np.ndarray, alpha_rad: float) -> np.ndarray:
    mu = (
        -np.cos(lon_rad) * np.cos(lat_rad) * np.sin(alpha_rad)
        + np.sin(lat_rad) * np.cos(alpha_rad)
    )
    return np.arcsin(np.clip(mu, -1.0, 1.0))


def tc2_eta(lon_rad: np.ndarray, lat_rad: np.ndarray, alpha_rad: float = 0.0) -> np.ndarray:
    mu = (
        -np.cos(lon_rad) * np.cos(lat_rad) * np.sin(alpha_rad)
        + np.sin(lat_rad) * np.cos(alpha_rad)
    )
    gh = GH0_SOLID - (R_EARTH * OMEGA * U0_SOLID + 0.5 * U0_SOLID**2) * mu**2
    return gh / G - H0_SOLID


def tc2_u(lon_rad: np.ndarray, lat_rad: np.ndarray, alpha_rad: float = 0.0) -> np.ndarray:
    return U0_SOLID * (
        np.cos(lat_rad) * np.cos(alpha_rad)
        + np.cos(lon_rad) * np.sin(lat_rad) * np.sin(alpha_rad)
    )


def tc2_v(lon_rad: np.ndarray, lat_rad: np.ndarray, alpha_rad: float = 0.0) -> np.ndarray:
    return -U0_SOLID * np.sin(lon_rad) * np.sin(alpha_rad)


def _compact_raw(lat_rad: np.ndarray) -> np.ndarray:
    lat = np.asarray(lat_rad, dtype=np.float64)
    out = np.zeros_like(lat)
    inside = (lat > TC3_LAT0) & (lat < TC3_LAT1)
    denom = (lat[inside] - TC3_LAT0) * (lat[inside] - TC3_LAT1)
    out[inside] = np.exp(1.0 / denom)
    midpoint = 0.5 * (TC3_LAT0 + TC3_LAT1)
    norm = np.exp(1.0 / ((midpoint - TC3_LAT0) * (midpoint - TC3_LAT1)))
    return out / norm


def tc3_speed_profile(lat_rad: np.ndarray) -> np.ndarray:
    return TC3_U0 * _compact_raw(lat_rad)


def _cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def _tc3_height_grid() -> np.ndarray:
    lat = TC3_LAT_GRID
    u = tc3_speed_profile(lat)
    safe_lat = np.clip(lat, -0.5 * np.pi + 1.0e-12, 0.5 * np.pi - 1.0e-12)
    integrand = -R_EARTH * u * (2.0 * OMEGA * np.sin(lat) + u * np.tan(safe_lat) / R_EARTH) / G
    primitive = _cumulative_trapezoid(integrand, lat)
    weights = np.cos(lat)
    primitive -= np.trapz(primitive * weights, lat) / np.trapz(weights, lat)
    return H0_SOLID + primitive


def _tc3_streamfunction_grid() -> np.ndarray:
    return _cumulative_trapezoid(R_EARTH * tc3_speed_profile(TC3_LAT_GRID), TC3_LAT_GRID)


TC3_HEIGHT_GRID = _tc3_height_grid()
TC3_STREAMFUNCTION_GRID = _tc3_streamfunction_grid()

TC4_U0 = 20.0
TC4_GH0 = 1.0e5
TC4_H0 = TC4_GH0 / G
TC4_F0 = 2.0 * OMEGA * np.sin(np.pi / 4.0)
TC4_PSI0 = -0.03 * TC4_GH0 / TC4_F0
TC4_SIGMA = 12.74244**2
TC4_LON0 = -7.0 * np.pi / 12.0
TC4_LAT0 = np.pi / 4.0
TC4_LAT_GRID = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 20001)


def _interp_lat(values: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    flat = np.interp(np.ravel(lat_rad), TC3_LAT_GRID, values)
    return flat.reshape(np.shape(lat_rad))


def tc3_eta(lon_rad: np.ndarray, lat_rad: np.ndarray, alpha_rad: float = 0.0) -> np.ndarray:
    lat_prime = rotated_latitude(lon_rad, lat_rad, alpha_rad)
    return _interp_lat(TC3_HEIGHT_GRID, lat_prime) - H0_SOLID


def _tc3_streamfunction(lon_rad: np.ndarray, lat_rad: np.ndarray, alpha_rad: float) -> np.ndarray:
    lat_prime = rotated_latitude(lon_rad, lat_rad, alpha_rad)
    return _interp_lat(TC3_STREAMFUNCTION_GRID, lat_prime)


def tc3_velocity_fields(alpha_rad: float = 0.0, nx: int = NX, ny: int = NY) -> tuple[np.ndarray, np.ndarray]:
    dlon = 2.0 * np.pi / nx
    dlat = np.pi / ny

    lon_u = np.arange(nx, dtype=np.float64) * dlon
    lat_s = -0.5 * np.pi + np.arange(ny, dtype=np.float64) * dlat
    lat_n = lat_s + dlat
    psi_n = _tc3_streamfunction(lon_u[None, :], lat_n[:, None], alpha_rad)
    psi_s = _tc3_streamfunction(lon_u[None, :], lat_s[:, None], alpha_rad)
    u = (psi_n - psi_s) / (R_EARTH * dlat)

    lon_w = np.arange(nx, dtype=np.float64) * dlon
    lon_e = lon_w + dlon
    lat_v = -0.5 * np.pi + np.arange(ny, dtype=np.float64) * dlat
    psi_w = _tc3_streamfunction(lon_w[None, :], lat_v[:, None], alpha_rad)
    psi_e = _tc3_streamfunction(lon_e[None, :], lat_v[:, None], alpha_rad)
    dx = R_EARTH * np.cos(lat_v) * dlon
    v = np.zeros((ny, nx), dtype=np.float64)
    valid = np.abs(dx) > 1.0e-14
    v[valid, :] = (psi_w[valid, :] - psi_e[valid, :]) / dx[valid, None]
    v[0, :] = 0.0
    v[-1, :] = 0.0
    return u, v


def tc4_background_u(lat_rad: np.ndarray, u0: float = TC4_U0) -> np.ndarray:
    return u0 * np.sin(2.0 * lat_rad) ** 14


def _tc4_background_geopotential_grid(u0: float = TC4_U0) -> np.ndarray:
    lat = TC4_LAT_GRID
    safe_lat = np.clip(lat, -0.5 * np.pi + 1.0e-12, 0.5 * np.pi - 1.0e-12)
    u = tc4_background_u(lat, u0=u0)
    f = 2.0 * OMEGA * np.sin(lat)
    integrand = (R_EARTH * f + u * np.tan(safe_lat)) * u
    return TC4_GH0 - _cumulative_trapezoid(integrand, lat)


TC4_BACKGROUND_GH_GRID = _tc4_background_geopotential_grid()


def _interp_tc4_background_gh(lat_rad: np.ndarray, u0: float = TC4_U0) -> np.ndarray:
    if u0 == TC4_U0:
        values = TC4_BACKGROUND_GH_GRID
    else:
        values = _tc4_background_geopotential_grid(u0=u0)
    flat = np.interp(np.ravel(lat_rad), TC4_LAT_GRID, values)
    return flat.reshape(np.shape(lat_rad))


def _tc4_psi_and_derivatives(
    lon_rad: np.ndarray,
    lat_rad: np.ndarray,
    time_s: float = 0.0,
    u0: float = TC4_U0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    phase = lon_rad - (u0 / R_EARTH) * time_s - TC4_LON0
    c = (
        np.sin(TC4_LAT0) * np.sin(lat_rad)
        + np.cos(TC4_LAT0) * np.cos(lat_rad) * np.cos(phase)
    )
    denom = np.maximum(1.0 + c, 1.0e-14)
    psi = TC4_PSI0 * np.exp(-TC4_SIGMA * (1.0 - c) / denom)
    scale = 2.0 * TC4_SIGMA / (denom * denom)
    c_lon = -np.cos(TC4_LAT0) * np.cos(lat_rad) * np.sin(phase)
    c_lat = (
        np.sin(TC4_LAT0) * np.cos(lat_rad)
        - np.cos(TC4_LAT0) * np.sin(lat_rad) * np.cos(phase)
    )
    psi_lon = psi * scale * c_lon
    psi_lat = psi * scale * c_lat
    return psi, psi_lon, psi_lat


def tc4_height(
    lon_rad: np.ndarray,
    lat_rad: np.ndarray,
    time_s: float = 0.0,
    u0: float = TC4_U0,
) -> np.ndarray:
    psi, _, _ = _tc4_psi_and_derivatives(lon_rad, lat_rad, time_s=time_s, u0=u0)
    gh = _interp_tc4_background_gh(lat_rad, u0=u0) + TC4_F0 * psi
    return gh / G


def tc4_eta(
    lon_rad: np.ndarray,
    lat_rad: np.ndarray,
    time_s: float = 0.0,
    u0: float = TC4_U0,
) -> np.ndarray:
    return tc4_height(lon_rad, lat_rad, time_s=time_s, u0=u0) - TC4_H0


def tc4_u(
    lon_rad: np.ndarray,
    lat_rad: np.ndarray,
    time_s: float = 0.0,
    u0: float = TC4_U0,
) -> np.ndarray:
    _, _, psi_lat = _tc4_psi_and_derivatives(lon_rad, lat_rad, time_s=time_s, u0=u0)
    return tc4_background_u(lat_rad, u0=u0) - psi_lat / R_EARTH


def tc4_v(
    lon_rad: np.ndarray,
    lat_rad: np.ndarray,
    time_s: float = 0.0,
    u0: float = TC4_U0,
) -> np.ndarray:
    _, psi_lon, _ = _tc4_psi_and_derivatives(lon_rad, lat_rad, time_s=time_s, u0=u0)
    cos_lat = np.cos(lat_rad)
    out = np.zeros_like(psi_lon)
    valid = np.abs(cos_lat) > 1.0e-14
    out[valid] = psi_lon[valid] / (R_EARTH * cos_lat[valid])
    return out


def tc5_mountain(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    dlon = (lon_rad - TC5_MOUNTAIN_LON + np.pi) % (2.0 * np.pi) - np.pi
    radius = np.sqrt(dlon * dlon + (lat_rad - TC5_MOUNTAIN_LAT) ** 2)
    mountain = np.zeros_like(radius)
    inside = radius < TC5_MOUNTAIN_RADIUS
    mountain[inside] = TC5_MOUNTAIN_HEIGHT * (1.0 - radius[inside] / TC5_MOUNTAIN_RADIUS)
    return mountain


def tc5_eta(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    geopotential_drop = (R_EARTH * OMEGA * TC5_U0 + 0.5 * TC5_U0**2) * np.sin(lat_rad) ** 2 / G
    return -geopotential_drop


def tc5_bathymetry(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    return -(TC5_H0 - tc5_mountain(lon_rad, lat_rad))


def tc5_u(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    return TC5_U0 * np.cos(lat_rad)


def tc5_v(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    return np.zeros_like(lon_rad)


def tc6_u(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    r = TC6_WAVE_NUMBER
    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)
    return (
        R_EARTH * TC6_WAVE_OMEGA * cos_lat
        + R_EARTH
        * TC6_K
        * cos_lat ** (r - 1)
        * (r * sin_lat * sin_lat - cos_lat * cos_lat)
        * np.cos(r * lon_rad)
    )


def tc6_v(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    r = TC6_WAVE_NUMBER
    return -R_EARTH * TC6_K * r * np.cos(lat_rad) ** (r - 1) * np.sin(lat_rad) * np.sin(r * lon_rad)


def tc6_height(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    r = TC6_WAVE_NUMBER
    cos_lat = np.cos(lat_rad)
    cos2 = cos_lat * cos_lat
    safe_cos2 = np.maximum(cos2, 1.0e-30)
    a = (
        0.5 * TC6_WAVE_OMEGA * (2.0 * OMEGA + TC6_WAVE_OMEGA) * cos2
        + 0.25
        * TC6_K**2
        * cos_lat ** (2 * r)
        * ((r + 1.0) * cos2 + (2.0 * r * r - r - 2.0) - 2.0 * r * r / safe_cos2)
    )
    b = (
        2.0
        * (OMEGA + TC6_WAVE_OMEGA)
        * TC6_K
        / ((r + 1.0) * (r + 2.0))
        * cos_lat**r
        * ((r * r + 2.0 * r + 2.0) - (r + 1.0) ** 2 * cos2)
    )
    c = 0.25 * TC6_K**2 * cos_lat ** (2 * r) * ((r + 1.0) * cos2 - (r + 2.0))
    return TC6_H0 + R_EARTH**2 * (a + b * np.cos(r * lon_rad) + c * np.cos(2 * r * lon_rad)) / G


def tc6_eta(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    return tc6_height(lon_rad, lat_rad) - TC6_H0


def flat_bathymetry(depth: float, nx: int = NX, ny: int = NY) -> np.ndarray:
    return np.full((ny, nx), -float(depth), dtype=np.float64)


def longitude_spectral_filter(field: np.ndarray, max_wavenumber: int) -> np.ndarray:
    field = np.asarray(field, dtype=np.float64)
    coeffs = np.fft.rfft(field, axis=1)
    cutoff = max(0, int(max_wavenumber))
    if cutoff + 1 < coeffs.shape[1]:
        coeffs[:, cutoff + 1 :] = 0.0
    return np.fft.irfft(coeffs, n=field.shape[1], axis=1)


def shapiro_smooth(field: np.ndarray, passes: int) -> np.ndarray:
    out = np.asarray(field, dtype=np.float64).copy()
    for _ in range(max(0, int(passes))):
        x_smooth = 0.25 * np.roll(out, 1, axis=1) + 0.5 * out + 0.25 * np.roll(out, -1, axis=1)
        padded = np.pad(x_smooth, ((1, 1), (0, 0)), mode="edge")
        out = 0.25 * padded[:-2, :] + 0.5 * x_smooth + 0.25 * padded[2:, :]
    return out


def preprocess_tc7_analysis_fields(
    eta: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    *,
    max_wavenumber: int = TC7_LONGITUDE_WAVENUMBER_CUTOFF,
    smoothing_passes: int = TC7_SHAPIRO_SMOOTHING_PASSES,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the large-scale TC7 preprocessing used by the local MITgcm setup."""
    eta_mean = float(np.mean(eta))
    eta_f = longitude_spectral_filter(eta - eta_mean, max_wavenumber) + eta_mean
    u_f = longitude_spectral_filter(u, max_wavenumber)
    v_f = longitude_spectral_filter(v, max_wavenumber)
    if smoothing_passes:
        eta_f = shapiro_smooth(eta_f, smoothing_passes)
        u_f = shapiro_smooth(u_f, smoothing_passes)
        v_f = shapiro_smooth(v_f, smoothing_passes)
    v_f[0, :] = 0.0
    v_f[-1, :] = 0.0
    return eta_f, u_f, v_f


def resolve_tc7_raw_file(base_dir: Path) -> Path:
    candidates = [
        base_dir / "raw" / "tc7_initial_conditions.npz",
        base_dir.parent / "input" / "raw" / "tc7_initial_conditions.npz",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        "TC7 requires input/raw/tc7_initial_conditions.npz with eta_m, u_m_s, and v_m_s arrays."
    )


def load_tc7_fields(base_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    path = resolve_tc7_raw_file(base_dir)
    with np.load(path) as data:
        missing = {"eta_m", "u_m_s", "v_m_s"} - set(data.files)
        if missing:
            raise ValueError(f"{path} is missing required arrays: {', '.join(sorted(missing))}")
        eta = np.array(data["eta_m"], dtype=np.float64)
        u = np.array(data["u_m_s"], dtype=np.float64)
        v = np.array(data["v_m_s"], dtype=np.float64)
        bathy = np.array(data["bathymetry_m"], dtype=np.float64) if "bathymetry_m" in data.files else flat_bathymetry(TC7_REFERENCE_DEPTH)

    for name, field in {"eta_m": eta, "u_m_s": u, "v_m_s": v, "bathymetry_m": bathy}.items():
        if field.shape != (NY, NX):
            raise ValueError(f"{path}:{name} must have shape {(NY, NX)}, got {field.shape}")
    v[0, :] = 0.0
    v[-1, :] = 0.0
    return eta, u, v, bathy
