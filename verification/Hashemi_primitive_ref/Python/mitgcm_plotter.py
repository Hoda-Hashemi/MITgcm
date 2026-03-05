#!/usr/bin/env python3
"""
MITgcm primitive-reference plotter.

Style and diagnostics follow the language used in:
verification/Hashemi/Python/plotter.py
"""
#%%
import os
import re
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.path.join(SCRIPT_DIR, "..", "run")
INPUT_BATHY_PATH = os.path.join(SCRIPT_DIR, "..", "input", "bathymetry.bin")
TARGET_IT = 0 #None  # 

#%%

def _parse_mds_meta(meta_path):
    """Parse MITgcm .meta file and return shape/precision info."""
    with open(meta_path, "r", encoding="utf-8") as f:
        txt = f.read()

    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", txt).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", txt, re.S).group(1)
    dim_vals = [int(v) for v in re.findall(r"[-+]?\d+", dim_block)]

    dims = []
    for d in range(n_dims):
        g, s, e = dim_vals[3 * d : 3 * d + 3]
        dims.append({"global": g, "start": s, "end": e, "n": e - s + 1})

    prec_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", txt)
    prec = prec_match.group(1).lower() if prec_match else "float32"

    nrecords_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", txt)
    nrecords = int(nrecords_match.group(1)) if nrecords_match else 1

    it_match = re.search(r"timeStepNumber\s*=\s*\[\s*(\d+)\s*\]", txt)
    it = int(it_match.group(1)) if it_match else 0

    return {"n_dims": n_dims, "dims": dims, "prec": prec, "nrecords": nrecords, "it": it}

def _discover_field_files(run_dir, var_name):
    files = os.listdir(run_dir)
    regex_global = re.compile(rf"^{re.escape(var_name)}\.(\d{{10}})\.meta$")
    regex_tiled = re.compile(rf"^{re.escape(var_name)}\.(\d{{10}})\.(\d{{3}})\.(\d{{3}})\.meta$")

    matches_global = []
    matches_tiled = []
    for fname in files:
        mg = regex_global.match(fname)
        if mg:
            matches_global.append({"it": int(mg.group(1)), "file": fname, "type": "global"})
            continue

        mt = regex_tiled.match(fname)
        if mt:
            matches_tiled.append({"it": int(mt.group(1)), "file": fname, "type": "tiled"})

    all_matches = matches_global + matches_tiled
    if not all_matches:
        raise FileNotFoundError(f"No output files found for variable '{var_name}' in {run_dir}")

    return all_matches

def read_mitgcm_field(run_dir, var_name, it=None, record=0):
    """
    Read MITgcm MDS output in either global or tiled format.
    Returns (field, iteration).
    """
    all_matches = _discover_field_files(run_dir, var_name)
    available_its = sorted({m["it"] for m in all_matches})

    if it is None:
        target_it = available_its[-1]
    else:
        target_it = int(it)
        if target_it not in available_its:
            raise FileNotFoundError(
                f"Timestep {target_it} not found for {var_name}. Available: {available_its}"
            )

    print(f"Reading {var_name} at iteration {target_it}...")
    target_files = [m for m in all_matches if m["it"] == target_it]
    file_type = target_files[0]["type"]

    if file_type == "global":
        meta_path = os.path.join(run_dir, target_files[0]["file"])
        data_path = meta_path.replace(".meta", ".data")
        meta = _parse_mds_meta(meta_path)
        dtype = ">f8" if meta["prec"] in ["float64", "real*8"] else ">f4"

        raw_data = np.fromfile(data_path, dtype=dtype)
        nx = meta["dims"][0]["global"]
        ny = meta["dims"][1]["global"]
        nrecords = meta["nrecords"]

        if meta["n_dims"] == 2:
            if nrecords > 1:
                field = raw_data.reshape((nrecords, ny, nx))[record]
            else:
                field = raw_data.reshape((ny, nx))
        elif meta["n_dims"] == 3:
            nz = meta["dims"][2]["global"]
            if nrecords > 1:
                field = raw_data.reshape((nrecords, nz, ny, nx))[record]
            else:
                field = raw_data.reshape((nz, ny, nx))
        else:
            raise ValueError(f"Unsupported n_dims={meta['n_dims']} for {var_name}")

        return field, target_it

    first_meta = os.path.join(run_dir, sorted([t["file"] for t in target_files])[0])
    meta_global = _parse_mds_meta(first_meta)
    dtype = ">f8" if meta_global["prec"] in ["float64", "real*8"] else ">f4"
    nx_glob = meta_global["dims"][0]["global"]
    ny_glob = meta_global["dims"][1]["global"]

    if meta_global["n_dims"] == 2:
        full_field = np.zeros((ny_glob, nx_glob), dtype=np.float64)
    elif meta_global["n_dims"] == 3:
        nz_glob = meta_global["dims"][2]["global"]
        full_field = np.zeros((nz_glob, ny_glob, nx_glob), dtype=np.float64)
    else:
        raise ValueError(f"Unsupported n_dims={meta_global['n_dims']} for {var_name}")

    for tile in target_files:
        meta_path = os.path.join(run_dir, tile["file"])
        data_path = meta_path.replace(".meta", ".data")
        meta = _parse_mds_meta(meta_path)
        tile_data = np.fromfile(data_path, dtype=dtype)

        x_s = meta["dims"][0]["start"] - 1
        x_e = meta["dims"][0]["end"]
        y_s = meta["dims"][1]["start"] - 1
        y_e = meta["dims"][1]["end"]

        nx_loc = x_e - x_s
        ny_loc = y_e - y_s
        if meta["n_dims"] == 2:
            tile_reshaped = tile_data.reshape((ny_loc, nx_loc))
            full_field[y_s:y_e, x_s:x_e] = tile_reshaped
        else:
            nz_loc = meta["dims"][2]["end"] - meta["dims"][2]["start"] + 1
            tile_reshaped = tile_data.reshape((nz_loc, ny_loc, nx_loc))
            full_field[:, y_s:y_e, x_s:x_e] = tile_reshaped

    return full_field, target_it

def _surface_layer(field):
    return field[0] if field.ndim == 3 else field

def _lon_center_deg(i_lon, n_lon):
    return (i_lon + 0.5) * (360.0 / n_lon)

def _lat_centers_deg(n_lat):
    return -90.0 + (np.arange(n_lat) + 0.5) * (180.0 / n_lat)

def report_ocean_extrema(field_ma, field_name, lat_deg, n_lon, depth_m):
    fld = ma.array(field_ma, copy=False)
    if fld.count() == 0:
        print(f"[{field_name}] no valid ocean points found.")
        return

    idx_min = np.unravel_index(ma.argmin(fld), fld.shape)
    idx_max = np.unravel_index(ma.argmax(fld), fld.shape)

    for tag, (ilat0, ilon0) in (("min", idx_min), ("max", idx_max)):
        val = float(fld[ilat0, ilon0])
        lat0 = float(lat_deg[ilat0])
        lon0 = float(_lon_center_deg(ilon0, n_lon))
        depth_h = max(0.0, float(depth_m[ilat0, ilon0]))
        print(
            f"[{field_name}] {tag} = {val:.8e} at "
            f"iLat={ilat0+1}, iLon={ilon0+1}, "
            f"lat={lat0:.3f} deg, lon={lon0:.3f} deg, "
            f"level=k=1 surface, bathymetry H={depth_h:.3f} m"
        )

def _safe_sym_limit(arr):
    arr_np = ma.filled(ma.array(arr, copy=False), np.nan)
    lim = np.nanmax(np.abs(arr_np))
    return lim if np.isfinite(lim) and lim > 0 else 1e-12

def _safe_max(arr):
    arr_np = ma.filled(ma.array(arr, copy=False), np.nan)
    lim = np.nanmax(arr_np)
    return lim if np.isfinite(lim) and lim > 0 else 1e-12

def _load_depth_and_ocean_mask(n_lat, n_lon):
    if os.path.exists(INPUT_BATHY_PATH):
        bathy = np.fromfile(INPUT_BATHY_PATH, dtype=">f4").reshape((n_lat, n_lon))
        if np.nanmin(bathy) < 0.0:
            ocean_mask = bathy < 0.0
            depth_m = np.maximum(0.0, -bathy)
            print(f"Using bathymetry mask from {INPUT_BATHY_PATH} (ocean where depth < 0).")
        else:
            ocean_mask = bathy > 0.0
            depth_m = np.maximum(0.0, bathy)
            print(f"Using bathymetry mask from {INPUT_BATHY_PATH} (ocean where depth > 0).")
        return depth_m, ocean_mask

    depth, _ = read_mitgcm_field(RUN_DIR, "Depth", it=0)
    depth_2d = _surface_layer(depth)
    ocean_mask = depth_2d > 0.0
    depth_m = np.maximum(0.0, depth_2d)
    print("bathymetry.bin not found; using Depth field from run/ (ocean where Depth > 0).")
    return depth_m, ocean_mask

def _load_optional_hfac_masks(it, n_lat, n_lon):
    tracer_mask = None
    u_mask = None
    v_mask = None
    try:
        hfac_c, _ = read_mitgcm_field(RUN_DIR, "hFacC", it=it)
        hfac_w, _ = read_mitgcm_field(RUN_DIR, "hFacW", it=it)
        hfac_s, _ = read_mitgcm_field(RUN_DIR, "hFacS", it=it)

        hfac_c_2d = _surface_layer(hfac_c)
        hfac_w_2d = _surface_layer(hfac_w)
        hfac_s_2d = _surface_layer(hfac_s)

        if hfac_c_2d.shape == (n_lat, n_lon):
            tracer_mask = hfac_c_2d > 0.0
        if hfac_w_2d.shape == (n_lat, n_lon):
            u_mask = hfac_w_2d > 0.0
        if hfac_s_2d.shape == (n_lat, n_lon):
            v_mask = hfac_s_2d > 0.0
        print("Using hFacC/hFacW/hFacS masks from MITgcm output.")
    except FileNotFoundError:
        pass
    return tracer_mask, u_mask, v_mask

#%%
U_raw, it_found = read_mitgcm_field(RUN_DIR, "U", it=TARGET_IT)
V_raw, _ = read_mitgcm_field(RUN_DIR, "V", it=it_found)
Eta_raw, _ = read_mitgcm_field(RUN_DIR, "Eta", it=it_found)

U_2d = _surface_layer(U_raw)
V_2d = _surface_layer(V_raw)
Eta_2d = _surface_layer(Eta_raw)

n_lat, n_lon = Eta_2d.shape
lat_deg = _lat_centers_deg(n_lat)

depth_m, ocean_mask = _load_depth_and_ocean_mask(n_lat, n_lon)
tracer_mask_hfac, u_mask_hfac, v_mask_hfac = _load_optional_hfac_masks(it_found, n_lat, n_lon)

tracer_mask = tracer_mask_hfac if tracer_mask_hfac is not None else ocean_mask
if u_mask_hfac is None:
    u_mask = tracer_mask & np.roll(tracer_mask, 1, axis=1)
else:
    u_mask = u_mask_hfac

if v_mask_hfac is None:
    v_mask = np.zeros_like(tracer_mask, dtype=bool)
    v_mask[1:, :] = tracer_mask[1:, :] & tracer_mask[:-1, :]
    v_mask[0, :] = tracer_mask[0, :]
else:
    v_mask = v_mask_hfac

Eta_ocean = ma.masked_where(~tracer_mask, Eta_2d)
U_west_face_ocean = ma.masked_where(~u_mask, U_2d)
V_south_face_ocean = ma.masked_where(~v_mask, V_2d)

u_wet = np.where(u_mask, U_2d, np.nan)
v_wet = np.where(v_mask, V_2d, np.nan)

u_center = np.nanmean(np.stack([u_wet, np.roll(u_wet, -1, axis=1)]), axis=0)
v_north_face = np.full_like(v_wet, np.nan)
v_north_face[:-1, :] = v_wet[1:, :]
v_center = np.nanmean(np.stack([v_wet, v_north_face]), axis=0)
v_center[-1, :] = v_wet[-1, :]

speed_center = np.sqrt(u_center * u_center + v_center * v_center)
u_center_ocean = ma.masked_where((~tracer_mask) | np.isnan(u_center), u_center)
v_center_ocean = ma.masked_where((~tracer_mask) | np.isnan(v_center), v_center)
speed_center_ocean = ma.masked_where((~tracer_mask) | np.isnan(speed_center), speed_center)

report_ocean_extrema(Eta_ocean, "eta(surface)", lat_deg, n_lon, depth_m)
report_ocean_extrema(u_center_ocean, "u(center from faces)", lat_deg, n_lon, depth_m)
report_ocean_extrema(v_center_ocean, "v(center from faces)", lat_deg, n_lon, depth_m)
report_ocean_extrema(speed_center_ocean, "speed(center from faces)", lat_deg, n_lon, depth_m)

lim_eta = _safe_sym_limit(Eta_ocean)
lim_u = _safe_sym_limit(U_west_face_ocean)
lim_v = _safe_sym_limit(V_south_face_ocean)
lim_speed = _safe_max(speed_center_ocean)

#%%
# --- Eta ---
plt.figure(figsize=(12,6))
im0 = plt.imshow(
    # Eta_ocean,
    Eta_2d,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=-lim_eta,
    vmax=lim_eta,
    aspect="auto",
)
plt.title(f"Global Ocean Free Surface Eta (MITgcm primitive, it={it_found})")
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.colorbar(im0, label="Eta [m]")
plt.tight_layout()
plt.show()

# --- U field ---
plt.figure(figsize=(12,6))
im1 = plt.imshow(
    # U_west_face_ocean,
    U_2d,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=-lim_u,
    vmax=lim_u,
    aspect="auto",
)
plt.title("C-grid U Field on West Faces")
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.colorbar(im1, label="u [m/s]")
plt.tight_layout()
plt.show()

# --- V field ---
plt.figure(figsize=(12,6))
im2 = plt.imshow(
    # V_south_face_ocean,
    V_2d,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=-lim_v,
    vmax=lim_v,
    aspect="auto",
)
plt.title("C-grid V Field on South Faces")
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.colorbar(im2, label="v [m/s]")
plt.tight_layout()
plt.show()

# --- Speed ---
plt.figure(figsize=(12,6))
im3 = plt.imshow(
    # speed_center_ocean,
    speed_center,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=0.0,
    vmax=lim_speed,
    aspect="auto",
)
plt.title("Velocity Magnitude at Tracer Centers (from C-grid faces)")
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.colorbar(im3, label="speed [m/s]")
plt.tight_layout()
plt.show()
#%%
fig, axs = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)

im0 = axs[0, 0].imshow(
    # Eta_ocean,
    Eta_2d,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=-lim_eta,
    vmax=lim_eta,
    aspect="auto",
)
axs[0, 0].set_title(f"Global Ocean Free Surface Eta (MITgcm primitive, it={it_found})")
axs[0, 0].set_xlabel("Longitude [deg]")
axs[0, 0].set_ylabel("Latitude [deg]")
plt.colorbar(im0, ax=axs[0, 0], label="Eta [m]")

im1 = axs[0, 1].imshow(
    # U_west_face_ocean,
    U_2d,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=-lim_u,
    vmax=lim_u,
    aspect="auto",
)
axs[0, 1].set_title("C-grid U Field on West Faces")
axs[0, 1].set_xlabel("Longitude [deg]")
axs[0, 1].set_ylabel("Latitude [deg]")
plt.colorbar(im1, ax=axs[0, 1], label="u [m/s]")

im2 = axs[1, 0].imshow(
    # V_south_face_ocean,
    V_2d,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=-lim_v,
    vmax=lim_v,
    aspect="auto",
)
axs[1, 0].set_title("C-grid V Field on South Faces")
axs[1, 0].set_xlabel("Longitude [deg]")
axs[1, 0].set_ylabel("Latitude [deg]")
plt.colorbar(im2, ax=axs[1, 0], label="v [m/s]")

im3 = axs[1, 1].imshow(
    # speed_center_ocean,
    speed_center,
    origin="lower",
    extent=[0, 360, -90, 90],
    cmap="seismic",
    vmin=0.0,
    vmax=lim_speed,
    aspect="auto",
)
axs[1, 1].set_title("Velocity Magnitude at Tracer Centers (from C-grid faces)")
axs[1, 1].set_xlabel("Longitude [deg]")
axs[1, 1].set_ylabel("Latitude [deg]")
plt.colorbar(im3, ax=axs[1, 1], label="speed [m/s]")

plt.show()

# %%

#!Surface anomaly generation check if it is fully on ocean:

# The above Python code is performing the following tasks:
import numpy as np
from pathlib import Path
NX, NY = 1440, 720
p = Path("../input/")
eta = np.fromfile(p/"eta_init.bin", dtype=">f4").reshape(NY, NX)
bathy = np.fromfile(p/"bathymetry.bin", dtype=">f4").reshape(NY, NX)
ocean = bathy < 0 if np.nanmin(bathy) < 0 else bathy > 0
land = ~ocean
tol = max(1e-10, 1e-6*np.nanmax(np.abs(eta)))
signal = np.abs(eta) > tol
north = np.zeros_like(land); north[1:,:] = land[:-1,:]
south = np.zeros_like(land); south[:-1,:] = land[1:,:]
coast_ocean = ocean & (north | south | np.roll(land,1,1) | np.roll(land,-1,1))
print("FAIL land leakage" if np.any(signal & land) else "PASS no land leakage")
print("FAIL coastline overlap" if np.any(signal & coast_ocean) else "PASS no coastline overlap")

# %%

