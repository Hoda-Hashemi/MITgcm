#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma

#! PARAMS
nLon, nLat = 1440, 720
dLat = np.pi/nLat  #0.00436 rad in degrees: 0.25
iLat = np.arange(nLat)

R = 6371000.0      

lat_rad = -np.pi/2 + dLat/2 + dLat*iLat
lat_deg = np.rad2deg(lat_rad)
#&Find and merge all 8 files
file_pattern = '../results/solutionPsi/PsiOcean.*.csv'
files = glob.glob(file_pattern)

#& Load all CSVs into one big DataFrame
df_list = [pd.read_csv(f) for f in files]
full_df = pd.concat(df_list, ignore_index=True)

#&Load bathymetry and build wet mask (MITgcm convention: ocean depth < 0)
bathy_path = '../input/bathymetry.bin'
bathy = np.fromfile(bathy_path, dtype='>f4').reshape((nLat, nLon))

ocean_mask = bathy < 0.0

#& psi_Global
psi_Global = np.zeros((nLat, nLon))

for _, row in full_df.iterrows():
    # Subtract 1 because Fortran indices were 1-based
    ilat = int(row['iLat']) - 1
    ilon = int(row['iLon']) - 1
    if 0 <= ilat < nLat and 0 <= ilon < nLon:
        psi_Global[ilat, ilon] = row['psi']

#& only ocean cells
psi_Ocean = ma.masked_where(~ocean_mask, psi_Global)

def _lon_center_deg(i_lon, n_lon):
    return (i_lon + 0.5) * (360.0 / n_lon)

def report_ocean_extrema(field_ma, field_name):
    fld = ma.array(field_ma, copy=False)
    if fld.count() == 0:
        print(f'[{field_name}] no valid ocean points found.')
        return

    idx_min = np.unravel_index(ma.argmin(fld), fld.shape)
    idx_max = np.unravel_index(ma.argmax(fld), fld.shape)

    for tag, (ilat0, ilon0) in (('min', idx_min), ('max', idx_max)):
        val = float(fld[ilat0, ilon0])
        lat0 = float(lat_deg[ilat0])
        lon0 = float(_lon_center_deg(ilon0, nLon))
        depth_h = max(0.0, float(-bathy[ilat0, ilon0]))
        print(
            f'[{field_name}] {tag} = {val:.8e} at '
            f'iLat={ilat0+1}, iLon={ilon0+1}, '
            f'lat={lat0:.3f} deg, lon={lon0:.3f} deg, '
            f'level=k=1 surface, bathymetry H={depth_h:.3f} m'
        )

report_ocean_extrema(psi_Ocean, 'psi')
#%% COMMENT: PLOTTING 2D 
plt.figure(figsize=(12, 6))
limit = np.nanmax(np.abs(psi_Global))
plt.imshow(psi_Global, 
           origin='lower', 
           extent=[0, 360, -90, 90], 
           cmap='seismic',
           vmin=-limit,
           vmax=limit)
plt.colorbar(label='Psi Value')
plt.title('Global MITgcm custom solver output')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
# plt.savefig('../results/Fortran/solution/global_psi_plot.png')
plt.show()

plt.figure(figsize=(14,6))

#! Use imshow with masked array
im = plt.imshow(psi_Ocean, 
                origin='lower', 
                extent=[0, 360, -90, 90], 
                cmap='seismic', 
                vmin=-np.nanmax(np.abs(psi_Ocean)), 
                vmax=np.nanmax(np.abs(psi_Ocean)))

plt.colorbar(im, label='Psi Value')

# Add labels and title
plt.xlabel('Longitude [deg]')
plt.ylabel('Latitude [deg]')
plt.title('Global Ocean Streamfunction (MITgcm)')
plt.show()

#%% COMMENT: VELOCITY FROM solutionVel CSV (2D + streamlines + extrema)
#! C-GRID AWARE INTERPOLATED VELOCITIES
vel_pattern = '../results/solutionVel/VelOcean.*.csv'
vel_files = sorted(glob.glob(vel_pattern))

vel_df_list = [pd.read_csv(f) for f in vel_files]
vel_df = pd.concat(vel_df_list, ignore_index=True)
print(f'Loaded {len(vel_files)} velocity CSV tiles from {vel_pattern}')

dLon_deg = 360.0 / nLon
dLat_deg = 180.0 / nLat
lon_center_deg = (np.arange(nLon) + 0.5) * dLon_deg
lat_center_deg = -90.0 + (np.arange(nLat) + 0.5) * dLat_deg

# C-grid storage from solver CSV:
# - u_west_face(j,i) is U at west face of tracer cell (i,j)
# - v_south_face(j,i) is V at south face of tracer cell (i,j)
u_west_face = np.full((nLat, nLon), np.nan, dtype=np.float64)
v_south_face = np.full((nLat, nLon), np.nan, dtype=np.float64)

ilat_idx = vel_df['iLat'].to_numpy(dtype=np.int64) - 1
ilon_idx = vel_df['iLon'].to_numpy(dtype=np.int64) - 1
valid = (
    (ilat_idx >= 0) & (ilat_idx < nLat) &
    (ilon_idx >= 0) & (ilon_idx < nLon)
)

u_vals = vel_df['uVel'].to_numpy(dtype=np.float64)
v_vals = vel_df['vVel'].to_numpy(dtype=np.float64)
u_west_face[ilat_idx[valid], ilon_idx[valid]] = u_vals[valid]
v_south_face[ilat_idx[valid], ilon_idx[valid]] = v_vals[valid]

# Build face wet masks from tracer wet mask.
# West face is between (j,i-1) and (j,i), periodic in longitude.
mask_w_face = ocean_mask & np.roll(ocean_mask, 1, axis=1)
# South face is between (j-1,i) and (j,i), non-periodic in latitude.
mask_s_face = np.zeros_like(ocean_mask, dtype=bool)
mask_s_face[1:, :] = ocean_mask[1:, :] & ocean_mask[:-1, :]
mask_s_face[0, :] = ocean_mask[0, :]

u_wet = np.where(mask_w_face, u_west_face, np.nan)
v_wet = np.where(mask_s_face, v_south_face, np.nan)

# Interpolate face velocities to tracer centers for center-based speed/streamplot.
u_center = np.nanmean(np.stack([u_wet, np.roll(u_wet, -1, axis=1)]), axis=0)
v_north_face = np.full_like(v_wet, np.nan)
v_north_face[:-1, :] = v_wet[1:, :]
v_center = np.nanmean(np.stack([v_wet, v_north_face]), axis=0)

# At the northernmost row, keep one-sided value from available south face.
v_center[-1, :] = v_wet[-1, :]

speed_center = np.sqrt(u_center*u_center + v_center*v_center)

u_center_ocean = ma.masked_where((~ocean_mask) | np.isnan(u_center), u_center)
v_center_ocean = ma.masked_where((~ocean_mask) | np.isnan(v_center), v_center)
speed_center_ocean = ma.masked_where((~ocean_mask) | np.isnan(speed_center), speed_center)

report_ocean_extrema(u_center_ocean, 'u(center from faces)')
report_ocean_extrema(v_center_ocean, 'v(center from faces)')
report_ocean_extrema(speed_center_ocean, 'speed(center from faces)')
#%%
#ALTERED: SPEED
plt.figure(figsize=(12, 6))
im_vel = plt.imshow(
    speed_center_ocean,
    origin='lower',
    extent=[0, 360, -90, 90],
    cmap='seismic',
    aspect='auto'
)
plt.colorbar(im_vel, label='Speed (m/s)')
plt.title('Velocity Magnitude at Tracer Centers (from C-grid faces)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

#ALTERED: U and V streamlines
u_stream = np.nan_to_num(ma.filled(u_center_ocean, np.nan), nan=0.0)
v_stream = np.nan_to_num(ma.filled(v_center_ocean, np.nan), nan=0.0)
speed_stream = np.sqrt(u_stream*u_stream + v_stream*v_stream)

plt.figure(figsize=(12, 6))
strm = plt.streamplot(
    lon_center_deg, lat_center_deg,
    u_stream, v_stream,
    color=speed_stream, cmap='seismic',
    density=5, linewidth=0.6, arrowsize=0.4
)
plt.colorbar(strm.lines, label='Speed (m/s)')
plt.title('Surface Streamlines at Tracer Centers (from C-grid faces)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

#ALTERED: U and V 
# Optional: show raw staggered fields at their native face locations.
u_west_face_ocean = ma.masked_where((~mask_w_face) | np.isnan(u_west_face), u_west_face)
v_south_face_ocean = ma.masked_where((~mask_s_face) | np.isnan(v_south_face), v_south_face)

plt.figure(figsize=(12, 6))
im_u = plt.imshow(
    u_west_face_ocean,
    origin='lower',
    extent=[-0.5*dLon_deg, 360.0 - 0.5*dLon_deg, -90.0, 90.0],
    cmap='seismic',
    aspect='auto'
)
plt.colorbar(im_u, label='u at west faces (m/s)')
plt.title('C-grid U Field on West Faces')
plt.xlabel('Longitude (west-face locations)')
plt.ylabel('Latitude (center locations)')
plt.show()

plt.figure(figsize=(12, 6))
im_v = plt.imshow(
    v_south_face_ocean,
    origin='lower',
    extent=[0.0, 360.0, -90.0 - 0.5*dLat_deg, 90.0 - 0.5*dLat_deg],
    cmap='seismic',
    aspect='auto'
)
plt.colorbar(im_v, label='v at south faces (m/s)')
plt.title('C-grid V Field on South Faces')
plt.xlabel('Longitude (center locations)')
plt.ylabel('Latitude (south-face locations)')
plt.show()
#%% ALTERED: norm
speed_norm = speed_center_ocean / np.nanmax(speed_center_ocean) 

plt.figure(figsize=(12, 6))
im_vel = plt.imshow(
    speed_norm,
    origin='lower',
    extent=[0, 360, -90, 90],
    cmap='seismic',
    aspect='auto'
)
plt.colorbar(im_vel, label='Speed (m/s)')
plt.title('Velocity Magnitude at Tracer Centers (from C-grid faces)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# %%
#!RESIDUALS
import pandas as pd
import re

df = pd.read_csv('../results/residuals.csv')

plt.figure(figsize=(10, 6))

# Use log scale for the Y-axis (Residuals)
plt.semilogy(df['iter'], df['residual'], color='blue', linewidth=1.5, label='CG2D Residual')

# Add labels and styling
plt.title('MITgcm CG2D Solver Convergence', fontsize=14)
plt.xlabel('Iteration Number', fontsize=12)
plt.ylabel('L2-Norm Residual (Log Scale)', fontsize=12)
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()

# Save the plot
# plt.savefig('../results/Plots/convergence_plot.png', dpi=300)
plt.show()

r = df['residual'].to_numpy()
it = df['iter'].to_numpy()

print('min residual ', min(r))

rel = r / r[0]
ratio = r[1:] / r[:-1]
imin = df['residual'].idxmin()
plt.figure(figsize=(10,6))
plt.semilogy(it, r, label='abs residual')
plt.semilogy(it, rel, '--', label='relative residual')
plt.scatter(df.loc[imin,'iter'], df.loc[imin,'residual'], 
            color='r', 
            label=f'Min: {min(r):.2e}')

plt.axhline(1e-8, color='k', ls=':', label='target abs=1e-8')
plt.grid(True, which='both', alpha=0.4); plt.legend()
# plt.savefig('../results/Plots/relative_residual.png', dpi=300)

plt.show()

plt.figure(figsize=(10,4))
plt.plot(it[1:], ratio)
plt.axhline(1.0, color='r', ls='--')
plt.ylabel('r_k / r_{k-1}'); plt.xlabel('iter')
plt.title('Per-iteration reduction factor')
# plt.savefig('../results/Plots/Per_iterationa_residual.png', dpi=300)
plt.grid(True, alpha=0.3); plt.show()

# %%
#%%
#!HERE! they are zero since mitgcm is unaware of the forcings that will be used to solve for the velocities inside the momentum equations
import os
import re

# MITgcm velocity notes (Arakawa C-grid):
# - U is the zonal velocity at the west/east cell faces (u-point).
# - V is the meridional velocity at the south/north cell faces (v-point).
# - W is the vertical velocity at vertical faces (w-point), diagnosed from
#   continuity and boundary conditions (surface/bottom constraints).
# So U,V,W are staggered by design; for a speed map at tracer centers we
# do a simple face-to-center average of U and V.

run_dir = '../../Hashemi_primitive_ref/run/'
target_it = None  # None -> latest available timestep in run_dir
import os
import re

# Get all unique iteration numbers for variable U
files = os.listdir(run_dir)
iterations = sorted(list(set([re.search(r'\.(\d{10})\.', f).group(1) 
                             for f in files if f.startswith('U.') and '.meta' in f])))

print(f"Available Timesteps: {iterations}")

def _parse_mds_meta(meta_path):
    txt = open(meta_path, 'r').read()

    n_dims = int(re.search(r'nDims\s*=\s*\[\s*(\d+)\s*\]', txt).group(1))
    dim_block = re.search(r'dimList\s*=\s*\[(.*?)\];', txt, re.S).group(1)
    dim_vals = [int(v) for v in re.findall(r'[-+]?\d+', dim_block)]
    dims = []
    for d in range(n_dims):
        g, s, e = dim_vals[3*d:3*d+3]
        dims.append({'global': g, 'start': s, 'end': e, 'n': e - s + 1})

    prec = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", txt).group(1).lower()
    nrecords = int(re.search(r'nrecords\s*=\s*\[\s*(\d+)\s*\]', txt).group(1))
    it = int(re.search(r'timeStepNumber\s*=\s*\[\s*(\d+)\s*\]', txt).group(1))

    return {'n_dims': n_dims, 'dims': dims, 'prec': prec, 'nrecords': nrecords, 'it': it}

def _prec_to_dtype(prec):
    if prec in ('float32', 'real*4'):
        return '>f4'
    if prec in ('float64', 'real*8'):
        return '>f8'
    raise ValueError(f'Unsupported dataprec: {prec}')

def _read_tile_data(data_path, meta, record=0):
    dtype = _prec_to_dtype(meta['prec'])
    raw = np.fromfile(data_path, dtype=dtype)

    local_shape = [d['n'] for d in meta['dims']]  # [nx, ny] or [nx, ny, nz]
    n_local = int(np.prod(local_shape))
    nrecords = meta['nrecords']
    expected = n_local * nrecords
    if raw.size != expected:
        raise ValueError(
            f'Size mismatch for {data_path}: got {raw.size}, expected {expected}'
        )

    offset = record * n_local
    buf = raw[offset:offset + n_local]

    if meta['n_dims'] == 2:
        nx_loc, ny_loc = local_shape
        # Fortran write order: i-fastest then j. Convert to [j,i].
        tile = buf.reshape((nx_loc, ny_loc), order='F').T
        return tile

    if meta['n_dims'] == 3:
        nx_loc, ny_loc, nz_loc = local_shape
        # Convert [i,j,k] -> [k,j,i]
        tile = buf.reshape((nx_loc, ny_loc, nz_loc), order='F').transpose(2, 1, 0)
        return tile

    raise ValueError(f'Unsupported nDims: {meta["n_dims"]}')

def read_mitgcm_tiled_field(run_dir, var_name, it=None, record=0):
    # Expected files: U.0000000000.003.001.meta (and .data)
    pattern = re.compile(rf'^{re.escape(var_name)}\.(\d+)\.(\d+)\.(\d+)\.meta$')
    meta_files = []
    for fn in os.listdir(run_dir):
        m = pattern.match(fn)
        if m:
            meta_files.append((int(m.group(1)), int(m.group(2)), int(m.group(3)),
                               os.path.join(run_dir, fn)))

    if not meta_files:
        raise FileNotFoundError(f'No {var_name} tiled meta files found in {run_dir}')

    if it is None:
        target = max(x[0] for x in meta_files)
    else:
        target = int(it)

    meta_files = [x for x in meta_files if x[0] == target]
    if not meta_files:
        raise FileNotFoundError(f'No {var_name} files for timestep {target}')

    # Build global array from first tile metadata
    m0 = _parse_mds_meta(meta_files[0][3])
    if m0['n_dims'] == 2:
        nx_glob = m0['dims'][0]['global']
        ny_glob = m0['dims'][1]['global']
        field = np.full((ny_glob, nx_glob), np.nan, dtype=np.float64)
    elif m0['n_dims'] == 3:
        nx_glob = m0['dims'][0]['global']
        ny_glob = m0['dims'][1]['global']
        nz_glob = m0['dims'][2]['global']
        field = np.full((nz_glob, ny_glob, nx_glob), np.nan, dtype=np.float64)
    else:
        raise ValueError(f'Unsupported nDims in meta: {m0["n_dims"]}')

    for _, _, _, meta_path in sorted(meta_files):
        meta = _parse_mds_meta(meta_path)
        data_path = meta_path.replace('.meta', '.data')
        tile = _read_tile_data(data_path, meta, record=record)

        x0 = meta['dims'][0]['start'] - 1
        x1 = meta['dims'][0]['end']
        y0 = meta['dims'][1]['start'] - 1
        y1 = meta['dims'][1]['end']

        if meta['n_dims'] == 2:
            field[y0:y1, x0:x1] = tile
        else:
            z0 = meta['dims'][2]['start'] - 1
            z1 = meta['dims'][2]['end']
            field[z0:z1, y0:y1, x0:x1] = tile

    return field, target

# Read U, V, W from MITgcm run directory
U_raw, it_found = read_mitgcm_tiled_field(run_dir, 'U', it=target_it)
V_raw, _ = read_mitgcm_tiled_field(run_dir, 'V', it=it_found)
W_raw, _ = read_mitgcm_tiled_field(run_dir, 'W', it=it_found)

# Handle either 2D (Nr=1) or 3D output; plot top level if 3D
U_2d = U_raw[0] if U_raw.ndim == 3 else U_raw
V_2d = V_raw[0] if V_raw.ndim == 3 else V_raw
W_2d = W_raw[0] if W_raw.ndim == 3 else W_raw

# Mask land for cleaner maps
U_2d = ma.masked_where(~ocean_mask, U_2d)
V_2d = ma.masked_where(~ocean_mask, V_2d)
W_2d = ma.masked_where(~ocean_mask, W_2d)

# Simple face->center estimate of horizontal speed
Ue = np.roll(U_2d, -1, axis=1)                 # east neighbor (periodic in lon)
Vn = np.vstack([V_2d[1:, :], V_2d[-1:, :]])    # north neighbor (clamped at north edge)
Uc = 0.5 * (U_2d + Ue)
Vc = 0.5 * (V_2d + Vn)
UV_speed = np.sqrt(Uc**2 + Vc**2)

fig, axs = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)

im0 = axs[0, 0].imshow(U_2d, origin='lower', extent=[0, 360, -90, 90], cmap='RdBu_r')
axs[0, 0].set_title(f'U at it={it_found}')
axs[0, 0].set_xlabel('Longitude')
axs[0, 0].set_ylabel('Latitude')
plt.colorbar(im0, ax=axs[0, 0], label='m/s')

im1 = axs[0, 1].imshow(V_2d, origin='lower', extent=[0, 360, -90, 90], cmap='RdBu_r')
axs[0, 1].set_title(f'V at it={it_found}')
axs[0, 1].set_xlabel('Longitude')
axs[0, 1].set_ylabel('Latitude')
plt.colorbar(im1, ax=axs[0, 1], label='m/s')

im2 = axs[1, 0].imshow(W_2d, origin='lower', extent=[0, 360, -90, 90], cmap='PuOr')
axs[1, 0].set_title(f'W at it={it_found}')
axs[1, 0].set_xlabel('Longitude')
axs[1, 0].set_ylabel('Latitude')
plt.colorbar(im2, ax=axs[1, 0], label='m/s')

im3 = axs[1, 1].imshow(UV_speed, origin='lower', extent=[0, 360, -90, 90], cmap='viridis')
axs[1, 1].set_title('Estimated |U,V| at tracer centers')
axs[1, 1].set_xlabel('Longitude')
axs[1, 1].set_ylabel('Latitude')
plt.colorbar(im3, ax=axs[1, 1], label='m/s')

plt.show()

# %%

#! Velocity representations (line plots, contour map, streamlines, quiver)
# Options shown here:
# 1) Zonal-mean profiles vs latitude
# 2) Meridional slice at one longitude
# 3) 2D contour of |U,V|
# 4) Streamlines (dark theme)
# 5) Coarse quiver vectors

lon_deg_plot = np.linspace(0.0, 360.0, nLon, endpoint=False)
lat_deg_plot = np.rad2deg(lat_rad)

Uc_plot = np.ma.filled(Uc, np.nan)
Vc_plot = np.ma.filled(Vc, np.nan)
W_plot = np.ma.filled(W_2d, np.nan)
spd_plot = np.ma.filled(UV_speed, np.nan)

# 1) Zonal means
U_zm = np.nanmean(Uc_plot, axis=1)
V_zm = np.nanmean(Vc_plot, axis=1)
W_zm = np.nanmean(W_plot, axis=1)
S_zm = np.nanmean(spd_plot, axis=1)

plt.figure(figsize=(10, 6))
plt.plot(lat_deg_plot, U_zm, lw=1.8, label='U zonal mean')
plt.plot(lat_deg_plot, V_zm, lw=1.8, label='V zonal mean')
plt.plot(lat_deg_plot, W_zm, lw=1.8, label='W zonal mean')
plt.plot(lat_deg_plot, S_zm, lw=2.0, color='k', label='|U,V| zonal mean')
plt.xlabel('Latitude (deg)')
plt.ylabel('Velocity (m/s)')
plt.title('MITgcm Velocity Profiles')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# 2) Meridional slices at lon=0
ilon0 = 0
plt.figure(figsize=(10, 6))
plt.plot(lat_deg_plot, Uc_plot[:, ilon0], lw=1.5, label='U at 0°')
plt.plot(lat_deg_plot, Vc_plot[:, ilon0], lw=1.5, label='V at 0°')
plt.plot(lat_deg_plot, W_plot[:, ilon0], lw=1.5, label='W at 0°')
plt.plot(lat_deg_plot, spd_plot[:, ilon0], lw=2.0, color='k', label='|U,V| at 0°')
plt.xlabel('Latitude (deg)')
plt.ylabel('Velocity (m/s)')
plt.title('Meridional Slice at 0° Longitude')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# 3) 2D imshow map of horizontal speed
plt.figure(figsize=(13, 5))
im_spd = plt.imshow(
    spd_plot, origin='lower', extent=[0, 360, -90, 90],
    cmap='viridis', aspect='auto'
)
plt.colorbar(im_spd, label='|U,V| (m/s)')
plt.xlabel('Longitude (deg)')
plt.ylabel('Latitude (deg)')
plt.title('Horizontal Speed Magnitude (imshow)')
plt.show()

# 4) Streamlines (dark theme)
u_dark = np.nan_to_num(Uc_plot, nan=0.0, posinf=0.0, neginf=0.0)
v_dark = np.nan_to_num(Vc_plot, nan=0.0, posinf=0.0, neginf=0.0)
s_dark = np.nan_to_num(spd_plot, nan=0.0, posinf=0.0, neginf=0.0)

fig, ax = plt.subplots(figsize=(14, 6), facecolor='#0b1020')
ax.set_facecolor('#0b1020')
strm = ax.streamplot(
    lon_deg_plot, lat_deg_plot, u_dark, v_dark,
    color=s_dark, cmap='turbo', density=2.0, linewidth=1.0, arrowsize=0.8
)
cbar = plt.colorbar(strm.lines, ax=ax, pad=0.02)
cbar.set_label('|U,V| (m/s)', color='white')
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(cbar.ax.get_yticklabels(), color='white')
ax.set_title('MITgcm Surface Streamlines', color='white')
ax.set_xlabel('Longitude (deg)', color='white')
ax.set_ylabel('Latitude (deg)', color='white')
ax.tick_params(colors='white')
ax.grid(color='white', alpha=0.12)
plt.show()

# 5) 2D imshow map of vertical velocity W (instead of quiver)
plt.figure(figsize=(13, 5))
wlim = np.nanpercentile(np.abs(W_plot), 99)
im_w = plt.imshow(
    W_plot, origin='lower', extent=[0, 360, -90, 90],
    cmap='RdBu_r', vmin=-wlim, vmax=wlim, aspect='auto'
)
plt.colorbar(im_w, label='W (m/s)')
plt.xlabel('Longitude (deg)')
plt.ylabel('Latitude (deg)')
plt.title('Vertical Velocity W (imshow)')
plt.show()
# %%
#! Relative vorticity zeta from U,V (on sphere)
# zeta = (1/(R cos(phi))) * ( dV/dlambda - d(U cos(phi))/dphi )

dphi = lat_rad[1] - lat_rad[0]
dlambda = 2.0 * np.pi / nLon
cos_phi = np.cos(lat_rad)[:, None]
cos_safe = np.where(np.abs(cos_phi) < 1e-10, np.nan, cos_phi)

dV_dlambda = (np.roll(Vc_plot, -1, axis=1) - np.roll(Vc_plot, 1, axis=1)) / (2.0 * dlambda)
Ucos = Uc_plot * cos_phi
dUcos_dphi = np.gradient(Ucos, dphi, axis=0, edge_order=2)

zeta = (dV_dlambda - dUcos_dphi) / (R * cos_safe)
zeta = np.where(ocean_mask, zeta, np.nan)

plt.figure(figsize=(13, 5))
zlim = np.nanpercentile(np.abs(zeta), 99)
imz = plt.contourf(lon_deg_plot, lat_deg_plot, zeta, levels=40, cmap='RdBu_r',
                   vmin=-zlim, vmax=zlim)
plt.colorbar(imz, label=r'$\zeta$ (s$^{-1}$)')
plt.xlabel('Longitude (deg)')
plt.ylabel('Latitude (deg)')
plt.title('Relative Vorticity (zeta)')
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(lat_deg_plot, np.nanmean(zeta, axis=1), color='purple', lw=1.8)
plt.xlabel('Latitude (deg)')
plt.ylabel(r'Zonal-mean $\zeta$ (s$^{-1}$)')
plt.title('Zonal-Mean Relative Vorticity')
plt.grid(True, alpha=0.3)
plt.show()

