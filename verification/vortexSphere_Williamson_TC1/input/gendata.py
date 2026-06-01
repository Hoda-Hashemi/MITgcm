import numpy as np

# =========================
# USER SETTINGS
# =========================
Nx = 1440
Ny = 720
a  = 6.37122e6              # Earth radius [m]
alpha = 0.5*np.pi          # over-the-pole
Hflat = 1000.0             # flat layer depth [m]
h0 = 1000.0                # cosine-bell max height [m]
Rbell = a/3.0              # bell radius [m]
u0 = 2.0*np.pi*a/(12.0*24.0*3600.0)   # ~38.6 m/s
lonc = 1.5*np.pi           # 3pi/2
latc = 0.0                 # 0

# output dtype / endianness for MITgcm
dtype = ">f8"

# =========================
# GRID
# =========================
dlon = 360.0 / Nx
dlat = 180.0 / Ny

lon_deg = (np.arange(Nx) + 0.5) * dlon
lat_deg = -90.0 + (np.arange(Ny) + 0.5) * dlat

lon = np.deg2rad(lon_deg)[None, :]              # (1,Nx)
lat = np.deg2rad(lat_deg)[:, None]              # (Ny,1)

# U points: longitude faces, latitude centers
lon_u_deg = np.arange(Nx) * dlon
lat_u_deg = lat_deg
lon_u = np.deg2rad(lon_u_deg)[None, :]
lat_u = np.deg2rad(lat_u_deg)[:, None]

# V points: longitude centers, latitude faces
lon_v_deg = lon_deg
lat_v_deg = -90.0 + np.arange(Ny) * dlat
lon_v = np.deg2rad(lon_v_deg)[None, :]
lat_v = np.deg2rad(lat_v_deg)[:, None]

# =========================
# INITIAL HEIGHT (cosine bell)
# =========================
arg = np.sin(latc)*np.sin(lat) + np.cos(latc)*np.cos(lat)*np.cos(lon-lonc)
arg = np.clip(arg, -1.0, 1.0)
r = a * np.arccos(arg)

Eta = np.zeros((Ny, Nx))
mask = r < Rbell
Eta[mask] = 0.5*h0*(1.0 + np.cos(np.pi*r[mask]/Rbell))

# =========================
# PRESCRIBED SOLID-BODY VELOCITY
# Williamson TC1 formulas
# u = u0 (cos(theta) cos(alpha) + sin(theta) cos(lambda) sin(alpha))
# v = -u0 sin(lambda) sin(alpha)
# =========================
U = u0 * (np.cos(lat_u)*np.cos(alpha) + np.sin(lat_u)*np.cos(lon_u)*np.sin(alpha))
V = -u0 * np.sin(lon_v) * np.sin(alpha)

# =========================
# FLAT DEPTH
# MITgcm uses positive bathymetry depth in meters
# =========================
Depth = Hflat * np.ones((Ny, Nx))

# =========================
# WRITE BINARIES
# =========================
Depth.astype(dtype).tofile("bathymetry_tc1.bin")
Eta.astype(dtype).tofile("eta_tc1.bin")
U.astype(dtype).tofile("u_init_tc1.bin")
V.astype(dtype).tofile("v_init_tc1.bin")

print("Wrote: bathymetry_tc1.bin, eta_tc1.bin, u_init_tc1.bin, v_init_tc1.bin")
print(f"u0 = {u0:.6f} m/s")

