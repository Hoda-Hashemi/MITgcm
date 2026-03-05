#%%
import matplotlib.pyplot as plt 
import numpy as np 
from MITgcmutils import rdmds
import re,os,glob
#%% 
#!COMMENT: reading the data of the experiment: verification/ajustment_32x32x1

path ='/home/hoda/Documents/Scripts/MITgcm/verification/adjustment.cs-32x32x1/run'

files = glob.glob(path + "/U.*.meta")

its = [int(re.findall(r'\.(\d+)\.meta$', f)[0]) for f in files]
its = sorted(its)

print("Found iterations:", its)
U_all = [rdmds(f"{path}/U", it) for it in its]
V_all = [rdmds(f"{path}/V", it) for it in its]
W_all = [rdmds(f"{path}/W", it) for it in its]

U_all = np.array(U_all)
V_all = np.array(V_all)
W_all = np.array(W_all)
XC = rdmds(f"{path}/XC")
YC = rdmds(f"{path}/YC")

#%%
print("U range:", U_all.min(), U_all.max())
print("V range:", V_all.min(), V_all.max())
print("W range:", W_all.min(), W_all.max())

vel_mag = np.sqrt(U_all**2 + V_all**2 + W_all**2)
print("Velocity magnitude range:", vel_mag.min(), vel_mag.max())



idx = 2 # iteration


vel_mag = np.sqrt(U_all[idx]**2 + V_all[idx]**2 + W_all[idx]**2)

hFacC = rdmds(f"{path}/hFacC")  # shape (Nx, Ny)
vel_mag_masked =vel_mag# np.ma.masked_where(hFacC==0, vel_mag)

plt.figure(figsize=(10,5))
plt.pcolormesh(XC, YC, vel_mag_masked, shading='auto', cmap='viridis')
plt.colorbar(label='Velocity magnitude (m/s)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title(f'Velocity Magnitude (iteration {its[idx]})')
plt.show()


#%% Overlay of multiple iterations
idxs = [0, 1, 2]
colors = ['r', 'g', 'b']

plt.figure(figsize=(10,5))
for i, idx in enumerate(idxs):
    vel = np.sqrt(U_all[idx]**2 + V_all[idx]**2 + W_all[idx]**2)
    vel_masked = np.ma.masked_where(hFacC==0, vel)
    # For overlay, take the "center face" (Face 3) as representative
    plt.plot(vel_masked[Nx:2*Nx, :Ny].mean(axis=0), color=colors[i], label=f'iter {its[idx]}')

plt.xlabel('Longitude index (Face 3)')
plt.ylabel('Velocity magnitude (m/s)')
plt.title('Velocity magnitude overlay (Face 3)')
plt.legend()
plt.show()

#%%
import matplotlib.pyplot as plt
import numpy as np
Nx, Ny = 32, 32 
idxs = [0,1,2] 
colors = ['r', 'g', 'b']
label_iters = [its[i] for i in idxs]


plt.figure(figsize=(10,5))
for i, idx in enumerate(idxs):
    vel = np.sqrt(U_all[idx]**2 + V_all[idx]**2 + W_all[idx]**2)
    vel_masked = vel# np.ma.masked_where(hFacC==0, vel)
    plt.plot(XC, vel_masked, color=colors[i])

plt.xlabel('Longitude')
plt.ylabel('Velocity magnitude (m/s)')
plt.title('Velocity magnitude overlay')
plt.legend()
plt.show()


# %%

#%% Utility function for cubed-sphere plotting
def plot_cubedsphere_faces(faces, title=''):
    """
    faces: list of 6 arrays (Nx, Ny)
    T-layout: 
         2
    1 3 5 6
         4
    """
    Nx, Ny = faces[0].shape
    canvas = np.full((3*Nx, 4*Ny), np.nan)

    canvas[0:Nx, Ny:2*Ny] = faces[1]      # Face 2 (top)
    canvas[Nx:2*Nx, 0:Ny] = faces[0]      # Face 1 (left)
    canvas[Nx:2*Nx, Ny:2*Ny] = faces[2]   # Face 3 (center)
    canvas[Nx:2*Nx, 2*Ny:3*Ny] = faces[4] # Face 5 (right)
    canvas[Nx:2*Nx, 3*Ny:4*Ny] = faces[5] # Face 6 (far right)
    canvas[2*Nx:3*Nx, Ny:2*Ny] = faces[3] # Face 4 (bottom)

    plt.figure(figsize=(12,8))
    plt.imshow(canvas, origin='upper', cmap='viridis')
    plt.colorbar(label='Velocity magnitude (m/s)')
    plt.title(title)
    plt.axis('off')
    plt.show()

#%% Single iteration example
idx = 2
vel_mag = np.sqrt(U_all[idx]**2 + V_all[idx]**2 + W_all[idx]**2)
vel_mag_masked = np.ma.masked_where(hFacC==0, vel_mag)  # mask land

# Assuming vel_mag shape is (6*Nx, Ny), split into 6 faces
Nx, Ny = 32, 32  # adjust to tiles*tile_size
faces =[vel_mag_masked[:, i*Ny:(i+1)*Ny] for i in range(6)]

# Plot using T-layout
plot_cubedsphere_faces(faces, title=f'Velocity magnitude (iteration {its[idx]})')



# %%
#%%
import numpy as np
import glob, re
from MITgcmutils import rdmds
import pyvista as pv

#%% --- Paths ---
path = '/home/hoda/Documents/Scripts/MITgcm/verification/adjustment.cs-32x32x1/run'

# Find iterations
files = glob.glob(path + "/U.*.meta")
its = sorted([int(re.findall(r'\.(\d+)\.meta$', f)[0]) for f in files])
print("Found iterations:", its)

# Read all fields
U_all = np.array([rdmds(f"{path}/U", it) for it in its])
V_all = np.array([rdmds(f"{path}/V", it) for it in its])
W_all = np.array([rdmds(f"{path}/W", it) for it in its])
hFacC = rdmds(f"{path}/hFacC")  # land mask

# Select iteration to export
idx = 2
U = U_all[idx]
V = V_all[idx]
W = W_all[idx]

# Compute velocity magnitude and mask land
vel_mag = np.sqrt(U**2 + V**2 + W**2)
vel_mag_masked = np.ma.masked_where(hFacC == 0, vel_mag)

#%% --- Split into 6 faces ---
Nx, Ny = 32, 32  # points per face
faces = [vel_mag_masked[:, i*Ny:(i+1)*Ny] for i in range(6)]

# --- Generate coordinates for cube faces ---
# Earth radius in meters
R = 6371e3
coords_list = []

# Simple cube-to-sphere mapping for visualization
def cube_to_sphere(ix, iy, Nx, Ny, face_id, R=6371e3):
    # Normalized coordinates [-1,1]
    a = 2*ix/(Nx-1) - 1
    b = 2*iy/(Ny-1) - 1
    if face_id == 0:   # Front
        x, y, z = 1, b, -a
    elif face_id == 1: # Right
        x, y, z = -b, 1, -a
    elif face_id == 2: # Back
        x, y, z = -1, -b, -a
    elif face_id == 3: # Left
        x, y, z = b, -1, -a
    elif face_id == 4: # Top
        x, y, z = a, b, 1
    elif face_id == 5: # Bottom
        x, y, z = a, b, -1
    # Normalize to sphere
    norm = np.sqrt(x**2 + y**2 + z**2)
    return R*x/norm, R*y/norm, R*z/norm

for fid, f in enumerate(faces):
    X = np.zeros_like(f)
    Y = np.zeros_like(f)
    Z = np.zeros_like(f)
    for ix in range(Nx):
        for iy in range(Ny):
            X[ix, iy], Y[ix, iy], Z[ix, iy] = cube_to_sphere(ix, iy, Nx, Ny, fid, R)
    coords_list.append((X, Y, Z, f))

#%% --- Combine all faces into a single structured grid ---
# Flatten faces along one axis to make a big structured grid
total_Nx = Nx * 2  # can adjust spacing
total_Ny = Ny * 3  # approximate for 6 faces T-layout
# We'll append faces using pyvista's Append filter
combined = pv.MultiBlock()
for fid, (X, Y, Z, F) in enumerate(coords_list):
    grid = pv.StructuredGrid(X, Y, Z)
    grid["Velocity"] = F.flatten(order="F")
    combined[f"face_{fid}"] = grid

#%% --- Save to VTK file ---
combined.save("cubedsphere_velocity.vtm")  # .vtm = MultiBlock for ParaView

print("Saved cubed-sphere velocity to cubedsphere_velocity.vtm")
