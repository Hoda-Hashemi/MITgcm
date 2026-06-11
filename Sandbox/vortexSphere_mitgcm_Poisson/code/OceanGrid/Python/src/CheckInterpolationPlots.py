#%%
#!COMMENT: THIS SCRIPT IS TO DOUBLE CHECK WHETHER PYTHON AND FORTRAN ARE INTERPOLATIONS ARE CORRECT ON RESOLUTION OF 720x1440 Ocean grid
#!COMMENT: IF YOU WISH TO CHANGE THE RESOLUTION, EXTRACT DEPTHS INTERPOLATED FROM FORTRAN-PARAVIEW AND THEN SAVE THEM IN DATA/
#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from src.gridOperations import OceanGrid  

ocean_grid = OceanGrid()
ocean_grid.readBathymetry()

csv_path = 'data/interpolatedDepthParaview.csv' 

df_fortran = pd.read_csv(csv_path)

df_fortran['depth'] = df_fortran['d₀']  
df_fortran['lat_deg'] = np.rad2deg(df_fortran['lat']) 
df_fortran['long_deg'] = np.rad2deg(df_fortran['long'])

depth_matrix = df_fortran['depth'].values.reshape(ocean_grid.nLat, ocean_grid.nPhi)
depth_matrix[depth_matrix <= 0] = 0  

depth_py_matrix = ocean_grid.d

lat_deg = np.rad2deg(np.linspace(-np.pi/2 + ocean_grid.dLat/2, np.pi/2 - ocean_grid.dLat/2, ocean_grid.nLat))
lon_deg = np.rad2deg(np.linspace(0 + ocean_grid.dPhi/2, 2*np.pi - ocean_grid.dPhi/2, ocean_grid.nPhi))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
#!COMMENT: FORTRAN DEPTHS
im1 = ax1.imshow(
    depth_matrix, 
    extent=[lon_deg.min(), lon_deg.max(), lat_deg.min(), lat_deg.max()], 
    origin='lower', cmap='Blues_r', vmin=0, vmax=5000
)
ax1.set_title('Fortran Depths (Exact Matrix)')
ax1.set_xlabel('Longitude (°E)')
ax1.set_ylabel('Latitude (°N)')
ax1.grid(True, alpha=0.3)
plt.colorbar(im1, ax=ax1, shrink=0.8, label='Depth (m)')

im2 = ax2.imshow(
    depth_py_matrix, 
    extent=[lon_deg.min(), lon_deg.max(), lat_deg.min(), lat_deg.max()], 
    origin='lower', cmap='Blues_r', vmin=0, vmax=5000
)
ax2.set_title('Python Depths (Exact Matrix)')
ax2.set_xlabel('Longitude (°E)')
ax2.set_ylabel('Latitude (°N)')
ax2.grid(True, alpha=0.3)
plt.colorbar(im2, ax=ax2, shrink=0.8, label='Depth (m)')
plt.tight_layout()
plt.show()

#!COMMENT: EXACT OVERLAY
# fig, ax = plt.subplots(figsize=(14, 8))
# ax.imshow(
#     depth_py_matrix, 
#     extent=[lon_deg.min(), lon_deg.max(), lat_deg.min(), lat_deg.max()], 
#     origin='lower', cmap='Blues_r', vmin=0, vmax=5000, alpha=0.6
# )
# ax.imshow(
#     depth_matrix, 
#     extent=[lon_deg.min(), lon_deg.max(), lat_deg.min(), lat_deg.max()], 
#     origin='lower', cmap='Reds_r', vmin=0, vmax=5000, alpha=0.6
# )
# ax.set_title('Exact Overlay: Fortran (Red) + Python (Blue)\n(Alpha Blend for Mismatch Spots)')
# ax.set_xlabel('Longitude (°E)')
# ax.set_ylabel('Latitude (°N)')
# ax.grid(True, alpha=0.3)
# plt.colorbar(im2, ax=ax, shrink=0.8, label='Depth (m)')  # Shared
# ax.invert_yaxis()  # Standard map: N at top
# plt.tight_layout()
# plt.show()

#!COMMENT: EXACT DIFFERENCE - IF BLANK THIS IS GOOD :)
diff_matrix = depth_matrix - depth_py_matrix
fig, ax = plt.subplots(figsize=(12, 6))
im_diff = ax.imshow(
    diff_matrix, 
    extent=[lon_deg.min(), lon_deg.max(), lat_deg.min(), lat_deg.max()], 
    origin='lower', cmap='RdBu_r', vmin=-100, vmax=100  
)
ax.set_title('Exact Depth Difference: Fortran - Python (m)')
ax.set_xlabel('Longitude (°E)')
ax.set_ylabel('Latitude (°N)')
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(im_diff, ax=ax, shrink=0.8)
cbar.set_label('Difference (m)')
plt.tight_layout()
plt.show()

#!COMMENT: VALIDATION STATS
print(f"Mean abs diff: {np.mean(np.abs(diff_matrix)):.2f} m")
print(f"Max abs diff: {np.max(np.abs(diff_matrix)):.2f} m")
print(f"Global mean Fortran: {np.mean(depth_matrix):.2f} m")
print(f"Global mean Python: {np.mean(depth_py_matrix):.2f} m")

# %%
#!COMMENT: SCATTER PLOTTING AND EVALUATING THE DEPTHS.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from src.gridOperations import OceanGrid  
ocean_grid = OceanGrid()
ocean_grid.readBathymetry()
d = ocean_grid.d
depths = d[d > 0].flatten()

plt.figure(figsize=(8,6))
plt.scatter(np.arange(len(depths)), depths, c=depths, s=20)
plt.colorbar(label='Depth')
plt.xlabel('Point Index')
plt.ylabel('Depth')
plt.title('Scatter Plot of Ocean Depths (d > 0)')
plt.show()
# %%
