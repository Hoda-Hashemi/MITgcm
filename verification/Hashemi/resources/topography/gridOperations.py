#%%
import numpy as np
import xarray as xr
import os

# --- CONFIGURATION ---
# Path to your ETOPO5 netCDF file
bathymetry_file_path = "resources/topography/etopo5.csv" 

# Target MITgcm Grid Dimensions
nLat, nPhi = 720, 1440
dLat = np.pi / nLat
dPhi = 2 * np.pi / nPhi

def create_mitgcm_bathymetry(input_path, output_filename="bathymetry.bin"):
    print(f"READING BATHYMETRY FROM: {input_path}")
    
    # 1. Load Dataset
    with xr.open_dataset(input_path) as ds:
        # Use ETOPO5 specific variable names
        lon_bathy = ds["topo_lon"].values
        lat_bathy = ds["topo_lat"].values
        topo = ds["topo"].values  # Elevation (meters)

    nLonB, nLatB_orig = topo.shape[1], topo.shape[0]
    
    # Padding for interpolation logic
    topo_padded = np.vstack([topo, topo[-1:, :]])
    topo_data = topo_padded.T  # Shape (lon, lat_padded)
    
    nLonB, nLatB = topo_data.shape
    dLatB = np.pi / (nLatB_orig - 1)
    dPhiB = 2 * np.pi / nLonB

    # 2. Interpolate to Model Grid
    depth_grid = np.zeros((nLat, nPhi))

    for iLat in range(nLat):
        # Model latitude center
        Lat0 = -0.5 * np.pi + 0.5 * dLat + dLat * iLat

        iLatB = 1 + int((Lat0 + 0.5 * np.pi - 0.5 * dLatB) / dLatB)
        iLatT = min(iLatB + 1, nLatB)

        LatB = -0.5 * np.pi + 0.5 * dLatB + dLatB * (iLatB - 1)
        LatT = -0.5 * np.pi + 0.5 * dLatB + dLatB * (iLatT - 1)
        alpha_lat = (Lat0 - LatB) / (LatT - LatB)

        for iPhi in range(nPhi):
            # Model longitude center
            phi0 = 0.5 * dPhi + dPhi * iPhi

            iPhiL = 1 + int((phi0 - 0.5 * dPhiB) / dPhiB)
            iPhiR = iPhiL + 1

            # Periodic wrap for longitude
            if iPhiL > nLonB: iPhiL -= nLonB
            if iPhiR > nLonB: iPhiR -= nLonB
            
            # Bilinear interpolation
            #!!!!!! FOR MITGCM THEY READ THE BATHYMETRY OF THE OCEAN AS NEGATIVE
            #!!!! Note: Removed the '-' sign here so ocean remains negative (MITgcm standard)
            zB_L = topo_data[iPhiL-1, iLatB-1]
            zT_L = topo_data[iPhiL-1, iLatT-1]
            zphiL = zB_L + alpha_lat * (zT_L - zB_L)

            zB_R = topo_data[iPhiR-1, iLatB-1]
            zT_R = topo_data[iPhiR-1, iLatT-1]
            zphiR = zB_R + alpha_lat * (zT_R - zB_R)

            phiL = 0.5 * dPhiB + dPhiB * (iPhiL - 1)
            phiR = 0.5 * dPhiB + dPhiB * (iPhiR - 1)
            alpha_lon = (phi0 - phiL) / (phiR - phiL)

            depth_grid[iLat, iPhi] = zphiL + alpha_lon * (zphiR - zphiL)

    # 3. MITgcm Specific Processing
    # Ensure land (positive elevation) is exactly 0
    depth_grid[depth_grid > 0] = 0.0
    
    # 4. Save as Big-Endian Float32 Binary
    # MITgcm usually expects Big-Endian ('>f4')
    depth_grid.astype('>f4').tofile(output_filename)
    
    print(f"DONE: {output_filename} created.")
    print(f"Min depth: {depth_grid.min()}, Max depth: {depth_grid.max()}")

# Execute
if __name__ == "__main__":
    create_mitgcm_bathymetry(bathymetry_file_path)

# %%

