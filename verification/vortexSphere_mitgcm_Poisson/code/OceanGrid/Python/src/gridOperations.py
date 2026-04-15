#%% !ROCKET
import numpy as np
import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
import math
from src.config import (
    nLat, nPhi, dLat, dPhi, D,
    readbathymetryFlag, bathymetryfilePath,bathymetryCSVfilePath
)

#!COMMENT: THIS CLASS READS OR SETS THE BATHYMETRY AND INTERPOLATES IT ONTO THE OCEAN GRID
class OceanGrid:
    def __init__(self):
        self.nLat = nLat
        self.nPhi = nPhi
        self.dLat = dLat
        self.dPhi = dPhi
        self.d  = np.zeros((nLat, nPhi))   

        self.readbathymetryFlag = readbathymetryFlag 
        self.bathymetryfilePath = bathymetryfilePath
        self.bathymetryCSVfilePath = bathymetryCSVfilePath

    def readBathymetry(self):
        if not self.readbathymetryFlag:
            print("FILLING GRID WITH DEFAULT DEPTH D")
            self.d.fill(D)  

        else:
            print("READING BATHYMETRY VALUES FROM FILE")

            ds = xr.open_dataset(self.bathymetryfilePath)
            print("DATASET STATS:")
            print(ds)
            lon = ds["topo_lon"].values
            lat = ds["topo_lat"].values
            elev = ds["topo"].values

            print(f"Longitude shape: {lon.shape}")
            print(f"Latitude shape: {lat.shape}")
            print(f"Elevation shape: {elev.shape}")

            print(f"Elevation min: {elev.min()}, max: {elev.max()}, mean: {elev.mean()}")
            # Run interpolation
            self.interpolateBathymetry(topo=elev)

            print("OCEAN GRID IS CONSTRUCTED")

            nOcean = np.sum(self.d > 0)
            nLand = np.sum(self.d <= 0)
            print(f"TOTAL GRID POINTS = {self.nLat * self.nPhi}")
            print(f"OCEAN POINTS (d>0) = {nOcean}")
            print(f"LAND POINTS (d<=0) = {nLand}")

            ds.close()
            # sys.exit()
    
    #!ROCKET: INTERPOLATE BATHYMETRY TO MODEL GRID (Manual Bilinear, Fortran-Matching)
    def saveBathymetry(self):
        with xr.open_dataset(bathymetryfilePath) as ds:
            topo = ds["topo"].values
            lat_bathy = ds["topo_lat"].values
            lon_bathy = ds["topo_lon"].values 

        nLatB, nLonB = topo.shape

        # Get dimensions
        nLatB, nLonB = topo.shape
        print(f"Bathymetry data loaded: {nLatB} x {nLonB}")

        # Ensure output directory exists
        save_dir = "../data"
        os.makedirs(save_dir, exist_ok=True)

        save_path = os.path.join(save_dir, "bathymetry_data.npz")
        np.savez_compressed(save_path, topo=topo, lat=lat_bathy, lon=lon_bathy)

        print(f"BATHYMETRY DATA IS SAVED, FIND THEM IN: {save_path}")
        return save_path

    def interpolateBathymetry(self, topo):
        orig_nLatB = topo.shape[0]
        topo_padded = np.vstack([topo, topo[-1:,:]])
        topo = topo_padded.T  # (lon, lat_padded)
        nLonB, nLatB = topo.shape
        dLatB = np.pi / (orig_nLatB - 1)  # Correct: π / 2160
        dPhiB = 2 * np.pi / nLonB

        self.d = np.zeros((self.nLat, self.nPhi))

        for iLat in range(self.nLat):
            Lat0 = -0.5 * np.pi + 0.5 * self.dLat + self.dLat * iLat

            iLatB = 1 + int((Lat0 + 0.5 * np.pi - 0.5 * dLatB) / dLatB)
            iLatT = iLatB + 1
            if iLatT > nLatB:
                iLatT = nLatB

            LatB = -0.5 * np.pi + 0.5 * dLatB + dLatB * (iLatB - 1)  # Fixed
            LatT = -0.5 * np.pi + 0.5 * dLatB + dLatB * (iLatT - 1)  # Fixed
            alpha_lat = (Lat0 - LatB) / (LatT - LatB)

            for iPhi in range(self.nPhi):
                phi0 = 0.5 * self.dPhi + self.dPhi * iPhi

                iPhiL = 1 + int((phi0 - 0.5 * dPhiB) / dPhiB)
                iPhiR = iPhiL + 1

                # Periodic wrap
                if iPhiL > nLonB:
                    iPhiL -= nLonB
                if iPhiR > nLonB:
                    iPhiR -= nLonB
                if iPhiL < 1:
                    iPhiL += nLonB
                if iPhiR < 1:
                    iPhiR += nLonB

                # 0-based indices
                iLatB0, iLatT0 = iLatB - 1, iLatT - 1
                iPhiL0, iPhiR0 = iPhiL - 1, iPhiR - 1

                phiL = 0.5 * dPhiB + dPhiB * (iPhiL - 1)  # Fixed
                phiR = 0.5 * dPhiB + dPhiB * (iPhiR - 1)  # Fixed
                alpha_lon = (phi0 - phiL) / (phiR - phiL)

                # Bilinear interp
                zB = topo[iPhiL0, iLatB0]
                zT = topo[iPhiL0, iLatT0]
                zphiL = zB + alpha_lat * (zT - zB)

                zB = topo[iPhiR0, iLatB0]
                zT = topo[iPhiR0, iLatT0]
                zphiR = zB + alpha_lat * (zT - zB)

                self.d[iLat, iPhi] = - (zphiL + alpha_lon * (zphiR - zphiL))

        # # Clamp land to 0
        # self.d[self.d < 0] = 0

        print("INTERPOLATION BILINEAR DONE")