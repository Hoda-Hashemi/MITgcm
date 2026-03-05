#%% !COMMENT: Global definitions
import numpy as np
import os

## !COMMENT: GLOBAL CONSTANTS & PHYSICAL PARAMETERS
GRAVITY= 9.81                                # (m/s²)
OMEGA = (2 * np.pi) / (24 * 3600)       # (rad/s)
R = 6371e3                              # (m)

## !COMMENT:  GRID PARAMETERS

#ROCKET!!
nLat = 720
nPhi = 2 * nLat  
dLat = np.pi / nLat
dPhi = dLat    
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PAY ATTENTION HODA LOOK HERE LOOK HERE LOOK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#& GridFortran:
Lat = -np.pi/2 + dLat/2 + dLat * np.arange(nLat)
#&GridPoles: 
# Lat = -np.pi/2 + dLat * np.arange(nLat)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PAY ATTENTION HODA LOOK HERE LOOK HERE LOOK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# #&Grid middle section:
# # min and max latitude in radians
# lat_min = -np.pi/4
# lat_max = np.pi/4

# # compute start and end indices
# i_start = int((lat_min + np.pi/2 - dLat/2) / dLat)
# i_end   = int((lat_max + np.pi/2 - dLat/2) / dLat) + 1  # +1 to include last point

# # full latitude array
# Lat = -np.pi/2 + dLat/2 + dLat * np.arange(nLat)

# # select middle section
# Lat_slice = Lat[i_start:i_end]

# # Lat = -(4*np.pi/9) + dLat/2 + dLat * np.arange(nLat) #start from -80 to 80

sinLat = np.sin(Lat)
cosLat = np.cos(Lat)

## !COMMENT: DYNAMICAL PARAMETERS
Lat0 = np.pi / 4.0
sinLat0 = np.sin(Lat0)

D = 1000.0 
f2OMEGA = 2 * OMEGA * sinLat0

## !COMMENT: Ld0 - Rossby radius of deformation
# Original placeholder values from Fortran:
# Ld0 = np.inf  
# Ld0 = R           
# Ld0 = 1e6             
# R/10, R/100, R/1000 
# Ld0 = sqrt(GRAVITY * D) / f2OMEGA

# Compute Coriolis parameter at reference latitude

# Dynamically computed Ld0
Ld0 = R/1000
# np.sqrt(GRAVITY * D) / f2OMEGA

#!ALTERED: DOUBLE CHECK FOR EACH RUN
includeFreeSurface = True

## !COMMENT: neighbor handling
neighborhandling = 'zero'#'zero'#'include'

##! COMMENT:
ScaleEquation = True  #False 
readbathymetryFlag= True #True #False 

bathymetryfilePath = '../resources/topography/etopo5.nc'
bathymetryCSVfilePath = '../resources/topography/etopo5.csv'

# OUTPUT_PATH='output/ComparingSolversPython/iterationConstraint/scaled/iter10000'
# OUTPUT_PATH='output/ComparingSolversPython/unscaled/'
# OUTPUT_PATH='output/ComparingSolversPython/Preconditioning/unscaled/'
# OUTPUT_PATH='output/FluxForm/ILU_preconditioning_rtol_overS/'
# OUTPUT_PATH='output/FluxForm/Focus/Scaled_PolesScaled/'
# OUTPUT_PATH='output/FluxForm/Focus/Scaled_NoPhi/'
# OUTPUT_PATH='output/FluxForm/Focus/Scaled_Neumann/'
# OUTPUT_PATH='output/FluxForm/Focus/Unscaled_WithNeumann/'
# OUTPUT_PATH='output/FluxForm/Focus/Final_OnlyPoles_Scaled_WithNeumann/'
# OUTPUT_PATH='output/FluxForm/Focus/Mod_Neum_Ld_theta/'
# OUTPUT_PATH='output/FluxForm/Wings/Coarserby2'
# OUTPUT_PATH='output/FluxForm/Wings/Coarserby4'
# OUTPUT_PATH='output/FluxForm/Wings/delta2/'
# OUTPUT_PATH='output/FluxForm/Wings/Finerby2/'
# OUTPUT_PATH='output/FluxForm/BRICK/' #-80 to 80

# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_varying/'
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_varying_scaled/' #R2 scaling
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_varying_scaled_stabilization/' #R2 scaling+stabilization
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_varying_fullscaled_stabilization/'

# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_constant/'
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_constant_phi/'

# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_varying_stab/'
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/RHS_varying_stab_multiply/' #best unscaled with stabilization

# OUTPUT_PATH = 'output/1_14AM/PolesGrid/RHS_varying/'
# OUTPUT_PATH = 'output/1_14AM/PolesGrid/RHS_varying_FortranGrid/'
# OUTPUT_PATH = 'output/1_14AM/PolesGrid/RHS_varying_stabilization/'

# OUTPUT_PATH = 'output/1_14AM/PolesGrid/RHS_varying_More_stabilization/'

# OUTPUT_PATH = 'output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled/' #R2 coslat dlat2 scaling fortran grid  poisson (no stabilizer)
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled_stabilized/' #R2 coslat dlat2 scaling fortran grid  poisson ( stabilizer)

# OUTPUT_PATH = 'output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_1000/' 
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_100000/' 
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover10/' 
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover100/' 
# OUTPUT_PATH = 'output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover1000/'

# OUTPUT_PATH = 'output/1_14AM/FortranGrid/TestingSphericalHarmonics/VariableScreening/'
# OUTPUT_PATH = 'output/Fortran/LambdaPython/unscreened/'

#! FortranFilling + without alpha on diagonals
# OUTPUT_PATH = 'output/Fortran/EquationOne/Python/unnormalized_Noaddition_diagonal/'

# #! FortranFilling + with alpha on diagonals
# OUTPUT_PATH = 'output/Fortran/EquationOne/Python/unnormalized_addition_diagonal/'

#! FortranFilling + normalization
# OUTPUT_PATH = 'output/Fortran/EquationOne/Python/normalized_Noaddition_diagonal/'

# #! FortranFilling + with alpha on diagonals + normalization
# OUTPUT_PATH = 'output/Fortran/EquationOne/Python/normalized_addition_diagonal/'

# #! MY way of flagging of land neighors + adding
# OUTPUT_PATH = 'output/Fortran/EquationOne/PythonHoda/unnormalized_addition_diagonal/'

#! MY way of flagging of land neighors +normalized+adding
# OUTPUT_PATH = 'output/Fortran/EquationOne/PythonHoda/normalized_addition_diagonal/'

OUTPUT_PATH = 'output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/'

os.makedirs(OUTPUT_PATH, exist_ok=True)

#%%
# import matplotlib.pyplot as plt
# scaling = R**2 * cosLat**2 * dLat**2

# # ---------------- PLOT ----------------
# plt.figure(figsize=(8,5))
# plt.plot(Lat * 180/np.pi, scaling, color='blue')
# plt.xlabel("Latitude (degrees)")
# plt.ylabel(r"Scaling term $R^2 \cos^2(\phi) \Delta\phi^2$")
# plt.title("Scaling factor vs Latitude")
# plt.grid(True)

# scaling = R**2 * cosLat**2 * dLat**2

# print('Min, Max = ', min(scaling), max(scaling))

# %%

