#%%
#COMMENT: Libraries imported
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import warnings
warnings.filterwarnings("ignore")
from src.config import nLat, nPhi, R, OMEGA, GRAVITY, dLat, Lat, sinLat,D , sinLat0

#%% 
#!ADDED 1D Solver
class SphereSolver1D:
   
    def __init__(self):
        self.nLat = nLat
        self.nPhi = nPhi
        self.R = R
        self.OMEGA = OMEGA
        self.GRAVITY = GRAVITY

        self.dLat = dLat
        self.Lat = Lat
        self.sinLat = sinLat
        self.sinLat0 = sinLat0

        self.CORIOLIS = 2.0*OMEGA*self.sinLat
    
    def build_matrix(self):
        nLat = self.nLat
        Lat = self.Lat 
        dLat = self.dLat

        cosLat = np.cos(Lat) 

        # main_diag = -(self.Ld**(-2)) * np.ones(nLat)
        main_diag =np.zeros(nLat)
        upper_diag = np.zeros(nLat-1)
        lower_diag = np.zeros(nLat-1)

        for i in range(1, nLat-1):
            cosLatPHalf = np.cos(Lat[i] + dLat/2)
            cosLatMHalf = np.cos(Lat[i] - dLat/2)
            #!ROCKET: UNSCALED CASE 1
            lower_diag[i-1] = cosLatMHalf / (dLat**2 * self.R**2 * cosLat[i])
            upper_diag[i] = cosLatPHalf / (dLat**2 * self.R**2 * cosLat[i])
            main_diag[i] = -(cosLatPHalf +cosLatMHalf) / (dLat**2 * self.R**2 * cosLat[i]) - ( (2.0*OMEGA*self.sinLat0)**2 /(GRAVITY * D) )

            #!ROCKET: SCALED CASE 1
            # lower_diag[i-1] = cosLatMHalf 
            # upper_diag[i] = cosLatPHalf 
            # main_diag[i] = -(cosLatPHalf +cosLatMHalf)  - ( (2.0*OMEGA*self.sinLat0)**2 /(GRAVITY * D) )*(dLat**2 * self.R**2 * cosLat[i])
            
        main_diag[0] = -1.0
        upper_diag[0] = 1.0
        main_diag[-1] = 1.0
        lower_diag[-1] = -1.0

        diagonals = [main_diag, upper_diag, lower_diag]
        A = diags(diagonals, [0, 1, -1], format='csr')
        
        # print(A)

        return A
    
    def build_rhs(self):
        cosLat = np.cos(Lat) 
        #!ROCKET: UNSCALED CASE 1
        rhs = - 2.0*OMEGA*self.sinLat

        #!ROCKET SCALED CASE 2
        # rhs= -(2.0*OMEGA*self.sinLat)*(dLat**2 * self.R**2 * cosLat)

        #!COMMENT: Boundary conditions (Neumann: zero gradient)
        rhs[0] = rhs[-1] = 0.0

        return rhs

    def solve(self):
        A = self.build_matrix()
        # print("COEFFICIENT MATRIX:\n", A)
        rhs = self.build_rhs()
        # print("RHS : -f :\n", rhs)
        psi = spsolve(A, rhs)
        print(f"ψ₁ range: min = {np.min(psi):.4e}, max = {np.max(psi):.4e}")

        return psi

    def compute_velocity(self, psi):
        dLat = self.dLat

        dpsi_dLat = np.zeros_like(psi)
        dpsi_dLat[1:-1] = (psi[2:] - psi[:-2]) / (2 * dLat)

        dpsi_dLat[0] =(psi[1] - psi[0]) / (dLat)
        dpsi_dLat[-1] = (psi[-1] - psi[-2]) / (dLat)

        u = - dpsi_dLat/self.R
        v = np.zeros_like(u)
        print(f"u_ψ₁ range: min = {np.min(u):.4e}, max = {np.max(u):.4e}")

        return u, v

