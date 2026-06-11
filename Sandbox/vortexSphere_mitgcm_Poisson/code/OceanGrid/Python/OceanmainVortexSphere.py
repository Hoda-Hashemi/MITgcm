#%%
import numpy as np
import warnings
import os

import sys
warnings.filterwarnings("ignore")

from src.gridOperations import OceanGrid 
from src.config import nLat,nPhi, ScaleEquation, OUTPUT_PATH
from src.SphereSolver_2D import SphereSolver2D
from src.SphereSolver_1D import SphereSolver1D
from src.utils import ParaviewSavedData, setup_logging
from src.plotting import plot_solution, plot_comparison_1d_paraview, plot_comparison_2d_paraview

# def main():

    # run_dir = setup_logging("output")
    # os.makedirs(run_dir, exist_ok=True)

print("FRAMEWORK:")
print(f"θ , ƛ  = {nLat} , {nPhi}")

N = nLat*nPhi

#!COMMENT: GRID OPERATIONS 
oceanGrid= OceanGrid()
oceanGrid.readBathymetry()

# # #!COMMENT: PSI Solver 
solver2D = SphereSolver2D(Ocean_Grid=oceanGrid)

PSI_direct, PSI_cg, PSI_bicg, PSI_adi = solver2D.PsiSolve()
# PSI_direct,PSI_cg, PSI_bicg = solver2D.PsiSolve()

#!Running some analysis on coefficient matrix:
# A = solver2D.Psi_build_matrix()
# RHS = solver2D.Psi_build_rhs()
# solver2D.interactive_coefficient_matrix_plot(A)
# solver2D.plot_A_1d_plotly(A)
# solver2D.plot_rhs_1d_plotly(RHS)
# #!ROCKET:NICE WAY TO SEE SMALL COEFFICIENT MATRIX IN INTERACTIVE WAY TRY NLAT=3.
# solver2D.interactive_plot_coefficient_small_n_full(A)

#%%
# ?SPARSITY CHECK: if small output means sparse
# sparse symmetric check
# diff = A - A.T
# norm_diff = np.linalg.norm(diff.data)
# print("Symmetric (sparse)? norm(A - A.T) =", norm_diff)

# #? SPD or not, if the smallest three eigenvalues are + -> SPD otherwise A is not PD
# from sksparse.cholmod import cholesky
# try:
#     factor = cholesky(A)   # sparse
#     print("Matrix is SPD (Cholesky succeeded).")
# except Exception as e:
#     print("Matrix is NOT SPD (Cholesky failed):", e)
# #? Smallest eigenvalues
# from scipy.sparse.linalg import eigsh
# vals, vecs = eigsh(A, k=5, sigma=0.0)
# print("smallest approx eigs:", vals)

# #?SPARSITY:
# print("nnz:", A.nnz)
# print("density:", A.nnz / (A.shape[0]**2))

# #?RHS norm:
# print("||RHS||2 =", np.linalg.norm(RHS))
# print("min RHS, max RHS =", RHS.min(), RHS.max())

# #? took 25 min
# from scipy.sparse.linalg import eigsh
# import numpy as np

# vals_small, _ = eigsh(A, k=3, which='SM') 
# vals_big, _   = eigsh(A, k=1, which='LM')   

# lambda_min = vals_small.min()
# lambda_max = vals_big[0]
# cond_est = lambda_max / max(lambda_min, np.finfo(float).eps)

# cond_Kappa = lambda_max / lambda_min

# print("lambda_min", lambda_min)
# print("lambda_max", lambda_max)
# print("cond_est", cond_est)
# print('Kappa = ', cond_Kappa)

#%%
import numpy as np
from scipy.sparse import issparse
from scipy.sparse.linalg import norm as sparse_norm, eigsh, svds
import matplotlib.pyplot as plt

# 1. Symmetry check: Compute norm(A - A.T) to verify exact/approximate symmetry
# Should be exactly 0 for scaled A; if >0 but small (e.g., <1e-10), it's numerical artifact
sym_diff = A - A.T
sym_norm = sparse_norm(sym_diff)  # Use sparse norm for efficiency
print(f"Norm of (A - A.T): {sym_norm:.2e}")
if sym_norm < 1e-12:  # Threshold for floating-point tolerance
    print("Matrix A is symmetric (within numerical tolerance).")
else:
    print("Matrix A is not symmetric.")

# 2. Sparsity metrics: nnz and density (you may already have this, but for completeness)
if issparse(A):
    nnz = A.nnz
    density = nnz / (N * N)
    print(f"Number of non-zeros (nnz): {nnz}")
    print(f"Density: {density:.2e}")
else:
    print("A is not sparse; convert to sparse format for large N.")

# 3. Eigenvalues: Compute k smallest and largest (in magnitude, since A is neg. semi-def.)
# Use eigsh for symmetric matrices; k=20 is a good balance for clustering without long runtime
# tol=1e-12 for better accuracy near zero; which='LM' for largest mag, 'SM' for smallest mag
k = 20  # Adjust if needed; higher k takes longer (~minutes for N=1e6)
eig_large, vec_large = eigsh(A, k=k, which='LM', tol=1e-12, return_eigenvectors=True)
eig_small, vec_small = eigsh(A, k=k, which='SM', tol=1e-12, return_eigenvectors=True)

# Print min/max (smallest/largest in value; since eigenvalues <=0, max is least negative)
lambda_min = np.min(eig_small)  # Most negative in cluster near zero? Wait, no: for SM, smallest algebraic
lambda_max = np.max(eig_large)  # Largest algebraic (least negative if all <=0)
print(f"Estimated lambda_min (smallest eigenvalue): {lambda_min:.2e}")
print(f"Estimated lambda_max (largest eigenvalue): {lambda_max:.2e}")
print(f"Smallest {k} eigenvalues (sorted): {np.sort(eig_small)}")
print(f"Largest {k} eigenvalues (sorted descending): {np.sort(eig_large)[::-1]}")

# 5. Nullspace checker: Test if constants are in nullspace (||A * ones|| ≈ 0)
# Also compute smallest singular value/vector via svds (handles non-symmetric but works here)
ones_vec = np.ones(N)
null_check = A @ ones_vec
null_norm = np.linalg.norm(null_check)
rel_null_norm = null_norm / np.sqrt(N)  # Normalize by vector length
print(f"Norm of A * ones(N): {null_norm:.2e} (relative: {rel_null_norm:.2e})")
if rel_null_norm < 1e-10:
    print("Constants are in the approximate nullspace (expected for Laplacian).")

# More thorough: Smallest singular value/vector (sigma_min ≈ 0 indicates singularity)
k_sv = 1  # Just the smallest one
sigma_small, u_small, vt_small = svds(A, k=k_sv, which='SM', tol=1e-12)
print(f"Smallest singular value (sigma_min): {sigma_small[0]:.2e}")
if sigma_small[0] < 1e-10:
    print("Matrix is numerically singular (nullspace dimension >=1).")
null_vec_approx = vt_small[0, :]  # Right singular vector for smallest sigma

# Bonus: Condition number estimate (from eigenvalues; avoid full cond for large N)
if lambda_min != 0:
    cond_est = abs(lambda_max / lambda_min)
else:
    cond_est = np.inf  # Singular
print(f"Estimated condition number (kappa): {cond_est:.2e}")
#%%
# cond = solver2D.analyze_and_save_conditioning(A, OUTPUT_PATH)

# PSI_direct, PSI_bicgstab, PSI_ILU = solver2D.PsiSolve()

# Psi_u2D, Psi_v2D, Psi_umag2d= solver2D.Psi_compute_velocity(PSI2D)

    # # #!COMMENT: PHI Solver 
    # PHI2D = solver2D.PhiSolve(Psi=PSI2D)
    # Phi_u2D, Phi_v2D, Phi_umag2d= solver2D.Phi_compute_velocity(PHI2D)

    #!COMMENT: TOTAL VELOCITY
    # Compute total velocity
    # U_total = Psi_u2D + Phi_u2D
    # V_total = Psi_v2D + Phi_v2D

    # Umag_total = np.sqrt(U_total**2 + V_total**2)

    # print(f"uTOT (θ) range:  [{U_total .min(): .4e}, {U_total .max(): .4e}]")
    # print(f"uTOT (λ) range: [{V_total.min(): .4e}, {V_total.max(): .4e}]")
    # print(f"Velocity TOT magnitude range: [{Umag_total.min(): .4e}, {Umag_total.max(): .4e}]")

    # plot_solution(solver2D, PSI2D, u=u2D, v=v2D, run_dir = run_dir)

    #!COMMENT: 1D SOLVER WORKING
    # print("\n Solving 1D System ")
    # solver1D = SphereSolver1D()
    # PSI1D = solver1D.PsiSolve()
    # u1D, v1D = solver1D.compute_velocity(PSI1D)
    # plot_solution(solver1D, PSI1D, u=u1D, v=v1D,run_dir = run_dir)

    #!COMMENT:PARAVIEW COMPARISON

    #COMMENT: CASE1: values from paraview WITHOUT Phi
    # DATA_DIR = "data/CASE1_Paraview.csv" 
   
    #COMMENT: CASE2: values form paraview WITH Phi
    # DATA_DIR = "data/CASE2_Paraview.csv"

    #COMMENT: CASE3: values form paraview WITH:Phi and f(theta)
    # DATA_DIR = "data/CASE3_Paraview.csv"
    
    # DATA_DIR = "data/ConstantLd.csv"          #!ADDED: Constant Ld Scaled
    # DATA_DIR = "data/LdTheta.csv"             #!ADDED: Ld(theta) Scaled
    # DATA_DIR = "data/ThetaLambdaLd.csv"       #!ADDED: Ld(theta,lambda) Scaled
    # DATA_DIR = "data/UnscaledLdConstant.csv"  #!ADDED: Constant Ld Unscaled
    # DATA_DIR = "data/UnscaledLdTheta.csv"     #!ADDED: Constant Ld(theta) Unscaled
    # DATA_DIR = "data/UnscaledLdThetaLambda.csv" #!ADDED: Constant Ld(theta,lambda) Unscaled

    # df_para = ParaviewSavedData(os.path.join(DATA_DIR))
    # print(f"Reading {os.path.join(DATA_DIR)}")

    # print("\n=== Paraview vs 2D Comparison ===")
    # plot_comparison_2d_paraview(df_para, solver2D, PSI2D,Psi_umag2d, run_dir = run_dir)

    # # # print("\n=== Paraview vs 1D Comparison ===")
    # # # plot_comparison_1d_paraview(df_para, solver1D, PSI1D,run_dir = run_dir)

    # print("\nALL PROCESSES ARE COMPLETED")
    # print(f"TERMINAL SCRIPT IS SAVED IN: output/'")
# # #!COMMENT:  EXECUTION GUARD
# if __name__ == "__main__":
#     main()

