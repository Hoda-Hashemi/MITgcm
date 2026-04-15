#%% 
#ROCKET: Reading PSI
#ROCKET: =============================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.config import (nLat, nPhi,Lat, sinLat0, dLat, GRAVITY ,D, sinLat, OMEGA,R, ScaleEquation, neighborhandling, includeFreeSurface,Ld0, OUTPUT_PATH)
#%% 
#=============================================
#ROCKET: ../MeetingWed03/ ~ Fortran code
#=============================================
psi_ocean_scaled = pd.read_csv("../MeetingWed03/8/scaled/PsiOcean.csv")

#=============================================
#ROCKET: iterationConstraint/
#=============================================
#DIRECT SOLVER: UNSCALED, iter 1000
psi_python_sp = pd.read_csv("output/ComparingSolversPython/iterationConstraint/iter1000/unscaled/PythonPsi.csv")

#BICGSTAB SOLVER: SCALED,  iter 1000
psi_python_bicgstab_SCALED_iter1000 = pd.read_csv("output/ComparingSolversPython/iterationConstraint/iter1000/scaled/PythonPsi_Bicgstab.csv")

#BICGSTAB SOLVER: SCALED,  iter 5000
psi_python_bicgstab_iter5000 = pd.read_csv("output/ComparingSolversPython/iterationConstraint/scaled/iter5000/PythonPsi_Bicgstab.csv")

#BICGSTAB SOLVER: UNSCALED,  iter 1000
psi_python_bicgstab_UNSCALED_iter1000 = pd.read_csv("output/ComparingSolversPython/iterationConstraint/iter1000/unscaled/PythonPsi_Bicgstab.csv")

#=============================================
#ROCKET: Preconditioning/
#=============================================
#BICGSTAB SOLVER: Preconditioning ILU: SCALED
psi_python_bicgstab_iter1000_preconditioning= pd.read_csv("output/ComparingSolversPython/Preconditioning/scaled/PythonPsi_Bicgstab_preconditioned.csv")

#BICGSTAB SOLVER: Preconditioning ILU: UNSCALED
psi_python_bicgstab_iter1000_preconditioning_unscaled= pd.read_csv("output/ComparingSolversPython/Preconditioning/unscaled/PythonPsi_Bicgstab_preconditioned.csv")

#=============================================
#ROCKET: JACOBI/
#=============================================
#BICGSTAB SOLVER: Preconditioning JACOBI: SCALED
psi_python_bicgstab_JACOBI=pd.read_csv("output/ComparingSolversPython/Jacobi/scaled/PythonPsi_Bicgstab_preconditioned.csv")

#=============================================
#ROCKET: FluxForm/
#=============================================
#With Neumann BC: cos theta pm 1/2 == 0 on the Boundaries
python_poles_bicgstab = pd.read_csv("output/FluxForm/scaled/PythonPsi_BICGSTAB.csv")
python_poles_direct = pd.read_csv("output/FluxForm/scaled/PythonPsi_direct.csv")

#No Neumann BC: θ , ƛ  = 722 , 1444
python_NoNeumann_poles_bicgstab = pd.read_csv("output/FluxForm/NoNeumannScaled/PythonPsi_BICGSTAB.csv")
python_NoNeumann_poles_direct = pd.read_csv("output/FluxForm/NoNeumannScaled/PythonPsi_direct.csv")

python_rtol_scaled = pd.read_csv("output/FluxForm/rtol/PythonPsi_BICGSTAB.csv")

python_rtol_scaled = pd.read_csv("output/FluxForm/ILU_preconditioning_rtol_same/PythonPsi_BICGSTAB.csv")
python_rtol_scaled_rtolOverS= pd.read_csv("output/FluxForm/ILU_preconditioning_rtol_overS/PythonPsi_BICGSTAB.csv")

#=============================================
#ROCKET: FluxForm/FOCUS/
#=============================================
#this is scaled no neumann
ScaledInterior_UnscaledPoles= pd.read_csv("output/FluxForm/Focus/Scaled_PolesUnscaled/PythonPsi_BICGSTAB.csv")
ScaledInterior_UnscaledPoles_sp= pd.read_csv("output/FluxForm/Focus/Scaled_PolesUnscaled/PythonPsi_direct.csv")

ScaledInterior_ScaledPoles= pd.read_csv("output/FluxForm/Focus/Scaled_PolesScaled/PythonPsi_BICGSTAB.csv")
ScaledInterior_ScaledPoles_sp= pd.read_csv("output/FluxForm/Focus/Scaled_PolesScaled/PythonPsi_direct.csv")

NoPhi_Scaled = pd.read_csv("output/FluxForm/Focus/Scaled_NoPhi/PythonPsi_BICGSTAB.csv")
NoPhi_Scaled_sp= pd.read_csv("output/FluxForm/Focus/Scaled_NoPhi/PythonPsi_direct.csv")

Neumann_Scaled = pd.read_csv("output/FluxForm/Focus/Scaled_Neumann/PythonPsi_BICGSTAB.csv")
Neumann_Scaled_sp= pd.read_csv("output/FluxForm/Focus/Scaled_Neumann/PythonPsi_direct.csv")

Unscaled_NoNeumann = pd.read_csv("output/FluxForm/Focus/Unscaled/PythonPsi_BICGSTAB.csv")
Unscaled_NoNeumann_sp= pd.read_csv("output/FluxForm/Focus/Unscaled/PythonPsi_direct.csv")

Unscaled_Neumann = pd.read_csv("output/FluxForm/Focus/Unscaled_WithNeumann/PythonPsi_BICGSTAB.csv")
Unscaled_Neumann_sp= pd.read_csv("output/FluxForm/Focus/Unscaled_WithNeumann/PythonPsi_direct.csv")

Scaled_NoNeuman_NoPhi= pd.read_csv("output/FluxForm/Focus/Scaled_NoNeumann_noPhi/PythonPsi_BICGSTAB.csv")
Scaled_NoNeuman_NoPhi_sp= pd.read_csv("output/FluxForm/Focus/Scaled_NoNeumann_noPhi/PythonPsi_direct.csv")

simple_modification_PolesOnly = pd.read_csv("output/FluxForm/Focus/Final_OnlyPoles_Scaled/PythonPsi_BICGSTAB.csv") 

simple_modification_Poles_Neumann =  pd.read_csv("output/FluxForm/Focus/Final_OnlyPoles_Scaled_WithNeumann/PythonPsi_BICGSTAB.csv") 

Ld_theta_modified_neumann =  pd.read_csv("output/FluxForm/Focus/Mod_Neum_Ld_theta/PythonPsi_BICGSTAB.csv")
Ld_theta_modified_neumann_sp = pd.read_csv("output/FluxForm/Focus/Mod_Neum_Ld_theta/PythonPsi_direct.csv")

#=============================================
#ROCKET: FluxForm/Wings/
#=============================================
big_normal = pd.read_csv("output/FluxForm/Wings/Normal/PythonPsi_BICGSTAB.csv")
big_pre_normal= pd.read_csv("output/FluxForm/Wings/Normal/PythonPsi_Preconditioned.csv")
direct_normal= pd.read_csv("output/FluxForm/Wings/Normal/PythonPsi_direct.csv")

big = pd.read_csv("output/FluxForm/Wings/Coarserby2/PythonPsi_BICGSTAB.csv")
big_pre= pd.read_csv("output/FluxForm/Wings/Coarserby2/PythonPsi_Preconditioned.csv")
direct = pd.read_csv("output/FluxForm/Wings/Coarserby2/PythonPsi_direct.csv")

big_by4 = pd.read_csv("output/FluxForm/Wings/Coarserby4/PythonPsi_BICGSTAB.csv")
big_pre_by4= pd.read_csv("output/FluxForm/Wings/Coarserby4/PythonPsi_Preconditioned.csv")
direct_by4 = pd.read_csv("output/FluxForm/Wings/Coarserby4/PythonPsi_direct.csv")

#! Delta/2 on the poles:
big_Deltaby2 = pd.read_csv("output/FluxForm/Wings/delta2/PythonPsi_BICGSTAB.csv")
big_pre_Deltaby2= pd.read_csv("output/FluxForm/Wings/delta2/PythonPsi_Preconditioned.csv")
direct_Deltaby2 = pd.read_csv("output/FluxForm/Wings/delta2/PythonPsi_direct.csv")

#! 722*2:
big_Finerby2 = pd.read_csv("output/FluxForm/Wings/Finerby2/PythonPsi_BICGSTAB.csv")
big_pre_Finerby2= pd.read_csv("output/FluxForm/Wings/Finerby2/PythonPsi_Preconditioned.csv")
direct_Finerby2 = pd.read_csv("output/FluxForm/Wings/Finerby2/PythonPsi_direct.csv")

#=============================================
#ROCKET: FluxForm/BRICK/
#=============================================
big_80 = pd.read_csv("output/FluxForm/BRICK/PythonPsi_BICGSTAB.csv")
big_pre_80= pd.read_csv("output/FluxForm/BRICK/PythonPsi_Preconditioned.csv")
direct_80= pd.read_csv("output/FluxForm/BRICK/PythonPsi_direct.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_varying
#=============================================
psi_RHSv_Fortran_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_varying/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying/PythonPsi_ADI.csv")
psi_RHSv_Fortran_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying/PythonPsi_direct.csv")
psi_RHSv_Fortran_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_varying_stab
#=============================================
psi_RHSv_stab_Fortran_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab/PythonPsi_BICGSTAB.csv")
psi_RHSv_stab_Fortran_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab/PythonPsi_ADI.csv")
psi_RHSv_stab_Fortran_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab/PythonPsi_direct.csv")
psi_RHSv_stab_Fortran_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_varying_stab_multiply
#=============================================
psi_RHSv_stab_multiply_Fortran_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab_multiply/PythonPsi_BICGSTAB.csv")
psi_RHSv_stab_multiply_Fortran_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab_multiply/PythonPsi_ADI.csv")
psi_RHSv_stab_multiply_Fortran_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab_multiply/PythonPsi_direct.csv")
psi_RHSv_stab_multiply_Fortran_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_stab_multiply/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_varying_scaled/
#=============================================
psi_RHSv_Fortran_Scaled_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_Scaled_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled/PythonPsi_ADI.csv")
psi_RHSv_Fortran_Scaled_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled/PythonPsi_direct.csv")
psi_RHSv_Fortran_Scaled_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_varying_scaled_stabilization/
#=============================================
psi_RHSv_Fortran_Scaled_stabilization_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled_stabilization/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_Scaled_stabilization_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled_stabilization/PythonPsi_ADI.csv")
psi_RHSv_Fortran_Scaled_stabilization_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled_stabilization/PythonPsi_direct.csv")
psi_RHSv_Fortran_Scaled_stabilization_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_scaled_stabilization/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_varying_fullscaled_stabilization/
#=============================================
psi_RHSv_Fortran_fullScaled_stabilization_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_fullscaled_stabilization/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_fullScaled_stabilization_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_fullscaled_stabilization/PythonPsi_ADI.csv")
psi_RHSv_Fortran_fullScaled_stabilization_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_fullscaled_stabilization/PythonPsi_direct.csv")
psi_RHSv_Fortran_fullScaled_stabilization_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_varying_fullscaled_stabilization/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_constant
#=============================================
# psi_RHSc_Fortran_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_constant/PythonPsi_BICGSTAB.csv")
# psi_RHSc_Fortran_adi= pd.read_csv("output/1_14AM/FortranGrid/RHS_constant/PythonPsi_ADI.csv")
psi_RHSc_Fortran_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_constant/PythonPsi_direct.csv")
# psi_RHSc_Fortran_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_constant/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/FortranGrid/RHS_constant_phi
#=============================================
psi_RHSc_phi_Fortran_bicg = pd.read_csv("output/1_14AM/FortranGrid/RHS_constant_phi/PythonPsi_BICGSTAB.csv")
psi_RHSc_phi_Fortran_direct= pd.read_csv("output/1_14AM/FortranGrid/RHS_constant_phi/PythonPsi_direct.csv")
psi_RHSc_phi_Fortran_cg= pd.read_csv("output/1_14AM/FortranGrid/RHS_constant_phi/PythonPsi_CG.csv")

#!Ld=100 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_bicg = pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_adi= pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_direct= pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_cg= pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled/PythonPsi_CG.csv")
#!Ld=100 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled_stabilized/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_stabilized_bicg = pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled_stabilized/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_stabilized_adi= pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled_stabilized/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_stabilized_direct= pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled_stabilized/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_stabilized_cg= pd.read_csv("output/1_14AM/FortranGrid/ScreenedPoisson_FG_All_scaled_stabilized/PythonPsi_CG.csv")
#%%
#!Testing another Ld :

#!Ld=1000 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_1000/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld1000_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_1000/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld1000_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_1000/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld1000_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_1000/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld1000_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_1000/PythonPsi_CG.csv")

#!Ld=100000 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_100000/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld100000_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_100000/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld100000_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_100000/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld100000_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_100000/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld100000_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_100000/PythonPsi_CG.csv")

#!Ld= R m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_R/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_R_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_R/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_R_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_R/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_R_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_R/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_R_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_R/PythonPsi_CG.csv")

#!Ld= R/10 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover10/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover10_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover10/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover10_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover10/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover10_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover10/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover10_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover10/PythonPsi_CG.csv")

#!Ld= R/100 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover100/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover100_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover100/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover100_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover100/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover100_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover100/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover100_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover100/PythonPsi_CG.csv")

#!Ld= R/1000 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover1000/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover1000_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover1000/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover1000_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover1000/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover1000_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover1000/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Rover1000_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/Ld_Rover1000/PythonPsi_CG.csv")

#%%
#!VArying Ld Screened- scaled sabulized:
#!Ld= R/1000 m
#=============================================
#ROCKET: 1_14AM/FortranGrid/TestingSphericalHarmonics/VariableScreening/
#=============================================
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Varying_bicg = pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/VariableScreening/PythonPsi_BICGSTAB.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Varying_adi= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/VariableScreening/PythonPsi_ADI.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Varying_direct= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/VariableScreening/PythonPsi_direct.csv")
psi_RHSv_Fortran_ScreenedPoisson_Scaled_Stabilized_Ld_Varying_cg= pd.read_csv("output/1_14AM/FortranGrid/TestingSphericalHarmonics/VariableScreening/PythonPsi_CG.csv")

#%%
#=============================================
#ROCKET: 1_14AM/PolesGrid/RHS_varying
#=============================================
psi_RHSv_Poles_bicg = pd.read_csv("output/1_14AM/PolesGrid/RHS_varying/PythonPsi_BICGSTAB.csv")
psi_RHSv_Poles_adi= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying/PythonPsi_ADI.csv")
psi_RHSv_Poles_direct= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying/PythonPsi_direct.csv")
psi_RHSv_Poles_cg= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/PolesGrid/RHS_varying_FortranGrid
#=============================================
psi_RHSv_Poles_FortranGrid_bicg = pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_FortranGrid/PythonPsi_BICGSTAB.csv")
psi_RHSv_Poles_FortranGrid_adi= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_FortranGrid/PythonPsi_ADI.csv")
psi_RHSv_Poles_FortranGrid_direct= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_FortranGrid/PythonPsi_direct.csv")
psi_RHSv_Poles_FortranGrid_cg= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_FortranGrid/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/PolesGrid/RHS_varying_stabilization
#=============================================
psi_RHSv_Poles_stabilization_bicg = pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_stabilization/PythonPsi_BICGSTAB.csv")
psi_RHSv_Poles_stabilization_adi= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_stabilization/PythonPsi_ADI.csv")
psi_RHSv_Poles_stabilization_direct= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_stabilization/PythonPsi_direct.csv")
psi_RHSv_Poles_stabilization_cg= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_stabilization/PythonPsi_CG.csv")

#=============================================
#ROCKET: 1_14AM/PolesGrid/RHS_varying_More_stabilization
#=============================================
psi_RHSv_Poles_More_stabilization_bicg = pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_More_stabilization/PythonPsi_BICGSTAB.csv")
psi_RHSv_Poles_More_stabilization_adi= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_More_stabilization/PythonPsi_ADI.csv")
psi_RHSv_Poles_More_stabilization_direct= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_More_stabilization/PythonPsi_direct.csv")
psi_RHSv_Poles_More_stabilization_cg= pd.read_csv("output/1_14AM/PolesGrid/RHS_varying_More_stabilization/PythonPsi_CG.csv")
#%%
#=============================================
#ROCKET: Fortran/
#=============================================
#going back to fortran and adding the updates: this is the raw case scaled no preconditioning no stabilizers
Raw_Fortran_Psi= pd.read_csv("output/Fortran/CaseZero/Unscreened_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Raw_Unscaled_Fortran_Psi= pd.read_csv("output/Fortran/CaseZero/Unscaled_Unscreened_NoPreconditioner_NoStabilizer/PsiOcean.csv")

#!CASEONE- stabilizer:
Scaled_Stab_Fortran_Psi= pd.read_csv("output/Fortran/CaseOne/Unscreened_NoPreconditioner_Stabilizer/PsiOcean.csv")

#!CaseTwo - stabilizer with screening
#!Ld=100

Scaled_Stab_Ld_100_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_100_NoPreconditioner_Stabilizer/PsiOcean.csv")
#! no stabilizer for this case:
Scaled_NoStab_Ld_100_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_100_NoPreconditioner_NoStabilizer/PsiOcean.csv")

#!Ld=1000
Scaled_NoStab_Ld_1000_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_1000_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Scaled_Stab_Ld_1000_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_1000_NoPreconditioner_Stabilizer/PsiOcean.csv")

#!Ld=100,000
Scaled_NoStab_Ld_100000_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_100000_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Scaled_Stab_Ld_100000_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_100000_NoPreconditioner_Stabilizer/PsiOcean.csv")

#!Ld=R
Scaled_NoStab_Ld_R_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_R_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Scaled_Stab_Ld_R_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Screened_Ld_R_NoPreconditioner_Stabilizer/PsiOcean.csv")

Unscaled_NoStab_Ld_R_Fortran_Psi = pd.read_csv("output/Fortran/CaseTwo/Unscaled_Screened_Ld_R_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Unscaled_Stab_Ld_R_Fortran_Psi = pd.read_csv("output/Fortran/CaseTwo/Unscaled_Screened_Ld_R_NoPreconditioner_Stabilizer/PsiOcean.csv")

#!Ld=R/10
Scaled_Stab_Ld_Rover10_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Scaled_Screened_Ld_Rover10_NoPreconditioner_Stabilizer/PsiOcean.csv")
Scaled_NoStab_Ld_Rover10_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Scaled_Screened_Ld_Rover10_NoPreconditioner_NoStabilizer/PsiOcean.csv")

# #!Ld=R/100
Scaled_Stab_Ld_Rover100_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Scaled_Screened_Ld_Rover100_NoPreconditioner_Stabilizer/PsiOcean.csv")
Scaled_NoStab_Ld_Rover100_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Scaled_Screened_Ld_Rover100_NoPreconditioner_NoStabilizer/PsiOcean.csv")

#!Ld=R/1000
Scaled_Stab_Ld_Rover1000_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Scaled_Screened_Ld_Rover1000_NoPreconditioner_Stabilizer/PsiOcean.csv")
Scaled_NoStab_Ld_Rover1000_Fortran_Psi= pd.read_csv("output/Fortran/CaseTwo/Scaled_Screened_Ld_Rover1000_NoPreconditioner_NoStabilizer/PsiOcean.csv")

#!CaseStupid
Stupid_Fortran_Psi= pd.read_csv("output/Fortran/CaseStupidMaybe/PhiDependence/PsiOcean.csv")
Stupid_stab_Fortran_Psi= pd.read_csv("output/Fortran/CaseStupidMaybe/PhiDependence_stab/PsiOcean.csv")

#!Case 3: varying screening term
Scaled_NoStab_Ld_theta_Fortran_Psi= pd.read_csv("output/Fortran/CaseThree/Scaled_Screened_Ld_theta_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Scaled_Stab_Ld_theta_Fortran_Psi= pd.read_csv("output/Fortran/CaseThree/Scaled_Screened_Ld_theta_NoPreconditioner_Stabilizer/PsiOcean.csv")

Unscaled_NoStab_Ld_theta_Fortran_Psi= pd.read_csv("output/Fortran/CaseThree/Unscaled_Screened_Ld_theta_NoPreconditioner_NoStabilizer/PsiOcean.csv")
Unscaled_Stab_Ld_theta_Fortran_Psi= pd.read_csv("output/Fortran/CaseThree/Unscaled_Screened_Ld_theta_NoPreconditioner_Stabilizer/PsiOcean.csv")
#!Case four scaled normalized ld = R/1000
CaseFour_Fortran_Psi= pd.read_csv("output/Fortran/CaseFour/PsiOcean.csv")
CaseFourA_Fortran_Psi= pd.read_csv("output/Fortran/CaseFourA/PsiOcean.csv")
#%%
#!Correcte Section:
#! the stabilizer for unscreened -fortran
Updated_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Unscreened/PsiOcean.csv")

#! Ld= 100
Updated_Ld_100_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_100/PsiOcean.csv")

Updated_Ld_1000_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_1000/PsiOcean.csv")
Updated_Ld_100000_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_100000/PsiOcean.csv")
Updated_Ld_R_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_R/PsiOcean.csv")
Updated_Ld_Rover10_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_Rover10/PsiOcean.csv")
Updated_Ld_Rover100_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_Rover100/PsiOcean.csv")

Updated_Ld_Rover1000_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_Rover1000/PsiOcean.csv")
Updated_Ld_theta_Fortran_Psi= pd.read_csv("output/Fortran/ConditionUpdated_Ld_theta/PsiOcean.csv")
#%%
#=============================================
#ROCKET: Fortran/Lambda/unscreened/PsiOcean.csv
#=============================================
Unscreened_Lambda_Fortran_Psi= pd.read_csv("output/Fortran/Lambda/unscreened/PsiOcean.csv")

#!PYthon unscreened - lambda dependence:
psi_Python_bicg = pd.read_csv("output/Fortran/LambdaPython/unscreened/PythonPsi_BICGSTAB.csv")
psi_Python_adi= pd.read_csv("output/Fortran/LambdaPython/unscreened/PythonPsi_ADI.csv")
psi_Python_direct= pd.read_csv("output/Fortran/LambdaPython/unscreened/PythonPsi_direct.csv")
psi_Python_cg= pd.read_csv("output/Fortran/LambdaPython/unscreened/PythonPsi_CG.csv")

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.config import (nLat, nPhi,Lat, sinLat0, dLat, GRAVITY ,D, sinLat, OMEGA,R, ScaleEquation, neighborhandling, includeFreeSurface,Ld0, OUTPUT_PATH)

#!Full equation with depth:
# Psi_Paraview_Fortran_unnormalized= pd.read_csv("output/Fortran/EquationOne/Unnormalized/PsiOcean.csv")
# Psi_Paraview_Fortran_normalized= pd.read_csv("output/Fortran/EquationOne/Normalized/PsiOcean.csv")

# # #!PYTHON CODE SAME AS FORTRAN - ADDING 1E-12
# Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_Noaddition_diagonal/PythonPsi_ADI.csv")
# Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_Noaddition_diagonal/PythonPsi_BICGSTAB.csv")
# Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_Noaddition_diagonal/PythonPsi_CG.csv")
# Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_Noaddition_diagonal/PythonPsi_direct.csv")

# #!PYTHON CODE SAME AS FORTRAN - ADDING 1E-12
# Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_addition_diagonal/PythonPsi_ADI.csv")
# Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_addition_diagonal/PythonPsi_BICGSTAB.csv")
# Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_addition_diagonal/PythonPsi_CG.csv")
# Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/Python/unnormalized_addition_diagonal/PythonPsi_direct.csv")

# # #!PYTHON CODE + ADDING 1E-12  + normalize
# Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/Python/normalized_addition_diagonal/PythonPsi_ADI.csv")
# Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/Python/normalized_addition_diagonal/PythonPsi_BICGSTAB.csv")
# Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/Python/normalized_addition_diagonal/PythonPsi_CG.csv")
# Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/Python/normalized_addition_diagonal/PythonPsi_direct.csv")

# # #!PYTHON CODE + normalize
# Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/Python/normalized_Noaddition_diagonal/PythonPsi_ADI.csv")
# Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/Python/normalized_Noaddition_diagonal/PythonPsi_BICGSTAB.csv")
# Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/Python/normalized_Noaddition_diagonal/PythonPsi_CG.csv")
# Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/Python/normalized_Noaddition_diagonal/PythonPsi_direct.csv")

# #!HODA NEIGHBOR - ADDING 1E-12
# Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/PythonHoda/unnormalized_addition_diagonal/PythonPsi_ADI.csv")
# Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/PythonHoda/unnormalized_addition_diagonal/PythonPsi_BICGSTAB.csv")
# Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/PythonHoda/unnormalized_addition_diagonal/PythonPsi_CG.csv")
# Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/PythonHoda/unnormalized_addition_diagonal/PythonPsi_direct.csv")
#!HODA NEIGHBOR - ADDING 1E-12 + normalized
# Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/PythonHoda/normalized_addition_diagonal/PythonPsi_ADI.csv")
# Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/PythonHoda/normalized_addition_diagonal/PythonPsi_BICGSTAB.csv")
# Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/PythonHoda/normalized_addition_diagonal/PythonPsi_CG.csv")
# Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/PythonHoda/normalized_addition_diagonal/PythonPsi_direct.csv")




# #!PYTHON CODE + MITGCM 1 WHEN LAND 
Psi_Python_ADI= pd.read_csv("output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/PythonPsi_ADI.csv")
Psi_Python_BICGSTAB= pd.read_csv("output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/PythonPsi_BICGSTAB.csv")
Psi_Python_CG= pd.read_csv("output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/PythonPsi_CG.csv")
Psi_Python_direct= pd.read_csv("output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/PythonPsi_direct.csv")
#%%
#=============================================
#ROCKET: PLOTTING PSI
#=============================================
# iLat1, psi1 = Unscreened_Lambda_Fortran_Psi.iloc[:, 0], Unscreened_Lambda_Fortran_Psi.iloc[:, 2]

# iLat2, psi2 = Psi_Paraview_Fortran_normalized.iloc[:, 0], Psi_Paraview_Fortran_normalized.iloc[:, 2]
# iLat3, psi3 = Psi_Paraview_Fortran_unnormalized.iloc[:, 0], Psi_Paraview_Fortran_unnormalized.iloc[:, 2]

iLat3, psi3 = Psi_Python_BICGSTAB.iloc[:, 0], Psi_Python_BICGSTAB.iloc[:, 2]
iLat4, psi4 = Psi_Python_ADI.iloc[:,0], Psi_Python_ADI.iloc[:, 2]
iLat5, psi5 = Psi_Python_direct.iloc[:,0], Psi_Python_direct.iloc[:, 2]
iLat6, psi6 = Psi_Python_CG.iloc[:,0], Psi_Python_CG.iloc[:, 2]

# plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Direct') 

nLat = len(iLat3)  # or your grid size
lat_deg = np.linspace(-90, 90, nLat) 
# plt.plot(lat_deg, psi2, color='#D62728', linewidth=2, label='BICGSTAB OMP Paraview')
# plt.plot(lat_deg, psi3, color='#7F7F7F', linewidth=2, label='BICGSTAB OMP Paraview')

# plt.plot(lat_deg, psi3, color='#D62728', linewidth=2, label='BICGSTAB')
# plt.plot(lat_deg, psi4, color='#2CA02C', linestyle='--',  linewidth=3, label='ADI') 
plt.plot(lat_deg, psi5, color='#7F7F7F', linestyle='--', linewidth=3, label='DIRECT')   
# plt.plot(lat_deg, psi6, color='#1F77B4', linestyle='--', linewidth=3, label='CG')  

plt.xlabel("Latitude (°)")
plt.ylabel("Psi")
plt.title("PCG") #PCG
plt.legend()
plt.grid(True)
plt.show()
#%%
#!for comaprsion:
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

# Top panel
ax1.plot(lat_deg, psi3, color='#D62728', linewidth=3, label='BICGSTAB')
ax1.plot(lat_deg, psi4, color='#2CA02C', linewidth=3, label='ADI')
ax1.plot(lat_deg, psi5, color='#7F7F7F', linewidth=3, label='DIRECT')
ax1.plot(lat_deg, psi6, color='#1F77B4', linewidth=3, linestyle='--', label='CG')

# ax1.plot(lat_deg, psi2, color='#D62728', linewidth=3, label='Unnormalized')
# ax1.plot(lat_deg, psi3, color='#2CA02C', linewidth=3, label='Normalized')
# ax1.plot(lat_deg, psi6, color='#FF8C00', linewidth=3, markeredgewidth=3,linestyle='--',label='CG')

ax1.set_ylabel('Ψ (Streamfunction)', fontsize=14)
ax1.set_title('Solver Comparison: Streamfunction Ψ (top)', fontsize=16, pad=20)

# ax1.set_title('Solver Comparison: Streamfunction Ψ (top) and Errors vs DIRECT (bottom)', fontsize=16, pad=20)
ax1.legend(fontsize=12, loc='upper right', frameon=True, fancybox=True)
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# === Bottom panel: Differences relative to DIRECT (magnifies tiny errors) ===
diff_bicg = psi3 - psi5
diff_adi  = psi4 - psi5
diff_cg   = psi6 - psi5

# diff_bicg = psi3 - psi2
# ax2.plot(lat_deg, diff_bicg, color='#D62728', linewidth=3, label='Normalized-Unnormalized')

# Bottom panel (same colors!)
ax2.plot(lat_deg, diff_bicg, color='#D62728', linewidth=3, label='BICGSTAB - DIRECT')
ax2.plot(lat_deg, diff_adi,  color='#2CA02C', linewidth=3, label='ADI - DIRECT')
ax2.plot(lat_deg, diff_cg,   color='#1F77B4', linewidth=3, linestyle='--', label='CG - DIRECT')
ax2.axhline(0, color='gray', linewidth=1, linestyle='-') #zero line for reference 
ax2.set_xlabel('Latitude (°)', fontsize=14) 
ax2.set_ylabel('ΔΨ (Error)', fontsize=14) 
ax2.legend(fontsize=11, loc='upper right') 
ax2.grid(True, linestyle=':', alpha=0.6) 
ax2.spines['top'].set_visible(False) 
ax2.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
#%%
# MITGCM
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
# 1. Define your global grid size (matches SIZE.h)
nx, ny = 1440, 720

# 2. Find and merge all 8 files
file_pattern ="output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/PythonPsi_BICGSTAB.csv"
# "output/Fortran/EquationOne/Normalized/PsiOcean.csv"
# "output/Fortran/EquationOne/PythonHoda/MITGCM_one_diag_Land/PythonPsi_BICGSTAB.csv"
# 'output/Fortran/EquationOne/Normalized/PsiOcean.csv'
files = glob.glob(file_pattern)

# Load all CSVs into one big DataFrame
df_list = [pd.read_csv(f) for f in files]
full_df = pd.concat(df_list, ignore_index=True)
psi_global = np.zeros((ny, nx))

for _, row in full_df.iterrows():
    # Subtract 1 because Fortran indices were 1-based
    ilat = int(row['iLat']) #-1
    ilon = int(row['iLon']) #-1
    # Boundary check to prevent crashes
    if 0 <= ilat < ny and 0 <= ilon < nx:
        psi_global[ilat, ilon] = row['psi']

# 4. Plot the results
plt.figure(figsize=(12, 6))
# Use pcolormesh or imshow to show the global map
plt.imshow(psi_global, origin='lower', extent=[0, 360, -90, 90], cmap='RdBu_r')
plt.colorbar(label='Psi Value')
plt.title('Global Fortran solver output')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
# plt.savefig('global_psi_plot.png')
plt.show()

#!FOR PYTHON
psi_global = np.full((ny, nx), np.nan)

for _, row in full_df.iterrows():
    # MANDATORY: Subtract 1 for 0-based Python indexing
    ilat = int(row['iLat']) - 1
    ilon = int(row['iLon']) - 1
    
    if 0 <= ilat < ny and 0 <= ilon < nx:
        psi_global[ilat, ilon] = row['psi']

# 2. Plot with robust scaling
plt.figure(figsize=(12, 6))
# vmin/vmax ensure the colorbar is centered on 0 if using RdBu
mag = np.nanmax(np.abs(psi_global)) 
im = plt.imshow(psi_global, origin='lower', extent=[0, 360, -90, 90], 
                cmap='RdBu_r', vmin=-mag, vmax=mag)
plt.title('Global Python solver output')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.colorbar(im, label='Psi Value')


#%%
#BUGTRACK: READING RESIDUALS
#BUGTRACK: =============================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#=============================================
#BUGTRACK: iterationConstraint/
#=============================================

#BICGSTAB SOLVER: SCALED,  iter 1,000
res = pd.read_csv("output/ComparingSolversPython/iterationConstraint/iter1000/scaled/bicgstab_residuals.csv")
iterations,residuals = res['Iteration'],res['Residual']

##BICGSTAB SOLVER: SCALED,  iter 5,000
res_5000 = pd.read_csv("output/ComparingSolversPython/iterationConstraint/scaled/iter5000/bicgstab_residuals.csv")
iterations_5000,residuals_5000 = res_5000['Iteration'],res_5000['Residual']

#BICGSTAB SOLVER: SCALED,  iter 10,000
res_10000 = pd.read_csv("output/ComparingSolversPython/iterationConstraint/scaled/iter10000/bicgstab_residuals.csv")
iterations_10000,residuals_10000 = res_10000['Iteration'],res_10000['Residual']

#BICGSTAB SOLVER: UNSCALED,  iter 1000
res_unscaled = pd.read_csv("output/ComparingSolversPython/iterationConstraint/iter1000/unscaled/bicgstab_residuals.csv")
iterations_unscaled,residuals_unscaled = res_unscaled['Iteration'], res_unscaled['Residual']

#=============================================
#BUGTRACK: Preconditioning/
#=============================================
#BICGSTAB SOLVER: Preconditioning ILU: SCALED
res_preconditioned= pd.read_csv("output/ComparingSolversPython/Preconditioning/scaled/bicgstab_residuals.csv")
iterations_preconditioned,residuals_preconditioned = res_preconditioned['Iteration'],res_preconditioned['Residual']

#BICGSTAB SOLVER: Preconditioning ILU: UNSCALED
# res_preconditioned_unscaled= pd.read_csv("output/ComparingSolversPython/Preconditioning/unscaled/bicgstab_residuals.csv")
# iterations_preconditioned_unscaled,residuals_preconditioned_unscaled = res_preconditioned_unscaled['Iteration'],res_preconditioned_unscaled['Residual']
# print('Min res, Max res: ' ,min(residuals_preconditioned_unscaled),',', max(residuals_preconditioned_unscaled))

#=============================================
#BUGTRACK: JACOBI/
#=============================================

res_jacobi_scaled= pd.read_csv("output/ComparingSolversPython/Jacobi/scaled/bicgstab_residuals.csv")
iter_jacobi_scaled,residuals_jacobi_scaled = res_jacobi_scaled['Iteration'],res_jacobi_scaled['Residual']
print('Min res, Max res: ' ,min(residuals_jacobi_scaled),',', max(residuals_jacobi_scaled))

#=============================================
#BUGTRACK: FluxForm/
#=============================================

#With Neumann BC: cos theta pm 1/2 == 0 on the Boundaries
python_poles_residuals = pd.read_csv("output/FluxForm/scaled/BICGSTAB_residuals.csv")
iter_poles_bicgstab,residuals_poles_bicgstab = python_poles_residuals['Iteration'],python_poles_residuals['Residual']
b_poles = pd.read_csv("output/FluxForm/scaled/b.csv")

#No Neumann BC
python_NoNeumann_poles_residuals = pd.read_csv("output/FluxForm/NoNeumannScaled/BICGSTAB_residuals.csv")
iter_NoNeumann_poles_bicgstab,residuals_NoNeumann_poles_bicgstab = python_NoNeumann_poles_residuals['Iteration'],python_NoNeumann_poles_residuals['Residual']

#=============================================
#BUGTRACK: FluxForm/FOCUS/
#=============================================
ScaledInterior_UnscaledPoles= pd.read_csv("output/FluxForm/Focus/Scaled_PolesUnscaled/BICGSTAB_residuals.csv")
iter_ScaledInterior_UnscaledPoles, res_ScaledInterior_UnscaledPoles = ScaledInterior_UnscaledPoles['Iteration'], ScaledInterior_UnscaledPoles['Residual']

ScaledInterior_ScaledPoles= pd.read_csv("output/FluxForm/Focus/Scaled_PolesScaled/BICGSTAB_residuals.csv")
iter_ScaledInterior_ScaledPoles, res_ScaledInterior_ScaledPoles = ScaledInterior_ScaledPoles['Iteration'], ScaledInterior_ScaledPoles['Residual']

NoPhi= pd.read_csv("output/FluxForm/Focus/Scaled_NoPhi/BICGSTAB_residuals.csv")
iter_NoPhi, res_NoPhi = NoPhi['Iteration'], NoPhi['Residual']

Neumann= pd.read_csv("output/FluxForm/Focus/Scaled_Neumann/BICGSTAB_residuals.csv")
iter_Neumann, res_Neumann = Neumann['Iteration'], Neumann['Residual']

Unscaled_Neumann= pd.read_csv("output/FluxForm/Focus/Unscaled_WithNeumann/BICGSTAB_residuals.csv")
iter_Unscaled_Neumann, res_Unscaled_Neumann = Unscaled_Neumann['Iteration'], Unscaled_Neumann['Residual']

Unscaled_No_Neumann= pd.read_csv("output/FluxForm/Focus/Unscaled/BICGSTAB_residuals.csv")
iter_Unscaled_No_Neumann, res_Unscaled_No_Neumann = Unscaled_No_Neumann['Iteration'], Unscaled_No_Neumann['Residual']

#=============================================
#BUGTRACK: FluxForm/BRICK
#=============================================
brick_80= pd.read_csv("output/FluxForm/BRICK/BICGSTAB_residuals.csv")
iter_brick_80, res_brick_80 = brick_80['Iteration'], brick_80['Residual']

#%%
#=============================================
#BUGTRACK: PLOTTING RESIDUALS
#=============================================
plt.figure()
plt.semilogy(iterations, residuals,color='#00F5FF',lw=1,label='Scaled BICGSTAB solver')
# residuals_unscaled
# plt.semilogy(iterations_unscaled, residuals_unscaled,color='#00F5FF',lw=1,label='Unscaled BICGSTAB solver')

# plt.semilogy(iter_Unscaled_No_Neumann, res_Unscaled_No_Neumann,color='#FF5F1F',lw=1, label='Unscaled Python No Neumann')

plt.semilogy(iter_brick_80, res_brick_80,color='royalblue',lw=1, label='Modified Bicgstab ')

plt.xlabel('Iteration')
plt.ylabel('Residual norm')
plt.title('BiCGSTAB Residuals')
plt.legend()
plt.grid(True)
plt.show()

#%% COMMENT: conditioning analysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Scaled: iterations 1000
kappa = pd.read_csv("output/ComparingSolversPython/scaled/conditioning_vs_lat.csv")
lat,k = kappa['Latitude(deg)'],kappa['Conditioning']

# kappa_unscaled = pd.read_csv("output/ComparingSolversPython/unscaled/conditioning_vs_lat.csv")
# lat_unscaled,k_unscaled = kappa_unscaled['Latitude(deg)'],kappa_unscaled['Conditioning']

plt.figure()
plt.semilogy(lat, k,color='#00F5FF',lw=3,label='Scaled python 1000 iters')
# plt.semilogy(lat_unscaled, k_unscaled,color='#FF5F1F',lw=1,label='Unscaled python 1000 iters')

plt.xlabel('latitude in deg')
plt.ylabel('k')
plt.title("Conditioning Term k vs Latitude")
plt.legend()
plt.grid(True)
plt.show()

#%%
#=============================================
#ROCKET: ../MeetingWed03/Case 8 scaled
#=============================================
psi_python_8_scaled = pd.read_csv("../MeetingWed03/8/scaled/PythonPsi.csv")

# psi_python_8_scaled = pd.read_csv("Half_PythonPsi.csv")

psi_ocean_scaled = pd.read_csv("../MeetingWed03/8/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/8/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python_8_scaled.iloc[:, 0], psi_python_8_scaled.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#royalblue', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()
#=============================================
#ROCKET: ../MeetingWed03/Case 8 unscaled 
#=============================================
psi_python_8_unscaled = pd.read_csv("../MeetingWed03/8/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/8/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/8/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python_8_unscaled.iloc[:, 0], psi_python_8_unscaled.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()
#%%
import pandas as pd
import numpy as np
#!COMMENT: THIS FUNCTION IS TO COMPARE TWO CSV FILES AND MAKE SURE THEY PRODUCE THE SAME RESULTS.
def compare_csv(file1, file2, tol=1e-12):
# Read CSVs
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    
    if df1.shape != df2.shape:
        print("Warning: CSV files have different shapes")
    
    # Compare numeric values
    numeric_cols = df1.select_dtypes(include=[np.number]).columns
    diff_mask = (df1[numeric_cols] - df2[numeric_cols]).abs() > tol
    
    # Keep rows where there is at least one difference
    differences = df1[diff_mask.any(axis=1)].copy()
    
    if differences.empty:
        print("No differences found (within tolerance)")
    else:
        print(f"Differences found in {len(differences)} rows")
        # Optionally, add columns showing differences
        for col in numeric_cols:
            differences[col + "_diff"] = df1[col] - df2[col]
    
    return differences

# Example usage
diffs = compare_csv("../MeetingWed03/8/PythonPsi.csv", "Half_PythonPsi.csv")
print(diffs)

# %%
#=============================================
#ROCKET: ../MeetingWed03/Case 7/ Ld = 100
#=============================================
#! 7. Case Ld =100 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 8
psi_python = pd.read_csv("../MeetingWed03/7/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/7/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/7/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 7")
plt.legend()
plt.grid(True)
plt.show()

# psi_python = pd.read_csv("../MeetingWed03/7/unscaled/PythonPsi.csv")
# psi_ocean_unscaled = pd.read_csv("../MeetingWed03/7/unscaled/PsiOcean.csv")
# psi_spherical_unscaled = pd.read_csv("../MeetingWed03/7/unscaled/PsiSpherical.csv")

# iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
# iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
# iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# # Plot ψ vs iLat for each dataset
# plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
# plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
# plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

# plt.xlabel("iLat")
# plt.ylabel("Psi")
# plt.title("psi vs iLat UNscaled Fortran codes")
# plt.legend()
# plt.grid(True)
# plt.show()

# %%
#=============================================
#ROCKET: ../MeetingWed03/Case 6/ Ld=1000
#=============================================
#! 6. Case Ld =1000
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 
# psi_python = pd.read_csv("../MeetingWed03/6/scaled/PythonPsi.csv")
# psi_ocean_scaled = pd.read_csv("../MeetingWed03/6/scaled/PsiOcean.csv")
# psi_spherical_scaled = pd.read_csv("../MeetingWed03/6/scaled/PsiSpherical.csv")
psi_python = pd.read_csv("../MeetingWed03/6/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/6/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/6/PsiSpherical.csv")
iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 6")
plt.legend()
plt.grid(True)
plt.show()

psi_python = pd.read_csv("../MeetingWed03/6/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/6/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/6/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()

# %%
#=============================================
#ROCKET: ../MeetingWed03/Case 5 scaled
#=============================================
#! 5. Case Ld = R/1000
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 
psi_python = pd.read_csv("../MeetingWed03/5/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/5/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/5/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 5")
plt.legend()
plt.grid(True)
plt.show()

psi_python = pd.read_csv("../MeetingWed03/5/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/5/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/5/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()

# %%
#=============================================
#ROCKET: ../MeetingWed03/Case 4 scaled
#=============================================
#! Case 4 Ld = R/100
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 
psi_python = pd.read_csv("../MeetingWed03/4/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/4/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/4/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 4")
plt.legend()
plt.grid(True)
plt.show()

psi_python = pd.read_csv("../MeetingWed03/4/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/4/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/4/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()

#%%
#=============================================
#ROCKET: ../MeetingWed03/Case3
#=============================================
#! Case 3 Ld = R/10

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 
psi_python = pd.read_csv("../MeetingWed03/3/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/3/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/3/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 3")
plt.legend()
plt.grid(True)
plt.show()

psi_python = pd.read_csv("../MeetingWed03/3/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/3/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/3/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()
#%%
#=============================================
#ROCKET: ../MeetingWed03/Case2
#=============================================
#! Case 2 : Ld = R

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 
psi_python = pd.read_csv("../MeetingWed03/2/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/2/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/2/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 2")
plt.legend()
plt.grid(True)
plt.show()

psi_python = pd.read_csv("../MeetingWed03/2/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/2/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/2/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()

#%%
#=============================================
#ROCKET: ../MeetingWed03/Case1
#=============================================
#! Case 1 Ld = infinity

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 
psi_python = pd.read_csv("../MeetingWed03/1/scaled/PythonPsi.csv")
psi_ocean_scaled = pd.read_csv("../MeetingWed03/1/scaled/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/1/scaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='--',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='--', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran codes Case 1")
plt.legend()
plt.grid(True)
plt.show()

psi_python = pd.read_csv("../MeetingWed03/1/unscaled/PythonPsi.csv")
psi_ocean_unscaled = pd.read_csv("../MeetingWed03/1/unscaled/PsiOcean.csv")
psi_spherical_unscaled = pd.read_csv("../MeetingWed03/1/unscaled/PsiSpherical.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]
iLat2, psi2 = psi_ocean_unscaled.iloc[:, 0], psi_ocean_unscaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_unscaled.iloc[:, 0], psi_spherical_unscaled.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat UNscaled Fortran codes")
plt.legend()
plt.grid(True)
plt.show()

# %%
#=============================================
#ROCKET: ../MeetingWed03/tolBlueBerries
#=============================================
#!tolBlueBerries:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#ROCKET: SCALED CASE 8
psi_ocean_scaled = pd.read_csv("../MeetingWed03/tolBlueBerries/PsiOcean.csv")
psi_spherical_scaled = pd.read_csv("../MeetingWed03/tolBlueBerries/PsiSpherical.csv")

iLat2, psi2 = psi_ocean_scaled.iloc[:, 0], psi_ocean_scaled.iloc[:, 2]
iLat3, psi3 = psi_spherical_scaled.iloc[:, 0], psi_spherical_scaled.iloc[:, 2]

psi_python = pd.read_csv("../MeetingWed03/8/PythonPsi.csv")

iLat1, psi1 = psi_python.iloc[:, 0], psi_python.iloc[:, 2]

# Plot ψ vs iLat for each dataset
plt.plot(iLat1, psi1, color='royalblue', linestyle='solid', linewidth=2, label='Python') 
plt.plot(iLat2, psi2, color='#FF5F1F', linestyle='-',  linewidth=2, label='Fortran OceanGrid')
plt.plot(iLat3, psi3, color='#00F5FF', linestyle='-', linewidth=1, label='Fortran SphericalGrid')

plt.xlabel("iLat")
plt.ylabel("Psi")
plt.title("psi vs iLat Scaled Fortran - old unadjusted tol - ocean adjusted")
plt.legend()
plt.grid(True)
plt.show()

# %%
#=============================================
#ROCKET: COEFFICIENT MATRIX
#=============================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#!COMMENT:  i am plotting AC, b, solution in the three codes for Case8 - scaled before changing the solver to make sure everything is implemented correctly.

#!scaled python 
AC_python = pd.read_csv("output/scaled_py/AC.csv")
iLat_AC_python = AC_python.iloc[:, 0]
iLon_AC_python = AC_python.iloc[:, 1]
value_AC_python= AC_python.iloc[:, 2]

b_python = pd.read_csv("output/scaled_py/b.csv")
iLat_b_python = b_python.iloc[:, 0]
iLon_b_python = b_python.iloc[:, 1]
value_b_python= b_python.iloc[:, 2]

PythonPsi_python = pd.read_csv("output/scaled_py/PythonPsi.csv")
iLat_PythonPsi_python = PythonPsi_python.iloc[:, 0]
iLon_PythonPsi_python = PythonPsi_python.iloc[:, 1]
value_PythonPsi_python= PythonPsi_python.iloc[:, 2]

#!scaled Fortran Ocean 
AC_ocean_fort = pd.read_csv("output/scaled_ocean_fort/AC.csv")
iLat_AC_ocean_fort = AC_ocean_fort.iloc[:, 0]
iLon_AC_ocean_fort = AC_ocean_fort.iloc[:, 1]
value_AC_ocean_fort= AC_ocean_fort.iloc[:, 2]

b_ocean_fort = pd.read_csv("output/scaled_ocean_fort/b.csv")
iLat_b_ocean_fort = b_ocean_fort.iloc[:, 0]
iLon_b_ocean_fort = b_ocean_fort.iloc[:, 1]
value_b_ocean_fort= b_ocean_fort.iloc[:, 2]

psi_ocean_fort = pd.read_csv("output/scaled_ocean_fort/PsiOcean.csv")
iLat_psi_ocean_fort = psi_ocean_fort.iloc[:, 0]
iLon_psi_ocean_fort = psi_ocean_fort.iloc[:, 1]
value_psi_ocean_fort= psi_ocean_fort.iloc[:, 2]

#!scaled Fortran Spherical 
AC_spherical_fort = pd.read_csv("output/scaled_spherical_fort/AC.csv")
iLat_AC_spherical_fort = AC_spherical_fort.iloc[:, 0]
iLon_AC_spherical_fort = AC_spherical_fort.iloc[:, 1]
value_AC_spherical_fort= AC_spherical_fort.iloc[:, 2]

b_spherical_fort = pd.read_csv("output/scaled_spherical_fort/b.csv")
iLat_b_spherical_fort = b_spherical_fort.iloc[:, 0]
iLon_b_spherical_fort = b_spherical_fort.iloc[:, 1]
value_b_spherical_fort= b_spherical_fort.iloc[:, 2]

psi_spherical_fort = pd.read_csv("output/scaled_spherical_fort/PsiSpherical.csv")
iLat_psi_spherical_fort = psi_spherical_fort.iloc[:, 0]
iLon_psi_spherical_fort = psi_spherical_fort.iloc[:, 1]
value_psi_spherical_fort= psi_spherical_fort.iloc[:, 2]

AC_FLUXFORM_NEUMANN_SCALED = pd.read_csv("output/FluxForm/scaled/AC.csv")
iLat_AC_FLUXFORM_NEUMANN_SCALED = AC_FLUXFORM_NEUMANN_SCALED.iloc[:, 0]
iLon_AC_FLUXFORM_NEUMANN_SCALED = AC_FLUXFORM_NEUMANN_SCALED.iloc[:, 1]
value_AC_FLUXFORM_NEUMANN_SCALED= AC_FLUXFORM_NEUMANN_SCALED.iloc[:, 2]

#! MODIFIED POLES:
#%%
AC_Python = pd.read_csv("output/Fortran/LambdaPython/unscreened/AC.csv")
iLat_AC_Python = AC_Python.iloc[:, 0]
iLon_AC_Python = AC_Python.iloc[:, 1]
value_AC_Python= AC_Python.iloc[:, 2]

b_Python = pd.read_csv("output/Fortran/LambdaPython/unscreened/b.csv")
iLat_b_Python = b_Python.iloc[:, 0]
iLon_b_Python = b_Python.iloc[:, 1]
value_b_Python= b_Python.iloc[:, 2]

# %% !COMMENT: PLOTS USED
c_py   = 'royalblue'
c_oc   = '#FF5F1F'
c_sph  = '#00F5FF'

# ----------------------- AC -----------------------
plt.figure()
plt.plot(iLat_AC_python,     value_AC_python,     color=c_py,  linewidth=2, label='Python')

plt.plot(iLat_AC_Python, value_AC_Python, color=c_oc, linestyle='-', linewidth=1, label='Modified No Neumann Scaled interior')

# plt.plot(iLat_AC_spherical_fort, value_AC_spherical_fort, color=c_sph,linestyle='--', linewidth=1,label='Fortran Spherical')
plt.xlabel("iLat"); plt.ylabel("AC value")
plt.title("AC comparison")
plt.legend(); plt.grid(True)
plt.show()

# ----------------------- b ------------------------
plt.figure()
plt.plot(iLat_b_python,     value_b_python,     color=c_py, linewidth=2,  label='Python')

plt.plot(iLat_b_Python, value_b_Python, color=c_oc, linestyle='-', linewidth=1, label='Modified No Neumann Scaled interior')

# plt.plot(iLat_b_spherical_fort, value_b_spherical_fort, color=c_sph,linestyle='--', linewidth=1, label='Fortran Spherical')
plt.xlabel("iLat"); plt.ylabel("b value")
plt.title("b comparison")
plt.legend(); plt.grid(True)
plt.show()

# ----------------------- psi -----------------------
# plt.figure()
# plt.plot(iLat_PythonPsi_python, value_PythonPsi_python, color=c_py,  linewidth=2, label='Python')

# plt.plot(iLat_psi_ocean_fort,   value_psi_ocean_fort,   color=c_oc, linestyle='-', linewidth=1, label='Fortran Ocean')

# plt.plot(iLat_psi_spherical_fort,   value_psi_spherical_fort,   color=c_sph,  linestyle='--', linewidth=1,label='Fortran Spherical')

# plt.xlabel("iLat"); plt.ylabel("Psi")
# plt.title("ψ comparison")
# plt.legend(); plt.grid(True)
# plt.show()

# %%
AC_python        = AC_python.values
AC_ocean_fort    = AC_ocean_fort.values
AC_spherical_fort = AC_spherical_fort.values
AC_Python =AC_Python.values

b_python        = b_python.values
b_ocean_fort    = b_ocean_fort.values
b_spherical_fort = b_spherical_fort.values
#%%
import numpy as np
import matplotlib.pyplot as plt

def row_sum(lat, val):
    u = np.unique(lat)
    idx = np.searchsorted(u, lat)
    return u, np.bincount(idx, weights=val)

def plot(lat, py, oc, sph, name):
    plt.figure()
    # plt.plot(lat, py-oc,  color=c_oc,  linestyle='--', linewidth=1,
            #  label='Fortran Ocean - Python')
    plt.plot(lat, py-sph, color=c_sph, linestyle='--', linewidth=1,
             label='Fortran Spherical - Python')
    plt.title(f"{name} diff vs lat")
    plt.grid()
    plt.legend()   # <-- THIS WAS MISSING
    plt.show()

# AC
lat_AC, AC_py = AC_python[:,0], AC_python[:,2]

l_AC, AC_py   = row_sum(lat_AC, AC_py)

_,   AC_modified   = row_sum(lat_AC, AC_Python[:,2])

# _,   AC_sph   = row_sum(lat_AC, AC_spherical_fort[:,2])
# plot(l_AC, AC_py, AC_oc, AC_sph, "AC")

# b
lat_b, b_py = b_python[:,0], b_python[:,2]
l_b, b_py   = row_sum(lat_b, b_py)
_,   b_oc   = row_sum(lat_b, b_ocean_fort[:,2])
_,   b_sph  = row_sum(lat_b, b_spherical_fort[:,2])
plot(l_b, b_py, b_oc, b_sph, "b")

print(AC_python.shape, AC_ocean_fort.shape, AC_spherical_fort.shape)
print(b_python.shape, b_ocean_fort.shape, b_spherical_fort.shape)

# %%

import numpy as np
from scipy.sparse import csr_matrix

def per_lat_diff(A_old, A_new, nLat_old, nPhi_old, nLat_new, nPhi_new):
    """
    Compute per-latitude row-block difference norms between two sparse matrices.
    Returns (lats_old, diff_old) and (lats_new, diff_new) where:
      - if shapes align, diffs are for the same lat indices (single array printed).
      - if new has +2 lat (poles) and same nPhi, compares new[1:-1] vs old and also shows pole rows.
      - otherwise returns separate per-lat norms for each matrix (no alignment).
    """
    A_old = csr_matrix(A_old)
    A_new = csr_matrix(A_new)

    # convenience
    rows_old = A_old.shape[0]
    rows_new = A_new.shape[0]

    # case 1: identical shape -> direct difference per lat
    if A_old.shape == A_new.shape and nLat_old * nPhi_old == rows_old:
        D = (A_new - A_old).tocsr()
        nLat = nLat_old
        nPhi = nPhi_old
        diffs = np.zeros(nLat)
        for i in range(nLat):
            r0 = i * nPhi
            r1 = (i + 1) * nPhi
            block = D[r0:r1, :].toarray()
            diffs[i] = np.linalg.norm(block, ord='fro')
        return {'mode':'aligned','nLat':nLat,'nPhi':nPhi,'diffs':diffs}

    # case 2: same nPhi, new has +2 lat (poles added)
    if nPhi_old == nPhi_new and nLat_new == nLat_old + 2:
        nPhi = nPhi_old
        # compute diffs for interior comparison (map old i -> new i+1)
        diffs_interior = np.zeros(nLat_old)
        for i in range(nLat_old):
            r_old0 = i * nPhi
            r_old1 = (i + 1) * nPhi
            r_new0 = (i + 1) * nPhi
            r_new1 = (i + 2) * nPhi
            block_diff = (A_new[r_new0:r_new1, :] - A_old[r_old0:r_old1, :]).toarray()
            diffs_interior[i] = np.linalg.norm(block_diff, ord='fro')
        # compute pole norms (how big are the pole rows)
        south_block = A_new[0:nPhi, :].toarray()
        north_block = A_new[(nLat_new-1)*nPhi: nLat_new*nPhi, :].toarray()
        pole_norms = {'south': np.linalg.norm(south_block, ord='fro'),
                      'north': np.linalg.norm(north_block, ord='fro')}
        return {'mode':'new_has_poles','nPhi':nPhi,
                'diffs_interior':diffs_interior,'pole_norms':pole_norms}

    # fallback: return per-lat norms separately (shapes incompatible)
    def per_lat_norms(A, nLat, nPhi):
        norms = np.zeros(nLat)
        for i in range(nLat):
            r0 = i * nPhi
            r1 = (i + 1) * nPhi
            block = A[r0:r1, :].toarray()
            norms[i] = np.linalg.norm(block, ord='fro')
        return norms

    res = {'mode':'separate'}
    if rows_old == nLat_old * nPhi_old:
        res['old_norms'] = per_lat_norms(A_old, nLat_old, nPhi_old)
    if rows_new == nLat_new * nPhi_new:
        res['new_norms'] = per_lat_norms(A_new, nLat_new, nPhi_new)
    return res

# ---------- Example usage ----------
# A_old, A_new are your sparse matrices (scipy.sparse)
# supply nLat_old, nPhi_old, nLat_new, nPhi_new accordingly
out = per_lat_diff(AC_python[:,2], AC_Python[:,2], 720, 1440, 722, 1444)

# %%

