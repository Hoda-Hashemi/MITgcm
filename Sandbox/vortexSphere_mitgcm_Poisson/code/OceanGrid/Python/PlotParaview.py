
#%% 
import os
from src.plotting import plot_paraview , plot_paraview_Ocean
from src.utils import ParaviewSavedData 

#%%
# DATA_DIR = "../results/Ld50000.csv"
# DATA_DIR = "../results/Ld1000000.csv"
# DATA_DIR = "../results/depth1000.csv"
# DATA_DIR = "../results/theta0.csv"

# DATA_DIR = "../results/Old_unscaled.csv"  #!COMMENT: old code unscaled 
# DATA_DIR = "../results/Old_unscaledLdTheta.csv" #!COMMENT: old code unscaled Ld(theta)

DATA_DIR = "../results/Old_unscaledCorrected.csv" #!COMMENT: old code unscaled - corrected run

print("Reading", DATA_DIR)

plot_paraview(DATA_DIR)

# %%
# DATA_DIR = "../results/OceanCorrected.csv"  #!COMMENT: Ocean code - corrected scaled
# DATA_DIR = "../results/OceanCorrected_unscaled.csv" 
# DATA_DIR = "../results/RHS0.csv" 
DATA_DIR = "../results/RHS0unscaled.csv" 
print("Reading", DATA_DIR)
df_para = ParaviewSavedData(os.path.join(DATA_DIR))

plot_paraview_Ocean(df_para,DATA_DIR)

# %%

