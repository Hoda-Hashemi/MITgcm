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

#%%
import matplotlib.pyplot as plt
import numpy as np

# Pick three iterations
idxs = [0,1,2]  # adjust based on your available iterations
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
