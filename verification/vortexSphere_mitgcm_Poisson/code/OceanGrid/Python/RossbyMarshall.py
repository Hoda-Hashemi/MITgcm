#%%
#!BUGTRACKING: 
import numpy as np
import matplotlib.pyplot as plt

#PARAMETERS: 
L = 8* 10**6      
k = 4 * np.pi / L  
l = k            
a = 5.1* 10**5        

#beta-plane
R = 6371* 10**3
OMEGA = (2 * np.pi)/(24 * 3600)
theta0 = np.deg2rad(0)
beta = (2*OMEGA*np.cos(theta0))/R


#Dispersion relation: 
omega = -beta * k / (k**2 + l**2)  


#Grid parameters: resolution of the mesh: 
nx, ny = 400, 400  
x = np.linspace(0, L, nx)  
y = np.linspace(0, L, ny) 
X, Y = np.meshgrid(x, y)   # 2D meshgrid 

def PSI(X, Y, t):
    return a*np.sin(k*X - omega*t) * np.sin(l*Y) 

def vorticity(X, Y, t):
    return -(k**2+l**2)* PSI(X, Y, t)

#!FIG2 LEFT
psi_t0 = PSI(X, Y, t=0)

plt.figure(figsize=(6,5))
im = plt.imshow(psi_t0/1e5,
                # extent=[0,8,0,8],
                origin='lower',
                cmap='seismic',
                aspect='equal')
                # vmin=-1.8, vmax=1.8)

cbar = plt.colorbar(im)
cbar.set_label(r'$\psi\;(\times10^5\;\mathrm{m^2/s})$')
plt.xlabel('x ($10^6$ m)')
plt.ylabel('y ($10^6$ m)')
plt.title('Stream function at t = 0')

plt.show()


#! Fig2 right parameters:
y0 = 5e6  

nt = 400

t_max = 230*86400        #days - > [s]
times = np.linspace(0, t_max, nt) 

days = times / 86400.0   #[s] -> days  

# index of y closest to y0
j = np.argmin(np.abs(y - y0))

# allocate matrix: rows=time, columns=x
hov = np.zeros((nt, nx))

for it, t in enumerate(times):
    psi_t = PSI(X, Y, t)
    hov[it, :] = psi_t[j, :]   # slice at y = y0


#!FIG2 RIGHT
plt.figure(figsize=(9,5))

extent = [0, L/1e6, days.max(), days.min()]  # x=0..8, y=0..t_max in days
plt.imshow(hov/1e5,
           cmap='seismic',
        #    vmin=-1.8, vmax=1.8,
           extent=extent,
           aspect='auto')

plt.colorbar(label=r'$\psi\;(\times10^5\;m^2/s)$')
plt.xlabel('x ($10^6$ m)')
plt.ylabel('time (days)')
plt.title(r'Hovmöller diagram at $y=5\times10^6$ m')

# Dashed analytical line
c = omega/k  
c_plot = c * 86400 / 1e6   # 10^6 m/day
n = 1
x0 = (np.pi/2 + 2*np.pi*n) / k / 1e6
x_line = x0 + c_plot * days
# plt.plot(x_line, days, 'k--', lw=2)

# L_millions = L / 1e6
# x_line_wrapped = np.mod(x_line, L_millions)
#!mask it for the sake of plotting
mask = (x_line >= 0) & (x_line <= L/1e6)
plt.plot(x_line[mask], days[mask], 'k--',lw=3)
# plt.plot(x_line, days, 'k--', lw=2)

plt.gca().invert_yaxis()
plt.show()



#%%
#!FIGURES 3:

#PARAMETERS: 
L = 8* 10**6      
k = 4 * np.pi / L  
l = k            
a = 5.1* 10**5        

#beta-plane
R = 6371* 10**3
OMEGA = (2 * np.pi)/(24 * 3600)
theta0 = np.deg2rad(0)
beta = (2*OMEGA*np.cos(theta0))/R


#Dispersion relation: 
omega = -beta * k / (k**2 + l**2)  


#Grid parameters: resolution of the mesh: 
nx, ny = 400,400
x = np.linspace(0, L, nx)  
y = np.linspace(0, L, ny) 
X, Y = np.meshgrid(x, y)   # 2D meshgrid 

T_period = (2*np.pi)/np.abs(omega)


psi_period = PSI(x,y0,t=T_period)

zeta_period = vorticity(x,y0,T_period)
# --- Streamfunction ---
plt.figure()
plt.plot(x/1e6,psi_period/1e5,color='blue'
, label='Analytical')
plt.xlabel('x ($10^6$ m)')
plt.ylabel('Stream Function ($10^5$ m$^2$/s)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()


# --- Vorticity ---
plt.figure()
plt.plot(x/1e6, zeta_period,
         color='darkorange',
         lw=2,
         label='Analytical')
plt.xlabel('x ($10^6$ m)')
plt.ylabel('Vorticity  s$^{-1}$)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()


# %%
#!Atmmosphere:

#!FIGURE4: a
#parameters: 
Omega = (2*np.pi)/(24*3600)

def Coriolis(thetas):
    return 2*Omega*np.cos(thetas)

def Rossby_Hauwartiz():
    return 
#!FIGURE4: b

#!------------------------------------------


#!OCEAN: 
#!FIGURE4: c
#!Using Linear scheme
Ld=100 #km
a = 1e6



#!FIGURE4: d

#%%
#!FIGURE 5: Comparison between the numerical and analytical solution
#! Ld=100 km , period t=2pi/omega , cubic interpolation (npt linear)

#parameters: 
A = 1e6

# Parameters
R = 6.371e6
Omega = 7.292e-5
m = 4
# A = 1.0
omega = -2 * Omega * m / ((m + 1) * (m + 2))

# Grid
nlat, nlon = 181, 360
TH = np.deg2rad(30) #np.linspace(-np.deg2rad(80),np.deg2rad(80), nlat)
LA = np.linspace(0, 2*np.pi, nlon)
# TH, LA = np.meshgrid(theta, lam, indexing="ij")

# Exact solution
t = 0.0
psi = -A * np.sin(TH) * np.cos(TH)**m * np.cos(m*LA - omega*t)

# zeta = (-m**2 / np.cos(theta)**2)*psi + ((-a*np.cos())/)

# Plot
# plt.figure()
# plt.contourf(LA, TH, psi, levels=30)
# plt.xlabel(" λ")
# plt.ylabel(" θ")
# plt.title("Rossby–Haurwitz Streamfunction ψ at t = 0")
# plt.colorbar(label="ψ")
# plt.show()

# lon_deg = np.degrees(lam)

# plt.figure()
# plt.plot(lon_deg, zeta/1e8, 'blue', lw=1, label='Exact')
# plt.xlabel('Longitude')
# plt.ylabel(r'$\zeta$')
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.legend()

plt.figure()
plt.plot(LA,psi/1e5, 'orange', lw=1, label='Exact')
plt.xlabel('Longitude')
plt.ylabel(r'$\psi$ *e5' )
plt.grid(True, linestyle='--', alpha=0.5)
# plt.legend()

plt.show()
# %%




# %%
