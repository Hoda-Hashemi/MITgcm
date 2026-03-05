# %%! COMMENT: Solving the ODE on python 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
R = 6371000
Omega = 0.0000727

# Analytic solution (regular at poles)
def psi_analytic(theta, C2=0.0):
    return R**2 * Omega * np.sin(theta) + C2

# ODE system: convert 2nd-order ODE to 1st-order system
# Original: d/dθ (cosθ * dψ/dθ) = -2 R² Ω sinθ cosθ
# Let y[0] = ψ, y[1] = dψ/dθ
def ode_system(theta, y):
    # Avoid singularities at poles by masking (not needed for interior)
    cos_t = np.cos(theta)
    if np.abs(cos_t) < 1e-12:
        # Use limit: RHS → 0 as θ→±π/2 (since sinθ cosθ → 0)
        d2psi = 0.0
    else:
        d2psi = (-2 * R**2 * Omega * np.sin(theta) * np.cos(theta) 
                 - (-np.sin(theta)) * y[1]) / cos_t
    return [y[1], d2psi]

# But better: use the first integral directly for numerical stability:
def ode_first_order(theta, psi_prime):
    # From: cosθ * ψ' = (R²Ω/2) cos(2θ) + C1
    # For regular solution, C1 = R²Ω/2 → RHS = (R²Ω/2)(cos2θ + 1) = R²Ω cos²θ
    # So ψ' = R²Ω cosθ
    return R**2 * Omega * np.cos(theta)

# Numerical integration of ψ' = R²Ω cosθ
theta_vals = np.linspace(-np.pi/2 + 1e-6, np.pi/2 - 1e-6, 500)
psi_num = R**2 * Omega * np.sin(theta_vals)  # Exact integral

# Plot
plt.figure(figsize=(8,5))
plt.plot(theta_vals, psi_analytic(theta_vals), 'r--', label='Analytic: $R^2\\Omega\\sin\\theta$')
plt.plot(theta_vals, psi_num, 'b:', label='Numerical (integrated)')
plt.xlabel(r'$\theta$ (radians)')
plt.ylabel(r'$\psi(\theta)$')
plt.title('Validation of Analytic Solution')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Print max difference (should be ~0)
print("Max difference:", np.max(np.abs(psi_analytic(theta_vals) - psi_num)))

# %% !COMMENT: only Direct Integral
import numpy as np
import matplotlib.pyplot as plt

Omega = 7.292e-5        
R = 6371000        

# === Create latitude grid ===
nlat = 720
dtheta = np.pi / nlat
theta = np.linspace(-np.pi/2 + dtheta/2, np.pi/2 - dtheta/2, nlat)  # your interior grid

# === Analytical solution ===
psi_analytical = Omega * R**2 * np.sin(theta)

# Optional: normalize / scale for nicer plotting (comment out if you want physical units)
# psi_analytical /= (Omega * R**2)   # → max value = 1

# === Plotting ===
plt.figure(figsize=(10, 6))

# 1. Main plot: ψ vs latitude
plt.plot(np.degrees(theta), psi_analytical,
         color='royalblue', linewidth=2.5, label='Analytical: ψ = Ω R² sin(θ)')

# 2. Add zero line and equator reference
plt.axhline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
plt.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

# Formatting
plt.xlabel('Latitude [degrees]', fontsize=12)
plt.ylabel('Streamfunction ψ  [m²/s]', fontsize=12)
plt.title('Analytical Solution\n∇²ψ = -2Ω sin(θ)    (axisymmetric case)', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11)

# Optional: add some annotation
plt.text(-80, psi_analytical.min()*0.9,
         f'Ω = {Omega:.2e} rad/s\nR = {R/1000:.0f} km',
         fontsize=10, bbox=dict(facecolor='white', alpha=0.7))

plt.tight_layout()
plt.show()

# %%
#!BUGTRACKING: Spherical harmonics and direct solution:                     
import numpy as np
from scipy.special import lpmv
import matplotlib.pyplot as plt

# Parameters
R = 6371e3         
Omega = (2*np.pi)/(24*3600) #7.292e-5   
nlat = 720   

# theta_lat = np.linspace(-np.pi/2, np.pi/2, nlat)  # latitude in radians

dlat = np.pi/nlat
delta2 = dlat**2

theta_lat = np.linspace(
    -np.pi/2 + dlat/2,
    +np.pi/2 - dlat/2,
    nlat
)

#? RHS
f = -2 * Omega * np.sin(theta_lat)

#? Mode: l=1, m=0 check my notes for full derivation
l=1
m=0
P1 = lpmv(m, l, np.sin(theta_lat)) 

# Normalization of Y_1^0
norm = np.sqrt((2*l + 1)/(4*np.pi))          # ≈ 0.4886

#! f(theta) = -2 Omega sin(theta) only (so unscaled)
#! f(theta) = sum f_l Yl = f_l ( norm * sin(theta)) 

f_l1 =  (-8/3)*np.pi *Omega*np.sqrt(3/(4*np.pi))
# ( -2 * Omega *np.sin(theta_lat)) /(norm*P1)

# theta0=np.pi/6
# f_l1 = (-2*Omega*np.sin(theta0))/(norm*P1)

Ld = 100
# np.sqrt(9.8*1000)/(2*Omega*np.sin(np.pi/4))
alpha = 0#1/Ld**2

# Solution coefficient: Ld= R m
psi_l1 = (f_l1) / (-l*(l + 1)/R**2 - alpha)     

psi_harmonic = psi_l1 * norm * P1

#!Compare with known analytical solution
psi_analytical = Omega * R**2 * np.sin(theta_lat)
# ────────────────────────────────────────────────
# print(f"Max difference between methods: {np.max(np.abs(psi_harmonic - psi_analytical)):.2e}")

plt.figure(figsize=(10,6))
# plt.plot(np.degrees(theta_lat), psi_analytical,'royalblue', lw=2.5, label='Direct analytical: Ω R² sin(θ)')
plt.plot(np.degrees(theta_lat), psi_harmonic, 'orangered', lw=2.2, ls='--', label=f'From spherical harmonic projection: Ld = {Ld} m')
        #   Ld = {Ld} m'
plt.xlabel('Latitude [°]')
plt.ylabel('ψ  [m²/s]')
plt.title('Analytical solution via spherical harmonic (ℓ=1 only)')
plt.legend()
plt.grid(alpha=0.3)
plt.show()
#%%
#!BUGTRACKING: Screened Poisson with full spherical Laplacian (m != 0)
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv
from math import factorial

# Equation solved:
# (1/(R^2 cos(theta))) d/dtheta(cos(theta) dpsi/dtheta)
# + (1/(R^2 cos(theta)^2)) d2psi/dlambda2 - alpha*psi = f(theta,lambda)

# Parameters
R = 6371e3
Omega = 2 * np.pi / (24 * 3600)
nlat = 720
nlon = 2*nlat
Lmax = 12

dlat = np.pi / nlat
dlon = 2 * np.pi / nlon
theta_lat = np.linspace(-np.pi / 2 + dlat / 2, np.pi / 2 - dlat / 2, nlat)
lambda_lon = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
Theta, Lambda = np.meshgrid(theta_lat, lambda_lon, indexing="ij")

# Choose a non-zonal wavenumber to include longitudinal structure.
# m_forcing = number of waves around a full 360-degree longitude circle.
m_forcing = 1
eps_lon = 0.25

# If eps_lon = 0, forcing is purely axisymmetric and only m=0 exists.
f = -2 * Omega * np.sin(Theta) * (1.0 + eps_lon * np.cos(m_forcing * Lambda))

Ld = np.inf
alpha = 0.0 if np.isinf(Ld) else 1.0 / Ld**2

def Ylm(l, m, theta, lon):
    mm = abs(m)
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - mm) / factorial(l + mm))
    P_lm = lpmv(mm, l, np.sin(theta))
    if m >= 0:
        return norm * P_lm * np.exp(1j * m * lon)
    return ((-1) ** mm) * norm * P_lm * np.exp(-1j * mm * lon)

# Project RHS onto spherical harmonics
cell_area_weight = np.cos(Theta) * dlat * dlon
f_lm = {}
for l in range(Lmax + 1):
    for m in range(-l, l + 1):
        Y = Ylm(l, m, Theta, Lambda)
        f_lm[(l, m)] = np.sum(f * np.conj(Y) * cell_area_weight)

# Solve mode-by-mode in spectral space
psi_lm = {}
for (l, m), coeff in f_lm.items():
    denom = -l * (l + 1) / R**2 - alpha
    if l == 0 and np.isclose(denom, 0.0):
        psi_lm[(l, m)] = 0.0
    else:
        psi_lm[(l, m)] = coeff / denom

# Reconstruct psi(theta, lambda)
psi = np.zeros((nlat, nlon), dtype=complex)
for (l, m), coeff in psi_lm.items():
    psi += coeff * Ylm(l, m, Theta, Lambda)
psi = np.real(psi)

# Plot latitude slice and 2D map
plt.figure(figsize=(10, 5))
plt.plot(np.degrees(theta_lat), psi[:, 0], lw=2.0, color="orangered")
plt.xlabel("Latitude [deg]")
plt.ylabel("psi [m^2/s]")
plt.title(f"Screened Poisson solution at lambda=0 (m_forcing={m_forcing})")
plt.grid(alpha=0.3)
plt.show()

plt.figure(figsize=(11, 4))
plt.pcolormesh(np.degrees(lambda_lon), np.degrees(theta_lat), psi, shading="auto", cmap="coolwarm")
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.title("psi(theta, lambda) with longitudinal structure")
plt.colorbar(label="psi [m^2/s]")
plt.tight_layout()
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv, factorial

# -------------------------------
# PARAMETERS
# -------------------------------
R      = 6371e3
Omega  = 2*np.pi/(24*3600)
nlat   = 720
nlon   = 2*nlat

theta_lat = np.linspace(-np.pi/2, np.pi/2, nlat)       # physical latitude
lambda_lon = np.linspace(0,2*np.pi,nlon,endpoint=False) # physical longitude

# full 2D forcing on physical grid
Theta, Lambda = np.meshgrid(theta_lat, lambda_lon, indexing='ij')
f_vals = -2*Omega*np.sin(Theta)

# -------------------------------
# SPHERICAL HARMONICS SETUP
# -------------------------------
Lmax = 10  # max spherical degree

# --- Gauss–Legendre quadrature on mu = sin(theta)
nquad = 2*Lmax + 10
mu, wmu = np.polynomial.legendre.leggauss(nquad)
theta_q = np.arcsin(mu)                      # quadrature latitudes
lambda_q = lambda_lon                        # keep same lon grid
dphi = 2*np.pi/len(lambda_q)

# forcing on quad grid
Theta_q, Lambda_q = np.meshgrid(theta_q, lambda_q, indexing='ij')
f_q = -2*Omega*np.sin(Theta_q)

# ----------------------------------------------------------------
# Project forcing onto spherical harmonics f_lm
# ----------------------------------------------------------------
f_lm = {}
for l in range(Lmax+1):
    for m in range(-l, l+1):
        # associated Legendre on mu
        Plm = lpmv(abs(m), l, mu)
        # normalization + phase
        norm = np.sqrt((2*l+1)/(4*np.pi) * factorial(l-abs(m))/factorial(l+abs(m)))
        phase = (-1)**m
        # construct Y_l^m on quad grid
        Y_lm_q = norm * phase * np.outer(Plm, np.exp(1j*m*lambda_q))
        # inner product ∫ Y* f dΩ
        integrand = np.conj(Y_lm_q) * f_q
        # integrate over mu (wmu), then over lambda (uniform)
        f_lm[(l,m)] = np.sum(integrand * wmu[:,None]) * dphi

# ----------------------------------------------------------------
# Solve spectral Poisson: psi_lm = f_lm / (-l(l+1)/R^2)
# ----------------------------------------------------------------
psi_lm = {}
for (l,m), val in f_lm.items():
    if l == 0:
        psi_lm[(l,m)] = 0.0  # no mean mode from forcing
    else:
        psi_lm[(l,m)] = val / (-l*(l+1)/R**2)

# ----------------------------------------------------------------
# Reconstruct ψ(θ,λ) on physical grid
# ----------------------------------------------------------------
psi_vals = np.zeros((nlat,nlon),dtype=complex)

for (l,m), coeff in psi_lm.items():
    Plm_phys = lpmv(abs(m), l, np.sin(theta_lat))      # physical lat
    norm = np.sqrt((2*l+1)/(4*np.pi)*factorial(l-abs(m))/factorial(l+abs(m)))
    phase = (-1)**m
    Y_lm_phys = np.outer(Plm_phys*norm*phase, np.exp(1j*m*lambda_lon))
    psi_vals += coeff * Y_lm_phys

psi_vals = np.real(psi_vals)

# ----------------------------------------------------------------
# PLOT (zonal slice at λ=0)
# ----------------------------------------------------------------
plt.figure(figsize=(8,5))
plt.plot(np.degrees(theta_lat), psi_vals[:,0], 'r', lw=2)
plt.xlabel('Latitude [°]')
plt.ylabel('ψ [m²/s]')
plt.title('Unscreened Poisson: ψ vs latitude (λ=0)')
plt.grid(alpha=0.3)
plt.show()

#%%
# #!Velocity:
# eq_text = r"$\frac{\partial \psi}{\partial \theta} \approx \frac{\psi_{i+1} - \psi_{i-1}}{2\,\Delta \theta}$"

# # ─── Finite-difference derivative ────────────────
# dpsi_dtheta = np.zeros_like(psi_harmonic)

# dpsi_dtheta[1:-1] = (psi_harmonic[2:] - psi_harmonic[:-2]) / (2*dlat)
# dpsi_dtheta[0]    = (psi_harmonic[1] - psi_harmonic[0]) / dlat
# dpsi_dtheta[-1]   = (psi_harmonic[-1] - psi_harmonic[-2]) / dlat

# u_lambda = -(1/R) * dpsi_dtheta

# plt.figure(figsize=(10,6))
# plt.plot(np.degrees(theta_lat), u_lambda, lw=2)
# plt.xlabel('Latitude [°]')
# plt.ylabel(r'$u^\lambda$  [m/s]')
# plt.title(f"Zonal velocity $u^\\lambda$ using central finite difference\n{eq_text}")
# plt.grid(alpha=0.3)
# plt.show()

# %% #! varying LD:

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval

nLat   = 720
R      = 6371e3
Omega  = 2*np.pi/(24*3600)
g      = 9.8
d0     = 1000.0
Lmax   = nLat // 2

#! Gauss–Legendre quadrature
# mu = sin(theta) in [-1,1]
nquad = 2*Lmax + 10
mu, wmu = np.polynomial.legendre.leggauss(nquad)
theta_lat = np.arcsin(mu)     # latitude

#! Zonal basis: Y_{l0}(theta)
#! Y_{l0} = sqrt((2l+1)/(4pi)) P_l(mu)
Y = {}
for ell in range(Lmax+1):
    P = legval(mu, [0]*ell + [1])   # P_ell(mu)
    Y[(ell,0)] = np.sqrt((2*ell+1)/(4*np.pi)) * P

#!Forcing: f(theta) = 2 Omega sin(theta)
#! RHS projection
f_vals = -2 * Omega * mu
f_vec  = np.zeros(Lmax+1)

#!only l=1 survives analytically, but keep general
for ell in range(Lmax+1):
    f_vec[ell] = 2*np.pi * np.sum(f_vals * Y[(ell,0)] * wmu)

#! Variable screening: Ld^{-2}(theta)
Ld2_vals = (2*Omega*mu)**2 / (g*d0)

a_k0 = np.zeros(Lmax+1)
for k in range(Lmax+1):
    a_k0[k] = 2*np.pi * np.sum(Ld2_vals * Y[(k,0)] * wmu)

#!Triple-product coupling C(L,k,l) = ∫ Y_L Y_k Y_l dΩ
def C(L, k, l):
    return 2*np.pi * np.sum(
        Y[(L,0)] * Y[(k,0)] * Y[(l,0)] * wmu
    )
# Build linear system
A = np.zeros((Lmax+1, Lmax+1))

for L in range(Lmax+1):
    for l in range(Lmax+1):
        #! Laplacian (diagonal)
        if L == l:
            A[L,l] += -L*(L+1)/R**2
        # Variable screening (mode coupling)
        for k in range(Lmax+1):
            if abs(a_k0[k]) > 1e-14:
                A[L,l] -= a_k0[k] * C(L, k, l)

# #!Stabilization - nothing changed
# s = np.max(np.abs(A.diagonal()))
# A     /= s
# f_vec /= s

#! Solve spectral system
psi_coeffs = np.linalg.solve(A, f_vec)

#! Reconstruct psi(theta)
psi_vals = np.zeros_like(mu)
for ell in range(Lmax+1):
    psi_vals += psi_coeffs[ell] * Y[(ell,0)]

# Plots
plt.figure(figsize=(8,5))
plt.plot(np.degrees(theta_lat), psi_vals, 'b', lw=2)
plt.axhline(0, color='gray', lw=0.8)
plt.xlabel('Latitude (degrees)')
plt.ylabel(r'$\psi(\theta)$')
plt.title(r'Solution $\psi(\theta)$ with variable screening')
plt.grid(alpha=0.3)
plt.show()

# plt.figure(figsize=(8,5))
# for ell in range(0, min(6, Lmax+1)):
#     plt.plot(
#         np.degrees(theta_lat),
#         Y[(ell,0)],
#         label=fr'$\ell={ell}$'
#     )
# plt.xlabel('Latitude (degrees)')
# plt.ylabel(r'$Y_{\ell 0}(\theta)$')
# plt.title(r'Zonal Legendre basis $Y_{\ell 0}$')
# plt.legend()
# plt.grid(alpha=0.3)
# plt.show()

# %%
#!BUGTRACKING MARSHALL'S Fig5:
import numpy as np
import matplotlib.pyplot as plt

a = 1e6
m = 4
OMEGA = 7.292e-5

lat_deg = np.linspace(-90, 90, 720)
lat = np.deg2rad(lat_deg)  # radians

phi = 0.0  # longitude
omega = -2 * OMEGA * m / ((m+1)*(m+2))
t = 2*np.pi / abs(omega)

psi_analytical = -a * np.sin(lat) * (np.cos(lat)**m) * np.cos(m*phi - omega*t)

plt.figure(figsize=(8,4))
plt.plot(lat_deg, psi_analytical/1e5, 'b-', lw=2, label='Analytical')
plt.xlabel('Latitude (deg)')
plt.ylabel(r'Streamfunction $\psi$ ($10^5$ m$^2$/s)')
plt.title('Fig.5 Right: Streamfunction')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.show()

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv  
# -------------------------------
# PARAMETERS
# -------------------------------
R = 6371e3            # Earth radius [m]
Omega = 7.2921e-5     # Rotation rate [rad/s]
L_d = R
ell = 2       # Total spherical harmonic degree
m = 0                # Zonal wavenumber
A = 1.0               # Amplitude

# latitude array (theta = latitude in radians)
nLat = 180
theta = np.linspace(-np.pi/2, np.pi/2, nLat)

# -------------------------------
# ASSOCIATED LEGENDRE POLYNOMIAL
# -------------------------------
# P_l^m(sin(theta))
P_lm = lpmv(m, ell, np.sin(theta))

# -------------------------------
# PHASE SPEED (with constant L_d)
# -------------------------------
lambda_factor = R**2 / L_d**2
omega = -2 * Omega * m / (ell*(ell+1) + lambda_factor)

# -------------------------------
# STREAMFUNCTION psi_RH
# -------------------------------
t=0
# choose lambda=0 (axisymmetric slice)
psi_RH = A * P_lm * np.cos(-omega*0)  # t=0

# -------------------------------
# OPTIONAL: COMPARISON WITH NUMERICAL psi
# Here we use a dummy numerical psi for example
# Replace psi_num with your actual numerical solution
# -------------------------------
psi_num = psi_RH + 0.1*np.random.randn(nLat)  # example perturbation

# -------------------------------
# PLOT AND ERROR
# -------------------------------
plt.figure(figsize=(8,5))
plt.plot(theta*180/np.pi, psi_RH, label='Analytical RH')
# plt.plot(theta*180/np.pi, psi_num, '--', label='Numerical psi')
plt.xlabel('Latitude [deg]')
plt.ylabel('Streamfunction ψ')
plt.title(f'Rossby-Haurwitz Wave l={ell}, m={m}, Ld={L_d/1e3} km')
plt.legend()
plt.grid(True)
plt.show()

# # ERROR NORMS
# L2_error = np.sqrt(np.sum((psi_num - psi_RH)**2)/nLat)
# max_error = np.max(np.abs(psi_num - psi_RH))
# print(f"L2 error = {L2_error:.4e}, Max error = {max_error:.4e}")

# %%

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv

# -------------------------------
# PARAMETERS
# -------------------------------
R = 6371e3        # Earth radius [m]
Omega = 2*np.pi/(24*3600)  # rotation rate

Ld = R         # deformation radius (for RH)
nlat = 720
nlon = 2*nlat

t = 0              # initial time
# t = 3600*6    # time in seconds (6 hours)

# Latitude and longitude arrays
theta = np.linspace(-np.pi/2, np.pi/2, nlat)
lambda_lon = np.linspace(0, 2*np.pi, nlon)

# Meshgrid for 2D
Theta, Lambda = np.meshgrid(theta, lambda_lon, indexing='ij')

# -------------------------------
# 1) ROSSBY WAVE ANALYTICAL
# -------------------------------
a = Omega* R**2
m_r = 4 # zonal wavenumber for Rossby
l=4

omega_r = -2 * Omega * m_r / (l*(l+1)+ (R**2)/(Ld**2))
psi_rossby = a * np.sin(Theta) * np.cos(m_r * Lambda - omega_r * t)

# -------------------------------
# 2) ROSSBY-HAURWITZ WAVE ANALYTICAL
# -------------------------------

m = 0
mu = np.sin(theta)
dmu = mu[1] - mu[0]

# Forcing
f_theta = -2 * Omega * np.sin(theta)

# Choose degree l
l = 1

# Axisymmetric spherical harmonic Y_l0
Y_l0 = np.sqrt((2*l+1)/(4*np.pi)) * lpmv(0, l, mu)

# Project forcing onto Y_l0
f_l = 2 * np.pi * np.sum(f_theta * Y_l0 * dmu)  # integral over sphere

# Compute amplitude A
A = f_l / (-l*(l+1)/R**2 - 1/Ld**2)

print("Amplitude A =", A)

P_lm = lpmv(m, l, np.sin(Theta)) 

omega_rh = -2 * Omega * m / (l*(l+1) + (R**2)/(Ld**2))

psi_rh =  A * P_lm * np.cos(m*Lambda - omega_rh * t)

# -------------------------------
# 3) PLOT LATITUDINAL SLICE (axisymmetric)
# -------------------------------
plt.figure(figsize=(8,5))
plt.plot(np.degrees(theta), psi_rossby[:,0], label='Rossby wave (axisymmetric)')
# plt.plot(np.degrees(theta), psi_rh[:,0], '--', label='Rossby-Haurwitz (axisymmetric)')
plt.xlabel('Latitude [deg]')
plt.ylabel('Streamfunction ψ')
plt.title('Analytical waves (axisymmetric slice at λ=0)')
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# 4) OPTIONAL: Compare with your numerical psi
# Replace psi_num with your numerical solution
# -------------------------------
# psi_num = your_numerical_solution_array
# L2_error = np.sqrt(np.mean((psi_num - psi_rh[:,0])**2))
# print("L2 error:", L2_error)

# %%

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import lpmv

# -------------------------------
# PARAMETERS
# -------------------------------
R = 6371e3
Omega = 2*np.pi/(24*3600)
Ld = 1e6  # constant deformation radius for screened
nlat = 180
nlon = 360
Lmax = 10  # max spherical harmonic degree

# latitude and longitude grids
theta = np.linspace(-np.pi/2, np.pi/2, nlat)        # latitude
lambda_lon = np.linspace(0, 2*np.pi, nlon)          # longitude
Theta, Lambda = np.meshgrid(theta, lambda_lon, indexing='ij')

# forcing f(θ,λ) on full sphere
f2d = -2 * Omega * np.sin(Theta)

# -------------------------------
# Spherical harmonic transform (real)
# -------------------------------

# storage for spectral coefficients
f_lm = {}
psi_lm_unscreened = {}
psi_lm_screened = {}

# normalization factor
def N_lm(ell, m):
    return np.sqrt((2*ell+1)/(4*np.pi) * math.factorial(ell-abs(m)) / math.factorial(ell+abs(m)))

# forward transform: project f onto Y_lm
for ell in range(Lmax+1):
    for m in range(-ell, ell+1):
        # real spherical harmonic basis
        if m > 0:
            Ylm = np.sqrt(2)*N_lm(ell,m)*lpmv(m,ell,np.sin(Theta))*np.cos(m*Lambda)
        elif m < 0:
            Ylm = np.sqrt(2)*N_lm(ell,abs(m))*lpmv(abs(m),ell,np.sin(Theta))*np.sin(abs(m)*Lambda)
        else:
            Ylm = N_lm(ell,0)*lpmv(0,ell,np.sin(Theta))

        # inner product ∫ f * Ylm dΩ = ∑ f*Ylm*sinθΔθΔλ
        integrand = f2d * Ylm * np.cos(Theta)
        f_lm[(ell,m)] = np.sum(integrand) * (np.pi/(nlat-1)) * (2*np.pi/(nlon-1))

        # unscreened solve
        if ell*(ell+1) != 0:
            psi_lm_unscreened[(ell,m)] = f_lm[(ell,m)] / (-ell*(ell+1)/R**2)
        else:
            psi_lm_unscreened[(ell,m)] = 0.0  # ℓ=0 mode trivial

        # screened solve
        denom = -ell*(ell+1)/R**2 - 1/(Ld**2)
        psi_lm_screened[(ell,m)] = f_lm[(ell,m)] / denom if abs(denom)>0 else 0.0

# -------------------------------
# Inverse transform: reconstruct ψ
# -------------------------------

psi_unscreened = np.zeros_like(Theta)
psi_screened   = np.zeros_like(Theta)

for ell in range(Lmax+1):
    for m in range(-ell, ell+1):
        if m > 0:
            Ylm = np.sqrt(2)*N_lm(ell,m)*lpmv(m,ell,np.sin(Theta))*np.cos(m*Lambda)
        elif m < 0:
            Ylm = np.sqrt(2)*N_lm(ell,abs(m))*lpmv(abs(m),ell,np.sin(Theta))*np.sin(abs(m)*Lambda)
        else:
            Ylm = N_lm(ell,0)*lpmv(0,ell,np.sin(Theta))

        psi_unscreened += psi_lm_unscreened[(ell,m)] * Ylm
        psi_screened   += psi_lm_screened[(ell,m)] * Ylm

# -------------------------------
# PLOTS
# -------------------------------

# Unscreeened
plt.figure(figsize=(8,5))
plt.contourf(lambda_lon*180/np.pi, theta*180/np.pi, psi_unscreened, 20, cmap='viridis')
plt.colorbar(label='ψ (unscreened)')
plt.xlabel('Longitude [deg]')
plt.ylabel('Latitude [deg]')
plt.title('Unscreened Poisson spectral solution')
plt.show()

# Screened
plt.figure(figsize=(8,5))
plt.contourf(lambda_lon*180/np.pi, theta*180/np.pi, psi_screened, 20, cmap='viridis')
plt.colorbar(label='ψ (screened)')
plt.xlabel('Longitude [deg]')
plt.ylabel('Latitude [deg]')
plt.title(f'Screened Poisson spectral solution (Ld={Ld/1000:.1f} km)')
plt.show()

# %%

