import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np

# Ensure output folder exists
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# #!ADDED: 2D Comparsion
def plot_comparison_2d_paraview(df, solver, psi_model,u_model, run_dir=None):
    #!COMMENT: FIND COLUMNS 
    lat_col = next((c for c in df.columns if 'lat' in c.lower()), None)
    psi_col   = next((c for c in df.columns if 'psi' in c.lower()), None)

    u0_col = next((c for c in df.columns if 'u0' in c.lower()), None)
    u1_col =  next((c for c in df.columns if 'u1' in c.lower()), None)
    u2_col =  next((c for c in df.columns if 'u2' in c.lower()), None)

    if not all([lat_col, psi_col, u0_col, u1_col, u2_col]):
        raise ValueError(
            f"Need 'lat', 'psi', 'u0', 'u1', 'u2'. Found: {list(df.columns)}"
    )

    df['lat_deg'] = np.degrees(df[lat_col])
    lat_deg = df['lat_deg']

    psi_csv  = df[psi_col]
    psi_py   = psi_model

    umag_csv = np.sqrt(df[u0_col]**2 +df[u1_col]**2 + df[u2_col]**2)
    umag_py  = u_model
    
    Lat_deg = np.degrees(solver.Lat)
    solver_type = "2D"
    run_dir = run_dir or "output"

    #!COMMENT: --- Side-by-side ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    lines1 = ax1.plot(lat_deg, psi_csv, color='#00FFFF', lw=2)
    lines2 = ax2.plot(Lat_deg, psi_py, color='#FFD300', lw=2)

    ax1.legend([lines1[0]], ['Paraview'], loc='best')
    ax2.legend([lines2[0]], ['Python Model'], loc='best')

    for ax in (ax1, ax2):
        ax.set_xlabel('Latitude [°]')
        ax.set_ylabel('ψ')
        ax.grid(alpha=0.3)

    ax1.set_title('Paraview Output')
    ax2.set_title(f'{solver_type} Python Model')

    fig.suptitle(f'{solver_type} Solver: Streamfunction ψ vs Latitude — Side-by-Side', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    path = os.path.join(run_dir, f"compare_side_{solver_type}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {path}")

      #!COMMENT: --- Overlay ---
    fig, ax = plt.subplots(figsize=(9, 5))
    line1 = ax.plot(lat_deg, psi_csv.iloc[:, 0] if psi_csv.ndim > 1 else psi_csv, color='#00FFFF', lw=3, label='Paraview')[0]

    line2 = ax.plot(Lat_deg, psi_py, color='#FFD300', lw=2, label='Python Model')[0]

    ax.set_xlabel('Latitude [°]'); ax.set_ylabel('ψ')
    ax.set_title(f'{solver_type} Solver: Streamfunction ψ — Overlay Comparison')
    ax.grid(alpha=0.3)
    ax.legend([line1, line2], ['Paraview', 'Python Model'], loc='best')

    path = os.path.join(run_dir, f"compare_overlay_{solver_type}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {path}")

    #!COMMENT: adding plots for velocity
    #!COMMENT: Side by Side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    lines1 = ax1.plot(lat_deg, umag_csv, color='#00FFFF', lw=2)
    lines2 = ax2.plot(Lat_deg, umag_py, color='#FFD300', lw=2)

    ax1.legend([lines1[0]], ['Paraview'], loc='best')
    ax2.legend([lines2[0]], ['Python Model'], loc='best')

    for ax in (ax1, ax2):
        ax.set_xlabel('Latitude [°]')
        ax.set_ylabel('uψ')
        ax.grid(alpha=0.3)

    ax1.set_title('Paraview Output')
    ax2.set_title(f'{solver_type} Python Model')

    fig.suptitle(f'{solver_type} Solver: uψ magnitude vs Latitude — Side-by-Side', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    path = os.path.join(run_dir, f"compare_side_{solver_type}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {path}")

      #!COMMENT: --- Overlay ---
    fig, ax = plt.subplots(figsize=(9, 5))
    line1 = ax.plot(lat_deg, umag_csv.iloc[:, 0] if umag_csv.ndim > 1 else umag_csv, color='#00FFFF', lw=3, label='Paraview')[0]

    line2 = ax.plot(Lat_deg, umag_py, color='#FFD300', lw=2, label='Python Model')[0]

    ax.set_xlabel('Latitude [°]'); ax.set_ylabel('uψ magnitude')
    ax.set_title(f'{solver_type} Solver: uψ magnitude — Overlay Comparison')
    ax.grid(alpha=0.3)
    ax.legend([line1, line2], ['Paraview', 'Python Model'], loc='best')

    path = os.path.join(run_dir, f"compare_overlay_{solver_type}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {path}")

#!ADDED: 1D
def plot_comparison_1d_paraview(df, solver, psi_model, run_dir=None):
    #!COMMENT: --- FIND COLUMNS ---
    lat_col = next((c for c in df.columns if 'lat' in c.lower()), None)
    psi_col   = next((c for c in df.columns if 'psi' in c.lower()), None)
    if not lat_col or not psi_col:
        raise ValueError(f"Need 'lat' and 'psi'. Found: {list(df.columns)}")

    df['lat_deg'] = np.degrees(df[lat_col])
    lat_deg, psi_csv = df['lat_deg'], df[psi_col]
    psi_py = psi_model if psi_model.ndim == 1 else psi_model.flatten()
    Lat_deg = np.degrees(solver.Lat)
    solver_type = "1D"
    run_dir = run_dir or "output"

    #!COMMENT: --- Side-by-side ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    ax1.plot(lat_deg, psi_csv, color='#0096FF', lw=2, label='Paraview')
    ax2.plot(Lat_deg, psi_py, color='#FF4500', lw=2, label=f'{solver_type} Python')

    for ax in (ax1, ax2):
        ax.set_xlabel('Latitude [°]')
        ax.set_ylabel('ψ')
        ax.grid(alpha=0.3)
        ax.legend(loc='best')

    ax1.set_title('Paraview Output')
    ax2.set_title(f'{solver_type} Python Model')

    fig.suptitle(f'{solver_type} Solver: Streamfunction ψ vs Latitude — Side-by-Side', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(run_dir, f"compare_side_{solver_type}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {path}")

    #!COMMENT: --- Overlay ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(lat_deg, psi_csv, color='#0096FF', lw=3, label='Paraview')
    ax.plot(Lat_deg, psi_py, color='#FF4500', lw=2, label=f'{solver_type} Python')

    ax.set_xlabel('Latitude [°]')
    ax.set_ylabel('ψ')
    ax.set_title(f'{solver_type} Solver: Streamfunction ψ — Overlay Comparison')
    ax.grid(alpha=0.3)
    ax.legend(loc='best')

    plt.tight_layout()
    path = os.path.join(run_dir, f"compare_overlay_{solver_type}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {path}")

#!ADDED Plotting (2,2) plots for both solvers 
def plot_solution(solver, psi, u=None, v=None, run_dir=None):
    print(f"Plotting {solver} results...")

    Lat_deg = np.degrees(solver.Lat)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    ax1, ax2, ax3, ax4 = axes.ravel()

    ax1.plot(Lat_deg, psi ,linewidth=2, color='royalblue')
    ax1.set_xlabel("Latitude [°]")
    ax1.set_ylabel("Streamfunction ψ₁")
    ax1.set_title("ψ₁ vs Latitude")
    ax1.grid(True, alpha=0.3)

    ax2.plot(Lat_deg, u ,linewidth=2, color='crimson')
    ax2.set_xlabel("Latitude [°]")
    ax2.set_ylabel("uθ")
    ax2.set_title("uθ vs Latitude")
    ax2.grid(True, alpha=0.3)

    ax3.plot(Lat_deg,v, linewidth=2, color='forestgreen')
    ax3.set_xlabel("Latitude [°]")
    ax3.set_ylabel("uλ")
    ax3.set_title("uλ vs Latitude")
    ax3.grid(True, alpha=0.3)

    ax4.plot(Lat_deg, solver.CORIOLIS, linewidth=2, color='darkorchid')
    ax4.set_xlabel("Latitude [°]")
    ax4.set_ylabel("CORIOLIS f")
    ax4.set_title("CORIOLIS vs Latitude")
    ax4.grid(True, alpha=0.3)

    #!ADDED: Title for the entire figure
    fig.suptitle(f"Solver: {solver}", fontsize=16, fontweight='bold')

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    plot_path = os.path.join(run_dir or "output", f"Plots_{solver}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot_path}")
    plt.show()

