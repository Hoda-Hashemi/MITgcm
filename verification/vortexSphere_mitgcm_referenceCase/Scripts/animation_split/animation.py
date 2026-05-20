#%%
from pathlib import Path
import sys
import importlib
import pandas as pd

SCRIPT_DIR = Path.cwd().resolve()

ANIM_DIR = SCRIPT_DIR / "animation_split"
RUN_DIR  = (SCRIPT_DIR / "../../run").resolve()
OUT_DIR  = (SCRIPT_DIR / "../docs").resolve()

sys.path.insert(0, str(ANIM_DIR))

import animationComponents as ac
importlib.reload(ac)

print(ac.__file__)
#%%
# Cell 2: Settings.
FIELD_NAME = "Eta" 
RECORD_NAME = None
RECORD_INDEX = None

CMAP = "RdBu_r"
COLORBAR_LABEL = ""

FIELD_UNITS = {
    "Eta": "m",
    "ETAN": "m",
    "PH": "m^2/s^2",
    "PhiVEL": "m^2/s",
    "PsiVEL": "m^3/s",
    "UVELMASS": "m/s",
    "VVELMASS": "m/s",
    "WVELMASS": "m/s",
}

FIELD_DESCRIPTIONS = {
    "Eta": "Free surface",
    "ETAN": "Free surface",
    "PH": "Geopotential",
    "PhiVEL": "Velocity potential",
    "PsiVEL": "Streamfunction",
}

LATITUDE_ORDER = "south_to_north"
DX_DEG = 0.25
DY_DEG = 0.25

#no skips set it to 1
HTML_SKIP = 1 #2
VIDEO_SKIP = 1 # 2

COLOR_LIMIT = None

WRITE_HTML = True
SHOW_FIGURE = True
OPEN_BROWSER = False
MAKE_VIDEO = True

HTML_FILE = "mitgcm_animation.html"
VIDEO_FILE = "mitgcm_animation.gif"
VIDEO_FALLBACK_FILE = "mitgcm_animation.gif"

FIG_WIDTH = 1600
FIG_HEIGHT = 850
VIDEO_WIDTH = 1800
VIDEO_HEIGHT = 950
VIDEO_FRAME_DURATION_MS = 180
#%%#%%
# Cell 3: Inventory table.
inventory = ac.collect_inventory(RUN_DIR)

static_df = pd.DataFrame({
    "type": "static",
    "variable": inventory["static_variables"],
    "iterations": "",
    "records": "",
})

timed_df = pd.DataFrame([
    {
        "type": "time-dependent",
        "variable": name,
        "iterations": values,
        "records": "",
    }
    for name, values in inventory["timed_variables"].items()
])

dyn_df = pd.DataFrame([
    {
        "type": "dynamic bundle",
        "variable": name,
        "iterations": data["iterations"],
        "records": ", ".join([f"{i}: {rec}" for i, rec in data["records"]]),
    }
    for name, data in inventory["dynamic_records"].items()
])

summary_df = pd.concat([static_df, timed_df, dyn_df], ignore_index=True)

print("RUN_DIR =", RUN_DIR)
display(summary_df)

# Example choices for next cell:
# FIELD_NAME = "Eta"
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PsiVEL"
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PhiVEL"
# FIELD_NAME = "dynDiag"; RECORD_NAME = "UVELMASS"

#%%
# Cell 3: Choose one variable/record from the printed list.
FIELD_NAME = "Eta"
RECORD_NAME = None
RECORD_INDEX = None

# For dynamic bundles, use this style:
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PsiVEL"
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PhiVEL"
# FIELD_NAME = "dynDiag"; RECORD_NAME = "UVELMASS"

#%%
# Cell 4: Load selected field series.
timed = ac.discover_timed_variables(RUN_DIR)
if FIELD_NAME not in timed:
    raise FileNotFoundError(f"{FIELD_NAME} not found. Run Cell 2 and choose one printed variable.")

iterations = timed[FIELD_NAME]
first_field, _ = ac.read_selected_field(RUN_DIR, FIELD_NAME, iterations[0], RECORD_NAME, RECORD_INDEX)
NY, NX = first_field.shape
land_mask = ac.load_land_mask(RUN_DIR, NY, NX)

fields = ac.read_series(
    RUN_DIR,
    FIELD_NAME,
    iterations,
    land_mask,
    HTML_SKIP,
    record_name=RECORD_NAME,
    record_index=RECORD_INDEX,
)

plot_name = RECORD_NAME or FIELD_NAME
unit = FIELD_UNITS.get(plot_name, FIELD_UNITS.get(FIELD_NAME, ""))
field_label = ac.title_name(plot_name, FIELD_DESCRIPTIONS)
colorbar_label = COLORBAR_LABEL or ac.label_with_unit(plot_name, FIELD_UNITS)
zmax = ac.symmetric_limit(fields, COLOR_LIMIT)

lon, lat = ac.grid_1d(NX, NY, DX_DEG, DY_DEG, LATITUDE_ORDER)
lon_p = lon[::HTML_SKIP]
lat_p = lat[::HTML_SKIP]

print(f"Selected: FIELD_NAME={FIELD_NAME}, RECORD_NAME={RECORD_NAME}, RECORD_INDEX={RECORD_INDEX}")
print(f"Iterations: {iterations}")
print(f"Grid: {NY} x {NX}; HTML_SKIP={HTML_SKIP}; VIDEO_SKIP={VIDEO_SKIP}")
print(f"Color scale: +/- {zmax:.6g} {unit}")

#%%
# Cell 5: High-quality slider HTML.
fig = ac.build_slider_figure(
    fields=fields,
    iterations=iterations,
    lon=lon_p,
    lat=lat_p,
    land_mask=land_mask,
    run_dir=RUN_DIR,
    field_label=field_label,
    unit=unit,
    zmax=zmax,
    cmap=CMAP,
    colorbar_label=colorbar_label,
    html_skip=HTML_SKIP,
    width=FIG_WIDTH,
    height=FIG_HEIGHT,
    frame_duration_ms=VIDEO_FRAME_DURATION_MS,
)

html_path = ac.output_path(SCRIPT_DIR, HTML_FILE)
ac.show_or_write_figure(fig, html_path, WRITE_HTML, OPEN_BROWSER, SHOW_FIGURE)

#%%
# Cell 6: High-quality video or GIF.
if MAKE_VIDEO:
    ac.write_video(
        run_dir=RUN_DIR,
        field_name=FIELD_NAME,
        iterations=iterations,
        land_mask=land_mask,
        zmax=zmax,
        cmap=CMAP,
        colorbar_label=colorbar_label,
        title=field_label,
        mp4_path=ac.output_path(SCRIPT_DIR, VIDEO_FILE),
        gif_path=ac.output_path(SCRIPT_DIR, VIDEO_FALLBACK_FILE),
        record_name=RECORD_NAME,
        record_index=RECORD_INDEX,
        video_skip=VIDEO_SKIP,
        frame_duration_ms=VIDEO_FRAME_DURATION_MS,
        width=VIDEO_WIDTH,
        height=VIDEO_HEIGHT,
        latitude_order=LATITUDE_ORDER,
    )

# %%
#%%
# Cell 6: Plotly GIF.

gif_path = ac.output_path(SCRIPT_DIR, "mitgcm_animation_plotly.gif")

ac.write_plotly_gif(
    fig=fig,
    iterations=iterations,
    gif_path=gif_path,
    width=FIG_WIDTH,
    height=FIG_HEIGHT,
    frame_duration_ms=VIDEO_FRAME_DURATION_MS,
)

# %%
#!!!!ADDED: 
#%%
# ============================================================
# MITgcm GIF GENERATOR
#
# Creates:
#   1) prognostic.gif
#   2) diagnostic.gif
#
# Features:
# - publication-style layout
# - massive spacing
# - rectangular panels
# - seismic colormap everywhere
# - fixed global colorbars
# - high 600
# - horizontal separator line
# - large headers
# - per-frame min/max
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

# ============================================================
# SETTINGS
# ============================================================

FPS = 4
CMAP = "seismic"

# ============================================================
# ITERATIONS
# ============================================================

timed = ac.discover_timed_variables(RUN_DIR)

eta_iters = timed["Eta"]
dyn_iters = timed["dyn_Aux"]

# ============================================================
# LOAD PROGNOSTICS
# ============================================================

eta_prog = ac.read_series(
    RUN_DIR,
    "Eta",
    eta_iters,
    land_mask,
    VIDEO_SKIP,
)

u_prog = ac.read_series(
    RUN_DIR,
    "U",
    eta_iters,
    land_mask,
    VIDEO_SKIP,
)

v_prog = ac.read_series(
    RUN_DIR,
    "V",
    eta_iters,
    land_mask,
    VIDEO_SKIP,
)

w_prog = ac.read_series(
    RUN_DIR,
    "W",
    eta_iters,
    land_mask,
    VIDEO_SKIP,
)

mag_prog = [
    np.sqrt(u**2 + v**2)
    for u, v in zip(u_prog, v_prog)
]

# ============================================================
# LOAD DIAGNOSTICS
# ============================================================

eta_diag = ac.read_series(
    RUN_DIR,
    "Eta",
    dyn_iters,
    land_mask,
    VIDEO_SKIP,
)

psi_diag = ac.read_series(
    RUN_DIR,
    "dyn_Aux",
    dyn_iters,
    land_mask,
    VIDEO_SKIP,
    record_name="PsiVEL",
)

phi_diag = ac.read_series(
    RUN_DIR,
    "dyn_Aux",
    dyn_iters,
    land_mask,
    VIDEO_SKIP,
    record_name="PhiVEL",
)

uvel_diag = ac.read_series(
    RUN_DIR,
    "dynDiag",
    dyn_iters,
    land_mask,
    VIDEO_SKIP,
    record_name="UVELMASS",
)

vvel_diag = ac.read_series(
    RUN_DIR,
    "dynDiag",
    dyn_iters,
    land_mask,
    VIDEO_SKIP,
    record_name="VVELMASS",
)

mag_diag = [
    np.sqrt(u**2 + v**2)
    for u, v in zip(uvel_diag, vvel_diag)
]

# ============================================================
# GLOBAL LIMITS
# ============================================================

eta_lim = max(
    np.nanmax(np.abs(f))
    for f in eta_prog + eta_diag
)

u_lim = max(
    np.nanmax(np.abs(f))
    for f in u_prog
)

v_lim = max(
    np.nanmax(np.abs(f))
    for f in v_prog
)

w_lim = max(
    np.nanmax(np.abs(f))
    for f in w_prog
)

mag_prog_lim = max(
    np.nanmax(np.abs(f))
    for f in mag_prog
)

psi_lim = max(
    np.nanmax(np.abs(f))
    for f in psi_diag
)

phi_lim = max(
    np.nanmax(np.abs(f))
    for f in phi_diag
)

mag_diag_lim = max(
    np.nanmax(np.abs(f))
    for f in mag_diag
)

# ============================================================
# ============================================================
# PROGNOSTIC GIF
# ============================================================
# ============================================================

prog_gif_path = ac.output_path(
    SCRIPT_DIR,
    "prognostic.gif"
)

fig, axes = plt.subplots(
    2,
    3,
    figsize=(60,32),
    dpi=300
)

# ============================================================
# SPACING
# ============================================================

fig.subplots_adjust(
    left=0.05,
    right=0.97,
    top=0.78,
    bottom=0.06,
    wspace=0.22,
    hspace=0.40
)

# ============================================================
# HEADER
# ============================================================

fig.text(
    0.5,
    0.975,
    "PROGNOSTIC VARIABLES",
    ha="center",
    va="top",
    fontsize=34,
    fontweight="bold"
)

fig.text(
    0.5,
    0.935,
    "Fields: Eta, U, V, W, |U,V|",
    ha="center",
    va="top",
    fontsize=34
)
# ============================================================
# PRINT ITERATION LIST
# ============================================================

iter_string = ", ".join(str(it) for it in eta_iters)

fig.text(
    0.5,
    0.895,
    f"Iterations: [{iter_string}]",
    ha="center",
    va="top",
    fontsize=34
)

# ============================================================
# HORIZONTAL LINE
# ============================================================

fig.add_artist(
    plt.Line2D(
        [0.08, 0.92],
        [0.875, 0.875],
        transform=fig.transFigure,
        color="black",
        linewidth=2.5
    )
)

# ============================================================
# AXES
# ============================================================

ax1, ax2, ax3 = axes[0]
ax4, ax5, ax6 = axes[1]

# ============================================================
# PLOTS
# ============================================================

im1 = ax1.imshow(
    eta_prog[0],
    cmap=CMAP,
    vmin=-eta_lim,
    vmax=eta_lim,
    origin="lower",
    aspect="auto"
)

im2 = ax2.imshow(
    u_prog[0],
    cmap=CMAP,
    vmin=-u_lim,
    vmax=u_lim,
    origin="lower",
    aspect="auto"
)

im3 = ax3.imshow(
    v_prog[0],
    cmap=CMAP,
    vmin=-v_lim,
    vmax=v_lim,
    origin="lower",
    aspect="auto"
)

im4 = ax4.imshow(
    w_prog[0],
    cmap=CMAP,
    vmin=-w_lim,
    vmax=w_lim,
    origin="lower",
    aspect="auto"
)

im5 = ax5.imshow(
    mag_prog[0],
    cmap=CMAP,
    vmin=0,
    vmax=mag_prog_lim,
    origin="lower",
    aspect="auto"
)

ax6.axis("off")

# ============================================================
# COLORBARS
# ============================================================

for im, ax in zip(
    [im1, im2, im3, im4, im5],
    [ax1, ax2, ax3, ax4, ax5]
):

    fig.colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        fraction=0.045,
        pad=0.08
    )

# ============================================================
# PANEL TITLES
# ============================================================

t1 = ax1.set_title("", fontsize=34, pad=16)
t2 = ax2.set_title("", fontsize=34, pad=16)
t3 = ax3.set_title("", fontsize=34, pad=16)
t4 = ax4.set_title("", fontsize=34, pad=16)
t5 = ax5.set_title("", fontsize=34, pad=16)

# ============================================================
# WRITE GIF
# ============================================================

writer = PillowWriter(fps=FPS)

with writer.saving(
    fig,
    prog_gif_path,
    dpi=300
):

    for i, it in enumerate(eta_iters):

        eta = eta_prog[i]
        u   = u_prog[i]
        v   = v_prog[i]
        w   = w_prog[i]
        mag = mag_prog[i]

        im1.set_data(eta)
        im2.set_data(u)
        im3.set_data(v)
        im4.set_data(w)
        im5.set_data(mag)

        t1.set_text(
            f"Eta | iteration = {it}\n"
            f"min = {np.nanmin(eta):.3e}    "
            f"max = {np.nanmax(eta):.3e}"
        )

        t2.set_text(
            f"U | iteration = {it}\n"
            f"min = {np.nanmin(u):.3e}    "
            f"max = {np.nanmax(u):.3e}"
        )

        t3.set_text(
            f"V | iteration = {it}\n"
            f"min = {np.nanmin(v):.3e}    "
            f"max = {np.nanmax(v):.3e}"
        )

        t4.set_text(
            f"W | iteration = {it}\n"
            f"min = {np.nanmin(w):.3e}    "
            f"max = {np.nanmax(w):.3e}"
        )

        t5.set_text(
            f"|U,V| | iteration = {it}\n"
            f"min = {np.nanmin(mag):.3e}    "
            f"max = {np.nanmax(mag):.3e}"
        )

        writer.grab_frame(facecolor="white")

plt.close()

print("Saved:", prog_gif_path)

# ============================================================
# ============================================================
# DIAGNOSTIC GIF
# ============================================================
# ============================================================

diag_gif_path = ac.output_path(
    SCRIPT_DIR,
    "diagnostic.gif"
)

fig, axes = plt.subplots(
    2,
    2,
    figsize=(60,32),
    dpi=300
)

# ============================================================
# SPACING
# ============================================================

fig.subplots_adjust(
    left=0.05,
    right=0.97,
    top=0.78,
    bottom=0.06,
    wspace=0.22,
    hspace=0.40
)

# ============================================================
# HEADER
# ============================================================

fig.text(
    0.5,
    0.975,
    "DIAGNOSTIC VARIABLES",
    ha="center",
    va="top",
    fontsize=34,
    fontweight="bold"
)

fig.text(
    0.5,
    0.935,
    "Fields: Eta, |UVELMASS,VVELMASS|, PsiVEL, PhiVEL",
    ha="center",
    va="top",
    fontsize=34
)

iter_string = ", ".join(str(it) for it in eta_iters)

fig.text(
    0.5,
    0.895,
    f"Iterations: [{iter_string}]",
    ha="center",
    va="top",
    fontsize=34
)

# ============================================================
# HORIZONTAL LINE
# ============================================================

fig.add_artist(
    plt.Line2D(
        [0.08, 0.92],
        [0.875, 0.875],
        transform=fig.transFigure,
        color="black",
        linewidth=2.5
    )
)

# ============================================================
# AXES
# ============================================================

ax1, ax2 = axes[0]
ax3, ax4 = axes[1]

# ============================================================
# PLOTS
# ============================================================

im1 = ax1.imshow(
    eta_diag[0],
    cmap=CMAP,
    vmin=-eta_lim,
    vmax=eta_lim,
    origin="lower",
    aspect="auto"
)

im2 = ax2.imshow(
    mag_diag[0],
    cmap=CMAP,
    vmin=0,
    vmax=mag_diag_lim,
    origin="lower",
    aspect="auto"
)

im3 = ax3.imshow(
    psi_diag[0],
    cmap=CMAP,
    vmin=-psi_lim,
    vmax=psi_lim,
    origin="lower",
    aspect="auto"
)

im4 = ax4.imshow(
    phi_diag[0],
    cmap=CMAP,
    vmin=-phi_lim,
    vmax=phi_lim,
    origin="lower",
    aspect="auto"
)

# ============================================================
# COLORBARS
# ============================================================

for im, ax in zip(
    [im1, im2, im3, im4],
    [ax1, ax2, ax3, ax4]
):

    fig.colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        fraction=0.045,
        pad=0.08
    )

# ============================================================
# PANEL TITLES
# ============================================================

t1 = ax1.set_title("", fontsize=34, pad=16)
t2 = ax2.set_title("", fontsize=34, pad=16)
t3 = ax3.set_title("", fontsize=34, pad=16)
t4 = ax4.set_title("", fontsize=34, pad=16)

# ============================================================
# WRITE GIF
# ============================================================

writer = PillowWriter(fps=FPS)

with writer.saving(
    fig,
    diag_gif_path,
    dpi=300
):

    for i, it in enumerate(dyn_iters):

        eta = eta_diag[i]
        mag = mag_diag[i]
        psi = psi_diag[i]
        phi = phi_diag[i]

        im1.set_data(eta)
        im2.set_data(mag)
        im3.set_data(psi)
        im4.set_data(phi)

        t1.set_text(
            f"Eta | iteration = {it}\n"
            f"min = {np.nanmin(eta):.3e}    "
            f"max = {np.nanmax(eta):.3e}"
        )

        t2.set_text(
            f"|UVELMASS,VVELMASS| | iteration = {it}\n"
            f"min = {np.nanmin(mag):.3e}    "
            f"max = {np.nanmax(mag):.3e}"
        )

        t3.set_text(
            f"PsiVEL | iteration = {it}\n"
            f"min = {np.nanmin(psi):.3e}    "
            f"max = {np.nanmax(psi):.3e}"
        )

        t4.set_text(
            f"PhiVEL | iteration = {it}\n"
            f"min = {np.nanmin(phi):.3e}    "
            f"max = {np.nanmax(phi):.3e}"
        )

        writer.grab_frame(facecolor="white")

plt.close()

print("Saved:", diag_gif_path)

# %%
#%%
# ============================================================
# HIGH-QUALITY WEBM VIDEO
#
# Output:
#   prognostic.webm
#
# Requirements:
#   brew install ffmpeg
#
# ============================================================

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FFMpegWriter

# ============================================================
# SETTINGS
# ============================================================

FPS = 8
CMAP = "seismic"

VIDEO_PATH = ac.output_path(
    SCRIPT_DIR,
    "prognostic.webm"
)

# ============================================================
# GLOBAL LIMITS
# ============================================================

eta_lim = max(
    np.nanmax(np.abs(f))
    for f in eta_prog
)

u_lim = max(
    np.nanmax(np.abs(f))
    for f in u_prog
)

v_lim = max(
    np.nanmax(np.abs(f))
    for f in v_prog
)

w_lim = max(
    np.nanmax(np.abs(f))
    for f in w_prog
)

mag_prog_lim = max(
    np.nanmax(np.abs(f))
    for f in mag_prog
)

# ============================================================
# FIGURE
# ============================================================

fig, axes = plt.subplots(
    2,
    3,
    figsize=(24, 16),
    dpi=300
)

# ============================================================
# SPACING
# ============================================================

fig.subplots_adjust(
    left=0.05,
    right=0.97,
    top=0.78,
    bottom=0.06,
    wspace=0.22,
    hspace=0.40
)

# ============================================================
# HEADER
# ============================================================

fig.text(
    0.5,
    0.975,
    "PROGNOSTIC VARIABLES",
    ha="center",
    va="top",
    fontsize=34,
    fontweight="bold"
)

fig.text(
    0.5,
    0.940,
    "Fields: Eta, U, V, W, |U,V|",
    ha="center",
    va="top",
    fontsize=18
)

# ============================================================
# ITERATIONS
# ============================================================

iter_string = ", ".join(str(it) for it in eta_iters)

fig.text(
    0.5,
    0.915,
    f"Iterations: [{iter_string}]",
    ha="center",
    va="top",
    fontsize=13
)

# ============================================================
# HORIZONTAL LINE
# ============================================================

fig.add_artist(
    plt.Line2D(
        [0.08, 0.92],
        [0.885, 0.885],
        transform=fig.transFigure,
        color="black",
        linewidth=2.0
    )
)
# ============================================================
# AXES
# ============================================================

ax1, ax2, ax3 = axes[0]
ax4, ax5, ax6 = axes[1]

# ============================================================
# PLOTS
# ============================================================

im1 = ax1.imshow(
    eta_prog[0],
    cmap=CMAP,
    vmin=-eta_lim,
    vmax=eta_lim,
    origin="lower",
    aspect="auto"
)

im2 = ax2.imshow(
    u_prog[0],
    cmap=CMAP,
    vmin=-u_lim,
    vmax=u_lim,
    origin="lower",
    aspect="auto"
)

im3 = ax3.imshow(
    v_prog[0],
    cmap=CMAP,
    vmin=-v_lim,
    vmax=v_lim,
    origin="lower",
    aspect="auto"
)

im4 = ax4.imshow(
    w_prog[0],
    cmap=CMAP,
    vmin=-w_lim,
    vmax=w_lim,
    origin="lower",
    aspect="auto"
)

im5 = ax5.imshow(
    mag_prog[0],
    cmap=CMAP,
    vmin=0,
    vmax=mag_prog_lim,
    origin="lower",
    aspect="auto"
)

# ============================================================
# INFO PANEL
# ============================================================

# ax6.text(
#     0.5,
#     0.5,
#     (
#         "PROGNOSTIC VARIABLES\n\n"
#         "Eta\n"
#         "U\n"
#         "V\n"
#         "W\n"
#         "|U,V|\n\n"
#         f"{len(eta_iters)} iterations"
#     ),
#     ha="center",
#     va="center",
#     fontsize=18,
#     fontweight="bold",
#     linespacing=1.6,
#     bbox=dict(
#         facecolor="white",
#         edgecolor="black",
#         boxstyle="round,pad=1.2"
#     )
# )

# ax6.set_xticks([])
# ax6.set_yticks([])
ax6.axis("off")
# ============================================================
# COLORBARS
# ============================================================

for im, ax in zip(
    [im1, im2, im3, im4, im5],
    [ax1, ax2, ax3, ax4, ax5]
):

    fig.colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        fraction=0.045,
        pad=0.08
    )

# ============================================================
# PANEL TITLES
# ============================================================

t1 = ax1.set_title("", fontsize=15, pad=16)
t2 = ax2.set_title("", fontsize=15, pad=16)
t3 = ax3.set_title("", fontsize=15, pad=16)
t4 = ax4.set_title("", fontsize=15, pad=16)
t5 = ax5.set_title("", fontsize=15, pad=16)

# ============================================================
# VIDEO WRITER
# ============================================================

writer = FFMpegWriter(
    fps=FPS,
    codec="libvpx-vp9",
    bitrate=20000
)

# ============================================================
# WRITE VIDEO
# ============================================================

with writer.saving(
    fig,
    VIDEO_PATH,
    dpi=300
):

    for i, it in enumerate(eta_iters):

        eta = eta_prog[i]
        u   = u_prog[i]
        v   = v_prog[i]
        w   = w_prog[i]
        mag = mag_prog[i]

        # ----------------------------------------------------
        # UPDATE DATA
        # ----------------------------------------------------

        im1.set_data(eta)
        im2.set_data(u)
        im3.set_data(v)
        im4.set_data(w)
        im5.set_data(mag)

        # ----------------------------------------------------
        # TITLES
        # ----------------------------------------------------

        t1.set_text(
            f"Eta | iteration = {it}\n"
            f"min = {np.nanmin(eta):.3e}    "
            f"max = {np.nanmax(eta):.3e}"
        )

        t2.set_text(
            f"U | iteration = {it}\n"
            f"min = {np.nanmin(u):.3e}    "
            f"max = {np.nanmax(u):.3e}"
        )

        t3.set_text(
            f"V | iteration = {it}\n"
            f"min = {np.nanmin(v):.3e}    "
            f"max = {np.nanmax(v):.3e}"
        )

        t4.set_text(
            f"W | iteration = {it}\n"
            f"min = {np.nanmin(w):.3e}    "
            f"max = {np.nanmax(w):.3e}"
        )

        t5.set_text(
            f"|U,V| | iteration = {it}\n"
            f"min = {np.nanmin(mag):.3e}    "
            f"max = {np.nanmax(mag):.3e}"
        )

        writer.grab_frame()

plt.close()

print("Saved:", VIDEO_PATH)

# %%
#!Diagnostic:
#%%
# ============================================================
# HIGH-QUALITY DIAGNOSTIC WEBM
#
# Output:
#   diagnostic.webm
#
# Requirements:
#   brew install ffmpeg
# ============================================================

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FFMpegWriter

# ============================================================
# SETTINGS
# ============================================================

FPS = 4
CMAP = "seismic"

VIDEO_PATH = ac.output_path(
    SCRIPT_DIR,
    "diagnostic.webm"
)

# ============================================================
# GLOBAL LIMITS
# ============================================================

eta_lim = max(
    np.nanmax(np.abs(f))
    for f in eta_diag
)

psi_lim = max(
    np.nanmax(np.abs(f))
    for f in psi_diag
)

phi_lim = max(
    np.nanmax(np.abs(f))
    for f in phi_diag
)

mag_diag_lim = max(
    np.nanmax(np.abs(f))
    for f in mag_diag
)

# ============================================================
# FIGURE
# ============================================================

fig, axes = plt.subplots(
    2,
    3,
    figsize=(24, 16),
    dpi=300
)

# ============================================================
# SPACING
# ============================================================

fig.subplots_adjust(
    left=0.05,
    right=0.97,
    top=0.78,
    bottom=0.06,
    wspace=0.22,
    hspace=0.40
)

# ============================================================
# HEADER
# ============================================================

fig.text(
    0.5,
    0.975,
    "DIAGNOSTIC VARIABLES",
    ha="center",
    va="top",
    fontsize=34,
    fontweight="bold"
)

fig.text(
    0.5,
    0.940,
    "Fields: Eta, |UVELMASS,VVELMASS|, PsiVEL, PhiVEL",
    ha="center",
    va="top",
    fontsize=18
)

# ============================================================
# ITERATIONS
# ============================================================

iter_string = ", ".join(str(it) for it in dyn_iters)

fig.text(
    0.5,
    0.915,
    f"Iterations: [{iter_string}]",
    ha="center",
    va="top",
    fontsize=13
)

# ============================================================
# HORIZONTAL LINE
# ============================================================

fig.add_artist(
    plt.Line2D(
        [0.08, 0.92],
        [0.885, 0.885],
        transform=fig.transFigure,
        color="black",
        linewidth=2.0
    )
)

# ============================================================
# AXES
# ============================================================

ax1, ax2, ax3 = axes[0]
ax4, ax5, ax6 = axes[1]

# ============================================================
# PLOTS
# ============================================================

im1 = ax1.imshow(
    eta_diag[0],
    cmap=CMAP,
    vmin=-eta_lim,
    vmax=eta_lim,
    origin="lower",
    aspect="auto"
)

im2 = ax2.imshow(
    mag_diag[0],
    cmap=CMAP,
    vmin=0,
    vmax=mag_diag_lim,
    origin="lower",
    aspect="auto"
)

im3 = ax3.imshow(
    psi_diag[0],
    cmap=CMAP,
    vmin=-psi_lim,
    vmax=psi_lim,
    origin="lower",
    aspect="auto"
)

im4 = ax4.imshow(
    phi_diag[0],
    cmap=CMAP,
    vmin=-phi_lim,
    vmax=phi_lim,
    origin="lower",
    aspect="auto"
)

# ============================================================
# EMPTY PANELS
# ============================================================

ax5.axis("off")
ax6.axis("off")

# ============================================================
# COLORBARS
# ============================================================

for im, ax in zip(
    [im1, im2, im3, im4],
    [ax1, ax2, ax3, ax4]
):

    fig.colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        fraction=0.045,
        pad=0.08
    )

# ============================================================
# PANEL TITLES
# ============================================================

t1 = ax1.set_title(
    "",
    fontsize=15,
    pad=16
)

t2 = ax2.set_title(
    "",
    fontsize=15,
    pad=16
)

t3 = ax3.set_title(
    "",
    fontsize=15,
    pad=16
)

t4 = ax4.set_title(
    "",
    fontsize=15,
    pad=16
)

# ============================================================
# VIDEO WRITER
# ============================================================

writer = FFMpegWriter(
    fps=FPS,
    codec="libvpx-vp9",
    bitrate=20000
)

# ============================================================
# WRITE VIDEO
# ============================================================

with writer.saving(
    fig,
    VIDEO_PATH,
    dpi=300
):

    for i, it in enumerate(dyn_iters):

        eta = eta_diag[i]
        mag = mag_diag[i]
        psi = psi_diag[i]
        phi = phi_diag[i]

        # ----------------------------------------------------
        # UPDATE DATA
        # ----------------------------------------------------

        im1.set_data(eta)
        im2.set_data(mag)
        im3.set_data(psi)
        im4.set_data(phi)

        # ----------------------------------------------------
        # TITLES
        # ----------------------------------------------------

        t1.set_text(
            f"Eta | iteration = {it}\n"
            f"min = {np.nanmin(eta):.3e}    "
            f"max = {np.nanmax(eta):.3e}"
        )

        t2.set_text(
            f"|UVELMASS,VVELMASS| | iteration = {it}\n"
            f"min = {np.nanmin(mag):.3e}    "
            f"max = {np.nanmax(mag):.3e}"
        )

        t3.set_text(
            f"PsiVEL | iteration = {it}\n"
            f"min = {np.nanmin(psi):.3e}    "
            f"max = {np.nanmax(psi):.3e}"
        )

        t4.set_text(
            f"PhiVEL | iteration = {it}\n"
            f"min = {np.nanmin(phi):.3e}    "
            f"max = {np.nanmax(phi):.3e}"
        )

        writer.grab_frame()

plt.close()

print("Saved:", VIDEO_PATH)

# %%
