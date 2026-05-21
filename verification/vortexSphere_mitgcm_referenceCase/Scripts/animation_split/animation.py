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
#%%
import os,re,glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

SCRIPT_DIR=os.getcwd()
RUN_DIR=os.path.abspath(os.path.join(SCRIPT_DIR,"..","run"))

Nx=1440
Ny=720
dtype=">f8"
CMAP="seismic"
EXPERIMENT_TITLE="Gaussian Patch_ Constant Bathymetry Depth = 4000 m"

if os.path.basename(SCRIPT_DIR) in ["animation_split","animationsplit"]:
    ANIMATION_ROOT=SCRIPT_DIR
else:
    ANIMATION_ROOT=os.path.join(SCRIPT_DIR,"animation_split")

safe_title=re.sub(r"[^A-Za-z0-9._-]+","_",EXPERIMENT_TITLE).strip("_")
ANIMATION_DIR=os.path.join(ANIMATION_ROOT,safe_title)
os.makedirs(ANIMATION_DIR,exist_ok=True)

def find_iters(name):
    files=glob.glob(os.path.join(RUN_DIR,f"{name}.*.data"))
    return sorted([int(os.path.basename(f).split(".")[1]) for f in files])

def common_iters(*names):
    sets=[set(find_iters(name)) for name in names]
    if any(len(s)==0 for s in sets):
        for name,s in zip(names,sets):
            print(name,sorted(s))
        raise ValueError("One variable has no files. Check RUN_DIR or file names.")
    return sorted(set.intersection(*sets))

def read_field(name,it):
    path=os.path.join(RUN_DIR,f"{name}.{it:010d}.data")
    return np.fromfile(path,dtype=dtype).reshape(Ny,Nx)

def read_dynDiag(it,nrecords=3):
    path=os.path.join(RUN_DIR,f"dynDiag.{it:010d}.data")
    return np.fromfile(path,dtype=dtype).reshape(nrecords,Ny,Nx)

eta_iters=common_iters("Eta","U","V","W")
dyn_iters=common_iters("dynDiag","Eta")

eta_prog=[read_field("Eta",it) for it in eta_iters]
u_prog=[read_field("U",it) for it in eta_iters]
v_prog=[read_field("V",it) for it in eta_iters]
w_prog=[read_field("W",it) for it in eta_iters]
mag_prog=[np.sqrt(u*u+v*v) for u,v in zip(u_prog,v_prog)]

eta_diag=[read_field("Eta",it) for it in dyn_iters]
dyn_raw=[read_dynDiag(it,3) for it in dyn_iters]
mag_diag=[f[0] for f in dyn_raw]
psi_diag=[f[1] for f in dyn_raw]
phi_diag=[f[2] for f in dyn_raw]

def lim_abs(fields):
    if len(fields)==0:
        raise ValueError("Empty field list.")
    return max(np.nanmax(np.abs(f)) for f in fields)

def make_webm(kind,fields,iters,fps,outname):
    names=list(fields.keys())
    limits={k:lim_abs(v) for k,v in fields.items()}
    VIDEO_PATH=os.path.join(ANIMATION_DIR,outname)
    fig,axes=plt.subplots(2,3,figsize=(24,16),dpi=300)
    fig.subplots_adjust(left=0.05,right=0.97,top=0.70,bottom=0.06,wspace=0.22,hspace=0.40)
    fig.text(0.5,0.975,f"{kind}\n{EXPERIMENT_TITLE}",ha="center",va="top",fontsize=34,linespacing=1.45,fontweight="bold")
    fig.text(0.5,0.865,"Fields: "+", ".join(names),ha="center",va="top",fontsize=18)
    fig.text(0.5,0.830,f"Iterations: [{', '.join(str(it) for it in iters)}]",ha="center",va="top",fontsize=13)
    fig.add_artist(plt.Line2D([0.08,0.92],[0.800,0.800],transform=fig.transFigure,color="black",linewidth=2.0))
    ax=axes.ravel()
    ims=[]
    titles=[]
    for j,name in enumerate(names):
        vmin=0 if "|" in name else -limits[name]
        vmax=limits[name]
        im=ax[j].imshow(fields[name][0],cmap=CMAP,vmin=vmin,vmax=vmax,origin="lower",aspect="auto")
        ims.append(im)
        titles.append(ax[j].set_title("",fontsize=15,pad=16))
        fig.colorbar(im,ax=ax[j],orientation="horizontal",fraction=0.045,pad=0.08)
    for j in range(len(names),6):
        ax[j].axis("off")
    writer=FFMpegWriter(fps=fps,codec="libvpx-vp9",bitrate=20000)
    with writer.saving(fig,VIDEO_PATH,dpi=300):
        for i,it in enumerate(iters):
            for name,im,t in zip(names,ims,titles):
                f=fields[name][i]
                im.set_data(f)
                t.set_text(f"{name} | iteration = {it}\nmin = {np.nanmin(f):.3e}    max = {np.nanmax(f):.3e}")
            writer.grab_frame()
    plt.close()
    print("Saved:",VIDEO_PATH)

print("RUN_DIR =",RUN_DIR)
print("ANIMATION_DIR =",ANIMATION_DIR)
print("Prognostic iterations =",eta_iters)
print("Diagnostic iterations =",dyn_iters)

make_webm("PROGNOSTIC VARIABLES",{"Eta":eta_prog,"U":u_prog,"V":v_prog,"W":w_prog,"|U,V|":mag_prog},eta_iters,8,"prognostic.webm")
make_webm("DIAGNOSTIC VARIABLES",{"Eta":eta_diag,"|UVELMASS,VVELMASS|":mag_diag,"PsiVEL":psi_diag,"PhiVEL":phi_diag},dyn_iters,4,"diagnostic.webm")
# %%

