#%%
#! Cell 1: Imports and paths

from pathlib import Path
import sys
import importlib
import re
import pandas as pd
import numpy as np
try:
    from IPython.display import display
except Exception:
    def display(x):
        print(x)

CWD = Path.cwd().resolve()
SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else CWD
ANIM_DIR = SCRIPT_DIR if (SCRIPT_DIR / "animationComponents.py").exists() else CWD / "animation_split"

sys.path.insert(0, str(ANIM_DIR))

import animationComponents as ac
importlib.reload(ac)

SCRIPT_DIR = ANIM_DIR.resolve()
RUN_DIR = (SCRIPT_DIR.parent.parent / "run_2").resolve()
OUT_DIR = (SCRIPT_DIR.parent / "docs").resolve()

print("animationComponents =", ac.__file__)
print("SCRIPT_DIR =", SCRIPT_DIR)
print("RUN_DIR =", RUN_DIR)
print("OUT_DIR =", OUT_DIR)

#%%
#! Cell 2: Settings

FIELD_NAME = "Eta"
RECORD_NAME = None
RECORD_INDEX = None

CMAP = "RdBu_r"
COLORBAR_LABEL = ""
EXPERIMENT_TITLE = "Gaussian Free-Surface Patch with Real Bathymetry"
# "Gaussian Free-Surface Patch with Constant Bathymetry d = -4000 m"
#! Simulation time settings.
NTIMESTEPS=14400
DELTA_T_SEC=60

#%%

FIELD_UNITS = {
    "Eta": "m",
    "ETAN": "m",
    "PH": "m^2/s^2",
    "PhiVEL": "m^2/s",
    "PsiVEL": "m^3/s",
    "UVELMASS": "m/s",
    "VVELMASS": "m/s",
    "WVELMASS": "m/s",
    "U": "m/s",
    "V": "m/s",
    "W": "m/s",
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

HTML_SKIP = 1
VIDEO_SKIP = 1
COLOR_LIMIT = None

WRITE_HTML = True
SHOW_FIGURE = True
OPEN_BROWSER = False
MAKE_VIDEO = True
MAKE_PLOTLY_GIF = False
MAKE_PANEL_WEBM = True

HTML_FILE = "mitgcm_animation.html"
VIDEO_FILE = "mitgcm_animation.gif"
VIDEO_FALLBACK_FILE = "mitgcm_animation.gif"

FIG_WIDTH = 1600
FIG_HEIGHT = 850
VIDEO_WIDTH = 1800
VIDEO_HEIGHT = 950
VIDEO_FRAME_DURATION_MS = 180

safe_title = re.sub(r"[^A-Za-z0-9._-]+", "_", EXPERIMENT_TITLE).strip("_")
ANIMATION_DIR = SCRIPT_DIR / safe_title
ANIMATION_DIR.mkdir(parents=True, exist_ok=True)

print("ANIMATION_DIR =", ANIMATION_DIR)

#%%
#! Cell 3: Inventory table

inventory = ac.collect_inventory(RUN_DIR)

static_df = pd.DataFrame({
    "type": "static",
    "variable": inventory["static_variables"],
    "iterations": "",
    "records": "",
})

timed_df = pd.DataFrame([
    {"type": "time-dependent", "variable": name, "iterations": values, "records": ""}
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

# Examples:
# FIELD_NAME = "Eta"
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PsiVEL"
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PhiVEL"
# FIELD_NAME = "dynDiag"; RECORD_NAME = "UVELMASS"

 #%%
#! Cell 4: Choose one variable/record

FIELD_NAME = "Eta"
RECORD_NAME = None
RECORD_INDEX = None

# For dynamic bundles:
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PsiVEL"
# FIELD_NAME = "dyn_Aux"; RECORD_NAME = "PhiVEL"
# FIELD_NAME = "dynDiag"; RECORD_NAME = "UVELMASS"

#%%
#! Cell 5: Load selected field series

timed = ac.discover_timed_variables(RUN_DIR)
if FIELD_NAME not in timed:
    raise FileNotFoundError(f"{FIELD_NAME} not found. Run Cell 3 and choose one printed variable.")

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
#! Cell 6: High-quality slider HTML

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

html_path = ac.output_path(OUT_DIR, HTML_FILE)
ac.show_or_write_figure(fig, html_path, WRITE_HTML, OPEN_BROWSER, SHOW_FIGURE)

#%%
#! Cell 7: High-quality video or GIF

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
        mp4_path=ac.output_path(ANIMATION_DIR, VIDEO_FILE),
        gif_path=ac.output_path(ANIMATION_DIR, VIDEO_FALLBACK_FILE),
        record_name=RECORD_NAME,
        record_index=RECORD_INDEX,
        video_skip=VIDEO_SKIP,
        frame_duration_ms=VIDEO_FRAME_DURATION_MS,
        width=VIDEO_WIDTH,
        height=VIDEO_HEIGHT,
        latitude_order=LATITUDE_ORDER,
    )

#%%
#! Cell 8: Plotly GIF

if MAKE_PLOTLY_GIF:
    gif_path = ac.output_path(ANIMATION_DIR, "mitgcm_animation_plotly.gif")
    ac.write_plotly_gif(
        fig=fig,
        iterations=iterations,
        gif_path=gif_path,
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        frame_duration_ms=VIDEO_FRAME_DURATION_MS,
    )

#%%
#! Cell 9: Multi-panel prognostic and diagnostic WebM

if MAKE_PANEL_WEBM:
    timed = ac.discover_timed_variables(RUN_DIR)

    eta_iters = ac.common_iterations(timed, ["Eta", "U", "V", "W"])
    if eta_iters:
        eta_prog = ac.read_variable_series(RUN_DIR, "Eta", eta_iters, land_mask)
        u_prog = ac.read_variable_series(RUN_DIR, "U", eta_iters, land_mask)
        v_prog = ac.read_variable_series(RUN_DIR, "V", eta_iters, land_mask)
        w_prog = ac.read_variable_series(RUN_DIR, "W", eta_iters, land_mask)
        mag_prog = [np.sqrt(u*u + v*v) for u, v in zip(u_prog, v_prog)]
        prognostic_fields = {"Eta": eta_prog, "U": u_prog, "V": v_prog, "W": w_prog, "|U,V|": mag_prog}
        ac.make_webm_adaptive("PROGNOSTIC VARIABLES",prognostic_fields,eta_iters,8,ANIMATION_DIR/"prognostic.webm",CMAP,EXPERIMENT_TITLE,NTIMESTEPS,DELTA_T_SEC)  
        #  ac.make_webm("PROGNOSTIC VARIABLES",prognostic_fields,eta_iters,8,ANIMATION_DIR/"prognostic.webm",CMAP,EXPERIMENT_TITLE,NTIMESTEPS,DELTA_T_SEC)    
    else:
        print("Skipped PROGNOSTIC VARIABLES: common Eta/U/V/W iterations not found.")

    diag_records = ["UVELMASS", "VVELMASS", "PsiVEL", "PhiVEL"]
    dyn_iters, sources = ac.common_iterations_for_records(RUN_DIR, diag_records, extra_vars=("Eta",))
    print("Diagnostic record sources =", sources)

    if dyn_iters:
        eta_diag = ac.read_variable_series(RUN_DIR, "Eta", dyn_iters, land_mask)
        diagnostic_fields = {"Eta": eta_diag}

        try:
            u_diag, _ = ac.read_record_series_any_bundle(RUN_DIR, "UVELMASS", dyn_iters, land_mask)
            v_diag, _ = ac.read_record_series_any_bundle(RUN_DIR, "VVELMASS", dyn_iters, land_mask)
            diagnostic_fields["|UVELMASS,VVELMASS|"] = [np.sqrt(u*u + v*v) for u, v in zip(u_diag, v_diag)]
        except Exception as exc:
            print("Skipped |UVELMASS,VVELMASS|:", exc)

        for rec_name in ["PsiVEL", "PhiVEL"]:
            try:
                diagnostic_fields[rec_name], bundle = ac.read_record_series_any_bundle(RUN_DIR, rec_name, dyn_iters, land_mask)
                print(f"{rec_name} source = {bundle}")
            except Exception as exc:
                print(f"Skipped {rec_name}:", exc)
        ac.make_webm_adaptive("DIAGNOSTIC VARIABLES",diagnostic_fields,dyn_iters,4,ANIMATION_DIR/"diagnostic.webm",CMAP,EXPERIMENT_TITLE,NTIMESTEPS,DELTA_T_SEC)

        # ac.make_webm("DIAGNOSTIC VARIABLES",diagnostic_fields,dyn_iters,4,ANIMATION_DIR/"diagnostic.webm",CMAP,EXPERIMENT_TITLE,NTIMESTEPS,DELTA_T_SEC)
    else:
        print("Skipped DIAGNOSTIC VARIABLES: no common diagnostic iterations found.")
# %%
