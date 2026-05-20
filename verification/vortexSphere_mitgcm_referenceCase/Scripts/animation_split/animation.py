#%%
# Cell 1: Imports and paths.
from pathlib import Path
import sys
import importlib
import pandas as pd

ANIM_DIR = Path("/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/animation_split").resolve()

print("ANIM_DIR =", ANIM_DIR)
print("exists =", ANIM_DIR.exists())
print("component exists =", (ANIM_DIR / "animationComponents.py").exists())

sys.path.insert(0, str(ANIM_DIR))

import importlib
import animationComponents as ac
importlib.reload(ac)

print("loaded from:", ac.__file__)

SCRIPT_DIR = ANIM_DIR
RUN_DIR = (SCRIPT_DIR / "../../run").resolve()
OUT_DIR = (SCRIPT_DIR / "../../docs").resolve()

print("SCRIPT_DIR =", SCRIPT_DIR)
print("RUN_DIR    =", RUN_DIR)
print("RUN exists =", RUN_DIR.exists())
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

