#! Cell 1: Imports and paths
#%%

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
OUT_DIR = (SCRIPT_DIR.parent / "docs").resolve()

RUN_CONFIGS = {
    "run_1": {
        "experiment_title": "Geostrophic Adjustment with Real Bathymetry",
        "n_timesteps": 960,
        "delta_t_sec": 900,
    },
    "run_2": {
        "experiment_title": "Gaussian Free-Surface Patch with Real Bathymetry",
        "n_timesteps": 14400,
        "delta_t_sec": 60,
    },
    "run_3": {
        "experiment_title": "Gaussian Free-Surface Patch with Constant Bathymetry d = -4000 m",
        "n_timesteps": 14400,
        "delta_t_sec": 60,
    },
}

RUN_NAME = "run_3"
# "run_2"
# "run_3"
RUN_DIR = (SCRIPT_DIR.parent.parent / RUN_NAME).resolve()

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
RUN_CONFIG = RUN_CONFIGS[RUN_NAME]
EXPERIMENT_TITLE = RUN_CONFIG["experiment_title"]
NTIMESTEPS = RUN_CONFIG["n_timesteps"]
DELTA_T_SEC = RUN_CONFIG["delta_t_sec"]



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
MAKE_SNAPSHOTS = True

HTML_FILE = "mitgcm_animation.html"
VIDEO_FILE = "mitgcm_animation.gif"
VIDEO_FALLBACK_FILE = "mitgcm_animation.gif"
SNAPSHOT_DAYS = [0, 1,  3, 6, 9, 10]

FIG_WIDTH = 1600
FIG_HEIGHT = 850
VIDEO_WIDTH = 1800
VIDEO_HEIGHT = 950
VIDEO_FRAME_DURATION_MS = 180
SNAPSHOT_DPI = 300
SINGLE_FIG_DPI = 320
SCALING_RATIO_THRESHOLD = 20.0
SINGLE_HEATMAP_SCALE_BY_FIELD = None
SCALING_ANALYSIS_DONE = False

safe_title = re.sub(r"[^A-Za-z0-9._-]+", "_", EXPERIMENT_TITLE).strip("_")
ANIMATION_DIR = SCRIPT_DIR / safe_title
ANIMATION_DIR.mkdir(parents=True, exist_ok=True)

print("ANIMATION_DIR =", ANIMATION_DIR)

#%%
#! Cell 2.5: Analyze whether single-heatmap colorbars should be constant or adaptive

print("!! Reminder: single-snapshot scaling can be misleading if you change runs and do not re-run this cell.")
print("How to decide:")
print(f"  constant: use when snapshot amplitudes stay within about {SCALING_RATIO_THRESHOLD:g}x of each other")
print(f"  adaptive: use when amplitudes spread by more than about {SCALING_RATIO_THRESHOLD:g}x and early snapshots get washed out")

def recommend_scale_modes(fields, threshold_ratio):
    recommendations = {}
    print("\nScale analysis:")
    for name, series in fields.items():
        peaks = [float(np.nanmax(np.abs(field))) for field in series]
        positive_peaks = [peak for peak in peaks if np.isfinite(peak) and peak > 0.0]
        if not positive_peaks:
            ratio = 1.0
        else:
            ratio = max(positive_peaks) / min(positive_peaks)
        mode = "adaptive" if ratio > threshold_ratio else "constant"
        recommendations[name] = mode
        print(
            f"  {name}: min_peak={min(positive_peaks) if positive_peaks else 0.0:.3e}, "
            f"max_peak={max(positive_peaks) if positive_peaks else 0.0:.3e}, "
            f"ratio={ratio:.3e} -> {mode}"
        )
    return recommendations

analysis_scale_by_field = {}

prog_iters_analysis = ac.common_iterations(ac.discover_timed_variables(RUN_DIR), ["Eta", "U", "V", "W"])
if prog_iters_analysis:
    prog_snapshots_analysis = ac.pick_snapshot_iterations(prog_iters_analysis, DELTA_T_SEC, SNAPSHOT_DAYS)
    prog_iter_list_analysis = [snap["iteration"] for snap in prog_snapshots_analysis]
    eta0_analysis, _ = ac.read_selected_field(RUN_DIR, "Eta", prog_iter_list_analysis[0])
    prog_land_mask_analysis = ac.load_land_mask(RUN_DIR, *eta0_analysis.shape)
    eta_analysis = ac.read_variable_series(RUN_DIR, "Eta", prog_iter_list_analysis, prog_land_mask_analysis)
    u_analysis = ac.read_variable_series(RUN_DIR, "U", prog_iter_list_analysis, prog_land_mask_analysis)
    v_analysis = ac.read_variable_series(RUN_DIR, "V", prog_iter_list_analysis, prog_land_mask_analysis)
    mag_analysis = [np.sqrt(u * u + v * v) for u, v in zip(u_analysis, v_analysis)]
    analysis_scale_by_field.update(
        recommend_scale_modes(
            {"Eta": eta_analysis, "U": u_analysis, "V": v_analysis, "|U,V|": mag_analysis},
            SCALING_RATIO_THRESHOLD,
        )
    )
else:
    print("No prognostic analysis fields found.")

diag_iters_analysis, diag_sources_analysis = ac.common_iterations_for_records(
    RUN_DIR,
    ["UVELMASS", "VVELMASS", "PsiVEL", "PhiVEL"],
)
print("\nDiagnostic analysis sources =", diag_sources_analysis)
if diag_iters_analysis:
    diag_snapshots_analysis = ac.pick_snapshot_iterations(diag_iters_analysis, DELTA_T_SEC, SNAPSHOT_DAYS)
    diag_iter_list_analysis = [snap["iteration"] for snap in diag_snapshots_analysis]
    diag0_analysis, _ = ac.read_record_any_bundle(RUN_DIR, "UVELMASS", diag_iter_list_analysis[0])
    diag_land_mask_analysis = ac.load_land_mask(RUN_DIR, *diag0_analysis.shape)
    uvelmass_analysis, _ = ac.read_record_series_any_bundle(RUN_DIR, "UVELMASS", diag_iter_list_analysis, diag_land_mask_analysis)
    vvelmass_analysis, _ = ac.read_record_series_any_bundle(RUN_DIR, "VVELMASS", diag_iter_list_analysis, diag_land_mask_analysis)
    psivel_analysis, _ = ac.read_record_series_any_bundle(RUN_DIR, "PsiVEL", diag_iter_list_analysis, diag_land_mask_analysis)
    phivel_analysis, _ = ac.read_record_series_any_bundle(RUN_DIR, "PhiVEL", diag_iter_list_analysis, diag_land_mask_analysis)
    transport_mag_analysis = [np.sqrt(u * u + v * v) for u, v in zip(uvelmass_analysis, vvelmass_analysis)]
    analysis_scale_by_field.update(
        recommend_scale_modes(
            {
                "UVELMASS": uvelmass_analysis,
                "VVELMASS": vvelmass_analysis,
                "|UVELMASS,VVELMASS|": transport_mag_analysis,
                "PsiVEL": psivel_analysis,
                "PhiVEL": phivel_analysis,
            },
            SCALING_RATIO_THRESHOLD,
        )
    )
else:
    print("No diagnostic analysis fields found.")

SINGLE_HEATMAP_SCALE_BY_FIELD = analysis_scale_by_field
SCALING_ANALYSIS_DONE = True

print("\nSelected scale modes:")
for name, mode in SINGLE_HEATMAP_SCALE_BY_FIELD.items():
    print(f"  {name}: {mode}")

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

#%%
#! Cell 10: Snapshot schedule for the current run

if MAKE_SNAPSHOTS:
    timed = ac.discover_timed_variables(RUN_DIR)
    prog_iters_for_snapshots = ac.common_iterations(timed, ["Eta", "U", "V", "W"])
    diag_iters_for_snapshots, _ = ac.common_iterations_for_records(
        RUN_DIR,
        ["UVELMASS", "VVELMASS", "PsiVEL", "PhiVEL"],
        extra_vars=("Eta",),
    )

    prognostic_snapshots = ac.pick_snapshot_iterations(prog_iters_for_snapshots, DELTA_T_SEC, SNAPSHOT_DAYS)
    diagnostic_snapshots = ac.pick_snapshot_iterations(diag_iters_for_snapshots, DELTA_T_SEC, SNAPSHOT_DAYS)

    print("Prognostic snapshot targets:")
    for snap in prognostic_snapshots:
        note = "" if snap["exact"] else " (nearest saved frame)"
        print(
            f"  {snap['requested_label']}: iteration {snap['iteration']} | "
            f"elapsed = {snap['actual_label']}{note}"
        )

    print("\nDiagnostic snapshot targets:")
    for snap in diagnostic_snapshots:
        note = "" if snap["exact"] else " (nearest saved frame)"
        print(
            f"  {snap['requested_label']}: iteration {snap['iteration']} | "
            f"elapsed = {snap['actual_label']}{note}"
        )

#%%
#! Cell 11: Save prognostic snapshot panels

if MAKE_SNAPSHOTS:
    timed = ac.discover_timed_variables(RUN_DIR)
    prog_iters = ac.common_iterations(timed, ["Eta", "U", "V", "W"])
    if prog_iters:
        eta0, _ = ac.read_selected_field(RUN_DIR, "Eta", prog_iters[0])
        prog_land_mask = ac.load_land_mask(RUN_DIR, *eta0.shape)
        eta_prog = ac.read_variable_series(RUN_DIR, "Eta", prog_iters, prog_land_mask)
        u_prog = ac.read_variable_series(RUN_DIR, "U", prog_iters, prog_land_mask)
        v_prog = ac.read_variable_series(RUN_DIR, "V", prog_iters, prog_land_mask)
        w_prog = ac.read_variable_series(RUN_DIR, "W", prog_iters, prog_land_mask)
        mag_prog = [np.sqrt(u * u + v * v) for u, v in zip(u_prog, v_prog)]
        prognostic_fields = {"Eta": eta_prog, "U": u_prog, "V": v_prog, "W": w_prog, "|U,V|": mag_prog}
        prognostic_snapshots = ac.pick_snapshot_iterations(prog_iters, DELTA_T_SEC, SNAPSHOT_DAYS)
        prognostic_snapshot_dir = ANIMATION_DIR / "snapshots" / "prognostic"
        ac.write_snapshot_panels(
            "PROGNOSTIC VARIABLES",
            prognostic_fields,
            prog_iters,
            prognostic_snapshots,
            prognostic_snapshot_dir,
            CMAP,
            EXPERIMENT_TITLE,
            NTIMESTEPS,
            DELTA_T_SEC,
            dpi=SNAPSHOT_DPI,
        )
    else:
        print("Skipped prognostic snapshots: common Eta/U/V/W iterations not found.")

#%%
#! Cell 12: Save diagnostic snapshot panels

if MAKE_SNAPSHOTS:
    diag_records = ["UVELMASS", "VVELMASS", "PsiVEL", "PhiVEL"]
    diag_iters, diag_sources = ac.common_iterations_for_records(RUN_DIR, diag_records, extra_vars=("Eta",))
    print("Diagnostic record sources =", diag_sources)
    if diag_iters:
        eta0, _ = ac.read_selected_field(RUN_DIR, "Eta", diag_iters[0])
        diag_land_mask = ac.load_land_mask(RUN_DIR, *eta0.shape)
        eta_diag = ac.read_variable_series(RUN_DIR, "Eta", diag_iters, diag_land_mask)
        diagnostic_fields = {"Eta": eta_diag}

        try:
            u_diag, _ = ac.read_record_series_any_bundle(RUN_DIR, "UVELMASS", diag_iters, diag_land_mask)
            v_diag, _ = ac.read_record_series_any_bundle(RUN_DIR, "VVELMASS", diag_iters, diag_land_mask)
            diagnostic_fields["|UVELMASS,VVELMASS|"] = [np.sqrt(u * u + v * v) for u, v in zip(u_diag, v_diag)]
        except Exception as exc:
            print("Skipped |UVELMASS,VVELMASS| snapshots:", exc)

        for rec_name in ["PsiVEL", "PhiVEL"]:
            try:
                diagnostic_fields[rec_name], bundle = ac.read_record_series_any_bundle(
                    RUN_DIR, rec_name, diag_iters, diag_land_mask
                )
                print(f"{rec_name} source = {bundle}")
            except Exception as exc:
                print(f"Skipped {rec_name} snapshots:", exc)

        diagnostic_snapshots = ac.pick_snapshot_iterations(diag_iters, DELTA_T_SEC, SNAPSHOT_DAYS)
        diagnostic_snapshot_dir = ANIMATION_DIR / "snapshots" / "diagnostic"
        ac.write_snapshot_panels(
            "DIAGNOSTIC VARIABLES",
            diagnostic_fields,
            diag_iters,
            diagnostic_snapshots,
            diagnostic_snapshot_dir,
            CMAP,
            EXPERIMENT_TITLE,
            NTIMESTEPS,
            DELTA_T_SEC,
            dpi=SNAPSHOT_DPI,
        )
    else:
        print("Skipped diagnostic snapshots: no common diagnostic iterations found.")

#%%
#! Cell 13: Clean single heatmap snapshots for requested days

if MAKE_SNAPSHOTS:
    if not SCALING_ANALYSIS_DONE or SINGLE_HEATMAP_SCALE_BY_FIELD is None:
        raise RuntimeError("!! Run Cell 2.5 first. It decides whether each field should use adaptive or constant scaling.")

    prog_iters = ac.common_iterations(ac.discover_timed_variables(RUN_DIR), ["Eta", "U", "V", "W"])
    if prog_iters:
        prog_snapshots = ac.pick_snapshot_iterations(prog_iters, DELTA_T_SEC, SNAPSHOT_DAYS)
        prog_iter_list = [snap["iteration"] for snap in prog_snapshots]

        eta0, _ = ac.read_selected_field(RUN_DIR, "Eta", prog_iter_list[0])
        prog_land_mask = ac.load_land_mask(RUN_DIR, *eta0.shape)

        eta_single = ac.read_variable_series(RUN_DIR, "Eta", prog_iter_list, prog_land_mask)
        u_single = ac.read_variable_series(RUN_DIR, "U", prog_iter_list, prog_land_mask)
        v_single = ac.read_variable_series(RUN_DIR, "V", prog_iter_list, prog_land_mask)
        mag_single = [np.sqrt(u * u + v * v) for u, v in zip(u_single, v_single)]

        prog_single_fields = {
            "Eta": eta_single,
            "U": u_single,
            "V": v_single,
            "|U,V|": mag_single,
        }
        single_snapshot_dir = ANIMATION_DIR / "snapshots" / "single_heatmaps"
        ac.write_single_heatmap_snapshots(
            prog_single_fields,
            prog_iter_list,
            prog_snapshots,
            single_snapshot_dir,
            CMAP,
            dpi=SINGLE_FIG_DPI,
            scale_mode_by_field=SINGLE_HEATMAP_SCALE_BY_FIELD,
        )
    else:
        print("Skipped prognostic single snapshots: common Eta/U/V/W iterations not found.")

    diag_records = ["UVELMASS", "VVELMASS", "PsiVEL", "PhiVEL"]
    diag_iters, diag_sources = ac.common_iterations_for_records(RUN_DIR, diag_records)
    print("Diagnostic single snapshot sources =", diag_sources)
    if diag_iters:
        diag_snapshots = ac.pick_snapshot_iterations(diag_iters, DELTA_T_SEC, SNAPSHOT_DAYS)
        diag_iter_list = [snap["iteration"] for snap in diag_snapshots]

        diag0_single, _ = ac.read_record_any_bundle(RUN_DIR, "UVELMASS", diag_iter_list[0])
        diag_land_mask = ac.load_land_mask(RUN_DIR, *diag0_single.shape)
        uvelmass_single, uvelmass_bundle = ac.read_record_series_any_bundle(RUN_DIR, "UVELMASS", diag_iter_list, diag_land_mask)
        vvelmass_single, vvelmass_bundle = ac.read_record_series_any_bundle(RUN_DIR, "VVELMASS", diag_iter_list, diag_land_mask)
        psivel_single, psi_bundle = ac.read_record_series_any_bundle(RUN_DIR, "PsiVEL", diag_iter_list, diag_land_mask)
        phivel_single, phi_bundle = ac.read_record_series_any_bundle(RUN_DIR, "PhiVEL", diag_iter_list, diag_land_mask)
        print(f"UVELMASS source = {uvelmass_bundle}")
        print(f"VVELMASS source = {vvelmass_bundle}")
        print(f"PsiVEL source = {psi_bundle}")
        print(f"PhiVEL source = {phi_bundle}")

        diag_single_fields = {
            "UVELMASS": uvelmass_single,
            "VVELMASS": vvelmass_single,
            "|UVELMASS,VVELMASS|": [np.sqrt(u * u + v * v) for u, v in zip(uvelmass_single, vvelmass_single)],
            "PsiVEL": psivel_single,
            "PhiVEL": phivel_single,
        }
        single_snapshot_dir = ANIMATION_DIR / "snapshots" / "single_heatmaps"
        ac.write_single_heatmap_snapshots(
            diag_single_fields,
            diag_iter_list,
            diag_snapshots,
            single_snapshot_dir,
            CMAP,
            dpi=SINGLE_FIG_DPI,
            scale_mode_by_field=SINGLE_HEATMAP_SCALE_BY_FIELD,
        )
    else:
        print("Skipped diagnostic single snapshots: common UVELMASS/VVELMASS/PsiVEL/PhiVEL iterations not found.")

# %%
