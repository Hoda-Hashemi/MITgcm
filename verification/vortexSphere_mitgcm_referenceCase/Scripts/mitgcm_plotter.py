#%%
# Cell 1/5: Imports and settings.
import numpy as np
import scripts as sp

TARGET_IT = 0
CASE_LABEL = "Vortex Sphere"
CASE_OUTPUT_PREFIX = "vortex_sphere"
CASE_TARGET_IT = TARGET_IT
CASE_DYN_TARGET_IT = None
SAVE_PUBLICATION_PLOTS = False #True
CASE_DIAGNOSTIC_VARS = {"dyn_Aux", "dynDiag"}
CASE_RESULTS_DIR = sp._results_target_dir(CASE_TARGET_IT)
PLOT_ALL_REGULAR_ITERATIONS = True
PLOT_ALL_DIAGNOSTIC_ITERATIONS = True

def _iterations_to_plot(it_list, requested_it, plot_all):
    if plot_all:
        return list(it_list)
    return [sp._select_reference_iteration(it_list, requested_it, "iteration", quiet=True)]

#%%
# Cell 2/5: vortexSphere_mitgcm_referenceCase simple run/ plotting.
with sp._capture_prints(
    CASE_TARGET_IT,
    heading=f"{CASE_LABEL} setup for {sp._iteration_folder_name(CASE_TARGET_IT)}",
):
    print(f"\n--- vortexSphere_mitgcm_referenceCase run plotting from: {sp.RUN_DIR}")
    print(
        f"{CASE_LABEL} plotting backend: matplotlib imshow "
        f"(MITgcmutils available={sp.HAVE_MITGCMUTILS}, used_for_vortex_case=False)"
    )

    case_timed_vars = sp._discover_timed_variables(sp.RUN_DIR)
    case_regular_timed_vars = {
        name: values for name, values in case_timed_vars.items() if name not in CASE_DIAGNOSTIC_VARS
    }
    case_diag_timed_vars = {
        name: values for name, values in case_timed_vars.items() if name in CASE_DIAGNOSTIC_VARS
    }

    print(f"Discovered {CASE_LABEL} regular time-dependent output variables:")
    for var_name, it_list in case_regular_timed_vars.items():
        print(f"  {var_name}: iterations {it_list}")
    if case_diag_timed_vars:
        print("Diagnostic bundles reserved for Cell 3/5:")
        for var_name, it_list in case_diag_timed_vars.items():
            print(f"  {var_name}: iterations {it_list}")

    depth_2d = sp.read_mitgcm_static_field(sp.RUN_DIR, "Depth")
    if np.ndim(depth_2d) == 3:
        depth_2d = sp._surface_layer(depth_2d)

    depth_m, ocean_mask = sp._load_depth_and_ocean_mask(*depth_2d.shape)
    land_mask = ~ocean_mask
    depth_plot = np.where(ocean_mask, depth_m, np.nan)
    lat_deg = sp._lat_centers_deg(depth_2d.shape[0])
    hfac_c_2d = sp.read_mitgcm_static_field(sp.RUN_DIR, "hFacC")
    if np.ndim(hfac_c_2d) == 3:
        hfac_c_2d = sp._surface_layer(hfac_c_2d)
    depth_consistency = sp._report_depth_consistency(depth_2d, depth_m, ocean_mask, hfac_c_2d=hfac_c_2d)

with sp._capture_prints(
    CASE_TARGET_IT,
    heading=f"{CASE_LABEL} depth plot for {sp._iteration_folder_name(CASE_TARGET_IT)}",
):
    sp._plot_global_latlon(
        depth_plot,
        land_mask,
        f"{CASE_LABEL} Depth (land masked)",
        "Depth [m]",
        cmap="seismic",
        positive_only=True,
        save_stem=f"{CASE_OUTPUT_PREFIX}_depth_masked",
        save_enabled=SAVE_PUBLICATION_PLOTS,
        save_dir=CASE_RESULTS_DIR,
    )
    sp._land_leak_check(depth_plot, land_mask, "Depth", allow_land_values=True)
    sp._print_ocean_min_max(depth_plot, ocean_mask, "Depth", lat_deg=lat_deg, depth_m=depth_m)

for var_name, it_list in case_regular_timed_vars.items():
    if (
        not PLOT_ALL_REGULAR_ITERATIONS
        and CASE_TARGET_IT is not None
        and int(CASE_TARGET_IT) not in it_list
    ):
        fallback_it = sp._select_reference_iteration(it_list, CASE_TARGET_IT, var_name, quiet=True)
        print(f"[{var_name}] requested it={int(CASE_TARGET_IT)} not available; using latest it={fallback_it}.")

    for it_plot in _iterations_to_plot(it_list, CASE_TARGET_IT, PLOT_ALL_REGULAR_ITERATIONS):
        with sp._capture_prints(
            it_plot,
            heading=f"{CASE_LABEL} regular field {var_name} for {sp._iteration_folder_name(it_plot)}",
        ):
            for rec_idx, rec_name in sp._iter_record_specs(sp.RUN_DIR, var_name, it_plot):
                fld_raw, it_used = sp.read_mitgcm_field(sp.RUN_DIR, var_name, it=it_plot, record=rec_idx)
                fld_2d = sp._surface_layer(fld_raw)
                cmap, symmetric = sp._auto_colormap_and_symmetry(fld_2d, ocean_mask)
                units = sp._unit_for_field(rec_name)
                cbar_label = f"{rec_name} [{units}]" if units else rec_name
                sp._land_leak_check(fld_2d, land_mask, f"{rec_name} (it={it_used})")
                sp._print_ocean_min_max(
                    fld_2d,
                    ocean_mask,
                    f"{rec_name} (it={it_used})",
                    lat_deg=lat_deg,
                    depth_m=depth_m,
                )
                sp._plot_global_latlon(
                    fld_2d,
                    land_mask,
                    f"{CASE_LABEL} {rec_name} (it={it_used})",
                    cbar_label,
                    cmap=cmap,
                    symmetric=symmetric,
                    save_stem=f"{CASE_OUTPUT_PREFIX}_{var_name}_{rec_name}_it{it_used:010d}",
                    save_enabled=SAVE_PUBLICATION_PLOTS,
                    save_dir=sp._results_target_dir(it_used),
                )

if "U" in case_regular_timed_vars and "V" in case_regular_timed_vars:
    common_its = sorted(set(case_regular_timed_vars["U"]).intersection(case_regular_timed_vars["V"]))
    speed_notice = None

    if common_its:
        if PLOT_ALL_REGULAR_ITERATIONS:
            speed_iterations = common_its
        elif CASE_TARGET_IT is not None and int(CASE_TARGET_IT) in common_its:
            speed_iterations = [int(CASE_TARGET_IT)]
        else:
            speed_iterations = [common_its[-1]]
            speed_notice = (
                f"[speed] requested it={CASE_TARGET_IT} not available for both U and V; "
                f"using latest common it={speed_iterations[0]}."
            )
    else:
        speed_iterations = []
        speed_notice = "[speed] no common U/V iterations found; skipping speed plot."

    if speed_iterations:
        for speed_it in speed_iterations:
            with sp._capture_prints(
                speed_it,
                heading=f"{CASE_LABEL} speed diagnostic for {sp._iteration_folder_name(speed_it)}",
            ):
                if speed_notice is not None:
                    print(speed_notice)
                    speed_notice = None
                u2d = sp._surface_layer(sp.read_mitgcm_field(sp.RUN_DIR, "U", it=speed_it)[0])
                v2d = sp._surface_layer(sp.read_mitgcm_field(sp.RUN_DIR, "V", it=speed_it)[0])
                speed2d = np.sqrt(u2d * u2d + v2d * v2d)
                sp._land_leak_check(speed2d, land_mask, f"speed (it={speed_it})")
                sp._print_ocean_min_max(
                    speed2d,
                    ocean_mask,
                    f"speed (it={speed_it})",
                    lat_deg=lat_deg,
                    depth_m=depth_m,
                )
                sp._plot_global_latlon(
                    speed2d,
                    land_mask,
                    f"{CASE_LABEL} Speed sqrt(U^2+V^2) (it={speed_it})",
                    "speed [m/s]",
                    cmap="seismic",
                    positive_only=True,
                    save_stem=f"{CASE_OUTPUT_PREFIX}_speed_it{speed_it:010d}",
                    save_enabled=SAVE_PUBLICATION_PLOTS,
                    save_dir=sp._results_target_dir(speed_it),
                )
    else:
        with sp._capture_prints(
            CASE_TARGET_IT,
            heading=f"{CASE_LABEL} speed diagnostic for {sp._iteration_folder_name(CASE_TARGET_IT)}",
        ):
            print(speed_notice)

#%%
# Cell 3/5: vortexSphere_mitgcm_referenceCase bundled diagnostics.
if case_diag_timed_vars:
    for var_name, it_list in case_diag_timed_vars.items():
        if (
            not PLOT_ALL_DIAGNOSTIC_ITERATIONS
            and CASE_DYN_TARGET_IT is not None
            and int(CASE_DYN_TARGET_IT) not in it_list
        ):
            fallback_it = sp._select_reference_iteration(it_list, CASE_DYN_TARGET_IT, var_name, quiet=True)
            print(
                f"[{var_name}] requested it={int(CASE_DYN_TARGET_IT)} not available; "
                f"using latest it={fallback_it}."
            )

        for it_plot in _iterations_to_plot(it_list, CASE_DYN_TARGET_IT, PLOT_ALL_DIAGNOSTIC_ITERATIONS):
            with sp._capture_prints(
                it_plot,
                dyn=True,
                heading=f"{CASE_LABEL} diagnostic bundle {var_name} for {sp._iteration_folder_name(it_plot)}",
            ):
                print(f"\n--- {CASE_LABEL} diagnostic bundle plotting ---")
                for rec_idx, rec_name in sp._iter_record_specs(sp.RUN_DIR, var_name, it_plot):
                    fld_raw, it_used = sp.read_mitgcm_field(sp.RUN_DIR, var_name, it=it_plot, record=rec_idx)
                    fld_2d = sp._surface_layer(fld_raw)
                    cmap, symmetric = sp._auto_colormap_and_symmetry(fld_2d, ocean_mask)
                    units = sp._unit_for_field(rec_name)
                    cbar_label = f"{rec_name} [{units}]" if units else rec_name
                    sp._land_leak_check(fld_2d, land_mask, f"{rec_name} (it={it_used})")
                    sp._print_ocean_min_max(
                        fld_2d,
                        ocean_mask,
                        f"{rec_name} (it={it_used})",
                        lat_deg=lat_deg,
                        depth_m=depth_m,
                    )
                    sp._plot_global_latlon(
                        fld_2d,
                        land_mask,
                        f"{CASE_LABEL} {rec_name} (it={it_used})",
                        cbar_label,
                        cmap=cmap,
                        symmetric=symmetric,
                        save_stem=f"{CASE_OUTPUT_PREFIX}_{var_name}_{rec_name}_it{it_used:010d}",
                        save_enabled=SAVE_PUBLICATION_PLOTS,
                        save_dir=sp._results_dyn_dir(it_used),
                    )
else:
    with sp._capture_prints(
        CASE_DYN_TARGET_IT,
        dyn=True,
        heading=f"{CASE_LABEL} diagnostics for {sp._iteration_folder_name(CASE_DYN_TARGET_IT)}",
    ):
        print("No dyn_* diagnostics found for Cell 3/5.")

with sp._capture_prints(
    CASE_TARGET_IT,
    heading=f"{CASE_LABEL} checks for {sp._iteration_folder_name(CASE_TARGET_IT)}",
):
    sp._print_eta_init_ocean_land_check()
    sp._suggest_additional_reference_plots(sp.RUN_DIR)

#%%
# Cell 4/5: Vortex Sphere bathymetry quicklook.
ny, nx = 720, 1440
bathy = np.fromfile(sp.INPUT_BATHY_PATH, dtype=">f4").reshape(ny, nx)
model_depth_2d = sp.read_mitgcm_static_field(sp.RUN_DIR, "Depth")
if np.ndim(model_depth_2d) == 3:
    model_depth_2d = sp._surface_layer(model_depth_2d)
hfac_c_bathy = sp.read_mitgcm_static_field(sp.RUN_DIR, "hFacC")
if np.ndim(hfac_c_bathy) == 3:
    hfac_c_bathy = sp._surface_layer(hfac_c_bathy)

input_depth_m, input_ocean_mask = sp._load_depth_and_ocean_mask(*model_depth_2d.shape)
depth_consistency_bathy = sp._report_depth_consistency(
    model_depth_2d,
    input_depth_m,
    input_ocean_mask,
    hfac_c_2d=hfac_c_bathy,
    verbose=False,
)

raw_depth = np.where(input_ocean_mask, input_depth_m, np.nan)
model_depth = np.where(model_depth_2d > 0.0, model_depth_2d, np.nan)
depth_diff = depth_consistency_bathy["depth_diff_map"]
hfac_plot = np.where(model_depth_2d > 0.0, hfac_c_bathy, np.nan)

with sp._capture_prints(
    CASE_TARGET_IT,
    heading=f"{CASE_LABEL} bathymetry quicklook for {sp._iteration_folder_name(CASE_TARGET_IT)}",
):
    sp._plot_global_latlon(
        raw_depth,
        ~input_ocean_mask,
        "Raw Bathymetry from bathymetry.bin",
        "Depth [m]",
        cmap="seismic",
        positive_only=True,
        save_stem=f"{CASE_OUTPUT_PREFIX}_raw_bathymetry",
        save_enabled=SAVE_PUBLICATION_PLOTS,
        save_dir=CASE_RESULTS_DIR,
    )
    sp._plot_global_latlon(
        model_depth,
        model_depth_2d <= 0.0,
        "Model-Resolved Depth from run/Depth",
        "Depth [m]",
        cmap="seismic",
        positive_only=True,
        save_stem=f"{CASE_OUTPUT_PREFIX}_model_depth",
        save_enabled=SAVE_PUBLICATION_PLOTS,
        save_dir=CASE_RESULTS_DIR,
    )
    sp._plot_global_latlon(
        depth_diff,
        ~depth_consistency_bathy["overlap_mask"],
        "Model Depth - Raw Bathymetry",
        "Depth difference [m]",
        cmap="seismic",
        symmetric=True,
        save_stem=f"{CASE_OUTPUT_PREFIX}_model_minus_raw_bathymetry",
        save_enabled=SAVE_PUBLICATION_PLOTS,
        save_dir=CASE_RESULTS_DIR,
    )
    sp._plot_global_latlon(
        hfac_plot,
        model_depth_2d <= 0.0,
        "One-Layer hFacC Carrying Bathymetry",
        "hFacC [-]",
        cmap="seismic",
        positive_only=True,
        save_stem=f"{CASE_OUTPUT_PREFIX}_hfacc",
        save_enabled=SAVE_PUBLICATION_PLOTS,
        save_dir=CASE_RESULTS_DIR,
    )
    sp._report_depth_consistency(
        model_depth_2d,
        input_depth_m,
        input_ocean_mask,
        hfac_c_2d=hfac_c_bathy,
        verbose=True,
    )

# %%
# !!!! SAVING SOME VARIABLES SO THEY WILL BE USED INSIDE FORTAN CODE
#%%
# dyn Cell: load + merge + save everything
#%%
import os
import numpy as np
import MITgcmutils as mit

OUT_DIR = os.path.join(sp.RUN_DIR, "csv_dyn")
os.makedirs(OUT_DIR, exist_ok=True)

it = 22

DXG = mit.rdmds(os.path.join(sp.RUN_DIR, "DXG"))
DYG = mit.rdmds(os.path.join(sp.RUN_DIR, "DYG"))
DRF = mit.rdmds(os.path.join(sp.RUN_DIR, "DRF"))[0]

fields = {}

for var_name in ["dynDiag", "dyn_Aux"]:

    for rec_idx, rec_name in sp._iter_record_specs(sp.RUN_DIR, var_name, it):

        fld_raw, _ = sp.read_mitgcm_field(
            sp.RUN_DIR, var_name, it=it, record=rec_idx
        )

        fields[rec_name] = sp._surface_layer(fld_raw)

if "UVELMASS" in fields:
    fields["uTrans"] = DYG * DRF * fields["UVELMASS"]

if "VVELMASS" in fields:
    fields["vTrans"] = DXG * DRF * fields["VVELMASS"]

for name, fld in fields.items():
    np.savetxt(
        os.path.join(OUT_DIR, f"{name}_it{it:010d}.csv"),
        fld,
        delimiter=","
    )

print("saved:", list(fields.keys()))
# %%
# --- Load mask fields
Depth = mit.rdmds(os.path.join(sp.RUN_DIR, "Depth"))
hFacC = mit.rdmds(os.path.join(sp.RUN_DIR, "hFacC"))

Depth = sp._surface_layer(Depth)
hFacC = sp._surface_layer(hFacC)

# --- Save Depth
np.savetxt(
    os.path.join(OUT_DIR, f"Depth_it{it:010d}.csv"),
    Depth,
    delimiter=","
)

# --- Save hFacC
np.savetxt(
    os.path.join(OUT_DIR, f"hFacC_it{it:010d}.csv"),
    hFacC,
    delimiter=","
)
#%%
# Print min/max of all dyn fields (iteration 22)

import numpy as np

it = 22

for var_name in ["dynDiag", "dyn_Aux"]:
    print(f"\n--- {var_name} ---")
    
    for rec_idx, rec_name in sp._iter_record_specs(sp.RUN_DIR, var_name, it):
        
        fld_raw, _ = sp.read_mitgcm_field(
            sp.RUN_DIR, var_name, it=it, record=rec_idx
        )
        
        fld = sp._surface_layer(fld_raw)
        
        print(f"{rec_name:15s}  min = {np.min(fld):.6e}   max = {np.max(fld):.6e}")

# %%
#!HTML
import glob
import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

Nx, Ny = 1440, 720
dx, dy = 0.25, 0.25

run_dir = "../run_Rossby_6hours"
# run_dir = "../run_34_MPI"
out_dir = "../docs"
skip = 4

os.makedirs(out_dir, exist_ok=True)

lon = (np.arange(Nx) + 0.5)*dx
lat = 90.0 - (np.arange(Ny) + 0.5)*dy
lon_map = lon - 180.0

lon_p = lon_map[::skip]
lat_p = lat[::skip]

colorscale_main ="balance" # "RdBu_r"

def read2d(path):
    return np.fromfile(path, dtype=">f8").reshape(Ny, Nx)

def shift(field):
    return np.roll(field, Nx//2, axis=1)

def down(field):
    return field[::skip, ::skip]

def sym_scale(fields):
    a = max(np.nanmax(np.abs(f)) for f in fields)
    return -a, a

def read_field_series(prefix):
    files = sorted(glob.glob(os.path.join(run_dir, f"{prefix}.*.data")))
    iters = [int(os.path.basename(f).split(".")[1]) for f in files]
    fields = [down(shift(read2d(f))) for f in files]
    print(prefix, "iterations:", iters)
    return iters, fields

Depth = read2d(os.path.join(run_dir, "Depth.data"))
Depth = np.where(Depth > 0, Depth, np.nan)
Depth_p = down(shift(Depth))
def write_map_html(iters, fields, name, title, outfile, symmetric=True, show_depth=True):
    if len(fields) == 0:
        print("No fields found for", name)
        return

    # Global color scale only
    global_min = min(np.nanmin(f) for f in fields)
    global_max = max(np.nanmax(f) for f in fields)

    if symmetric:
        zmax_abs = max(abs(global_min), abs(global_max))
        zmin, zmax = -zmax_abs, zmax_abs
    else:
        zmin, zmax = global_min, global_max

    fig = go.Figure()

    # Per-iteration min/max text
    def info_text(i):
        return (
            f"{name}<br>"
            f"iter = {iters[i]}<br>"
            f"min = {np.nanmin(fields[i]):.4e}<br>"
            f"max = {np.nanmax(fields[i]):.4e}<br>"
            f"skip = {skip}"
        )

    if show_depth:
        fig.add_trace(
            go.Heatmap(
                x=lon_p,
                y=lat_p,
                z=Depth_p,
                colorscale=colorscale_main,
                opacity=0.25,
                showscale=False,
                hovertemplate="lon=%{x}<br>lat=%{y}<br>Depth=%{z:.2f} m<extra></extra>"
            )
        )

    fig.add_trace(
        go.Heatmap(
            x=lon_p,
            y=lat_p,
            z=fields[0],
            colorscale=colorscale_main,
            zmin=zmin,
            zmax=zmax,
            opacity=0.95,
            colorbar=dict(title=name),
            hovertemplate=f"lon=%{{x}}<br>lat=%{{y}}<br>{name}=%{{z:.4e}}<extra></extra>"
        )
    )

    # Put box on right side, not title
    fig.add_annotation(
        text=info_text(0),
        x=1.22,
        y=0.92,
        xref="paper",
        yref="paper",
        showarrow=False,
        align="left",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1,
        font=dict(size=13)
    )

    frames = []

    for i, it in enumerate(iters):
        traces = []

        if show_depth:
            traces.append(
                go.Heatmap(
                    x=lon_p,
                    y=lat_p,
                    z=Depth_p,
                    colorscale=colorscale_main,
                    opacity=0.25,
                    showscale=False
                )
            )

        traces.append(
            go.Heatmap(
                x=lon_p,
                y=lat_p,
                z=fields[i],
                colorscale=colorscale_main,
                zmin=zmin,
                zmax=zmax,
                opacity=0.95,
                colorbar=dict(title=name)
            )
        )

        frames.append(
            go.Frame(
                data=traces,
                name=str(it),
                layout=go.Layout(
                    annotations=[
                        dict(
                            text=info_text(i),
                            x=1.22,
                            y=0.92,
                            xref="paper",
                            yref="paper",
                            showarrow=False,
                            align="left",
                            bgcolor="white",
                            bordercolor="black",
                            borderwidth=1,
                            font=dict(size=13)
                        )
                    ]
                )
            )
        )

    fig.frames = frames

    steps = [
        dict(
            method="animate",
            args=[
                [str(it)],
                dict(
                    mode="immediate",
                    frame=dict(duration=400, redraw=True),
                    transition=dict(duration=0),
                    redraw=True
                )
            ],
            label=f"Iter {it}"
        )
        for it in iters
    ]

    fig.update_layout(
        title=dict(text=title, x=0.5),
        width=1700,
        height=900,
        margin=dict(l=70, r=420, t=90, b=120),
        xaxis_title="Longitude [deg]",
        yaxis_title="Latitude [deg]",
        xaxis=dict(range=[-180, 180]),
        yaxis=dict(range=[90, -90]),
        sliders=[dict(active=0, steps=steps, x=0.08, y=0.02, len=0.85)]
    )

    fig.write_html(os.path.join(out_dir, outfile))
    fig.show()

#%%
# Read fields and write HTML maps

eta_iters, eta_fields = read_field_series("Eta")

write_map_html(
    eta_iters,
    eta_fields,
    name="h_s [m]",
    title="MITgcm Free-Surface Anomaly h_s over Depth",
    outfile="mitgcm_eta.html",
    symmetric=True
)

ph_iters, ph_fields = read_field_series("PH")

write_map_html(
    ph_iters,
    ph_fields,
    name="PH",
    title="MITgcm Geopotential PH over Depth",
    outfile="mitgcm_ph.html",
    symmetric=True
)

u_iters, u_fields = read_field_series("U")
v_iters, v_fields = read_field_series("V")

common_iters = sorted(set(u_iters).intersection(v_iters))

mag_fields = []

for it in common_iters:
    tag = f"{it:010d}"

    U = read2d(os.path.join(run_dir, f"U.{tag}.data"))
    V = read2d(os.path.join(run_dir, f"V.{tag}.data"))

    W_path = os.path.join(run_dir, f"W.{tag}.data")
    W = read2d(W_path) if os.path.exists(W_path) else np.zeros_like(U)

    Mag = np.sqrt(U**2 + V**2 + W**2)
    mag_fields.append(down(shift(Mag)))

write_map_html(
    common_iters,
    mag_fields,
    name="|u| [m/s]",
    title="MITgcm Velocity Magnitude over Depth",
    outfile="mitgcm_velocity_magnitude.html",
    symmetric=False
)
# %%

