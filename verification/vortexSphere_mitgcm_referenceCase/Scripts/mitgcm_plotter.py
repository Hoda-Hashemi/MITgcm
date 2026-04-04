#%%
# Cell 1/5: Imports and settings.
import numpy as np
import scripts as sp

TARGET_IT = 0
CASE_LABEL = "Vortex Sphere"
CASE_OUTPUT_PREFIX = "vortex_sphere"
CASE_TARGET_IT = TARGET_IT
CASE_DYN_TARGET_IT = None
SAVE_PUBLICATION_PLOTS = True
CASE_DIAGNOSTIC_VARS = {"dyn_Aux", "dynDiag"}
CASE_RESULTS_DIR = sp._results_target_dir(CASE_TARGET_IT)

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
    it_plot = sp._select_reference_iteration(it_list, CASE_TARGET_IT, var_name, quiet=True)
    with sp._capture_prints(
        it_plot,
        heading=f"{CASE_LABEL} regular field {var_name} for {sp._iteration_folder_name(it_plot)}",
    ):
        if CASE_TARGET_IT is not None and int(CASE_TARGET_IT) not in it_list:
            print(f"[{var_name}] requested it={int(CASE_TARGET_IT)} not available; using latest it={it_plot}.")
        for rec_idx, rec_name in sp._iter_record_specs(sp.RUN_DIR, var_name, it_plot):
            fld_raw, it_used = sp.read_mitgcm_field(sp.RUN_DIR, var_name, it=it_plot, record=rec_idx)
            fld_2d = sp._surface_layer(fld_raw)
            cmap, symmetric = sp._auto_colormap_and_symmetry(fld_2d, ocean_mask)
            units = sp._unit_for_field(rec_name)
            cbar_label = f"{rec_name} [{units}]" if units else rec_name
            sp._land_leak_check(fld_2d, land_mask, f"{rec_name} (it={it_used})")
            sp._print_ocean_min_max(fld_2d, ocean_mask, f"{rec_name} (it={it_used})", lat_deg=lat_deg, depth_m=depth_m)
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

    if CASE_TARGET_IT is not None and int(CASE_TARGET_IT) in common_its:
        speed_it = int(CASE_TARGET_IT)
    elif common_its:
        speed_it = common_its[-1]
        speed_notice = (
            f"[speed] requested it={CASE_TARGET_IT} not available for both U and V; "
            f"using latest common it={speed_it}."
        )
    else:
        speed_it = None
        speed_notice = "[speed] no common U/V iterations found; skipping speed plot."

    if speed_it is not None:
        with sp._capture_prints(
            speed_it,
            heading=f"{CASE_LABEL} speed diagnostic for {sp._iteration_folder_name(speed_it)}",
        ):
            if speed_notice is not None:
                print(speed_notice)
            u2d = sp._surface_layer(sp.read_mitgcm_field(sp.RUN_DIR, "U", it=speed_it)[0])
            v2d = sp._surface_layer(sp.read_mitgcm_field(sp.RUN_DIR, "V", it=speed_it)[0])
            speed2d = np.sqrt(u2d * u2d + v2d * v2d)
            sp._land_leak_check(speed2d, land_mask, f"speed (it={speed_it})")
            sp._print_ocean_min_max(speed2d, ocean_mask, f"speed (it={speed_it})", lat_deg=lat_deg, depth_m=depth_m)
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
        it_plot = sp._select_reference_iteration(it_list, CASE_DYN_TARGET_IT, var_name, quiet=True)
        with sp._capture_prints(
            it_plot,
            dyn=True,
            heading=f"{CASE_LABEL} diagnostic bundle {var_name} for {sp._iteration_folder_name(it_plot)}",
        ):
            print(f"\n--- {CASE_LABEL} diagnostic bundle plotting ---")
            if CASE_DYN_TARGET_IT is not None and int(CASE_DYN_TARGET_IT) not in it_list:
                print(
                    f"[{var_name}] requested it={int(CASE_DYN_TARGET_IT)} not available; "
                    f"using latest it={it_plot}."
                )
            for rec_idx, rec_name in sp._iter_record_specs(sp.RUN_DIR, var_name, it_plot):
                fld_raw, it_used = sp.read_mitgcm_field(sp.RUN_DIR, var_name, it=it_plot, record=rec_idx)
                fld_2d = sp._surface_layer(fld_raw)
                cmap, symmetric = sp._auto_colormap_and_symmetry(fld_2d, ocean_mask)
                units = sp._unit_for_field(rec_name)
                cbar_label = f"{rec_name} [{units}]" if units else rec_name
                sp._land_leak_check(fld_2d, land_mask, f"{rec_name} (it={it_used})")
                sp._print_ocean_min_max(fld_2d, ocean_mask, f"{rec_name} (it={it_used})", lat_deg=lat_deg, depth_m=depth_m)
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
