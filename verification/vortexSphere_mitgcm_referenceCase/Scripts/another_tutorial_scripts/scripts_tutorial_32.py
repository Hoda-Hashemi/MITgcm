import copy
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from MITgcmutils import cs as mitgcm_cs

import scripts as sp


def _split_cs_faces(field_2d, n_faces=6):
    ny, nx_total = field_2d.shape
    nx_face = nx_total // n_faces
    return [field_2d[:, index * nx_face : (index + 1) * nx_face] for index in range(n_faces)]


def _plot_cs_faces(field_2d, title, cmap="seismic", cbar_label="", symmetric=False, land_mask=None):
    fld = ma.array(field_2d, copy=False)
    if land_mask is not None:
        fld = ma.masked_where(land_mask, fld)
    faces = _split_cs_faces(fld, n_faces=6)
    fld_np = ma.filled(fld, np.nan)
    plot_kwargs = sp._plot_scale_kwargs(fld_np, cmap=cmap, symmetric=symmetric)
    cmap_obj = copy.copy(plt.get_cmap(cmap))
    cmap_obj.set_bad((1.0, 1.0, 1.0, 0.0))

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.2), constrained_layout=True)
    image = None
    for index, ax in enumerate(axes.flat):
        image = ax.imshow(faces[index], origin="lower", cmap=cmap_obj, aspect="equal", **plot_kwargs)
        ax.set_title(f"Face {index + 1}")
        ax.set_xlabel("i")
        ax.set_ylabel("j")

    fig.suptitle(title)
    sp._add_right_colorbar(fig, axes.ravel().tolist(), image, cbar_label, shrink=0.92)
    plt.show()


def _plot_cs_global_with_mitgcmutils(
    xc_2d,
    yc_2d,
    field_2d,
    title,
    cmap="seismic",
    cbar_label="",
    symmetric=False,
    land_mask=None,
):
    fld = ma.array(field_2d, copy=False)
    if land_mask is not None:
        fld = ma.masked_where(land_mask, fld)
    fld_np = ma.filled(fld, np.nan)
    plot_kwargs = sp._plot_scale_kwargs(fld_np, cmap=cmap, symmetric=symmetric)
    cmap_obj = copy.copy(plt.get_cmap(cmap))
    cmap_obj.set_bad("0.7")

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.sca(ax)
    image = mitgcm_cs.pcol(xc_2d, yc_2d, fld, cmap=cmap_obj, **plot_kwargs)
    ax.set_title(title)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    sp._add_right_colorbar(fig, ax, image[0], cbar_label)
    fig.tight_layout()
    plt.show()
    return True


def _print_wet_value_summary(field_2d, wet_mask, field_name):
    wet_vals = np.asarray(field_2d)[wet_mask]
    if wet_vals.size == 0:
        print(f"[{field_name}] no wet points available.")
        return
    print(
        f"[{field_name}] wet min={float(np.nanmin(wet_vals)):.6g}, "
        f"max={float(np.nanmax(wet_vals)):.6g}, "
        f"unique_wet_values={np.unique(wet_vals).size}"
    )


def _plot_compact_curve_diagnostics(field_2d, land_mask, title, y_label):
    fld = ma.masked_where(land_mask, ma.array(field_2d, copy=False))
    fld_np = ma.filled(fld, np.nan)
    ny, nx = fld_np.shape
    lat_idx = np.arange(ny)
    lon_idx = np.arange(nx)
    row_ids = sorted({ny // 4, ny // 2, (3 * ny) // 4})
    col_ids = sorted({nx // 8, nx // 2, (7 * nx) // 8})

    fig, axes = plt.subplots(1, 2, figsize=(14, 4), constrained_layout=True)

    for j0 in row_ids:
        axes[0].plot(lon_idx, fld_np[j0, :], lw=1.4, label=f"j={j0}")
    axes[0].plot(lon_idx, np.nanmean(fld_np, axis=0), "k--", lw=2.0, label="mean over j")
    axes[0].set_title(f"{title}: zonal transects")
    axes[0].set_xlabel("compact longitude index i")
    axes[0].set_ylabel(y_label)
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    for i0 in col_ids:
        axes[1].plot(lat_idx, fld_np[:, i0], lw=1.4, label=f"i={i0}")
    axes[1].plot(lat_idx, np.nanmean(fld_np, axis=1), "k--", lw=2.0, label="mean over i")
    axes[1].set_title(f"{title}: meridional transects")
    axes[1].set_xlabel("compact latitude index j")
    axes[1].set_ylabel(y_label)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    plt.show()


def _load_adjustment_cs_initial_eta(run_dir, ny, nx_total):
    ssh_path = os.path.join(run_dir, "ssh_eq.bin")
    if not os.path.exists(ssh_path):
        return None
    n_faces = 6
    nx_face = nx_total // n_faces
    raw = np.fromfile(ssh_path, dtype=">f8")
    arr = raw.reshape((nx_face, n_faces, ny), order="F")
    faces = [arr[:, face, :].T for face in range(n_faces)]
    return np.concatenate(faces, axis=1)


def _plot_cs_pcol(
    xg_2d,
    yg_2d,
    field_2d,
    land_mask,
    title,
    cmap,
    cbar_label,
    symmetric=False,
    projection=None,
    overlay_contours=False,
    contour_levels=12,
    contour_color="k",
    xc_2d=None,
    yc_2d=None,
):
    fld = ma.masked_where(land_mask, ma.array(field_2d, copy=False))
    fld_np = ma.filled(fld, np.nan)
    plot_kwargs = sp._plot_scale_kwargs(fld_np, cmap=cmap, symmetric=symmetric)
    cmap_obj = copy.copy(plt.get_cmap(cmap))
    cmap_obj.set_bad("0.7")

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    plt.sca(ax)
    hh = mitgcm_cs.pcol(xg_2d, yg_2d, fld, projection=projection, cmap=cmap_obj, **plot_kwargs)
    mappable = hh[0] if np.size(hh) > 0 else None

    if overlay_contours and projection is None and xc_2d is not None and yc_2d is not None:
        xc_flat = ma.filled(ma.masked_where(land_mask, xc_2d), np.nan).ravel()
        yc_flat = ma.filled(ma.masked_where(land_mask, yc_2d), np.nan).ravel()
        fld_flat = fld_np.ravel()
        good = np.isfinite(xc_flat) & np.isfinite(yc_flat) & np.isfinite(fld_flat)
        if np.count_nonzero(good) > 8:
            zc = fld_flat[good]
            if np.nanmax(zc) > np.nanmin(zc):
                ax.tricontour(
                    xc_flat[good],
                    yc_flat[good],
                    zc,
                    levels=contour_levels,
                    colors=contour_color,
                    linewidths=0.45,
                    alpha=0.6,
                )

    if mappable is not None:
        sp._add_right_colorbar(fig, ax, mappable, cbar_label)
    ax.set_title(title)
    if projection is None:
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
    fig.tight_layout()
    plt.show()


def run_adjustment_plots(experiment_dir, target_it=12, face_target_it=None):
    adjustment_run_dir = os.fspath(Path(experiment_dir).parent / "adjustment.cs-32x32x1" / "run")
    depth_adj, _ = sp.read_mitgcm_field_or_static(adjustment_run_dir, "Depth", it=0)
    hfac_c_adj, _ = sp.read_mitgcm_field_or_static(adjustment_run_dir, "hFacC", it=0)
    xg_adj, _ = sp.read_mitgcm_field_or_static(adjustment_run_dir, "XG")
    yg_adj, _ = sp.read_mitgcm_field_or_static(adjustment_run_dir, "YG")
    xc_adj, _ = sp.read_mitgcm_field_or_static(adjustment_run_dir, "XC")
    yc_adj, _ = sp.read_mitgcm_field_or_static(adjustment_run_dir, "YC")

    depth_adj_2d = sp._surface_layer(depth_adj)
    hfac_c_adj_2d = sp._surface_layer(hfac_c_adj)
    land_mask_adj = hfac_c_adj_2d <= 0.0
    wet_mask_adj = ~land_mask_adj

    print(f"\n--- Adjustment.cs-32x32x1 tutorial plots from: {adjustment_run_dir}")
    print("MITgcmutils options in this env: cs.pcol.")
    print(
        "Depth map note: projected global cs.pcol can look non-rectangular because cubed-sphere "
        "faces are warped to lon/lat."
    )
    _print_wet_value_summary(depth_adj_2d, wet_mask_adj, "Depth")
    if np.unique(depth_adj_2d[wet_mask_adj]).size == 1:
        print("Depth has one wet value in this experiment: ocean is flat-bottom, land is masked.")

    eta0_from_file = _load_adjustment_cs_initial_eta(adjustment_run_dir, depth_adj_2d.shape[0], depth_adj_2d.shape[1])
    if eta0_from_file is None:
        eta0_from_file = np.zeros_like(depth_adj_2d)

    _plot_cs_pcol(xg_adj, yg_adj, eta0_from_file, land_mask_adj, "CS-32x32x1 Initial Surface Anomaly (Global)", "seismic", "eta [m]", symmetric=True)
    _plot_cs_pcol(xg_adj, yg_adj, depth_adj_2d, land_mask_adj, "CS-32x32x1 Depth (Global)", "seismic", "Depth [m]")
    _plot_cs_pcol(
        xg_adj,
        yg_adj,
        depth_adj_2d,
        land_mask_adj,
        "CS-32x32x1 Depth (Global + Contour Overlay)",
        "seismic",
        "Depth [m]",
        overlay_contours=True,
        contour_levels=8,
        xc_2d=xc_adj,
        yc_2d=yc_adj,
    )

    var_to_its = sp._discover_timed_variables(adjustment_run_dir)
    print("\nDiscovered time-dependent output variables:")
    for var_name, it_list in var_to_its.items():
        print(f"  {var_name}: iterations {it_list}")

    for var_name, it_list in var_to_its.items():
        it_plot = sp._select_reference_iteration(it_list, target_it, var_name)
        rec_specs = sp._iter_record_specs(adjustment_run_dir, var_name, it_plot)

        if len(rec_specs) > 1:
            print(f"[{var_name}] bundled diagnostics with nrecords={len(rec_specs)}.")

        for rec_idx, rec_name in rec_specs:
            fld_raw, it_used = sp.read_mitgcm_field(adjustment_run_dir, var_name, it=it_plot, record=rec_idx)
            fld_2d = sp._surface_layer(fld_raw)
            if fld_2d.shape != land_mask_adj.shape:
                print(
                    f"Skipping {var_name} record={rec_idx}: shape {fld_2d.shape} "
                    f"does not match CS surface shape {land_mask_adj.shape}."
                )
                continue

            cmap, symmetric = sp._auto_colormap_and_symmetry(fld_2d, wet_mask_adj)
            units = sp._unit_for_field(rec_name)
            cbar_label = f"{rec_name} [{units}]" if units else rec_name

            _plot_cs_pcol(
                xg_adj,
                yg_adj,
                fld_2d,
                land_mask_adj,
                f"CS-32x32x1 {rec_name} (Global, it={it_used})",
                cmap,
                cbar_label,
                symmetric=symmetric,
            )

            if rec_name in ("Eta", "ETAN", "T", "S"):
                _plot_compact_curve_diagnostics(
                    fld_2d,
                    land_mask_adj,
                    f"CS-32x32x1 {rec_name} (it={it_used})",
                    cbar_label,
                )

    if "U" in var_to_its and "V" in var_to_its:
        speed_it = sp._select_reference_iteration(
            sorted(set(var_to_its["U"]).intersection(var_to_its["V"])),
            target_it,
            "Adjustment speed",
        )
        u2d = sp._surface_layer(sp.read_mitgcm_field(adjustment_run_dir, "U", it=speed_it)[0])
        v2d = sp._surface_layer(sp.read_mitgcm_field(adjustment_run_dir, "V", it=speed_it)[0])
        speed2d = np.sqrt(u2d * u2d + v2d * v2d)
        _plot_cs_pcol(
            xg_adj,
            yg_adj,
            speed2d,
            land_mask_adj,
            f"CS-32x32x1 Speed sqrt(U^2+V^2) (Global, it={speed_it})",
            "seismic",
            "speed [m/s]",
        )
        _plot_compact_curve_diagnostics(speed2d, land_mask_adj, "CS-32x32x1 Speed", "speed [m/s]")

    _plot_compact_curve_diagnostics(depth_adj_2d, land_mask_adj, "CS-32x32x1 Depth", "Depth [m]")

    print("Grid type: Cubed-Sphere (6 compact faces).")
    eta_adj_latest, it_adj = sp.read_mitgcm_field(adjustment_run_dir, "Eta", it=face_target_it)
    u_adj, _ = sp.read_mitgcm_field(adjustment_run_dir, "U", it=it_adj)
    v_adj, _ = sp.read_mitgcm_field(adjustment_run_dir, "V", it=it_adj)
    w_adj, _ = sp.read_mitgcm_field(adjustment_run_dir, "W", it=it_adj)
    phivel_adj, it_aux = sp.read_mitgcm_field(adjustment_run_dir, "dyn_Aux", it=None, record=1)
    psivel_adj, _ = sp.read_mitgcm_field(adjustment_run_dir, "dyn_Aux", it=it_aux, record=2)

    eta_adj_2d = sp._surface_layer(eta_adj_latest)
    u_adj_2d = sp._surface_layer(u_adj)
    v_adj_2d = sp._surface_layer(v_adj)
    w_adj_2d = sp._surface_layer(w_adj)
    phivel_adj_2d = sp._surface_layer(phivel_adj)
    psivel_adj_2d = sp._surface_layer(psivel_adj)
    speed_h_adj_2d = np.sqrt(u_adj_2d * u_adj_2d + v_adj_2d * v_adj_2d)

    print(
        "Global cs.pcol depth map is projected from cubed-sphere faces to lon/lat, "
        "so straight face edges can look curved."
    )
    _print_wet_value_summary(depth_adj_2d, wet_mask_adj, "Depth")
    if np.unique(depth_adj_2d[wet_mask_adj]).size == 1:
        print("Depth is intentionally uniform over ocean in this tutorial case; only land points differ.")

    eta0_from_file = _load_adjustment_cs_initial_eta(adjustment_run_dir, eta_adj_2d.shape[0], eta_adj_2d.shape[1])
    if eta0_from_file is None:
        print("ssh_eq.bin not found or incompatible; using zeros for initial eta plot.")
        eta0_from_file = np.zeros_like(eta_adj_2d)

    eta0_plot = ma.masked_where(land_mask_adj, eta0_from_file)
    depth_plot = ma.masked_where(land_mask_adj, depth_adj_2d)
    eta_plot = ma.masked_where(land_mask_adj, eta_adj_2d)
    u_plot = ma.masked_where(land_mask_adj, u_adj_2d)
    v_plot = ma.masked_where(land_mask_adj, v_adj_2d)
    w_plot = ma.masked_where(land_mask_adj, w_adj_2d)
    phivel_plot = ma.masked_where(land_mask_adj, phivel_adj_2d)
    psivel_plot = ma.masked_where(land_mask_adj, psivel_adj_2d)
    speed_h_plot = ma.masked_where(land_mask_adj, speed_h_adj_2d)

    _plot_cs_faces(xc_adj, "CS-32x32x1 Grid: XC", cmap="twilight", cbar_label="XC [deg]")
    _plot_cs_faces(yc_adj, "CS-32x32x1 Grid: YC", cmap="twilight", cbar_label="YC [deg]")
    _plot_cs_faces(wet_mask_adj.astype(float), "CS-32x32x1 Ocean Mask / Continent", cmap="Blues", cbar_label="wet mask")
    _plot_cs_faces(eta0_plot, "CS-32x32x1 Initial Surface Anomaly", cmap="seismic", cbar_label="eta [m]", symmetric=True)
    _plot_cs_faces(depth_plot, "CS-32x32x1 Depth", cmap="seismic", cbar_label="Depth [m]")
    _plot_cs_faces(eta_plot, f"CS-32x32x1 Eta (it={it_adj})", cmap="seismic", cbar_label="eta [m]", symmetric=True)
    _plot_cs_faces(u_plot, f"CS-32x32x1 U Velocity (it={it_adj})", cmap="seismic", cbar_label="u [m/s]", symmetric=True)
    _plot_cs_faces(v_plot, f"CS-32x32x1 V Velocity (it={it_adj})", cmap="seismic", cbar_label="v [m/s]", symmetric=True)
    _plot_cs_faces(w_plot, f"CS-32x32x1 W Velocity (it={it_adj})", cmap="seismic", cbar_label="w [m/s]", symmetric=True)
    _plot_cs_faces(phivel_plot, f"CS-32x32x1 PhiVEL (it={it_aux})", cmap="seismic", cbar_label="PhiVEL [m^2/s]", symmetric=True)
    _plot_cs_faces(psivel_plot, f"CS-32x32x1 PsiVEL (it={it_aux})", cmap="seismic", cbar_label="PsiVEL [m^3/s]", symmetric=True)
    _plot_cs_faces(speed_h_plot, f"CS-32x32x1 Speed sqrt(U^2+V^2) (it={it_adj})", cmap="seismic", cbar_label="speed [m/s]")
