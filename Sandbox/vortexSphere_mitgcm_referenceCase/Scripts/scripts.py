import copy
import os
import re
import sys
from contextlib import contextmanager, redirect_stdout
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib import colors as mcolors

sys.dont_write_bytecode = True
HAVE_MITGCMUTILS = True

EXPERIMENT_DIR = Path(__file__).resolve().parent.parent
SCRIPT_DIR = os.fspath(EXPERIMENT_DIR / "Scripts")
RUN_DIR = os.fspath(EXPERIMENT_DIR / "run")
INPUT_BATHY_PATH = os.fspath(EXPERIMENT_DIR / "input" / "bathymetry.bin")
SCRIPT_RESULTS_DIR = os.fspath(Path(SCRIPT_DIR) / "results")
RUN_PLOT_DIR = SCRIPT_RESULTS_DIR

def set_style():
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 400,
            "savefig.facecolor": "white",
            "font.size": 17,
        }
    )

set_style()

def _list_existing_dir(path, label):
    return os.listdir(os.fspath(path))

def _slugify_filename(text):
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    return stem.strip("._") or "figure"

def _add_right_colorbar(fig, ax, mappable, label, shrink=None):
    kwargs = {"ax": ax, "location": "right", "pad": 0.02, "fraction": 0.05}
    if shrink is not None:
        kwargs["shrink"] = shrink
    cbar = fig.colorbar(mappable, **kwargs)
    cbar.set_label(label)
    return cbar

def _save_publication_figure(fig, stem, enabled, out_dir=RUN_PLOT_DIR):
    if not enabled:
        return
    os.makedirs(out_dir, exist_ok=True)
    stem = _slugify_filename(stem)
    png_path = os.path.join(out_dir, f"{stem}.png")
    pdf_path = os.path.join(out_dir, f"{stem}.pdf")
    fig.savefig(png_path, dpi=400, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=400, bbox_inches="tight")
    print(f"Saved publication plots:\n  {png_path}\n  {pdf_path}")

def _iteration_folder_name(iteration):
    if iteration is None:
        return "latest"
    return f"it{int(iteration):010d}"

def _results_target_dir(iteration):
    return os.path.join(SCRIPT_RESULTS_DIR, _iteration_folder_name(iteration))

def _results_dyn_dir(iteration):
    return os.path.join(SCRIPT_RESULTS_DIR, "dyn", _iteration_folder_name(iteration))

class _TeeStream:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()

@contextmanager
def _capture_prints(iteration, dyn=False, heading=None):
    out_dir = _results_dyn_dir(iteration) if dyn else _results_target_dir(iteration)
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "plot_output.txt")
    with open(log_path, "a", encoding="utf-8") as log_file:
        tee = _TeeStream(sys.stdout, log_file)
        with redirect_stdout(tee):
            print("\n" + "=" * 80)
            if heading is not None:
                print(heading)
            print(f"Logging console output to: {log_path}")
            yield log_path

def _parse_mds_meta(meta_path):
    with open(meta_path, "r", encoding="utf-8") as meta_file:
        text = meta_file.read()

    n_dims = int(re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1))
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_vals = [int(value) for value in re.findall(r"[-+]?\d+", dim_block)]

    dims = []
    for index in range(n_dims):
        global_size, start, end = dim_vals[3 * index : 3 * index + 3]
        dims.append({"global": global_size, "start": start, "end": end, "n": end - start + 1})

    prec_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]", text)
    nrecords_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)\s*\]", text)
    it_match = re.search(r"timeStepNumber\s*=\s*\[\s*(\d+)\s*\]", text)

    return {
        "n_dims": n_dims,
        "dims": dims,
        "prec": prec_match.group(1).lower() if prec_match else "float32",
        "nrecords": int(nrecords_match.group(1)) if nrecords_match else 1,
        "it": int(it_match.group(1)) if it_match else 0,
    }

def _discover_field_files(run_dir, var_name):
    files = _list_existing_dir(run_dir, "MITgcm run directory")
    regex_global = re.compile(rf"^{re.escape(var_name)}\.(\d{{10}})\.meta$")
    regex_tiled = re.compile(rf"^{re.escape(var_name)}\.(\d{{10}})\.(\d{{3}})\.(\d{{3}})\.meta$")
    matches_global = []
    matches_tiled = []

    for fname in files:
        match_global = regex_global.match(fname)
        if match_global:
            matches_global.append({"it": int(match_global.group(1)), "file": fname, "type": "global"})
            continue

        match_tiled = regex_tiled.match(fname)
        if match_tiled:
            matches_tiled.append({"it": int(match_tiled.group(1)), "file": fname, "type": "tiled"})

    all_matches = matches_global + matches_tiled
    if not all_matches:
        raise FileNotFoundError(f"No output files found for variable '{var_name}' in {run_dir}")
    return all_matches

def _discover_timed_variables(run_dir):
    files = _list_existing_dir(run_dir, "MITgcm run directory")
    regex_global = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.(\d{10})\.meta$")
    regex_tiled = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.(\d{10})\.(\d{3})\.(\d{3})\.meta$")
    out = {}

    for fname in files:
        match = regex_global.match(fname)
        if match is None:
            match = regex_tiled.match(fname)
        if match is None:
            continue
        var_name = match.group(1)
        iteration = int(match.group(2))
        out.setdefault(var_name, set()).add(iteration)

    return {name: sorted(values) for name, values in sorted(out.items(), key=lambda item: item[0])}

def _discover_static_variables(run_dir):
    files = _list_existing_dir(run_dir, "MITgcm run directory")
    regex_static = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.meta$")
    regex_timed = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.(\d{10})(?:\.(\d{3})\.(\d{3}))?\.meta$")
    timed_names = set()
    static_names = set()

    for fname in files:
        match_timed = regex_timed.match(fname)
        if match_timed is not None:
            timed_names.add(match_timed.group(1))
            continue

        match_static = regex_static.match(fname)
        if match_static is not None:
            static_names.add(match_static.group(1))

    return sorted(static_names - timed_names)

def _meta_path_for_var_it(run_dir, var_name, iteration):
    matches = _discover_field_files(run_dir, var_name)
    files = sorted([match["file"] for match in matches if match["it"] == int(iteration)])
    if not files:
        raise FileNotFoundError(f"No meta file for {var_name} at iteration {iteration}")
    return os.path.join(run_dir, files[0])

def _parse_meta_fld_list(meta_path):
    with open(meta_path, "r", encoding="utf-8") as meta_file:
        text = meta_file.read()
    block = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    if block is None:
        return []
    raw = re.findall(r"'([^']+)'", block.group(1))
    return [name.strip() for name in raw]

def read_mitgcm_field(run_dir, var_name, it=None, record=0):
    all_matches = _discover_field_files(run_dir, var_name)
    available_its = sorted({match["it"] for match in all_matches})
    target_it = available_its[-1] if it is None else int(it)

    if target_it not in available_its:
        raise FileNotFoundError(
            f"Timestep {target_it} not found for {var_name}. Available: {available_its}"
        )

    print(f"Reading {var_name} at iteration {target_it}...")
    target_files = [match for match in all_matches if match["it"] == target_it]
    file_type = target_files[0]["type"]

    if file_type == "global":
        meta_path = os.path.join(run_dir, target_files[0]["file"])
        data_path = meta_path.replace(".meta", ".data")
        meta = _parse_mds_meta(meta_path)
        dtype = ">f8" if meta["prec"] in ["float64", "real*8"] else ">f4"
        raw_data = np.fromfile(data_path, dtype=dtype)
        nx = meta["dims"][0]["global"]
        ny = meta["dims"][1]["global"]

        if meta["n_dims"] == 2 and meta["nrecords"] == 1:
            return raw_data.reshape((ny, nx)), target_it
        if meta["n_dims"] == 2:
            return raw_data.reshape((meta["nrecords"], ny, nx))[record], target_it

        nz = meta["dims"][2]["global"]
        if meta["nrecords"] == 1:
            return raw_data.reshape((nz, ny, nx)), target_it
        return raw_data.reshape((meta["nrecords"], nz, ny, nx))[record], target_it

    first_meta = os.path.join(run_dir, sorted([tile["file"] for tile in target_files])[0])
    meta_global = _parse_mds_meta(first_meta)
    dtype = ">f8" if meta_global["prec"] in ["float64", "real*8"] else ">f4"
    nx_glob = meta_global["dims"][0]["global"]
    ny_glob = meta_global["dims"][1]["global"]

    if meta_global["n_dims"] == 2:
        full_field = np.zeros((ny_glob, nx_glob), dtype=np.float64)
    else:
        nz_glob = meta_global["dims"][2]["global"]
        full_field = np.zeros((nz_glob, ny_glob, nx_glob), dtype=np.float64)

    for tile in target_files:
        meta_path = os.path.join(run_dir, tile["file"])
        data_path = meta_path.replace(".meta", ".data")
        meta = _parse_mds_meta(meta_path)
        tile_data = np.fromfile(data_path, dtype=dtype)
        x_start = meta["dims"][0]["start"] - 1
        x_end = meta["dims"][0]["end"]
        y_start = meta["dims"][1]["start"] - 1
        y_end = meta["dims"][1]["end"]
        nx_local = x_end - x_start
        ny_local = y_end - y_start

        if meta["n_dims"] == 2:
            full_field[y_start:y_end, x_start:x_end] = tile_data.reshape((ny_local, nx_local))
            continue

        nz_local = meta["dims"][2]["end"] - meta["dims"][2]["start"] + 1
        full_field[:, y_start:y_end, x_start:x_end] = tile_data.reshape((nz_local, ny_local, nx_local))

    return full_field, target_it

def _surface_layer(field):
    return field[0] if field.ndim == 3 else field

def _lon_center_deg(i_lon, n_lon):
    return (i_lon + 0.5) * (360.0 / n_lon)

def _lon_centers_deg(n_lon):
    return (np.arange(n_lon) + 0.5) * (360.0 / n_lon)

def _lat_centers_deg(n_lat):
    return -90.0 + (np.arange(n_lat) + 0.5) * (180.0 / n_lat)

def report_ocean_extrema(field_ma, field_name, lat_deg, n_lon, depth_m):
    fld = ma.array(field_ma, copy=False)
    if fld.count() == 0:
        print(f"[{field_name}] no valid ocean points found.")
        return

    idx_min = np.unravel_index(ma.argmin(fld), fld.shape)
    idx_max = np.unravel_index(ma.argmax(fld), fld.shape)
    pairs = (("min", idx_min), ("max", idx_max))

    for tag, (ilat0, ilon0) in pairs:
        val = float(fld[ilat0, ilon0])
        lat0 = float(lat_deg[ilat0])
        lon0 = float(_lon_center_deg(ilon0, n_lon))
        depth_here = max(0.0, float(depth_m[ilat0, ilon0]))
        print(
            f"[{field_name}] {tag} = {val:.8e} at "
            f"iLat={ilat0+1}, iLon={ilon0+1}, "
            f"lat={lat0:.3f} deg, lon={lon0:.3f} deg, "
            f"level=k=1 surface, bathymetry H={depth_here:.3f} m"
        )

def _masked_min_max(field_ma):
    fld = ma.array(field_ma, copy=False)
    if fld.count() == 0:
        return np.nan, np.nan
    return float(ma.min(fld)), float(ma.max(fld))

def _print_ocean_min_max(field_2d, ocean_mask, field_name, lat_deg=None, depth_m=None):
    fld = ma.masked_where(~ocean_mask, ma.array(field_2d, copy=False))
    vmin, vmax = _masked_min_max(fld)
    print(f"[{field_name}] ocean min={vmin:.8e}, max={vmax:.8e}")
    if lat_deg is not None and depth_m is not None:
        report_ocean_extrema(fld, field_name, lat_deg, field_2d.shape[1], depth_m)

def _safe_sym_limit(arr):
    arr_np = ma.filled(ma.array(arr, copy=False), np.nan)
    limit = np.nanmax(np.abs(arr_np))
    return limit if np.isfinite(limit) and limit > 0 else 1e-12

def _safe_max(arr):
    arr_np = ma.filled(ma.array(arr, copy=False), np.nan)
    limit = np.nanmax(arr_np)
    return limit if np.isfinite(limit) and limit > 0 else 1e-12

def _load_depth_and_ocean_mask(n_lat, n_lon):
    bathy = np.fromfile(INPUT_BATHY_PATH, dtype=">f4").reshape((n_lat, n_lon))
    ocean_mask = bathy < 0.0
    depth_m = np.maximum(0.0, -bathy)
    print(f"Using bathymetry mask from {INPUT_BATHY_PATH} (ocean where depth < 0).")
    return depth_m, ocean_mask

def _load_optional_hfac_masks(it, n_lat, n_lon):
    hfac_c, _ = read_mitgcm_field(RUN_DIR, "hFacC", it=it)
    hfac_w, _ = read_mitgcm_field(RUN_DIR, "hFacW", it=it)
    hfac_s, _ = read_mitgcm_field(RUN_DIR, "hFacS", it=it)
    tracer_mask = _surface_layer(hfac_c) > 0.0
    u_mask = _surface_layer(hfac_w) > 0.0
    v_mask = _surface_layer(hfac_s) > 0.0
    print("Using hFacC/hFacW/hFacS masks from MITgcm output.")
    return tracer_mask, u_mask, v_mask

def _report_depth_consistency(model_depth_2d, input_depth_m, input_ocean_mask, hfac_c_2d=None, verbose=True):
    model_depth = np.asarray(model_depth_2d, dtype=float)
    input_depth = np.asarray(input_depth_m, dtype=float)
    input_ocean_mask = np.asarray(input_ocean_mask, dtype=bool)
    model_wet_mask = model_depth > 0.0
    overlap_mask = input_ocean_mask & model_wet_mask
    wet_mismatch_count = int(np.count_nonzero(input_ocean_mask ^ model_wet_mask))
    depth_diff_map = np.full(model_depth.shape, np.nan, dtype=float)

    if np.any(overlap_mask):
        depth_diff_map[overlap_mask] = model_depth[overlap_mask] - input_depth[overlap_mask]

    unique_depth_count = int(np.unique(model_depth[overlap_mask]).size) if np.any(overlap_mask) else 0
    max_abs_diff = float(np.nanmax(np.abs(depth_diff_map[overlap_mask]))) if np.any(overlap_mask) else np.nan
    match_fraction = (
        float(np.mean(np.isclose(model_depth[overlap_mask], input_depth[overlap_mask], rtol=0.0, atol=1e-9)))
        if np.any(overlap_mask)
        else np.nan
    )

    if verbose:
        print("Bathymetry consistency check:")
        print(f"  input ocean cells = {int(np.count_nonzero(input_ocean_mask))}")
        print(f"  model wet cells   = {int(np.count_nonzero(model_wet_mask))}")
        print(f"  overlap cells     = {int(np.count_nonzero(overlap_mask))}")
        print(f"  wet-mask mismatch = {wet_mismatch_count}")
        if np.any(overlap_mask):
            print(f"  max |Depth - input bathy| = {max_abs_diff:.3e} m")
            print(f"  bathy match within 1e-9 m = {100.0 * match_fraction:.2f}% of overlap cells")
            print(f"  unique model wet depths   = {unique_depth_count}")
            if unique_depth_count > 1 and wet_mismatch_count == 0 and max_abs_diff <= 1e-9:
                print("  CONFIRMED: wet-cell Depth matches the bathy file and is not uniform.")
        if hfac_c_2d is not None and np.any(overlap_mask):
            wet_hfac = np.asarray(hfac_c_2d, dtype=float)[overlap_mask]
            print(
                "  hFacC wet range           = "
                f"[{float(np.nanmin(wet_hfac)):.6g}, {float(np.nanmax(wet_hfac)):.6g}]"
            )
            print(f"  hFacC unique wet values   = {int(np.unique(wet_hfac).size)}")

    return {
        "model_wet_mask": model_wet_mask,
        "overlap_mask": overlap_mask,
        "depth_diff_map": depth_diff_map,
        "max_abs_diff": max_abs_diff,
        "wet_mismatch_count": wet_mismatch_count,
        "unique_depth_count": unique_depth_count,
        "match_fraction": match_fraction,
    }

def read_mitgcm_static_field(run_dir, var_name, record=0):
    meta_path = os.path.join(run_dir, f"{var_name}.meta")
    data_path = os.path.join(run_dir, f"{var_name}.data")
    if not os.path.exists(meta_path) or not os.path.exists(data_path):
        raise FileNotFoundError(f"Static field {var_name}.meta/.data not found in {run_dir}")

    meta = _parse_mds_meta(meta_path)
    dtype = ">f8" if meta["prec"] in ["float64", "real*8"] else ">f4"
    raw_data = np.fromfile(data_path, dtype=dtype)

    if meta["n_dims"] == 1 and meta["nrecords"] == 1:
        return raw_data.reshape((meta["dims"][0]["global"],))
    if meta["n_dims"] == 1:
        return raw_data.reshape((meta["nrecords"], meta["dims"][0]["global"]))[record]

    if meta["n_dims"] == 2 and meta["nrecords"] == 1:
        return raw_data.reshape((meta["dims"][1]["global"], meta["dims"][0]["global"]))
    if meta["n_dims"] == 2:
        return raw_data.reshape((meta["nrecords"], meta["dims"][1]["global"], meta["dims"][0]["global"]))[record]

    if meta["nrecords"] == 1:
        return raw_data.reshape((meta["dims"][2]["global"], meta["dims"][1]["global"], meta["dims"][0]["global"]))
    return raw_data.reshape(
        (meta["nrecords"], meta["dims"][2]["global"], meta["dims"][1]["global"], meta["dims"][0]["global"])
    )[record]

def read_mitgcm_field_or_static(run_dir, var_name, it=None, record=0):
    timed_meta = list(Path(run_dir).glob(f"{var_name}.*.meta"))
    if timed_meta:
        return read_mitgcm_field(run_dir, var_name, it=it, record=record)
    return read_mitgcm_static_field(run_dir, var_name, record=record), 0

def _finite_plot_limits(fld_np, symmetric=False):
    if symmetric:
        limit = np.nanmax(np.abs(fld_np))
        if not np.isfinite(limit) or limit <= 0:
            limit = 1e-12
        return -limit, limit

    vmin = np.nanmin(fld_np)
    vmax = np.nanmax(fld_np)
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return 0.0, 1.0
    if vmin == vmax:
        delta = max(abs(vmin) * 0.05, 1.0)
        return vmin - delta, vmax + delta
    return vmin, vmax

def _plot_scale_kwargs(fld_np, cmap="seismic", symmetric=False, positive_only=False):
    fld_np = np.asarray(fld_np, dtype=float)
    cmap_name = str(cmap).lower()

    if cmap_name == "seismic":
        limit = np.nanmax(np.abs(fld_np))
        if not np.isfinite(limit) or limit <= 0:
            limit = 1e-12
        return {"norm": mcolors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)}

    if positive_only:
        vmax = np.nanmax(fld_np)
        if not np.isfinite(vmax) or vmax <= 0.0:
            vmax = 1e-12
        return {"vmin": 0.0, "vmax": vmax}

    vmin, vmax = _finite_plot_limits(fld_np, symmetric=symmetric)
    return {"vmin": vmin, "vmax": vmax}

def _unit_for_field(field_name):
    units = {
        "Eta": "m",
        "ETAN": "m",
        "ETANSQ": "m^2",
        "DETADT2": "m^2/s^2",
        "U": "m/s",
        "V": "m/s",
        "W": "m/s",
        "UVELMASS": "m/s",
        "VVELMASS": "m/s",
        "WVELMASS": "m/s",
        "T": "degC",
        "S": "psu",
        "PH": "m^2/s^2",
        "PHL": "m^2/s^2",
        "PhiVEL": "m^2/s",
        "PsiVEL": "m^3/s",
        "Depth": "m",
    }
    return units.get(field_name, "")

def _auto_colormap_and_symmetry(field_2d, wet_mask):
    wet_vals = np.asarray(field_2d)[wet_mask]
    wet_vals = wet_vals[np.isfinite(wet_vals)]
    if wet_vals.size == 0:
        return "seismic", False
    if np.nanmin(wet_vals) < 0.0 and np.nanmax(wet_vals) > 0.0:
        return "seismic", True
    return "seismic", False

def _plot_global_latlon(
    field_2d,
    invalid_mask,
    title,
    cbar_label,
    cmap="seismic",
    symmetric=False,
    positive_only=False,
    save_stem=None,
    save_enabled=False,
    save_dir=None,
):
    fld = ma.masked_where(invalid_mask, ma.array(field_2d, copy=False))
    fld_np = ma.filled(fld, np.nan)
    plot_kwargs = _plot_scale_kwargs(fld_np, cmap=cmap, symmetric=symmetric, positive_only=positive_only)
    cmap_obj = copy.copy(plt.get_cmap(cmap))
    cmap_obj.set_bad("0.9")

    fig, ax = plt.subplots(figsize=(12.5, 6.2), constrained_layout=True)
    im = ax.imshow(
        fld,
        origin="lower",
        extent=[0, 360, -90, 90],
        cmap=cmap_obj,
        aspect="equal",
        interpolation="none",
        **plot_kwargs,
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.set_yticks(np.arange(-90.0, 91.0, 30.0))
    _add_right_colorbar(fig, ax, im, cbar_label)
    _save_publication_figure(fig, save_stem or title, save_enabled, out_dir=save_dir or RUN_PLOT_DIR)
    plt.show()

def _select_reference_iteration(it_list, requested_it, label, quiet=False):
    if requested_it is None:
        return it_list[-1]
    requested_it = int(requested_it)
    if requested_it in it_list:
        return requested_it
    fallback_it = it_list[-1]
    if not quiet:
        print(f"[{label}] requested it={requested_it} not available; using latest it={fallback_it}.")
    return fallback_it

def _iter_record_specs(run_dir, var_name, it_plot):
    meta_path = _meta_path_for_var_it(run_dir, var_name, it_plot)
    meta = _parse_mds_meta(meta_path)
    field_names = _parse_meta_fld_list(meta_path)

    if meta["nrecords"] == 1:
        return [(0, var_name)]
    if len(field_names) == meta["nrecords"]:
        return [(index, field_names[index]) for index in range(meta["nrecords"])]
    return [(index, f"{var_name}:record{index}") for index in range(meta["nrecords"])]

def _suggest_additional_reference_plots(run_dir):
    static_vars = set(_discover_static_variables(run_dir))
    suggestions = []
    extra_maps = [name for name in ("hFacC", "hFacW", "hFacS", "PHrefC", "PHrefF") if name in static_vars]
    extra_profiles = [name for name in ("RC", "RF", "DRF", "RhoRef") if name in static_vars]

    if extra_maps:
        suggestions.append("2-D maps: " + ", ".join(extra_maps))
    if extra_profiles:
        suggestions.append("1-D profiles: " + ", ".join(extra_profiles))
    if os.path.exists(os.path.join(run_dir, "dynStDiag.0000000000.txt")):
        suggestions.append("time-summary text file: dynStDiag.0000000000.txt")

    if suggestions:
        print("Other things worth plotting from this run:")
        for item in suggestions:
            print(f"  - {item}")

def _land_leak_check(field_2d, land_mask, field_name, allow_land_values=False):
    arr = np.asarray(field_2d)
    land_vals = arr[land_mask]
    land_vals = land_vals[np.isfinite(land_vals)]

    if land_vals.size == 0:
        print(f"[{field_name}] land check: no finite land points.")
        return

    ocean_vals = arr[~land_mask]
    ocean_vals = ocean_vals[np.isfinite(ocean_vals)]
    scale = np.nanmax(np.abs(ocean_vals)) if ocean_vals.size else 0.0
    tol = max(1e-12, 1e-6 * scale)
    bad = np.abs(land_vals) > tol

    if np.any(bad):
        max_bad = float(np.nanmax(np.abs(land_vals[bad])))
        msg = (
            f"[{field_name}] WARNING land leak: {bad.sum()} points "
            f"(|val| max={max_bad:.3e}, tol={tol:.3e})"
        )
        if allow_land_values:
            print(msg + " (allowed for this field)")
        else:
            print(msg)
        return

    print(f"[{field_name}] land check: PASS (tol={tol:.3e}).")

def _print_eta_init_ocean_land_check():
    nx, ny = 1440, 720
    input_dir = Path(SCRIPT_DIR).parent / "input"
    eta_init = np.fromfile(input_dir / "eta_init.bin", dtype=">f4").reshape(ny, nx)
    bathy = np.fromfile(input_dir / "bathymetry.bin", dtype=">f4").reshape(ny, nx)
    ocean = bathy < 0
    land = ~ocean
    tol = max(1e-10, 1e-6 * np.nanmax(np.abs(eta_init)))
    signal = np.abs(eta_init) > tol
    north = np.zeros_like(land)
    north[1:, :] = land[:-1, :]
    south = np.zeros_like(land)
    south[:-1, :] = land[1:, :]
    coast_ocean = ocean & (north | south | np.roll(land, 1, 1) | np.roll(land, -1, 1))

    print("eta_init location check:")
    print("  FAIL land leakage" if np.any(signal & land) else "  PASS no land leakage")
    print("  FAIL coastline overlap" if np.any(signal & coast_ocean) else "  PASS no coastline overlap")
