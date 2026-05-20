#%%
from pathlib import Path
import glob
import importlib.util
import os
import re
import shutil

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

def env_flag(name, default):
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}

def resolve_script_dir():
    if "__file__" in globals():
        return Path(__file__).resolve().parent

    cwd = Path.cwd().resolve()
    candidates = [
        cwd,
        cwd / "Scripts",
        cwd / "verification" / "vortexSphere_mitgcm_referenceCase" / "Scripts",
    ]

    for candidate in candidates:
        if (candidate.parent / "run").exists():
            return candidate

    return cwd

def output_path(script_dir, file_name):
    path = Path(file_name)
    if path.is_absolute():
        return path
    return Path(script_dir) / path

def set_plotly_renderer():
    if os.environ.get("PLOTLY_RENDERER"):
        pio.renderers.default = os.environ["PLOTLY_RENDERER"]
    elif os.environ.get("VSCODE_PID"):
        pio.renderers.default = "vscode"

def parse_mds_meta(meta_path):
    text = Path(meta_path).read_text(encoding="utf-8")
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

def parse_meta_fld_list(meta_path):
    text = Path(meta_path).read_text(encoding="utf-8")
    block = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    if block is None:
        return []
    return [name.strip() for name in re.findall(r"'([^']+)'", block.group(1))]

def discover_timed_variables(run_dir):
    run_dir = Path(run_dir)
    regex_global = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.(\d{10})\.meta$")
    regex_tiled = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.(\d{10})\.(\d{3})\.(\d{3})\.meta$")
    out = {}

    for path in run_dir.glob("*.meta"):
        match = regex_global.match(path.name) or regex_tiled.match(path.name)
        if match is None:
            continue
        out.setdefault(match.group(1), set()).add(int(match.group(2)))

    return {name: sorted(values) for name, values in sorted(out.items())}

def discover_static_variables(run_dir):
    run_dir = Path(run_dir)
    regex_static = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.meta$")
    regex_timed = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\.(\d{10})(?:\.(\d{3})\.(\d{3}))?\.meta$")
    static_names = set()
    timed_names = set()

    for path in run_dir.glob("*.meta"):
        timed = regex_timed.match(path.name)
        if timed is not None:
            timed_names.add(timed.group(1))
            continue
        static = regex_static.match(path.name)
        if static is not None:
            static_names.add(static.group(1))

    return sorted(static_names - timed_names)

def discover_binary_iterations(data_dir, field_name):
    files = sorted(glob.glob(str(Path(data_dir) / f"{field_name}.*.data")))
    iterations = []
    for file_path in files:
        parts = Path(file_path).name.split(".")
        if len(parts) >= 3 and parts[0] == field_name:
            try:
                iterations.append(int(parts[1]))
            except ValueError:
                pass
    return sorted(set(iterations))

def meta_path_for_var_it(run_dir, var_name, iteration):
    run_dir = Path(run_dir)
    patterns = [
        f"{var_name}.{int(iteration):010d}.meta",
        f"{var_name}.{int(iteration):010d}.*.*.meta",
    ]
    matches = []
    for pattern in patterns:
        matches.extend(sorted(run_dir.glob(pattern)))
    if not matches:
        raise FileNotFoundError(f"No meta file for {var_name} at iteration {iteration}")
    return matches[0]

def iter_record_specs(run_dir, var_name, iteration):
    meta_path = meta_path_for_var_it(run_dir, var_name, iteration)
    meta = parse_mds_meta(meta_path)
    field_names = parse_meta_fld_list(meta_path)

    if meta["nrecords"] == 1:
        return [(0, var_name)]
    if len(field_names) == meta["nrecords"]:
        return [(index, field_names[index]) for index in range(meta["nrecords"])]
    return [(index, f"{var_name}:record{index}") for index in range(meta["nrecords"])]

def print_inventory(run_dir):
    run_dir = Path(run_dir)
    timed = discover_timed_variables(run_dir)
    static = discover_static_variables(run_dir)

    print(f"RUN_DIR = {run_dir}")
    print("\nStatic variables:")
    for name in static:
        print(f"  {name}")

    print("\nTime-dependent variables:")
    for name, iterations in timed.items():
        print(f"  {name}: {iterations}")

    print("\nDynamic/bundled records:")
    for name, iterations in timed.items():
        if not iterations:
            continue
        try:
            records = iter_record_specs(run_dir, name, iterations[0])
        except Exception:
            continue
        if len(records) <= 1:
            continue
        print(f"  {name}: iterations {iterations}")
        for rec_idx, rec_name in records:
            print(f"    record {rec_idx}: {rec_name}")

    return {"static": static, "timed": timed}

def dtype_from_meta(meta):
    return ">f8" if meta["prec"] in {"float64", "real*8"} else ">f4"

def read_global_or_tiled_field(run_dir, var_name, iteration, record=0):
    run_dir = Path(run_dir)
    target = int(iteration)
    global_meta = run_dir / f"{var_name}.{target:010d}.meta"

    if global_meta.exists():
        meta = parse_mds_meta(global_meta)
        raw = np.fromfile(global_meta.with_suffix(".data"), dtype=dtype_from_meta(meta))
        nx = meta["dims"][0]["global"]
        ny = meta["dims"][1]["global"]

        if meta["n_dims"] == 2 and meta["nrecords"] == 1:
            return raw.reshape((ny, nx)), target
        if meta["n_dims"] == 2:
            return raw.reshape((meta["nrecords"], ny, nx))[record], target

        nz = meta["dims"][2]["global"]
        if meta["nrecords"] == 1:
            return raw.reshape((nz, ny, nx)), target
        return raw.reshape((meta["nrecords"], nz, ny, nx))[record], target

    tile_metas = sorted(run_dir.glob(f"{var_name}.{target:010d}.*.*.meta"))
    if not tile_metas:
        raise FileNotFoundError(f"No data found for {var_name} at iteration {target}")

    first_meta = parse_mds_meta(tile_metas[0])
    dtype = dtype_from_meta(first_meta)
    nx_glob = first_meta["dims"][0]["global"]
    ny_glob = first_meta["dims"][1]["global"]

    if first_meta["n_dims"] == 2:
        full = np.zeros((ny_glob, nx_glob), dtype=np.float64)
    else:
        nz_glob = first_meta["dims"][2]["global"]
        full = np.zeros((nz_glob, ny_glob, nx_glob), dtype=np.float64)

    for meta_path in tile_metas:
        meta = parse_mds_meta(meta_path)
        raw = np.fromfile(meta_path.with_suffix(".data"), dtype=dtype)
        x0 = meta["dims"][0]["start"] - 1
        x1 = meta["dims"][0]["end"]
        y0 = meta["dims"][1]["start"] - 1
        y1 = meta["dims"][1]["end"]
        nx = x1 - x0
        ny = y1 - y0

        if meta["n_dims"] == 2:
            tile = raw.reshape((meta["nrecords"], ny, nx))[record] if meta["nrecords"] > 1 else raw.reshape((ny, nx))
            full[y0:y1, x0:x1] = tile
            continue

        nz = meta["dims"][2]["end"] - meta["dims"][2]["start"] + 1
        tile = raw.reshape((meta["nrecords"], nz, ny, nx))[record] if meta["nrecords"] > 1 else raw.reshape((nz, ny, nx))
        full[:, y0:y1, x0:x1] = tile

    return full, target

def read_static_field(run_dir, var_name, record=0):
    meta_path = Path(run_dir) / f"{var_name}.meta"
    data_path = Path(run_dir) / f"{var_name}.data"
    if not meta_path.exists() or not data_path.exists():
        raise FileNotFoundError(f"Static field {var_name}.meta/.data not found in {run_dir}")

    meta = parse_mds_meta(meta_path)
    raw = np.fromfile(data_path, dtype=dtype_from_meta(meta))

    if meta["n_dims"] == 1 and meta["nrecords"] == 1:
        return raw.reshape((meta["dims"][0]["global"],))
    if meta["n_dims"] == 1:
        return raw.reshape((meta["nrecords"], meta["dims"][0]["global"]))[record]

    if meta["n_dims"] == 2 and meta["nrecords"] == 1:
        return raw.reshape((meta["dims"][1]["global"], meta["dims"][0]["global"]))
    if meta["n_dims"] == 2:
        return raw.reshape((meta["nrecords"], meta["dims"][1]["global"], meta["dims"][0]["global"]))[record]

    if meta["nrecords"] == 1:
        return raw.reshape((meta["dims"][2]["global"], meta["dims"][1]["global"], meta["dims"][0]["global"]))
    return raw.reshape((meta["nrecords"], meta["dims"][2]["global"], meta["dims"][1]["global"], meta["dims"][0]["global"]))[record]

def surface_layer(field):
    return field[0] if np.ndim(field) == 3 else field

def choose_record_index(run_dir, field_name, iteration, record_name=None, record_index=None):
    if record_index is not None:
        return int(record_index)
    if record_name is None:
        return 0

    for rec_idx, rec_name in iter_record_specs(run_dir, field_name, iteration):
        if rec_name == record_name:
            return rec_idx
    raise ValueError(f"Record {record_name!r} not found in {field_name}")

def read_selected_field(run_dir, field_name, iteration, record_name=None, record_index=None):
    rec_idx = choose_record_index(run_dir, field_name, iteration, record_name, record_index)
    field, it_used = read_global_or_tiled_field(run_dir, field_name, iteration, record=rec_idx)
    return surface_layer(field).astype(np.float32, copy=False), it_used

def load_land_mask(run_dir, ny=None, nx=None):
    run_dir = Path(run_dir)
    try:
        depth = surface_layer(read_static_field(run_dir, "Depth"))
        return depth <= 0.0
    except Exception:
        depth_path = run_dir / "Depth.data"
        if not depth_path.exists() or ny is None or nx is None:
            print("No usable Depth field found; land will not be masked.")
            return None
        depth = np.fromfile(depth_path, dtype=">f8")
        if depth.size != ny * nx:
            print("Depth.data size does not match grid; land will not be masked.")
            return None
        return depth.reshape(ny, nx) <= 0.0

def mask_land(field, land_mask):
    if land_mask is None:
        return field
    out = field.copy()
    out[land_mask] = np.nan
    return out

def grid_1d(nx, ny, dx_deg=None, dy_deg=None, latitude_order="south_to_north", lon_shift=False):
    dx = 360.0 / nx if dx_deg is None else dx_deg
    dy = 180.0 / ny if dy_deg is None else dy_deg
    lon = (np.arange(nx) + 0.5) * dx
    if lon_shift:
        lon = lon - 180.0

    if latitude_order == "north_to_south":
        lat = 90.0 - (np.arange(ny) + 0.5) * dy
    elif latitude_order == "south_to_north":
        lat = -90.0 + (np.arange(ny) + 0.5) * dy
    else:
        raise ValueError('latitude_order must be "north_to_south" or "south_to_north"')

    return lon, lat

def read_series(run_dir, field_name, iterations, land_mask, html_skip, record_name=None, record_index=None):
    fields = []
    for iteration in iterations:
        field, _ = read_selected_field(run_dir, field_name, iteration, record_name, record_index)
        fields.append(mask_land(field, land_mask)[::html_skip, ::html_skip])
    return fields

def finite_stats(field):
    return {
        "min": float(np.nanmin(field)),
        "max": float(np.nanmax(field)),
        "mean": float(np.nanmean(field)),
        "std": float(np.nanstd(field)),
    }

def nice_number(value):
    return f"{value:.6g}"

def value_text(value, unit=""):
    text = nice_number(value)
    return f"{text} {unit}" if unit else text

def symmetric_limit(fields, color_limit=None):
    if color_limit is not None:
        limit = abs(float(color_limit))
        if np.isfinite(limit) and limit > 0.0:
            return limit
    zmax = max(float(np.nanmax(np.abs(field))) for field in fields)
    return zmax if np.isfinite(zmax) and zmax > 0.0 else 1e-12

def label_with_unit(name, units):
    unit = units.get(name, "")
    return f"{name} [{unit}]" if unit else name

def title_name(name, descriptions):
    description = descriptions.get(name, "")
    return f"{description} {name}" if description else name

def stats_box_text(field_label, run_dir, iterations, iteration, field, all_stats, zmax, unit, html_skip):
    current = finite_stats(field)
    return (
        f"<b>{field_label}</b><br>"
        f"path: {run_dir}<br>"
        f"iterations: {len(iterations)} found<br>"
        f"first/last: {iterations[0]} / {iterations[-1]}<br>"
        f"current iteration: {iteration}<br>"
        f"grid skip: {html_skip}<br>"
        f"color scale: +/- {value_text(zmax, unit)}<br><br>"
        f"<b>Selected iteration</b><br>"
        f"min: {value_text(current['min'], unit)}<br>"
        f"max: {value_text(current['max'], unit)}<br>"
        f"mean: {value_text(current['mean'], unit)}<br>"
        f"std: {value_text(current['std'], unit)}<br><br>"
        f"<b>All loaded frames</b><br>"
        f"min: {value_text(all_stats['min'], unit)}<br>"
        f"max: {value_text(all_stats['max'], unit)}<br>"
        f"mean: {value_text(all_stats['mean'], unit)}"
    )

def info_annotation(text, info_box_x):
    return dict(
        text=text,
        x=info_box_x,
        y=0.95,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        align="left",
        showarrow=False,
        bordercolor="rgba(30, 30, 30, 0.35)",
        borderwidth=1,
        borderpad=10,
        bgcolor="rgba(255, 255, 255, 0.94)",
        font=dict(size=12),
    )

def plotly_colorscale(cmap):
    cmap = str(cmap).strip()

    matplotlib_to_plotly = {
        "seismic": "bluered",
        "seismic_r": "bluered_r",
        "RdBu": "rdbu",
        "RdBu_r": "rdbu_r",
        "coolwarm": "balance",
        "coolwarm_r": "balance_r",
        "Spectral": "spectral",
        "Spectral_r": "spectral_r",
        "PiYG": "piyg",
        "PiYG_r": "piyg_r",
        "PRGn": "prgn",
        "PRGn_r": "prgn_r",
        "BrBG": "brbg",
        "BrBG_r": "brbg_r",
        "viridis": "viridis",
        "plasma": "plasma",
        "inferno": "inferno",
        "magma": "magma",
        "cividis": "cividis",
        "turbo": "turbo",
    }

    return matplotlib_to_plotly.get(cmap, cmap.lower())

def make_heatmap(x, y, field, zmax, colorscale, colorbar_label, colorbar_x, show_colorbar=True):
    return go.Heatmap(
        x=x,
        y=y,
        z=field,
        zmin=-zmax,
        zmax=zmax,
        colorscale=colorscale,
        zsmooth="best",
        showscale=show_colorbar,
        colorbar=dict(
            title=colorbar_label,
            x=colorbar_x,
            xanchor="left",
            y=0.56,
            len=0.74,
            thickness=18,
        ),
        hovertemplate=(
            "lon: %{x:.2f} deg<br>"
            "lat: %{y:.2f} deg<br>"
            f"{colorbar_label}: "
            "%{z:.6g}<extra></extra>"
        ),
    )

def land_heatmap(x, y, land_mask, html_skip):
    if land_mask is None:
        return None
    land_z = np.where(land_mask[::html_skip, ::html_skip], 1.0, np.nan).astype(np.float32)
    return go.Heatmap(
        x=x,
        y=y,
        z=land_z,
        zmin=0.0,
        zmax=1.0,
        colorscale=[
    [0.0, "rgb(235,235,235)"],
    [1.0, "rgb(235,235,235)"],
],
opacity=0.65,
        showscale=False,
        hoverinfo="skip",
        name="land",
    )

def build_slider_figure(
    fields,
    iterations,
    lon,
    lat,
    land_mask,
    run_dir,
    field_label,
    unit,
    zmax,
    cmap,
    colorbar_label,
    html_skip,
    width=1600,
    height=850,
    map_domain_x=(0.0, 0.76),
    colorbar_x=0.80,
    info_box_x=0.86,
    frame_duration_ms=180,
):
    colorscale = plotly_colorscale(cmap)
    all_stats = finite_stats(np.asarray(fields))

    def annotation_for(index):
        text = stats_box_text(
            field_label,
            run_dir,
            iterations,
            iterations[index],
            fields[index],
            all_stats,
            zmax,
            unit,
            html_skip,
        )
        return info_annotation(text, info_box_x)

    data = [make_heatmap(lon, lat, fields[0], zmax, colorscale, colorbar_label, colorbar_x, True)]
    land = land_heatmap(lon, lat, land_mask, html_skip)
    if land is not None:
        data.append(land)

    fig = go.Figure(data=data)
    fig.frames = [
        go.Frame(
            data=[go.Heatmap(z=field, zsmooth="best")],
            traces=[0],
            layout=dict(annotations=[annotation_for(i)]),
            name=str(iteration),
        )
        for i, (field, iteration) in enumerate(zip(fields, iterations))
    ]

    steps = [
        dict(
            method="animate",
            args=[
                [str(iteration)],
                dict(mode="immediate", frame=dict(duration=350, redraw=True), transition=dict(duration=0)),
            ],
            label=str(iteration),
        )
        for iteration in iterations
    ]

    buttons = [
        dict(
            type="buttons",
            direction="left",
            x=0.08,
            y=0,
            showactive=False,
            buttons=[
                dict(
                    label="Play",
                    method="animate",
                    args=[
                        None,
                        dict(
                            frame=dict(duration=frame_duration_ms, redraw=True),
                            transition=dict(duration=0),
                            fromcurrent=True,
                        ),
                    ],
                ),
                dict(
                    label="Pause",
                    method="animate",
                    args=[
                        [None],
                        dict(frame=dict(duration=0, redraw=False), mode="immediate", transition=dict(duration=0)),
                    ],
                ),
            ],
        )
    ]

    fig.update_layout(
        title=dict(text=f"MITgcm {field_label} evolution", x=0.36, xanchor="center"),
        width=width,
        height=height,
        template="plotly_white",
        margin=dict(l=55, r=35, t=75, b=60),
        sliders=[dict(active=0, steps=steps, x=0.12, y=0, len=0.60, currentvalue=dict(prefix="Iteration: "))],
        updatemenus=buttons,
        xaxis=dict(domain=list(map_domain_x), title="Longitude [deg]", range=[float(np.nanmin(lon)), float(np.nanmax(lon))]),
        yaxis=dict(domain=[0.08, 1.0], title="Latitude [deg]", range=[float(np.nanmin(lat)), float(np.nanmax(lat))], scaleanchor="x", scaleratio=1.0),
        annotations=[annotation_for(0)],
        font=dict(size=14),
    )
    return fig

def show_or_write_figure(fig, html_path, write_html=True, open_browser=False, show_figure=True):
    html_path = Path(html_path)
    if write_html:
        fig.write_html(
            html_path,
            include_plotlyjs=True,
            auto_open=open_browser,
            auto_play=False,
            config={"responsive": True, "toImageButtonOptions": {"format": "png", "scale": 4}},
        )
        print(f"Interactive HTML written to: {html_path}")

    if show_figure:
        try:
            fig.show()
        except Exception as exc:
            print(f"fig.show() failed: {exc}")
            print(f"Open the written HTML file instead: {html_path}")

def ffmpeg_available():
    return shutil.which("ffmpeg") is not None or importlib.util.find_spec("imageio_ffmpeg") is not None

def choose_video_output(mp4_path, gif_path):
    mp4_path = Path(mp4_path)
    if mp4_path.suffix.lower() == ".mp4" and not ffmpeg_available():
        gif_path = Path(gif_path)
        print(f"ffmpeg is not available. Writing GIF instead: {gif_path}")
        return gif_path
    return mp4_path

def write_video(
    run_dir,
    field_name,
    iterations,
    land_mask,
    zmax,
    cmap,
    colorbar_label,
    title,
    mp4_path,
    gif_path,
    record_name=None,
    record_index=None,
    video_skip=2,
    frame_duration_ms=180,
    width=1800,
    height=950,
    latitude_order="south_to_north",
):
    import imageio.v2 as imageio

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors

    video_path = choose_video_output(mp4_path, gif_path)
    fps = max(1, int(round(1000 / frame_duration_ms)))

    first, _ = read_selected_field(run_dir, field_name, iterations[0], record_name, record_index)
    first = mask_land(first, land_mask)[::video_skip, ::video_skip]

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("0.78")
    norm = mcolors.TwoSlopeNorm(vmin=-zmax, vcenter=0.0, vmax=zmax)

    fig, ax = plt.subplots(figsize=(width / 200.0, height / 200.0), dpi=200, constrained_layout=True)
    image = ax.imshow(
        first,
        origin="lower" if latitude_order == "south_to_north" else "upper",
        extent=[0.0, 360.0, -90.0, 90.0],
        cmap=cmap_obj,
        norm=norm,
        interpolation="hanning",
        aspect="auto",
    )
    cbar = fig.colorbar(image, ax=ax, location="right", pad=0.02, fraction=0.04)
    cbar.set_label(colorbar_label)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.set_yticks(np.arange(-90.0, 91.0, 30.0))

    with imageio.get_writer(video_path, mode="I", fps=fps) as writer:
        for n, iteration in enumerate(iterations, start=1):
            field, _ = read_selected_field(run_dir, field_name, iteration, record_name, record_index)
            field = mask_land(field, land_mask)[::video_skip, ::video_skip]
            image.set_data(field)
            ax.set_title(f"{title} | iteration {iteration}")
            fig.canvas.draw()
            rgba = np.asarray(fig.canvas.buffer_rgba())
            writer.append_data(rgba[:, :, :3].copy())
            if n == 1 or n % 25 == 0 or n == len(iterations):
                print(f"  rendered video frame {n}/{len(iterations)}")

    plt.close(fig)
    print(f"Animation written to: {video_path}")

def collect_inventory(run_dir):
    static_variables = discover_static_variables(run_dir)
    timed_variables = discover_timed_variables(run_dir)

    dynamic_records = {}
    for name in timed_variables:
        try:
            it0 = timed_variables[name][0]
            records = iter_record_specs(run_dir, name, it0)
            if len(records) > 1:
                dynamic_records[name] = {
                    "iterations": timed_variables[name],
                    "records": records,
                }
        except Exception:
            pass

    return {
        "static_variables": static_variables,
        "timed_variables": timed_variables,
        "dynamic_records": dynamic_records,
    }
def write_plotly_gif(
    fig,
    iterations,
    gif_path,
    width=1600,
    height=850,
    frame_duration_ms=180,
):
    import copy
    import tempfile
    from pathlib import Path

    import imageio.v2 as imageio
    import plotly.graph_objects as go
    import plotly.io as pio

    gif_path = Path(gif_path)
    gif_path.parent.mkdir(parents=True, exist_ok=True)

    fps = max(1, int(round(1000 / frame_duration_ms)))

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        png_paths = []

        for i, frame in enumerate(fig.frames):
            frame_dict = frame.to_plotly_json()

            frame_fig = go.Figure(
                data=copy.deepcopy(fig.data),
                layout=copy.deepcopy(fig.layout),
            )

            for trace_i, trace_data in zip(frame_dict.get("traces", [0]), frame_dict["data"]):
                frame_fig.data[trace_i].update(trace_data)

            frame_fig.update_layout(frame_dict.get("layout", {}))

            png_path = tmp / f"frame_{i:05d}.png"
            pio.write_image(
                frame_fig,
                png_path,
                width=width,
                height=height,
                scale=2,
            )

            png_paths.append(png_path)
            print(f"rendered Plotly frame {i + 1}/{len(fig.frames)}")

        with imageio.get_writer(gif_path, mode="I", fps=fps) as writer:
            for png_path in png_paths:
                writer.append_data(imageio.imread(png_path))

    print(f"Plotly GIF written to: {gif_path}")
# %%

