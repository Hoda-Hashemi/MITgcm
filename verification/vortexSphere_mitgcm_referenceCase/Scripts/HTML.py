#%%
# Cell 1: Packages and parameters.
from pathlib import Path
import glob
import importlib.util
import os
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
        if (candidate / "HTML.py").exists() and (candidate.parent / "run").exists():
            return candidate

    return cwd


# Main knobs.
SCRIPT_DIR = resolve_script_dir()
RUN_DIR = (SCRIPT_DIR / "../run").resolve()
DYN_SUBDIR = "dyn"
FIELD_NAME = "Eta"
FIELD_KIND = "prognostic"  # "prognostic" or "diagnostic"
ASK_FIELD_KIND = False

# Grid and binary format.
NX = 1440
NY = 720
DX_DEG = 0.25
DY_DEG = 0.25
DT_SECONDS = 60
DATA_DTYPE = ">f8"
LATITUDE_ORDER = "south_to_north"  # MITgcm grid in this case starts at ygOrigin=-90.

# Plot controls.
CMAP = "RdBu_r"
COLORBAR_LABEL = ""  # blank means auto-fill from FIELD_NAME
FIELD_UNITS = {"Eta": "m", "ETAN": "m"}
FIELD_DESCRIPTIONS = {"Eta": "Free surface", "ETAN": "Free surface"}
DEFAULT_COLOR_LIMIT = "0.12" if FIELD_NAME in {"Eta", "ETAN"} else ""
COLOR_LIMIT_TEXT = os.environ.get("MITGCM_COLOR_LIMIT", DEFAULT_COLOR_LIMIT).strip()
COLOR_LIMIT = float(COLOR_LIMIT_TEXT) if COLOR_LIMIT_TEXT else None
WRITE_HTML = True
OPEN_BROWSER = env_flag("MITGCM_OPEN_BROWSER", False)
SHOW_FIGURE = env_flag("MITGCM_SHOW_FIGURE", True)
HTML_FILE = "mitgcm_surface.html"
FIG_WIDTH = 1200
FIG_HEIGHT = 720

# Browser responsiveness controls. All iterations are kept; the spatial grid is
# thinned only for the interactive HTML. Lower this for more detail if needed.
HTML_SKIP = int(os.environ.get("MITGCM_HTML_SKIP", "6"))
VIDEO_SKIP = int(os.environ.get("MITGCM_VIDEO_SKIP", "4"))

# Layout lanes: left is the map, right is colorbar + information box.
MAP_DOMAIN_X = [0.0, 0.76]
COLORBAR_X = 0.80
INFO_BOX_X = 0.86

# Animation/video controls. MP4 needs ffmpeg; without it this writes a GIF.
MAKE_VIDEO = env_flag("MITGCM_MAKE_VIDEO", True)
VIDEO_FILE = "mitgcm_surface.mp4"
VIDEO_FALLBACK_FILE = "mitgcm_surface.gif"
VIDEO_FRAME_DURATION_MS = 180
VIDEO_WIDTH = 1200
VIDEO_HEIGHT = 620

# Fancy 3-D PyVista controls.
EARTH_RADIUS_M = 6_371_000.0
PYVISTA_SKIP = int(os.environ.get("MITGCM_PYVISTA_SKIP", "4"))
PYVISTA_FRAME_STEP = int(os.environ.get("MITGCM_PYVISTA_FRAME_STEP", "3"))
PYVISTA_VERTICAL_SCALE = float(os.environ.get("MITGCM_PYVISTA_VERTICAL_SCALE", "250000.0"))
PYVISTA_LAND_OFFSET_M = float(os.environ.get("MITGCM_PYVISTA_LAND_OFFSET_M", "12000.0"))
PYVISTA_SHOW = env_flag("MITGCM_PYVISTA_SHOW", False)
PYVISTA_OFF_SCREEN = env_flag("MITGCM_PYVISTA_OFF_SCREEN", True)
PYVISTA_MAKE_GIF = env_flag("MITGCM_PYVISTA_MAKE_GIF", True)
PYVISTA_GIF_FILE = "mitgcm_surface_pyvista.gif"
PYVISTA_SCREENSHOT_FILE = "mitgcm_surface_pyvista.png"
PYVISTA_WINDOW_WIDTH = int(os.environ.get("MITGCM_PYVISTA_WIDTH", "1600"))
PYVISTA_WINDOW_HEIGHT = int(os.environ.get("MITGCM_PYVISTA_HEIGHT", "1000"))

# Let VS Code's interactive window render inline when it is available.
if os.environ.get("PLOTLY_RENDERER"):
    pio.renderers.default = os.environ["PLOTLY_RENDERER"]
elif os.environ.get("VSCODE_PID"):
    pio.renderers.default = "vscode"



#%%
# Cell 2: Helpers.
def output_path(file_name):
    path = Path(file_name)
    if path.is_absolute():
        return path
    return SCRIPT_DIR / path


def output_dir_for_field(run_dir, field_kind, dyn_subdir):
    if field_kind.lower() == "diagnostic":
        candidate = run_dir / dyn_subdir
        return candidate if candidate.exists() else run_dir
    if field_kind.lower() == "prognostic":
        return run_dir
    raise ValueError('FIELD_KIND must be "prognostic" or "diagnostic".')


def plotly_colorscale(cmap):
    if cmap.lower() == "coolwarm":
        return [
            [0.0, "rgb(59,76,192)"],
            [0.25, "rgb(130,165,253)"],
            [0.5, "rgb(221,221,221)"],
            [0.75, "rgb(244,152,122)"],
            [1.0, "rgb(180,4,38)"],
        ]

    return cmap


def choose_field_kind(default_kind):
    try:
        answer = input(
            f"Is {FIELD_NAME} prognostic or diagnostic? "
            f"[prognostic/diagnostic, default={default_kind}]: "
        ).strip().lower()
    except EOFError:
        answer = ""

    if not answer:
        return default_kind
    if answer in {"p", "prog", "prognostic"}:
        return "prognostic"
    if answer in {"d", "diag", "diagnostic"}:
        return "diagnostic"

    raise ValueError('Please answer "prognostic" or "diagnostic".')


def discover_iterations(data_dir, field_name):
    pattern = str(data_dir / f"{field_name}.*.data")
    files = sorted(glob.glob(pattern))
    iterations = []

    for file_path in files:
        parts = os.path.basename(file_path).split(".")
        if len(parts) >= 3 and parts[0] == field_name:
            try:
                iterations.append(int(parts[1]))
            except ValueError:
                pass

    return sorted(iterations)


def field_file(data_dir, field_name, iteration):
    return data_dir / f"{field_name}.{iteration:010d}.data"


def read_surface(data_dir, field_name, iteration):
    path = field_file(data_dir, field_name, iteration)
    field = np.fromfile(path, dtype=DATA_DTYPE)
    expected_size = NY * NX

    if field.size != expected_size:
        raise ValueError(
            f"{path} has {field.size} values, expected {expected_size} "
            f"for shape ({NY}, {NX})."
        )

    return field.reshape(NY, NX).astype(np.float32, copy=False)


def read_land_mask():
    depth_path = RUN_DIR / "Depth.data"
    if not depth_path.exists():
        print(f"No Depth.data found at {depth_path}; land will not be masked.")
        return None

    depth = np.fromfile(depth_path, dtype=DATA_DTYPE)
    expected_size = NY * NX
    if depth.size != expected_size:
        print(
            f"{depth_path} has {depth.size} values, expected {expected_size}; "
            "land will not be masked."
        )
        return None

    return depth.reshape(NY, NX) <= 0.0


def mask_land(field, land_mask):
    if land_mask is None:
        return field

    field = field.copy()
    field[land_mask] = np.nan
    return field


def read_plot_field(data_dir, field_name, iteration, land_mask, skip):
    return mask_land(read_surface(data_dir, field_name, iteration), land_mask)[::skip, ::skip]


def read_raw_field(data_dir, field_name, iteration, skip):
    return read_surface(data_dir, field_name, iteration)[::skip, ::skip]


def grid_1d():
    lon = (np.arange(NX) + 0.5) * DX_DEG

    if LATITUDE_ORDER == "north_to_south":
        lat = 90.0 - (np.arange(NY) + 0.5) * DY_DEG
    elif LATITUDE_ORDER == "south_to_north":
        lat = -90.0 + (np.arange(NY) + 0.5) * DY_DEG
    else:
        raise ValueError('LATITUDE_ORDER must be "north_to_south" or "south_to_north".')

    return lon, lat


def stats_for_field(field):
    return {
        "min": float(np.nanmin(field)),
        "max": float(np.nanmax(field)),
        "mean": float(np.nanmean(field)),
        "std": float(np.nanstd(field)),
    }


def safe_symmetric_limit(fields):
    zmax = max(float(np.nanmax(np.abs(field))) for field in fields)
    if not np.isfinite(zmax) or zmax <= 0.0:
        return 1e-12
    return zmax


def display_symmetric_limit(fields):
    if COLOR_LIMIT is not None:
        limit = abs(float(COLOR_LIMIT))
        if np.isfinite(limit) and limit > 0.0:
            return limit

    return safe_symmetric_limit(fields)


def nice_number(value):
    return f"{value:.6g}"


def label_with_unit(field_name):
    unit = FIELD_UNITS.get(field_name, "")
    if unit:
        return f"{field_name} [{unit}]"

    return field_name


def title_name(field_name):
    description = FIELD_DESCRIPTIONS.get(field_name, "")
    if description:
        return f"{description} {field_name}"

    return field_name


def value_text(value, unit):
    if unit:
        return f"{nice_number(value)} {unit}"

    return nice_number(value)


def stats_box_text(
    field_name,
    field_kind,
    data_dir,
    iterations,
    iteration,
    field,
    all_stats,
    zmax,
):
    current_stats = stats_for_field(field)
    unit = FIELD_UNITS.get(field_name, "")
    runtime_hours = iteration * DT_SECONDS / 3600.0

    return (
        f"<b>{title_name(field_name)}</b><br>"
        f"kind: {field_kind}<br>"
        f"path: {data_dir}<br>"
        f"iterations: {len(iterations)} found<br>"
        f"first/last: {iterations[0]} / {iterations[-1]}<br>"
        f"current iteration: {iteration}<br>"
        f"current time: {runtime_hours:.2f} hours<br>"
        f"grid: {NY} x {NX}, HTML skip={HTML_SKIP}<br>"
        f"color scale: +/- {value_text(zmax, unit)}"
        f"{' (fixed display)' if COLOR_LIMIT is not None else ''}<br><br>"
        f"<b>Selected iteration</b><br>"
        f"min: {value_text(current_stats['min'], unit)}<br>"
        f"max: {value_text(current_stats['max'], unit)}<br>"
        f"mean: {value_text(current_stats['mean'], unit)}<br>"
        f"std: {value_text(current_stats['std'], unit)}<br><br>"
        f"<b>All loaded frames</b><br>"
        f"min: {value_text(all_stats['min'], unit)}<br>"
        f"max: {value_text(all_stats['max'], unit)}<br>"
        f"mean: {value_text(all_stats['mean'], unit)}"
    )


def info_annotation(text):
    return dict(
        text=text,
        x=INFO_BOX_X,
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


def eta_heatmap(field, show_colorbar=True):
    return go.Heatmap(
        x=Lon_p,
        y=Lat_p,
        z=field,
        zmin=-zmax,
        zmax=zmax,
        colorscale=colorscale,
        showscale=show_colorbar,
        colorbar=dict(
            title=colorbar_label,
            x=COLORBAR_X,
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


def land_heatmap(land_mask):
    if land_mask is None:
        return None

    land_z = np.where(land_mask[::HTML_SKIP, ::HTML_SKIP], 1.0, np.nan).astype(np.float32)
    return go.Heatmap(
        x=Lon_p,
        y=Lat_p,
        z=land_z,
        zmin=0.0,
        zmax=1.0,
        colorscale=[[0.0, "rgb(198,198,198)"], [1.0, "rgb(198,198,198)"]],
        showscale=False,
        hoverinfo="skip",
        name="land",
    )


def frame_annotation(iteration, field):
    return info_annotation(
        stats_box_text(
            FIELD_NAME,
            FIELD_KIND,
            data_dir,
            iterations,
            iteration,
            field,
            all_stats,
            zmax,
        )
    )


def show_interactive_figure(fig, html_path):
    if WRITE_HTML:
        fig.write_html(
            html_path,
            include_plotlyjs=True,
            auto_open=OPEN_BROWSER,
            auto_play=False,
            config={"responsive": True},
        )
        print(f"Interactive HTML written to: {html_path}")

    if not SHOW_FIGURE:
        return

    try:
        fig.show()
    except Exception as exc:
        print(f"Plotly display was skipped because fig.show() failed: {exc}")
        print(f"Open the written HTML file instead: {html_path}")


def ffmpeg_available():
    return shutil.which("ffmpeg") is not None or importlib.util.find_spec("imageio_ffmpeg") is not None


def choose_video_output():
    requested = output_path(VIDEO_FILE)
    if requested.suffix.lower() == ".mp4" and not ffmpeg_available():
        fallback = output_path(VIDEO_FALLBACK_FILE)
        print(f"ffmpeg is not available, so writing GIF instead of MP4: {fallback}")
        return fallback

    return requested


def write_video(data_dir, field_name, iterations, land_mask, zmax):
    import imageio.v2 as imageio

    mpl_config_dir = Path(os.environ.get("MPLCONFIGDIR", "/tmp/matplotlib"))
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors

    video_path = choose_video_output()
    fps = max(1, int(round(1000 / VIDEO_FRAME_DURATION_MS)))
    first_field = read_plot_field(data_dir, field_name, iterations[0], land_mask, VIDEO_SKIP)

    cmap_obj = plt.get_cmap(CMAP).copy()
    cmap_obj.set_bad("0.78")
    norm = mcolors.TwoSlopeNorm(vmin=-zmax, vcenter=0.0, vmax=zmax)

    fig, ax = plt.subplots(
        figsize=(VIDEO_WIDTH / 100.0, VIDEO_HEIGHT / 100.0),
        dpi=100,
        constrained_layout=True,
    )
    image = ax.imshow(
        first_field,
        origin="lower" if LATITUDE_ORDER == "south_to_north" else "upper",
        extent=[0.0, 360.0, -90.0, 90.0],
        cmap=cmap_obj,
        norm=norm,
        interpolation="nearest",
        aspect="auto",
    )
    cbar = fig.colorbar(image, ax=ax, location="right", pad=0.02, fraction=0.04)
    cbar.set_label(colorbar_label)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    ax.set_yticks(np.arange(-90.0, 91.0, 30.0))

    writer_kwargs = {"mode": "I", "fps": fps}
    with imageio.get_writer(video_path, **writer_kwargs) as writer:
        for frame_number, iteration in enumerate(iterations, start=1):
            field = read_plot_field(data_dir, field_name, iteration, land_mask, VIDEO_SKIP)
            runtime_hours = iteration * DT_SECONDS / 3600.0

            image.set_data(field)
            ax.set_title(
                f"{title} | iteration {iteration} | time {runtime_hours:.2f} hours"
            )
            fig.canvas.draw()
            rgba = np.asarray(fig.canvas.buffer_rgba())
            writer.append_data(rgba[:, :, :3].copy())

            if frame_number == 1 or frame_number % 25 == 0 or frame_number == len(iterations):
                print(f"  rendered video frame {frame_number}/{len(iterations)}")

    plt.close(fig)
    print(f"Animation written to: {video_path}")


#%%
# Cell 3: Discover available files and read the selected variable.
if ASK_FIELD_KIND:
    FIELD_KIND = choose_field_kind(FIELD_KIND)

data_dir = output_dir_for_field(RUN_DIR, FIELD_KIND, DYN_SUBDIR)
iterations = discover_iterations(data_dir, FIELD_NAME)

print(f"Looking for {FIELD_NAME} as {FIELD_KIND} output in: {data_dir}")
print("Available iterations found:")
print(iterations)

if not iterations:
    raise FileNotFoundError(f"No files found matching {data_dir / (FIELD_NAME + '.*.data')}")

land_mask = read_land_mask()
Lon, Lat = grid_1d()
Lon_p = Lon[::HTML_SKIP]
Lat_p = Lat[::HTML_SKIP]

fields = [
    read_plot_field(data_dir, FIELD_NAME, iteration, land_mask, HTML_SKIP)
    for iteration in iterations
]

zmax = display_symmetric_limit(fields)
all_stats = stats_for_field(np.asarray(fields))
colorbar_label = COLORBAR_LABEL or label_with_unit(FIELD_NAME)
title = f"MITgcm {title_name(FIELD_NAME)} Evolution ({FIELD_KIND}, {len(iterations)} iterations)"
colorscale = plotly_colorscale(CMAP)


#%%
# Cell 4: Interactive HTML output.
# This view keeps all iterations and uses a lighter heatmap so the browser can
# render the full free-surface evolution instead of opening a blank 3D page.
data_traces = [eta_heatmap(fields[0], show_colorbar=True)]
land_trace = land_heatmap(land_mask)
if land_trace is not None:
    data_traces.append(land_trace)

fig = go.Figure(data=data_traces)

fig.frames = [
    go.Frame(
        data=[go.Heatmap(z=field)],
        traces=[0],
        layout=dict(annotations=[frame_annotation(iteration, field)]),
        name=str(iteration),
    )
    for field, iteration in zip(fields, iterations)
]

slider_steps = [
    dict(
        method="animate",
        args=[
            [str(iteration)],
            dict(
                mode="immediate",
                frame=dict(duration=350, redraw=True),
                transition=dict(duration=0),
            ),
        ],
        label=str(iteration),
    )
    for iteration in iterations
]

play_buttons = [
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
                        frame=dict(duration=VIDEO_FRAME_DURATION_MS, redraw=True),
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
                    dict(
                        frame=dict(duration=0, redraw=False),
                        mode="immediate",
                        transition=dict(duration=0),
                    ),
                ],
            ),
        ],
    )
]

fig.update_layout(
    title=dict(text=title, x=0.36, xanchor="center"),
    width=FIG_WIDTH,
    height=FIG_HEIGHT,
    template="plotly_white",
    margin=dict(l=40, r=35, t=75, b=50),
    sliders=[
        dict(
            active=0,
            steps=slider_steps,
            x=0.12,
            y=0,
            len=0.60,
            currentvalue=dict(prefix="Iteration: "),
        )
    ],
    updatemenus=play_buttons,
    xaxis=dict(
        domain=MAP_DOMAIN_X,
        title="Longitude [deg]",
        range=[0.0, 360.0],
        constrain="domain",
    ),
    yaxis=dict(
        domain=[0.08, 1.0],
        title="Latitude [deg]",
        range=[-90.0, 90.0],
        scaleanchor="x",
        scaleratio=1.0,
    ),
    annotations=[
        frame_annotation(iterations[0], fields[0])
    ],
)

html_path = output_path(HTML_FILE)
show_interactive_figure(fig, html_path)


#%%
# Cell 5: Optional animation export.
# This path renders directly with Matplotlib/imageio and does not use Kaleido.
if MAKE_VIDEO:
    write_video(data_dir, FIELD_NAME, iterations, land_mask, zmax)

