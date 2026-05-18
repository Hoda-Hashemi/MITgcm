#%%
# Cell 1: Packages and parameters.
from pathlib import Path
import glob
import os

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

# Main knobs.
SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
RUN_DIR = (SCRIPT_DIR / "../run").resolve()
DYN_SUBDIR = "dyn"
FIELD_NAME = "Eta"
FIELD_KIND = "prognostic"  # "prognostic" or "diagnostic"
ASK_FIELD_KIND = True

# Grid and binary format.
NX = 1440
NY = 720
DX_DEG = 0.25
DY_DEG = 0.25
DT_SECONDS = 900
DATA_DTYPE = ">f8"
SKIP = 4
LATITUDE_ORDER = "north_to_south"  # "north_to_south" or "south_to_north"

# Plot controls.
CMAP = "RdBu_r"  # good alternatives: "RdBu_r", "coolwarm"
COLORBAR_LABEL = ""  # blank means auto-fill from FIELD_NAME
FIELD_UNITS = {"Eta": "m"}
FIELD_DESCRIPTIONS = {"Eta": "Free surface"}
WRITE_HTML = True
OPEN_BROWSER = True
HTML_FILE = "mitgcm_surface.html"
FIG_WIDTH = 1500
FIG_HEIGHT = 850

# Layout lanes: left is the 3D surface, right is colorbar + information box.
SURFACE_DOMAIN_X = [0.0, 0.72]
COLORBAR_X = 0.78
INFO_BOX_X = 0.86

# Animation/video controls.
MAKE_VIDEO = False
VIDEO_FILE = "mitgcm_surface.mp4"
VIDEO_FRAME_DURATION_MS = 180


#%%
# Cell 2: Helpers.
def output_dir_for_field(run_dir, field_kind, dyn_subdir):
    if field_kind.lower() == "diagnostic":
        return run_dir / dyn_subdir
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

    return field.reshape(NY, NX)


def grid():
    lon = (np.arange(NX) + 0.5) * DX_DEG

    if LATITUDE_ORDER == "north_to_south":
        lat = 90.0 - (np.arange(NY) + 0.5) * DY_DEG
    elif LATITUDE_ORDER == "south_to_north":
        lat = -90.0 + (np.arange(NY) + 0.5) * DY_DEG
    else:
        raise ValueError('LATITUDE_ORDER must be "north_to_south" or "south_to_north".')

    return np.meshgrid(lon, lat)


def stats_for_field(field):
    return {
        "min": float(np.nanmin(field)),
        "max": float(np.nanmax(field)),
        "mean": float(np.nanmean(field)),
        "std": float(np.nanstd(field)),
    }


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
        f"grid: {NY} x {NX}, skip={SKIP}<br>"
        f"color scale: +/- {value_text(zmax, unit)}<br><br>"
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
        y=0.88,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        align="left",
        showarrow=False,
        bordercolor="rgba(30, 30, 30, 0.35)",
        borderwidth=1,
        borderpad=10,
        bgcolor="rgba(255, 255, 255, 0.92)",
        font=dict(size=13),
    )


def surface_trace(field, show_colorbar=False):
    return go.Surface(
        x=Lon_p,
        y=Lat_p,
        z=field,
        surfacecolor=field,
        colorscale=colorscale,
        cmin=-zmax,
        cmax=zmax,
        showscale=show_colorbar,
        colorbar=dict(
            title=colorbar_label,
            x=COLORBAR_X,
            xanchor="left",
            y=0.56,
            len=0.74,
            thickness=18,
        ),
        contours=dict(
            z=dict(show=True, usecolormap=True, highlightcolor="white", project_z=True)
        ),
        hovertemplate=(
            "lon: %{x:.2f} deg<br>"
            "lat: %{y:.2f} deg<br>"
            f"{colorbar_label}: "
            "%{z:.6g}<extra></extra>"
        ),
    )


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

Lon, Lat = grid()
Lon_p = Lon[::SKIP, ::SKIP]
Lat_p = Lat[::SKIP, ::SKIP]

fields = [
    read_surface(data_dir, FIELD_NAME, iteration)[::SKIP, ::SKIP]
    for iteration in iterations
]

zmax = max(np.nanmax(np.abs(field)) for field in fields)
all_stats = stats_for_field(np.asarray(fields))
colorbar_label = COLORBAR_LABEL or label_with_unit(FIELD_NAME)
title = f"MITgcm {title_name(FIELD_NAME)} Evolution ({FIELD_KIND}, {len(iterations)} iterations)"
colorscale = plotly_colorscale(CMAP)


#%%
# Cell 4: Surface plotter with browser HTML output.
# This is a real 3D surface: x=longitude, y=latitude, z=field value.
fig = go.Figure(
    data=[surface_trace(fields[0], show_colorbar=True)]
)

fig.frames = [
    go.Frame(
        data=[surface_trace(field, show_colorbar=True)],
        layout=dict(
            annotations=[
                info_annotation(
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
            ]
        ),
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
    margin=dict(l=20, r=40, t=80, b=45),
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
    scene=dict(
        domain=dict(x=SURFACE_DOMAIN_X, y=[0.08, 1.0]),
        xaxis_title="Longitude [deg]",
        yaxis_title="Latitude [deg]",
        zaxis_title=colorbar_label,
        zaxis=dict(range=[-zmax, zmax]),
        camera=dict(eye=dict(x=1.45, y=-1.55, z=0.85)),
        aspectmode="auto",
    ),
    annotations=[
        info_annotation(
            stats_box_text(
                FIELD_NAME,
                FIELD_KIND,
                data_dir,
                iterations,
                iterations[0],
                fields[0],
                all_stats,
                zmax,
            )
        )
    ],
)

if WRITE_HTML:
    fig.write_html(HTML_FILE, auto_open=OPEN_BROWSER)

fig.show()


#%%
# Cell 5: Optional video export.
# Requires kaleido and imageio:
#   pip install -U kaleido imageio
if MAKE_VIDEO:
    import imageio.v2 as imageio

    png_frames = []

    for iteration, field in zip(iterations, fields):
        fig.data[0].z = field
        fig.data[0].surfacecolor = field
        fig.update_layout(title=dict(text=f"{title} | iteration {iteration}"))

        png_name = f"_frame_{FIELD_NAME}_{iteration:010d}.png"
        fig.write_image(png_name, width=FIG_WIDTH, height=FIG_HEIGHT, scale=1)
        png_frames.append(png_name)

    with imageio.get_writer(VIDEO_FILE, fps=max(1, int(1000 / VIDEO_FRAME_DURATION_MS))) as writer:
        for png_name in png_frames:
            writer.append_data(imageio.imread(png_name))

    print(f"Video written to: {VIDEO_FILE}")
