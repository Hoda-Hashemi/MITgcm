# %%
"""Compare the passive tracer from two TC1 run folders.

Edit RUN_1_NAME and RUN_2_NAME, then press Shift+Enter on this cell.
The figures are displayed in the IPython kernel and also saved as PDF and PNG.
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_DIR = SCRIPT_DIR.parents[1]

# ================= EDIT ONLY THE RUN DIRECTORY NAMES =================
RUN_1_NAME = "run_alpha_0"
RUN_2_NAME = "run_alpha_1.52"
# =====================================================================

RUN_1 = CASE_DIR / RUN_1_NAME
RUN_2 = CASE_DIR / RUN_2_NAME

FIELD = "S"
SAVE_FIGURES = True

def read_gendata_value(run_folder, variable_name):
    gendata_files = sorted(run_folder.glob("gendata*.py"))
    if not gendata_files:
        raise FileNotFoundError(f"No gendata Python file found in {run_folder}")

    text = gendata_files[0].read_text()
    match = re.search(
        rf"^\s*{re.escape(variable_name)}\s*=\s*([^#\n]+)",
        text,
        re.MULTILINE,
    )
    if match is None:
        raise ValueError(
            f"{variable_name} was not found in {gendata_files[0].name}"
        )

    expression = match.group(1).strip()
    return float(
        eval(expression, {"__builtins__": {}}, {"np": np})
    )

def read_data_value(run_folder, variable_name):
    data_path = run_folder / "data"
    text = data_path.read_text()
    match = re.search(
        rf"\b{re.escape(variable_name)}\s*=\s*([+\-0-9.eEdD]+)",
        text,
        re.IGNORECASE,
    )
    if match is None:
        raise ValueError(f"{variable_name} was not found in {data_path}")
    return float(match.group(1).replace("D", "e").replace("d", "e"))

def read_meta(meta_path):
    text = meta_path.read_text()

    n_dims = int(
        re.search(r"nDims\s*=\s*\[\s*(\d+)\s*\]", text).group(1)
    )
    dim_text = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S).group(1)
    dim_values = [int(value) for value in re.findall(r"\d+", dim_text)]
    dimensions = [dim_values[3 * index] for index in range(n_dims)]

    precision_match = re.search(
        r"dataprec\s*=\s*\[\s*'([^']+)'", text
    )
    precision = precision_match.group(1).lower()
    dtype = ">f8" if precision in ("float64", "real*8") else ">f4"

    return dimensions, dtype

def read_field(run_folder, field_name, iteration=None):
    if iteration is None:
        file_stem = run_folder / field_name
    else:
        file_stem = run_folder / f"{field_name}.{iteration:010d}"

    meta_path = Path(f"{file_stem}.meta")
    data_path = Path(f"{file_stem}.data")
    dimensions, dtype = read_meta(meta_path)

    values = np.fromfile(data_path, dtype=dtype)
    nx, ny = dimensions[:2]
    return values.reshape(ny, nx).astype(float)

def available_iterations(run_folder, field_name):
    pattern = re.compile(rf"{re.escape(field_name)}\.(\d{{10}})\.meta")
    return sorted(
        int(pattern.fullmatch(path.name).group(1))
        for path in run_folder.glob(f"{field_name}.*.meta")
        if pattern.fullmatch(path.name)
    )

def matching_outputs(run_1, run_2, field_name, delta_t_1, delta_t_2):
    times_1 = {
        round(iteration * delta_t_1, 6): iteration
        for iteration in available_iterations(run_1, field_name)
    }
    times_2 = {
        round(iteration * delta_t_2, 6): iteration
        for iteration in available_iterations(run_2, field_name)
    }

    common_times = sorted(set(times_1) & set(times_2))
    return [
        (time_seconds, times_1[time_seconds], times_2[time_seconds])
        for time_seconds in common_times
    ]

def tracer_center(xc, yc, tracer):
    weights = np.maximum(tracer, 0.0) * np.cos(np.deg2rad(yc))

    lon = np.deg2rad(xc)
    lat = np.deg2rad(yc)
    x = np.sum(weights * np.cos(lat) * np.cos(lon))
    y = np.sum(weights * np.cos(lat) * np.sin(lon))
    z = np.sum(weights * np.sin(lat))

    center_lon = np.rad2deg(np.arctan2(y, x)) % 360.0
    center_lat = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
    return center_lon, center_lat

def center_vector(center):
    lon = np.deg2rad(center[0])
    lat = np.deg2rad(center[1])
    return np.array(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )

def comparison_center(center_1, center_2):
    vector = center_vector(center_1) + center_vector(center_2)
    vector /= np.linalg.norm(vector)

    lon = np.rad2deg(np.arctan2(vector[1], vector[0])) % 360.0
    lat = np.rad2deg(np.arctan2(vector[2], np.hypot(vector[0], vector[1])))
    return lon, lat

def center_separation(center_1, center_2):
    dot_product = np.dot(center_vector(center_1), center_vector(center_2))
    return np.rad2deg(np.arccos(np.clip(dot_product, -1.0, 1.0)))

def local_fields(xc, yc, tracer_1, tracer_2, center):
    lon_offset = ((xc[0] - center[0] + 180.0) % 360.0) - 180.0
    lat_offset = yc[:, 0] - center[1]

    lon_order = np.argsort(lon_offset)
    lat_order = np.argsort(lat_offset)
    lon_offset = lon_offset[lon_order]
    lat_offset = lat_offset[lat_order]
    tracer_1 = tracer_1[np.ix_(lat_order, lon_order)]
    tracer_2 = tracer_2[np.ix_(lat_order, lon_order)]

    support = (tracer_1 > 0.01) | (tracer_2 > 0.01)
    rows, columns = np.where(support)
    half_width = max(
        np.max(np.abs(lon_offset[columns])),
        np.max(np.abs(lat_offset[rows])),
    ) + 4.0

    lon_keep = np.abs(lon_offset) <= half_width
    lat_keep = np.abs(lat_offset) <= half_width
    lon, lat = np.meshgrid(lon_offset[lon_keep], lat_offset[lat_keep])

    return (
        lon,
        lat,
        tracer_1[np.ix_(lat_keep, lon_keep)],
        tracer_2[np.ix_(lat_keep, lon_keep)],
    )

def plot_comparison(
    lon,
    lat,
    tracer_1,
    tracer_2,
    day,
    separation,
    label_1,
    label_2,
    contour_levels,
):
    fig, ax = plt.subplots(figsize=(6.2, 6.2))

    levels = contour_levels[
        contour_levels <= min(tracer_1.max(), tracer_2.max())
    ]

    contour_1 = ax.contour(
        lon,
        lat,
        tracer_1,
        levels=levels,
        colors="#0072B2",
        linewidths=1.0,
        linestyles="-",
    )
    ax.contour(
        lon,
        lat,
        tracer_2,
        levels=levels,
        colors="#D55E00",
        linewidths=1.0,
        linestyles="--",
    )
    ax.clabel(contour_1, levels=levels[1::2], fmt="%g", fontsize=8)

    ax.legend(
        handles=[
            Line2D([0], [0], color="#0072B2", lw=1.5, label=label_1),
            Line2D(
                [0],
                [0],
                color="#D55E00",
                lw=1.5,
                linestyle="--",
                label=label_2,
            ),
        ],
        frameon=False,
        loc="upper right",
    )

    ax.set_title(
        f"TC1 passive-tracer overlay - day {day:g}\n"
        f"tracer-center separation = {separation:.3f} degrees"
    )
    ax.set_xlabel("Longitude offset from comparison center [deg]")
    ax.set_ylabel("Latitude offset from comparison center [deg]")
    ax.set_aspect("equal")
    ax.tick_params(top=True, right=True)
    fig.tight_layout()

    return fig

RUN_1 = RUN_1.expanduser().resolve()
RUN_2 = RUN_2.expanduser().resolve()

if not RUN_1.is_dir():
    raise FileNotFoundError(f"Run folder not found: {RUN_1}")
if not RUN_2.is_dir():
    raise FileNotFoundError(f"Run folder not found: {RUN_2}")

alpha_1 = read_gendata_value(RUN_1, "ALPHA_RAD")
alpha_2 = read_gendata_value(RUN_2, "ALPHA_RAD")
height_1 = read_gendata_value(RUN_1, "H0")
height_2 = read_gendata_value(RUN_2, "H0")
delta_t_1 = read_data_value(RUN_1, "deltaT")
delta_t_2 = read_data_value(RUN_2, "deltaT")

label_1 = (
    f"{RUN_1.name}: alpha = {alpha_1:g} rad "
    f"({np.rad2deg(alpha_1):.3f} deg)"
)
label_2 = (
    f"{RUN_2.name}: alpha = {alpha_2:g} rad "
    f"({np.rad2deg(alpha_2):.3f} deg)"
)

height_for_levels = min(height_1, height_2)
contour_levels = height_for_levels * np.arange(0.08, 0.81, 0.08)

outputs = matching_outputs(
    RUN_1, RUN_2, FIELD, delta_t_1, delta_t_2
)
if not outputs:
    raise FileNotFoundError(
        "The two folders have no tracer outputs at matching physical times."
    )

xc_1 = read_field(RUN_1, "XC")
yc_1 = read_field(RUN_1, "YC")
xc_2 = read_field(RUN_2, "XC")
yc_2 = read_field(RUN_2, "YC")

if not np.allclose(xc_1, xc_2) or not np.allclose(yc_1, yc_2):
    raise ValueError("The two runs do not use the same grid.")

output_folder = (
    RUN_1.parent
    / "Scripts"
    / "output"
    / f"compare_{RUN_1.name}_vs_{RUN_2.name}"
    / "PassiveTracerOverlay"
)
if SAVE_FIGURES:
    output_folder.mkdir(parents=True, exist_ok=True)

for time_seconds, iteration_1, iteration_2 in outputs:
    tracer_1 = read_field(RUN_1, FIELD, iteration_1)
    tracer_2 = read_field(RUN_2, FIELD, iteration_2)

    center_1 = tracer_center(xc_1, yc_1, tracer_1)
    center_2 = tracer_center(xc_1, yc_1, tracer_2)
    center = comparison_center(center_1, center_2)
    separation = center_separation(center_1, center_2)

    lon, lat, local_1, local_2 = local_fields(
        xc_1, yc_1, tracer_1, tracer_2, center
    )
    day = time_seconds / 86400.0
    figure = plot_comparison(
        lon,
        lat,
        local_1,
        local_2,
        day,
        separation,
        label_1,
        label_2,
        contour_levels,
    )

    if SAVE_FIGURES:
        name = (
            f"tc1_tracer_overlay_day_{day:05.2f}_"
            f"iter1_{iteration_1:010d}_iter2_{iteration_2:010d}"
        )
        figure.savefig(output_folder / f"{name}.pdf", bbox_inches="tight")
        figure.savefig(
            output_folder / f"{name}.png",
            dpi=300,
            bbox_inches="tight",
        )

    plt.show()

print(f"Run 1: alpha={alpha_1:g} rad, H0={height_1:g}, deltaT={delta_t_1:g} s")
print(f"Run 2: alpha={alpha_2:g} rad, H0={height_2:g}, deltaT={delta_t_2:g} s")
print(f"Compared {len(outputs)} outputs at matching physical times.")
if SAVE_FIGURES:
    print(f"Figures saved in: {output_folder}")

# %%

