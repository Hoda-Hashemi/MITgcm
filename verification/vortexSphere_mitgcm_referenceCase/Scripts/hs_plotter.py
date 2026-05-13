#%%
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

Nx, Ny = 1440, 720

lon = np.linspace(0.0, 360.0, Nx, endpoint=False)
lat = np.linspace(-90.0 + 0.125, 90.0 - 0.125, Ny)
Lon, Lat = np.meshgrid(lon, lat)

h_s = np.fromfile("../input/eta_init.bin", dtype=">f4").reshape(Ny, Nx)

skip = 4
Lon_p = Lon[::skip, ::skip]
Lat_p = Lat[::skip, ::skip]
h_s_p = h_s[::skip, ::skip]

cmax = np.nanmax(np.abs(h_s_p))
cmin = - cmax

fig = go.Figure(
    data=[
        go.Surface(
            x=Lon_p,
            y=Lat_p,
            z=h_s_p,
            surfacecolor=h_s_p,
            colorscale="RdBu_r",
            cmin=cmin,
            cmax=cmax,
            colorbar=dict(title="h_s [m]")
        )
    ]
)

fig.update_layout(
    title="Initial Free-Surface Perturbation h_s",
    width=1200,
    height=800,
    scene=dict(
        xaxis_title="Longitude [deg]",
        yaxis_title="Latitude [deg]",
        zaxis_title="h_s [m]",
        aspectmode="auto"
    )
)

zmax = np.max(h_s_p)

fig.update_layout(
    scene=dict(
        zaxis=dict(
            range=[-zmax, zmax]
        )
    )
)
fig.show()

# %%
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

# Grid

Nx = 1440
Ny = 720

lon = np.linspace(0.0, 360.0, Nx, endpoint=False)

lat = np.linspace(
    -90.0 + 0.125,
     90.0 - 0.125,
    Ny
)

Lon, Lat = np.meshgrid(lon, lat)

# Iterations

iters = [0, 12, 24]

# Read fields

fields = []

skip = 4

for it in iters:

    file_name = f"../run/Eta.{it:010d}.data"

    h_s = np.fromfile(
        file_name,
        dtype=">f8"
    ).reshape(Ny, Nx)

    h_s = h_s[::skip, ::skip]

    fields.append(h_s)

Lon_p = Lon[::skip, ::skip]
Lat_p = Lat[::skip, ::skip]

# Symmetric scale

zmax = max(
    np.nanmax(np.abs(f))
    for f in fields
)

# Initial figure

fig = go.Figure(

    data=[

        go.Surface(

            x=Lon_p,
            y=Lat_p,
            z=fields[0],

            surfacecolor=fields[0],

            colorscale="RdBu_r",

            cmin=-zmax,
            cmax= zmax,

            colorbar=dict(
                title="h_s [m]"
            )
        )
    ]
)

# Animation frames

frames = []

for i, it in enumerate(iters):

    frame = go.Frame(

        data=[

            go.Surface(

                x=Lon_p,
                y=Lat_p,
                z=fields[i],

                surfacecolor=fields[i],

                colorscale="RdBu_r",

                cmin=-zmax,
                cmax= zmax
            )
        ],

        name=str(it)
    )

    frames.append(frame)

fig.frames = frames

# Slider

steps = []

for it in iters:

    step = dict(

        method="animate",

        args=[
            [str(it)],
            dict(
                mode="immediate",
                frame=dict(duration=500, redraw=True),
                transition=dict(duration=0)
            )
        ],

        label=f"Iter {it}"
    )

    steps.append(step)

sliders = [dict(

    active=0,

    steps=steps,

    x=0.1,
    y=0,

    len=0.8
)]

# Layout

fig.update_layout(

    title="MITgcm Free-Surface Evolution",

    width=1200,
    height=800,

    sliders=sliders,

    scene=dict(

        xaxis_title="Longitude [deg]",
        yaxis_title="Latitude [deg]",
        zaxis_title="h_s [m]",

        zaxis=dict(
            range=[-zmax, zmax]
        ),

        aspectmode="auto"
    )
)
fig.write_html("mitgcm_surface.html")
fig.show()

# %%

import glob
import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

Nx = 1440
Ny = 720
dx = 0.25
dy = 0.25

run_dir = "../run_34_MPI"

lon = (np.arange(Nx) + 0.5)*dx
lat = 90.0 - (np.arange(Ny) + 0.5)*dy

Lon, Lat = np.meshgrid(lon, lat)

files = sorted(glob.glob(os.path.join(run_dir, "Eta.*.data")))

iters = [
    int(os.path.basename(f).split(".")[1])
    for f in files
]

print("Available Eta iterations:")
print(iters)
#%%
skip = 4
fields = []

for it in iters:
    file_name = f"../run/Eta.{it:010d}.data"

    h_s = np.fromfile(
        file_name,
        dtype=">f8"
    ).reshape(Ny, Nx)

    fields.append(h_s[::skip, ::skip])

Lon_p = Lon[::skip, ::skip]
Lat_p = Lat[::skip, ::skip]

zmax = max(np.nanmax(np.abs(f)) for f in fields)

fig = go.Figure(
    data=[
        go.Surface(
            x=Lon_p,
            y=Lat_p,
            z=fields[0],
            surfacecolor=fields[0],
            colorscale="RdBu_r",
            cmin=-zmax,
            cmax=zmax,
            colorbar=dict(title="h_s [m]")
        )
    ]
)

frames = []

for i, it in enumerate(iters):
    frames.append(
        go.Frame(
            data=[
                go.Surface(
                    x=Lon_p,
                    y=Lat_p,
                    z=fields[i],
                    surfacecolor=fields[i],
                    colorscale="RdBu_r",
                    cmin=-zmax,
                    cmax=zmax
                )
            ],
            name=str(it)
        )
    )

fig.frames = frames

steps = []

for it in iters:
    steps.append(
        dict(
            method="animate",
            args=[
                [str(it)],
                dict(
                    mode="immediate",
                    frame=dict(duration=500, redraw=True),
                    transition=dict(duration=0)
                )
            ],
            label=f"Iter {it}"
        )
    )

fig.update_layout(
    title="MITgcm Free-Surface Evolution on Latitude--Longitude Grid",
    width=1200,
    height=800,
    sliders=[
        dict(
            active=0,
            steps=steps,
            x=0.1,
            y=0,
            len=0.8
        )
    ],
    scene=dict(
        xaxis_title="Longitude λ [deg]",
        yaxis_title="Latitude θ [deg]",
        zaxis_title="h_s [m]",
        zaxis=dict(range=[-zmax, zmax]),
        aspectmode="auto"
    )
)
fig.add_annotation(
    text=
    "c = sqrt(gH) ≈ 328 m/s<br>"
    "T_circ ≈ 34 hours<br>"
    "nIter = 136<br>"
    "Δt = 900 s<br>"
    "T_run = 136 x 900 = 122400 s ≈ 34 hours<br>",

    x=0.5,
    y=1.08,

    xref="paper",
    yref="paper",

    showarrow=False,
    font=dict(size=16)
)
fig.write_html("mitgcm_surface_34hrs.html")
fig.show()

# %%

