#!/usr/bin/env python3
"""Generate the GitHub Pages homepage under docs/index.html."""
import html
import csv
import math
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

REPO_DIR = Path(__file__).resolve().parents[2]
SANDBOX_DIR = REPO_DIR / "Sandbox"
DOCS_DIR = REPO_DIR / "docs"
ASSET_ROOT = DOCS_DIR / "assets" / "williamson"
FRAGMENT_DIR = DOCS_DIR / "fragments"
SANDBOX_OUTPUT_ROOT = SANDBOX_DIR / "output"
HTML_TITLE = "MITgcm Verification"
SITE_TITLE = "WILLIAMSON TESTS"
SITE_SUBTITLE = (
    "MITgcm shallow-water verification pages for the Williamson standard spherical test suite."
)
AUTHOR_NAME = "Hoda Hashemi"
REPO_MAIN_URL = "https://github.com/Hoda-Hashemi/MITgcm/tree/main"
FAVICON_HREF = (
    "data:image/svg+xml,"
    "%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'%3E"
    "%3Crect width='100' height='100' fill='black'/%3E"
    "%3Ctext x='50' y='68' text-anchor='middle' font-family='Arial, sans-serif' font-size='58' fill='white'%3EM%3C/text%3E"
    "%3C/svg%3E"
)

SECTIONS = [
    {
        "slug": "case1_constant_bathymetry",
        "title": "Gaussian Patch Free Surface Initialization with Constant Bathymetry",
        "media": [
            {
                "label": "Prognostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Constant_Bathymetry_d_-4000_m" / "prognostic.webm",
            },
            {
                "label": "Diagnostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Constant_Bathymetry_d_-4000_m" / "diagnostic.webm",
            },
        ],
    },
    {
        "slug": "case2_real_bathymetry",
        "title": "Gaussian Patch Free Surface Initialization with Real Bathymetry",
        "media": [
            {
                "label": "Prognostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Real_Bathymetry" / "prognostic.webm",
            },
            {
                "label": "Diagnostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Real_Bathymetry" / "diagnostic.webm",
            },
        ],
    },
    {
        "slug": "case3_geostrophic_adjustment",
        "title": "Geostrophic Free Surface Initialization with Real Bathymetry",
        "media": [
            {
                "label": "Prognostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Geostrophic_Adjustment_with_Real_Bathymetry" / "prognostic.webm",
            },
            {
                "label": "Diagnostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Geostrophic_Adjustment_with_Real_Bathymetry" / "diagnostic.webm",
            },
        ],
    },
    {
        "slug": "testcase1",
        "title": "Williamson TC1: Advection of Cosine Bell",
        "summary": "",
        "cases": [
            {
                "label": "TC1 passive tracer",
                "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC1",
                "snapshot_fields": [
                    ("tracer", "Passive Tracer"),
                    ("etan", "ETAN"),
                    ("psi", "PsiVEL"),
                    ("phi", "PhiVEL"),
                    ("velocity_magnitude", "Velocity Magnitude"),
                ],
            },
            # {
            #     "label": "TC1 prime free surface",
            #     "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC1_prime",
            #     "snapshot_fields": [
            #         ("eta", "Eta"),
            #         ("etan", "ETAN"),
            #         ("psi", "PsiVEL"),
            #         ("phi", "PhiVEL"),
            #         ("velocity_magnitude", "Velocity Magnitude"),
            #     ],
            # },
        ],
        "media": [],
        "dynamic_gallery": "tc1",
    },
    {
        "slug": "testcase2",
        "title": "Williamson TC2: Global Steady-State Nonlinear Zonal Geostrophic Flow",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC2",
        "snapshot_fields": [
            ("eta", "Eta"),
            ("etan", "ETAN"),
            ("psi", "PsiVEL"),
            ("phi", "PhiVEL"),
            ("velocity_magnitude", "Velocity Magnitude"),
        ],
        "media": [],
        "dynamic_gallery": "tc2",
    },
    {
        "slug": "testcase3",
        "title": "Williamson TC3: Steady Geostrophic Flow with Compact Support",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC3",
        "snapshot_fields": [
            ("eta", "Eta"),
            ("etan", "ETAN"),
            ("psi", "PsiVEL"),
            ("phi", "PhiVEL"),
            ("velocity_magnitude", "Velocity Magnitude"),
        ],
        "media": [],
        "dynamic_gallery": "tc3",
    },
    {
        "slug": "testcase4",
        "title": "Williamson TC4: Forced Nonlinear Translating Low",
        "summary": "The U0=20 m s^-1 c1 HPC run is verified through day 5; snapshots, conservation assets, and analytic path checks are published.",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC4",
        "key_days": (0.0, 1.0, 2.0, 3.0, 4.0, 5.0),
        "snapshot_fields": [
            ("eta", "Eta"),
            ("etan", "ETAN"),
            ("psi", "PsiVEL"),
            ("phi", "PhiVEL"),
            ("velocity_magnitude", "Velocity Magnitude"),
        ],
        "media": [],
    },
    {
        "slug": "testcase5",
        "title": "Williamson TC5: Zonal Flow Over an Isolated Mountain",
        "summary": "The archived output becomes non-finite before the required day-10/day-15 checks; the corrected setup uses H0=5400 m, 15 s steps, and viscosity before rerun.",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC5",
        "key_days": (0.0, 3.0, 6.0, 9.0, 10.0, 12.0, 15.0),
        "snapshot_fields": [
            ("eta", "Eta"),
            ("etan", "ETAN"),
            ("psi", "PsiVEL"),
            ("phi", "PhiVEL"),
            ("velocity_magnitude", "Velocity Magnitude"),
        ],
        "media": [],
    },
    {
        "slug": "testcase6",
        "title": "Williamson TC6: Rossby-Haurwitz Wave",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC6",
        "key_days": (0.0, 3.0, 6.0, 9.0, 10.0, 12.0, 14.0),
        "snapshot_fields": [
            ("eta", "Eta"),
            ("etan", "ETAN"),
            ("psi", "PsiVEL"),
            ("phi", "PhiVEL"),
            ("velocity_magnitude", "Velocity Magnitude"),
        ],
        "media": [],
        "dynamic_gallery": "tc6",
    },
    {
        "slug": "testcase7",
        "title": "Williamson TC7: Analyzed 500 mb Initial State",
        "summary": "Analyzed-height and wind input is present; the corrected setup uses 25 s steps, large-scale filtering, polar tapering, and viscosity before rerun.",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC7",
        "snapshot_fields": [
            ("eta", "Eta"),
            ("etan", "ETAN"),
            ("psi", "PsiVEL"),
            ("phi", "PhiVEL"),
            ("velocity_magnitude", "Velocity Magnitude"),
        ],
        "media": [],
    },
]

DATA_COLUMNS = [
    "rhoConst",
    "gravity",
    "deltaT",
    "nTimeSteps",
    "delR",
    "delX",
    "delY",
]
SIZE_COLUMNS = ["sNx", "sNy", "nPx", "nPy", "Nx", "Ny", "mpi_ranks"]
DATA_COLUMN_META = {
    "rhoConst": ("rhoConst", "kg m<sup>-3</sup>"),
    "gravity": ("gravity", "m s<sup>-2</sup>"),
    "deltaT": ("deltaT", "s"),
    "nTimeSteps": ("nTimeSteps", "steps"),
    "delR": ("delR", "m"),
    "delX": ("delX", "deg"),
    "delY": ("delY", "deg"),
}
SIZE_COLUMN_META = {
    "sNx": ("sNx", "cells"),
    "sNy": ("sNy", "cells"),
    "nPx": ("nPx", "MPI tiles"),
    "nPy": ("nPy", "MPI tiles"),
    "Nx": ("Nx", "cells"),
    "Ny": ("Ny", "cells"),
    "mpi_ranks": ("mpi_ranks", "ranks"),
}
DAY_PATTERN = re.compile(r"_day_([0-9]+(?:\.[0-9]+)?)")
DEFAULT_KEY_SNAPSHOT_DAYS = (0.0, 3.0, 6.0, 9.0, 10.0, 12.0)
READY_ASSET_EXTENSIONS = {
    ".csv",
    ".json",
    ".md",
    ".png",
    ".jpg",
    ".jpeg",
    ".svg",
    ".txt",
    ".webp",
}
DIAGNOSIS_IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".svg", ".webp"}
DIAGNOSIS_DATA_EXTENSIONS = {".csv", ".json", ".md", ".txt"}
CASE_OUTPUT_NAMES = {
    "TC1_prime": "TestCase1Prime",
    "TC1": "TestCase1",
    "TC2": "TestCase2",
    "TC3": "TestCase3",
    "TC4": "TestCase4",
    "TC5": "TestCase5",
    "TC6": "TestCase6",
    "TC7": "TestCase7",
}
SNAPSHOT_FIELD_LABELS = {
    "tracer": "Passive Tracer",
    "eta": "Eta",
    "etan": "ETAN",
    "psi": "PsiVEL",
    "phi": "PhiVEL",
    "velocity_magnitude": "Velocity Magnitude",
}
SNAPSHOT_FIELD_UNITS = {
    "tracer": "psu",
    "eta": "m",
    "etan": "m",
    "psi": "m<sup>3</sup> s<sup>-1</sup>",
    "phi": "m<sup>2</sup> s<sup>-1</sup>",
    "velocity_magnitude": "m s<sup>-1</sup>",
}
SNAPSHOT_FIELD_ORDER = {
    field: index
    for index, field in enumerate(("tracer", "eta", "etan", "psi", "phi", "velocity_magnitude"))
}
FIELD_TAB_ORDER = (
    ("eta", "Eta"),
    ("etan", "ETAN"),
    ("tracer", "Tracer"),
    ("psi", "PsiVEL"),
    ("phi", "PhiVEL"),
    ("velocity_magnitude", "Velocity Magnitude"),
)

AUDIT_SECTIONS = (
    ("cfl-deltat-audit", "Time-step safety (CFL)", "Run checks", "Pending validation"),
    ("conservation-audit", "Conservation and health", "Run checks", "Pending validation"),
)

STATUS_META = {
    "testcase1": ("Verified", "verified"),
    "testcase2": ("Issues", "issues"),
    "testcase3": ("Issues", "issues"),
    "testcase4": ("Verified", "verified"),
    "testcase5": ("Rerun in progress", "pending"),
    "testcase6": ("Pending validation", "pending"),
    "testcase7": ("Submitted / preflight pass", "pending"),
    "case1_constant_bathymetry": ("Verified", "verified"),
    "case2_real_bathymetry": ("Verified", "verified"),
    "case3_geostrophic_adjustment": ("Verified", "verified"),
}

NAV_LABELS = {
    "testcase1": "TC1 - Cosine bell",
    "testcase2": "TC2 - Steady zonal",
    "testcase3": "TC3 - Compact jet",
    "testcase4": "TC4 - Translating low",
    "testcase5": "TC5 - Mountain flow",
    "testcase6": "TC6 - Rossby wave",
    "testcase7": "TC7 - 500 mb state",
    "case1_constant_bathymetry": "Patch - Constant depth",
    "case2_real_bathymetry": "Patch - Real bathymetry",
    "case3_geostrophic_adjustment": "Geostrophic adjustment",
}

DETAIL_ORDER = (
    ("definition", 0),
    ("what this test measures", 1),
    ("expected output", 2),
    ("initial conditions", 3),
    ("mitgcm equations", 3),
    ("parameters", 4),
)

def math_block(*equations: str) -> str:
    lines = "".join(f"<div class='math-line'>\\[{equation}\\]</div>" for equation in equations)
    return f"<div class='equation-box'>{lines}</div>"

COMMON_EQUATION_HTML = (
    "<p>The active Williamson cases are run as one-layer rotating shallow-water problems. "
    "The MITgcm fields reduce to horizontal velocity <code>u=(u,v)</code> and free surface "
    "<code>&eta;</code> over depth <code>H</code>.</p>"
    + math_block(
        "\\frac{D\\mathbf{u}}{Dt}+f\\hat{k}\\times\\mathbf{u}"
        "=-g\\nabla_s\\eta+\\mathcal{D}_{\\mathbf{u}}+\\mathbf{F}",
        "\\frac{\\partial \\eta}{\\partial t}+\\nabla_s\\cdot[(H+\\eta)\\mathbf{u}]=S_{\\eta}",
        "\\frac{Dq}{Dt}=\\mathcal{D}_{q}+S_q \\quad \\text{when passive tracers are enabled}",
    )
    +
    "<p>Each testcase changes the analytic initial state, forcing, bathymetry, or expected "
    "diagnostic behavior of this same system.</p>"
)
WILLIAMSON_DETAILS: Dict[str, List[Dict[str, str]]] = {
    "case1_constant_bathymetry": [
        {
            "title": "Definition",
            "body": (
                "<p>A Gaussian free-surface perturbation is initialized over constant "
                "bathymetry. The prognostic and diagnostic runs show how the same initial "
                "height anomaly adjusts under the MITgcm free-surface dynamics.</p>"
            ),
        },
        {
            "title": "MITgcm Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                + math_block(
                    "H(x,y)=4000\\,\\mathrm{m}",
                    "\\eta(x,y,0)=\\eta_0\\exp\\left[-\\frac{r^2}{2\\sigma^2}\\right]",
                )
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The Gaussian patch should radiate a clean gravity-wave adjustment pattern. "
                "The diagnostic initialization should look smoother, with less fast initial "
                "noise than the direct prognostic start.</p>"
            ),
        },
    ],
    "case2_real_bathymetry": [
        {
            "title": "Definition",
            "body": (
                "<p>The same Gaussian free-surface perturbation is placed over spatially varying "
                "bathymetry, so wave propagation and adjustment are shaped by the local depth field.</p>"
            ),
        },
        {
            "title": "MITgcm Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                + math_block(
                    "H=H(x,y)",
                    "c(x,y)=\\sqrt{gH(x,y)}",
                )
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The response should follow the bathymetric structure: faster propagation over "
                "deeper regions, slower propagation over shallow regions, and a diagnostic run "
                "that suppresses unbalanced high-frequency content.</p>"
            ),
        },
    ],
    "case3_geostrophic_adjustment": [
        {
            "title": "Definition",
            "body": (
                "<p>A free-surface field over real bathymetry is initialized near geostrophic "
                "balance, then compared between prognostic and diagnostic evolution.</p>"
            ),
        },
        {
            "title": "MITgcm Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                + math_block(
                    "f\\hat{k}\\times\\mathbf{u}\\approx -g\\nabla_s\\eta",
                    "\\nabla_s\\cdot\\left[(H+\\eta)\\mathbf{u}\\right]\\approx 0",
                )
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The balanced initialization should hold a coherent geostrophic pattern with "
                "limited inertial ringing. The diagnostic output should show less adjustment "
                "noise and cleaner large-scale structure.</p>"
            ),
        },
    ],
    "testcase1": [
        {
            "title": "Definition",
            "body": (
                "<p>Passive advection of a compact cosine bell by solid-body rotation on "
                "the sphere. The rotation axis can be tilted relative to the latitude-longitude "
                "grid, so the bell can cross the coordinate poles.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                "<p>The tracer is a bell centered at <code>(&lambda;c,&theta;c)=(3&pi;/2,0)</code> "
                "with compact radius <code>R=a/3</code>:</p>"
                + math_block(
                    "q=\\frac{q_0}{2}\\left[1+\\cos\\left(\\frac{\\pi r}{R}\\right)\\right],\\quad r<R",
                    "q=0,\\quad r\\ge R",
                    "\\frac{Dq}{Dt}=\\frac{\\partial q}{\\partial t}+\\mathbf{v}\\cdot\\nabla_s q=0",
                )
                +
                "<p>The prescribed velocity is the TC2 solid-body field, but momentum is not "
                "stepped in this local MITgcm case.</p>"
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li><code>a=6.371e6 m</code>, <code>g=9.81 m s^-2</code>, "
                "<code>&Omega;=7.292e-5 s^-1</code> in the local runs.</li>"
                "<li><code>q0=1000</code>, <code>R=a/3</code>, "
                "<code>U0=2&pi;a/(12 days)</code>.</li>"
                "<li>Rendered tilt angles: <code>&alpha;=0, 0.05, 1.52, 1.57 rad</code>.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>This is the clean transport test: conservation, phase accuracy, numerical "
                "diffusion, overshoots/undershoots, and the ability to move a compact structure "
                "through polar coordinate singularities without distortion.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>After one 12-day revolution the bell should return to its starting location. "
                "The exact and MITgcm contours should nearly overlap, the tracer integral should "
                "stay flat, and signed-error contours should remain small and symmetric rather "
                "than growing long noisy tails.</p>"
            ),
        },
    ],
    "testcase2": [
        {
            "title": "Definition",
            "body": (
                "<p>A full nonlinear shallow-water steady state: global solid-body zonal flow "
                "balanced by the free-surface height field. The flow is exact for every time "
                "when the continuous balance is preserved.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                "<p>With longitude <code>&lambda;</code>, latitude <code>&theta;</code>, and "
                "tilt <code>&alpha;</code>:</p>"
                + math_block(
                    "u=U_0\\left(\\cos\\theta\\cos\\alpha+\\cos\\lambda\\sin\\theta\\sin\\alpha\\right)",
                    "v=-U_0\\sin\\lambda\\sin\\alpha",
                    "\\mu=-\\cos\\lambda\\cos\\theta\\sin\\alpha+\\sin\\theta\\cos\\alpha",
                    "g h=g h_0-\\left(a\\Omega U_0+\\frac{1}{2}U_0^2\\right)\\mu^2",
                )
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li><code>U0=2&pi;a/(12 days)</code>.</li>"
                "<li><code>g h0=2.94e4 m^2 s^-2</code>, so local <code>H0=g h0/g=2996.94 m</code>.</li>"
                "<li>Standard tilts: <code>&alpha;=0, 0.05, &pi;/2-0.05, &pi;/2</code>; "
                "local outputs use <code>0, 0.05, 1.52, 1.57 rad</code>.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>Balanced nonlinear dynamics, metric terms, Coriolis terms, pole treatment, "
                "and preservation of a known steady state. The paper also identifies TC2 as the "
                "performance benchmark, typically using <code>&alpha;=&pi;/4</code> for a 5-day run.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>A good run stays almost motionless in the rotating balance: height and velocity "
                "errors should remain near discretization error, with little gravity-wave noise. "
                "Mesh refinement should reduce the 5-day error norms cleanly.</p>"
            ),
        },
    ],
    "testcase3": [
        {
            "title": "Definition",
            "body": (
                "<p>A steady nonlinear geostrophic flow like TC2, but the wind is a smooth compact "
                "jet. The active balanced flow occupies only a limited latitudinal band.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                "<p>The local implementation builds the velocity from a compact bump in rotated "
                "latitude <code>&theta;'</code> and computes the height by integrating the nonlinear "
                "geostrophic balance:</p>"
                + math_block(
                    "u'=U_0 B(\\theta'),\\quad v'=0",
                    "B(\\theta')=0\\quad \\text{outside }[-\\pi/6,\\pi/2]",
                    "h=H_0+\\int \\text{balanced pressure-gradient contribution}\\,d\\theta'",
                )
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li><code>U0=2&pi;a/(12 days)</code>, <code>H0=2996.94 m</code>.</li>"
                "<li>Compact support interval: <code>&theta;0=-&pi;/6</code> to "
                "<code>&theta;1=&pi;/2</code>.</li>"
                "<li>Runs can be rotated by <code>&alpha;</code> to move the jet away from grid lines.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>Whether a model can maintain a localized balanced flow without leakage, phase "
                "drift, ringing near smooth compact-support edges, or spurious gravity waves.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The compact jet and free-surface anomaly should stay steady. Errors should be "
                "localized near the jet edges if they appear at all, and conservation diagnostics "
                "should show little long-term drift.</p>"
            ),
        },
    ],
    "testcase4": [
        {
            "title": "Definition",
            "body": (
                "<p>A forced nonlinear shallow-water exact-solution test with a translating low. "
                "The prescribed forcing makes the unsteady analytic disturbance the target solution.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                "<p>The local TC4 code now injects the translating-low analytic source terms through "
                "<code>code/tc4_forcing.F</code> and <code>code/external_forcing_surf.F</code>. "
                "The initial fields use the same streamfunction and geopotential definitions as the forcing:</p>"
                + math_block(
                    "u_b(\\theta)=U_0\\sin^{14}(2\\theta)",
                    "\\psi=\\psi_0\\exp\\left[-\\sigma\\frac{1-c}{1+c}\\right]",
                    "c=\\sin\\theta_0\\sin\\theta+\\cos\\theta_0\\cos\\theta\\cos(\\lambda-U_0t/a-\\lambda_0)",
                    "h=\\frac{\\Phi_b(\\theta)+f_0\\psi}{g},\\quad \\eta=h-H_0",
                    "\\mathbf{F}=\\partial_t\\mathbf{u}+\\mathbf{u}\\cdot\\nabla_s\\mathbf{u}+f\\hat{k}\\times\\mathbf{u}+g\\nabla_s\\eta",
                )
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li><code>U0=20 m s^-1</code> for completed c1 run <code>run_u0_20</code>; the same code can use <code>TC4_U0_VALUE=40</code> for c2.</li>"
                "<li><code>H0=gh0/g=10193.67991845056 m</code>, flat bathymetry.</li>"
                "<li>Published snapshots and conservation assets use the verified day-0 through day-5 <code>run_u0_20</code> output.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>TC4 targets unsteady nonlinear accuracy with a known forced "
                "solution: advection of a strong localized low, source-term balance, and error "
                "growth in a moving feature.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The translating low should keep its shape as it moves eastward. The completed c1 "
                "run passes the first analytic check: the low center follows the expected path and the "
                "day-5 relative L2 errors are small for <code>eta</code>, <code>u</code>, and <code>v</code>. "
                "The published snapshots and conservation assets support marking this c1 run verified.</p>"
            ),
        },
    ],
    "testcase5": [
        {
            "title": "Definition",
            "body": (
                "<p>Zonal shallow-water flow impinging on an isolated conical mountain. Unlike "
                "TC2 and TC3, this is not an exact steady-state comparison; topography actively "
                "generates waves.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                "<p>The initialized flow is zonal and geostrophically adjusted around mountain "
                "topography:</p>"
                + math_block(
                    "u=U_0\\cos\\theta,\\quad v=0",
                    "h_s=h_{s0}\\left(1-\\frac{r}{R}\\right),\\quad r<R",
                    "h_s=0,\\quad r\\ge R",
                    "b=-(H_0-h_s)",
                )
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li><code>U0=20 m s^-1</code>, <code>H0=5400 m</code>.</li>"
                "<li>Mountain height <code>h_s0=2000 m</code>, radius <code>R=&pi;/9</code>.</li>"
                "<li>Mountain center <code>(&lambda;c,&theta;c)=(3&pi;/2,&pi;/6)</code>.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>Topographic source-term handling, conservation over variable bathymetry, wave "
                "generation, and whether the discretization creates unphysical noise near steep terrain.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The mountain should produce a coherent downstream wave response. There is no "
                "simple analytic truth field, so evaluate against reference solutions, global "
                "invariants, smoothness, and absence of terrain-locked artifacts.</p>"
            ),
        },
    ],
    "testcase6": [
        {
            "title": "Definition",
            "body": (
                "<p>A Rossby-Haurwitz wave initialized from a nondivergent planetary wave pattern. "
                "The shallow-water version is a demanding nonlinear wave propagation benchmark.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                "<p>The local fields follow the standard wave-number-4 streamfunction/height form:</p>"
                + math_block(
                    "r=4",
                    "u=a\\omega\\cos\\theta+aK\\cos^{r-1}\\theta\\left(r\\sin^2\\theta-\\cos^2\\theta\\right)\\cos(r\\lambda)",
                    "v=-aKr\\cos^{r-1}\\theta\\sin\\theta\\sin(r\\lambda)",
                    "h=H_0+\\frac{a^2}{g}\\left[A(\\theta)+B(\\theta)\\cos(r\\lambda)+C(\\theta)\\cos(2r\\lambda)\\right]",
                )
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li><code>r=4</code>, <code>K=7.848e-6 s^-1</code>, "
                "<code>&omega;=7.848e-6 s^-1</code>.</li>"
                "<li><code>H0=8000 m</code>, flat bathymetry.</li>"
                "<li>Snapshots include the common key days 0, 3, 6, 9, 10, and 12, plus day 14 when available.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>Longer-range nonlinear wave propagation, phase and amplitude errors, pole behavior, "
                "and conservation of geopotential, energy, vorticity, divergence, and potential enstrophy.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The wave should retain a coherent planetary-scale structure with controlled phase "
                "drift. Look for smooth contour evolution and slow, explainable drift in global invariants "
                "rather than noisy small-scale breakdown.</p>"
            ),
        },
    ],
    "testcase7": [
        {
            "title": "Definition",
            "body": (
                "<p>Realistic analyzed 500 mb height and wind states used as shallow-water initial "
                "conditions. This test gives the suite atmospheric texture rather than another exact solution.</p>"
            ),
        },
        {
            "title": "Initial Conditions And Equations",
            "body": (
                f"{COMMON_EQUATION_HTML}"
                "<p>The paper uses analyzed atmospheric states truncated to T42 spectral resolution, "
                "with optional nonlinear normal-mode initialization. This repository applies "
                "a large-scale filter to its local analysis file before staging MITgcm input. "
                "It has "
                "<code>input/raw/tc7_initial_conditions.npz</code> and generated binary "
                "height/wind fields staged for the submitted run.</p>"
            ),
        },
        {
            "title": "Parameters",
            "body": (
                "<ul>"
                "<li>Local reference depth <code>H0=8000 m</code>.</li>"
                "<li>Paper states include 0000 GMT 21 Dec 1978, 16 Jan 1979, and 9 Jan 1979 examples.</li>"
                "<li>Local grid: <code>1440 x 720</code> at <code>0.25 deg</code>; the submitted rerun uses a 25 s step, large-scale filtering, polar velocity tapering, and <code>viscAh=1e5 m2 s^-1</code>.</li>"
                "</ul>"
            ),
        },
        {
            "title": "What This Test Measures",
            "body": (
                "<p>Robustness on realistic large-scale atmospheric states: polar flow, blocking-like "
                "patterns, strong zonal flow, and sensitivity of a numerical method to the character "
                "of the initial state.</p>"
            ),
        },
        {
            "title": "Expected Output",
            "body": (
                "<p>The shallow-water model is not expected to reproduce the real atmosphere perfectly. "
                "Use high-resolution reference comparisons, day-1/day-5 forecast-error plots, and the "
                "same global invariants as TC6 to judge behavior.</p>"
            ),
        },
    ],
}
COPIED_ASSETS: Set[Path] = set()
MISSING_SNAPSHOT_ROOTS: Set[Path] = set()

def parse_data_file(path: Path) -> Dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values: dict[str, str] = {}
    for key in DATA_COLUMNS:
        match = re.search(rf"\b{re.escape(key)}\s*=\s*([^,\n/]+)", text, re.IGNORECASE)
        values[key] = match.group(1).strip() if match else "—"
    return values

def parse_size_file(path: Path) -> Dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values: dict[str, int] = {}
    for key in ("sNx", "sNy", "nSx", "nSy", "nPx", "nPy"):
        match = re.search(rf"\b{re.escape(key)}\s*=\s*([0-9]+)", text)
        values[key] = int(match.group(1)) if match else 0

    nx = values["sNx"] * values["nSx"] * values["nPx"]
    ny = values["sNy"] * values["nSy"] * values["nPy"]
    mpi_ranks = values["nPx"] * values["nPy"]
    return {
        "sNx": str(values["sNx"]),
        "sNy": str(values["sNy"]),
        "nPx": str(values["nPx"]),
        "nPy": str(values["nPy"]),
        "Nx": str(nx),
        "Ny": str(ny),
        "mpi_ranks": str(mpi_ranks),
    }

def relative_to_root(path: Path, root: Path) -> Optional[Path]:
    try:
        return path.relative_to(root)
    except ValueError:
        return None

def ensure_docs_asset(source: Path, slug: str) -> str:
    source = source.resolve()
    docs_relative = relative_to_root(source, DOCS_DIR.resolve())
    if docs_relative is not None:
        return docs_relative.as_posix()

    sandbox_output_relative = relative_to_root(source, SANDBOX_OUTPUT_ROOT.resolve())
    if sandbox_output_relative is not None:
        target = ASSET_ROOT / sandbox_output_relative
    else:
        sandbox_relative = relative_to_root(source, SANDBOX_DIR.resolve())
        if sandbox_relative is not None:
            target = ASSET_ROOT / sandbox_relative
        else:
            target = ASSET_ROOT / slug / source.name

    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    COPIED_ASSETS.add(target)
    return target.relative_to(DOCS_DIR).as_posix()

def sync_ready_output_assets() -> None:
    if not SANDBOX_OUTPUT_ROOT.exists():
        return

    for source in sorted(SANDBOX_OUTPUT_ROOT.glob("TestCase*/**/*"), key=lambda path: path.as_posix()):
        if not source.is_file() or source.suffix.lower() not in READY_ASSET_EXTENSIONS:
            continue
        if "Diagnosis" not in source.parts and "Snapshots" not in source.parts:
            continue
        ensure_docs_asset(source, "ready-output")

def natural_key(text: str) -> List[str]:
    return [part.zfill(12) if part.isdigit() else part for part in re.split(r"(\d+)", text)]

def path_sort_key(path: Path) -> Tuple[object, ...]:
    day = extract_day(path)
    day_key = f"{day:012.6f}" if day is not None else ""
    return (*natural_key(path.as_posix()), day_key)

def extract_day(path: Path) -> Optional[float]:
    match = DAY_PATTERN.search(path.name)
    return float(match.group(1)) if match else None

def format_number(value: float) -> str:
    return f"{value:.2f}".rstrip("0").rstrip(".")

def day_caption(path: Path) -> str:
    day = extract_day(path)
    return "Day unknown" if day is None else f"Day {format_number(day)}"

def path_preference(path: Path) -> Tuple[int, int]:
    return (1 if "_iter" in path.name else 0, len(path.name))

def alpha_from_name(name: str) -> str:
    if name.startswith("alpha_"):
        return name[len("alpha_"):]
    if name.startswith("run_alpha_"):
        return name[len("run_alpha_"):]
    return name

def label_from_token(token: str) -> str:
    return SNAPSHOT_FIELD_LABELS.get(token, token.replace("_", " ").title())

def unit_from_token(token: str) -> str:
    return SNAPSHOT_FIELD_UNITS.get(token, "")

def label_with_unit(label: str, unit: str) -> str:
    label_html = html.escape(label)
    if not unit:
        return label_html
    return f"{label_html} <span class='unit'>[{unit}]</span>"

def infer_case_output_name(case_dir: Path) -> Optional[str]:
    case_name = case_dir.name.lower()
    for case_code, output_name in sorted(CASE_OUTPUT_NAMES.items(), key=lambda item: len(item[0]), reverse=True):
        if case_code.lower() in case_name:
            return output_name
    return None

def section_case_configs(section: Dict[str, object]) -> List[Dict[str, object]]:
    configured_cases = section.get("cases")
    if configured_cases:
        return [dict(case) for case in configured_cases]  # type: ignore[arg-type]

    case_dir = section.get("case_dir")
    if case_dir is None:
        return []

    case: dict[str, object] = {"case_dir": case_dir}
    for key in ("label", "output_name", "snapshot_fields", "key_days"):
        if key in section:
            case[key] = section[key]
    return [case]

def case_display_label(case: Dict[str, object]) -> str:
    label = case.get("label")
    if label:
        return str(label)
    case_dir = case.get("case_dir")
    return Path(case_dir).name if case_dir is not None else ""

def case_snapshot_root(case: Dict[str, object]) -> Optional[Path]:
    case_dir = case.get("case_dir")
    if case_dir is None:
        return None
    case_path = Path(case_dir)
    output_name = case.get("output_name") or infer_case_output_name(case_path)
    if output_name is None:
        return None
    return SANDBOX_OUTPUT_ROOT / str(output_name) / "Snapshots"

def case_output_root(case: Dict[str, object]) -> Optional[Path]:
    case_dir = case.get("case_dir")
    if case_dir is None:
        return None
    case_path = Path(case_dir)
    output_name = case.get("output_name") or infer_case_output_name(case_path)
    if output_name is None:
        return None
    return SANDBOX_OUTPUT_ROOT / str(output_name)

def case_diagnosis_root(case: Dict[str, object]) -> Optional[Path]:
    output_root = case_output_root(case)
    if output_root is None:
        return None
    return output_root / "Diagnosis"

def normalize_snapshot_field(field: object) -> Tuple[str, str]:
    if isinstance(field, dict):
        folder = str(field["folder"])
        return folder, str(field.get("label", label_from_token(folder)))
    if isinstance(field, (list, tuple)):
        folder = str(field[0])
        label = str(field[1]) if len(field) > 1 else label_from_token(folder)
        return folder, label
    folder = str(field)
    return folder, label_from_token(folder)

def snapshot_field_sort_key(field: Tuple[str, str]) -> Tuple[int, List[str]]:
    folder, label = field
    return (SNAPSHOT_FIELD_ORDER.get(folder, len(SNAPSHOT_FIELD_ORDER)), natural_key(label))

def discover_snapshot_fields(snapshot_root: Path) -> List[Tuple[str, str]]:
    fields: set[str] = set()
    for alpha_dir in snapshot_root.glob("alpha_*"):
        if not alpha_dir.is_dir():
            continue
        for field_dir in alpha_dir.iterdir():
            if field_dir.is_dir() and any(field_dir.glob("*.png")):
                fields.add(field_dir.name)
    return sorted(((field, label_from_token(field)) for field in fields), key=snapshot_field_sort_key)

def case_snapshot_fields(case: Dict[str, object], snapshot_root: Path) -> List[Tuple[str, str]]:
    configured_fields = case.get("snapshot_fields")
    if configured_fields:
        return [normalize_snapshot_field(field) for field in configured_fields]  # type: ignore[union-attr]
    return discover_snapshot_fields(snapshot_root)

def reduce_to_one_image_per_day(paths: List[Path]) -> List[Path]:
    by_day: dict[float | str, Path] = {}
    for path in paths:
        key: float | str = extract_day(path)
        if key is None:
            key = path.name
        old_path = by_day.get(key)
        if old_path is None or path_preference(path) > path_preference(old_path):
            by_day[key] = path
    return sorted(by_day.values(), key=path_sort_key)

def snapshot_key_days(config: Dict[str, object] | None = None) -> Tuple[float, ...]:
    if config is None:
        return DEFAULT_KEY_SNAPSHOT_DAYS
    days = config.get("key_days")
    if not days:
        return DEFAULT_KEY_SNAPSHOT_DAYS
    return tuple(float(day) for day in days)  # type: ignore[union-attr]

def is_key_snapshot_day(path: Path, key_days: Tuple[float, ...]) -> bool:
    day = extract_day(path)
    return day is not None and any(abs(day - key_day) < 0.01 for key_day in key_days)

def case_snapshot_items(snapshot_root: Path, folder: str, field_label: str) -> List[Dict[str, object]]:
    items: list[dict[str, object]] = []
    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        field_dir = alpha_dir / folder
        if not field_dir.is_dir():
            continue
        alpha = alpha_from_name(alpha_dir.name)
        paths = sorted(field_dir.glob("*.png"), key=path_sort_key)
        items.extend(
            {
                "source": path,
                "caption": f"alpha {alpha}, {field_label}, {day_caption(path)}",
            }
            for path in reduce_to_one_image_per_day(paths)
        )
    return items

def case_alpha_snapshot_items(
    snapshot_root: Path,
    alpha_dir: Path,
    folder: str,
    field_label: str,
    key_days: Tuple[float, ...],
) -> List[Dict[str, object]]:
    field_dir = alpha_dir / folder
    if not field_dir.is_dir():
        return []

    alpha = alpha_from_name(alpha_dir.name)
    paths = [
        path
        for path in reduce_to_one_image_per_day(sorted(field_dir.glob("*.png"), key=path_sort_key))
        if is_key_snapshot_day(path, key_days)
    ]
    return [
        {
            "source": path,
            "caption": f"alpha {alpha}, {field_label}, {day_caption(path)}",
        }
        for path in paths
    ]

def render_plot_figure(item: Dict[str, object], slug: str, figure_class: str = "subfigure") -> str:
    source = Path(item["source"])
    if not source.exists():
        return ""
    relative_path = ensure_docs_asset(source, slug)
    caption = html.escape(str(item["caption"]))
    return (
        f"<figure class='{html.escape(figure_class)}'>"
        "<button class='plot-open' type='button' "
        f"data-modal-src='{html.escape(relative_path)}' "
        f"data-modal-caption='{caption}'>"
        f"<img src='{html.escape(relative_path)}' alt='{caption}' loading='lazy' />"
        "</button>"
        f"<figcaption>{caption}</figcaption>"
        "</figure>"
    )

def render_snapshot_grid(items: List[Dict[str, object]], slug: str) -> str:
    figures = "".join(render_plot_figure(item, slug) for item in items)
    if not figures:
        return "<p class='empty'>No key-day snapshots available.</p>"
    return f"<div class='snapshot-grid'>{figures}</div>"

def render_field_units_note(fields: List[Tuple[str, str]]) -> str:
    entries = [
        label_with_unit(label, unit_from_token(folder))
        for folder, label in fields
        if unit_from_token(folder)
    ]
    if not entries:
        return ""
    return f"<p class='field-units'>Units: {', '.join(entries)}.</p>"

def render_case_snapshot_browser(case: Dict[str, object], slug: str) -> str:
    snapshot_root = case_snapshot_root(case)
    if snapshot_root is None:
        return ""
    if not snapshot_root.exists():
        MISSING_SNAPSHOT_ROOTS.add(snapshot_root)
        return ""

    configured_fields = dict(case_snapshot_fields(case, snapshot_root))
    alpha_dirs = [
        alpha_dir
        for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name))
        if alpha_dir.is_dir()
    ]
    if not alpha_dirs:
        return ""

    key_days = snapshot_key_days(case)
    alpha_blocks: list[str] = []
    for alpha_index, alpha_dir in enumerate(alpha_dirs):
        alpha = alpha_from_name(alpha_dir.name)
        tab_base = f"{slug}-{alpha_dir.name.replace('.', '_')}"
        buttons: list[str] = []
        panels: list[str] = []
        first_panel = True

        for folder, button_label in FIELD_TAB_ORDER:
            field_label = configured_fields.get(folder, label_from_token(folder))
            items = case_alpha_snapshot_items(snapshot_root, alpha_dir, folder, field_label, key_days)
            if not items:
                continue

            tab_id = f"{tab_base}-{folder}"
            active_class = " is-active" if first_panel else ""
            selected = "true" if first_panel else "false"
            hidden = "" if first_panel else " hidden"
            button_label_html = label_with_unit(button_label, unit_from_token(folder))
            buttons.append(
                "<button class='field-tab"
                f"{active_class}' type='button' data-tab-target='{html.escape(tab_id)}' "
                f"aria-selected='{selected}'>{button_label_html}</button>"
            )
            panels.append(
                f"<div class='tab-panel{active_class}' data-tab-panel='{html.escape(tab_id)}'{hidden}>"
                f"{render_snapshot_grid(items, slug)}"
                "</div>"
            )
            first_panel = False

        if not panels:
            continue

        diagnostics_label = "Diagnostics" if slug == "testcase6" else "Error"
        buttons.append(
            f"<a class='field-tab field-tab-link' href='#{html.escape(slug)}-errors'>"
            f"{html.escape(diagnostics_label)}</a>"
        )
        alpha_blocks.append(
            f"<details class='alpha-panel' {'open' if alpha_index == 0 else ''}>"
            f"<summary>&alpha;={html.escape(alpha)}</summary>"
            "<div class='field-tabs'>"
            f"<div class='tab-list'>{''.join(buttons)}</div>"
            f"{''.join(panels)}"
            "</div>"
            "</details>"
        )

    if not alpha_blocks:
        return ""

    days = ", ".join(format_number(day) for day in key_days)
    label = html.escape(case_display_label(case))
    field_units_note = render_field_units_note(case_snapshot_fields(case, snapshot_root))
    return (
        "<div class='media-section snapshot-browser'>"
        f"<h3>{label}: key-day snapshots</h3>"
        f"<p class='section-note'>Rendered days: {html.escape(days)}.</p>"
        f"{field_units_note}"
        f"{''.join(alpha_blocks)}"
        "</div>"
    )

def render_case_snapshot_galleries(case: Dict[str, object], slug: str) -> str:
    return render_case_snapshot_browser(case, slug)

def comparison_caption(path: Path) -> str:
    parts = list(path.parts)
    for part in parts:
        match = re.match(r"compare_run_alpha_([^_]+)_vs_run_alpha_([^/]+)$", part)
        if match:
            return f"alpha {match.group(1)} vs alpha {match.group(2)}"

    if "comparison" in parts:
        index = parts.index("comparison")
        if index + 2 < len(parts):
            alpha_1 = alpha_from_name(parts[index + 1])
            alpha_2_name = parts[index + 2]
            if alpha_2_name.startswith("vs_"):
                alpha_2_name = alpha_2_name[len("vs_"):]
            alpha_2 = alpha_from_name(alpha_2_name)
            return f"alpha {alpha_1} vs alpha {alpha_2}"

    return "TC1 tracer comparison"

def final_day_from_error_table(error_dir: Path) -> Optional[float]:
    table = error_dir / "tc1_error_table.csv"
    if not table.exists():
        return None
    with table.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    value = rows[-1].get("day")
    return float(value) if value not in (None, "") else None

def tc1_snapshot_items() -> List[Dict[str, object]]:
    snapshot_root = SANDBOX_OUTPUT_ROOT / "TestCase1" / "Snapshots"
    items: list[dict[str, object]] = []
    if not snapshot_root.exists():
        return items

    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        alpha = alpha_from_name(alpha_dir.name)
        for folder in ("tracer", "eta"):
            paths = sorted((alpha_dir / folder).glob("*.png"), key=path_sort_key)
            if paths:
                by_day: dict[float | str, Path] = {}
                for path in paths:
                    key: float | str = extract_day(path)
                    if key is None:
                        key = path.name
                    old_path = by_day.get(key)
                    if old_path is None or path_preference(path) > path_preference(old_path):
                        by_day[key] = path
                items.extend(
                    {
                        "source": path,
                        "caption": f"alpha {alpha}, {day_caption(path)}",
                    }
                    for path in sorted(by_day.values(), key=path_sort_key)
                )
                break
    return items

def tc1_tracer_overlay_items() -> List[Dict[str, object]]:
    central = (
        SANDBOX_OUTPUT_ROOT
        / "TestCase1"
        / "Diagnosis"
        / "comparison"
    ).glob("alpha_*/vs_alpha_*/passive_tracer_overlay/*.png")
    legacy = (
        SANDBOX_DIR
        / "vortexSphere_Williamson_TC1"
        / "Scripts"
        / "output"
    ).glob("compare_*/PassiveTracerOverlay/tc1_tracer_overlay_day_*.png")
    paths = sorted([*central, *legacy], key=path_sort_key)
    by_case_day: dict[tuple[str, float | str], Path] = {}
    for path in paths:
        day_key: float | str = extract_day(path)
        if day_key is None:
            day_key = path.name
        key = (comparison_caption(path), day_key)
        old_path = by_case_day.get(key)
        if old_path is None or path_preference(path) > path_preference(old_path):
            by_case_day[key] = path
    return [
        {
            "source": path,
            "caption": f"{comparison_caption(path)}, {day_caption(path)}",
        }
        for path in sorted(by_case_day.values(), key=path_sort_key)
    ]

def tc1_error_contour_items() -> List[Dict[str, object]]:
    items: list[dict[str, object]] = []
    seen_alpha: set[str] = set()

    central_paths = sorted(
        (
            SANDBOX_OUTPUT_ROOT
            / "TestCase1"
            / "Diagnosis"
            / "error"
        ).glob("alpha_*/tc1_signed_error_contours.png"),
        key=path_sort_key,
    )
    legacy_paths = sorted(
        (
            SANDBOX_DIR
            / "vortexSphere_Williamson_TC1"
            / "Scripts"
            / "output"
        ).glob("run_alpha_*/TC1Error/tc1_signed_error_contours.png"),
        key=path_sort_key,
    )

    for path in [*central_paths, *legacy_paths]:
        alpha = alpha_from_name(path.parent.name)
        if path.parent.name == "TC1Error":
            alpha = alpha_from_name(path.parent.parent.name)
        if alpha in seen_alpha:
            continue
        seen_alpha.add(alpha)
        final_day = final_day_from_error_table(path.parent)
        day_text = f", Day {format_number(final_day)}" if final_day is not None else ""
        items.append(
            {
                "source": path,
                "caption": f"alpha {alpha}, signed error contours{day_text}",
            }
        )
    return items

def final_day_from_tc2_table(error_dir: Path) -> Optional[float]:
    table = error_dir / "TC2_error_norms.csv"
    if not table.exists():
        return None
    with table.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    value = rows[-1].get("day")
    return float(value) if value not in (None, "") else None

def tc2_snapshot_items() -> List[Dict[str, object]]:
    snapshot_root = SANDBOX_OUTPUT_ROOT / "TestCase2" / "Snapshots"
    items: list[dict[str, object]] = []
    if not snapshot_root.exists():
        return items

    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        alpha = alpha_from_name(alpha_dir.name)
        for folder, label in (("eta", "eta"), ("etan", "ETAN")):
            paths = sorted((alpha_dir / folder).glob("*.png"), key=path_sort_key)
            if not paths:
                continue
            by_day: dict[float | str, Path] = {}
            for path in paths:
                key: float | str = extract_day(path)
                if key is None:
                    key = path.name
                old_path = by_day.get(key)
                if old_path is None or path_preference(path) > path_preference(old_path):
                    by_day[key] = path
            items.extend(
                {
                    "source": path,
                    "caption": f"alpha {alpha}, {label}, {day_caption(path)}",
                }
                for path in sorted(by_day.values(), key=path_sort_key)
            )
            break
    return items

def tc2_error_norm_items() -> List[Dict[str, object]]:
    error_root = SANDBOX_OUTPUT_ROOT / "TestCase2" / "Diagnosis" / "error"
    items: list[dict[str, object]] = []
    if not error_root.exists():
        return items

    for path in sorted(error_root.glob("alpha_*/TC2_error_norms.png"), key=path_sort_key):
        alpha = alpha_from_name(path.parent.name)
        final_day = final_day_from_tc2_table(path.parent)
        day_text = f", final day {format_number(final_day)}" if final_day is not None else ""
        items.append(
            {
                "source": path,
                "caption": f"alpha {alpha}, error norms{day_text}",
            }
        )
    return items

def tc3_error_norm_items() -> List[Dict[str, object]]:
    error_root = SANDBOX_OUTPUT_ROOT / "TestCase3" / "Diagnosis" / "error"
    items: list[dict[str, object]] = []
    if not error_root.exists():
        return items

    for path in sorted(error_root.glob("alpha_*/TC3_error_norms.png"), key=path_sort_key):
        alpha = alpha_from_name(path.parent.name)
        table = path.parent / "TC3_error_norms.csv"
        final_day = None
        if table.exists():
            with table.open("r", newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            if rows and rows[-1].get("day"):
                final_day = float(rows[-1]["day"])
        day_text = f", final day {format_number(final_day)}" if final_day is not None else ""
        items.append(
            {
                "source": path,
                "caption": f"alpha {alpha}, error norms{day_text}",
            }
        )
    return items


def tc6_postprocessing_items() -> List[Dict[str, object]]:
    root = SANDBOX_OUTPUT_ROOT / "TestCase6" / "Diagnosis" / "postprocessing"
    labels = (
        ("postprocessing_vorticity_pv.png", "relative vorticity and Q diagnostics"),
        ("postprocessing_energies.png", "energy components"),
        ("postprocessing_mass.png", "mass anomaly"),
    )
    items: list[dict[str, object]] = []
    if not root.exists():
        return items

    for alpha_dir in sorted(root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        if not alpha_dir.is_dir():
            continue
        alpha = alpha_from_name(alpha_dir.name)
        for filename, label in labels:
            path = alpha_dir / filename
            if path.exists():
                items.append(
                    {
                        "source": path,
                        "caption": f"TC6 alpha {alpha}, {label}",
                    }
                )
    return items


def render_column_header(column: str, column_meta: Dict[str, Tuple[str, str]]) -> str:
    label, unit = column_meta.get(column, (column, ""))
    unit_html = f"<span class='table-unit'>[{unit}]</span>" if unit else ""
    return f"<th><span class='table-symbol'>{html.escape(label)}</span>{unit_html}</th>"

def render_table(
    title: str,
    columns: List[str],
    values: Dict[str, str],
    column_meta: Dict[str, Tuple[str, str]],
) -> str:
    header = "".join(render_column_header(column, column_meta) for column in columns)
    row = "".join(f"<td>{html.escape(values.get(column, '—'))}</td>" for column in columns)
    return (
        f"<div class='table-block data-panel'><h3>{html.escape(title)}</h3>"
        f"<div class='table-scroll'><table><thead><tr>{header}</tr></thead>"
        f"<tbody><tr>{row}</tr></tbody></table></div></div>"
    )

def render_case_tables(case: Dict[str, object], show_label: bool) -> str:
    case_dir = case.get("case_dir")
    if case_dir is None:
        return ""

    case_path = Path(case_dir)
    data_path = case_path / "input" / "data"
    size_path = case_path / "code" / "SIZE.h"
    if not data_path.exists() or not size_path.exists():
        return ""

    data_values = parse_data_file(data_path)
    size_values = parse_size_file(size_path)
    heading = (
        f"<h3 class='case-heading'>{html.escape(case_display_label(case))}</h3>"
        if show_label
        else ""
    )
    return (
        "<details class='details-panel case-block'>"
        "<summary>Numerical settings</summary>"
        f"{heading}"
        "<div class='tables'>"
        f"{render_table('Numerical Settings', DATA_COLUMNS, data_values, DATA_COLUMN_META)}"
        f"{render_table('Grid And MPI Layout', SIZE_COLUMNS, size_values, SIZE_COLUMN_META)}"
        "</div>"
        "</details>"
    )

def count_case_key_snapshots(case: Dict[str, object]) -> Tuple[int, int, int]:
    snapshot_root = case_snapshot_root(case)
    if snapshot_root is None or not snapshot_root.exists():
        return 0, 0, 0

    alphas: set[str] = set()
    fields: set[str] = set()
    count = 0
    key_days = snapshot_key_days(case)
    configured_fields = dict(case_snapshot_fields(case, snapshot_root))
    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        if not alpha_dir.is_dir():
            continue
        alpha_has_items = False
        for folder in configured_fields:
            paths = [
                path
                for path in reduce_to_one_image_per_day(sorted((alpha_dir / folder).glob("*.png"), key=path_sort_key))
                if is_key_snapshot_day(path, key_days)
            ]
            if paths:
                fields.add(folder)
                count += len(paths)
                alpha_has_items = True
        if alpha_has_items:
            alphas.add(alpha_from_name(alpha_dir.name))
    return len(alphas), len(fields), count

def slug_token(text: str) -> str:
    token = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return token or "case"

def diagnosis_asset_label(path: Path, root: Path) -> str:
    relative_path = path.relative_to(root)
    stem_parts = relative_path.with_suffix("").parts
    return " / ".join(part.replace("_", " ") for part in stem_parts)

def discover_case_diagnosis_assets(case: Dict[str, object]) -> Tuple[List[Path], List[Path]]:
    root = case_diagnosis_root(case)
    if root is None or not root.exists():
        return [], []

    images: list[Path] = []
    data: list[Path] = []
    for path in sorted(root.rglob("*"), key=path_sort_key):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix in DIAGNOSIS_IMAGE_EXTENSIONS:
            images.append(path)
        elif suffix in DIAGNOSIS_DATA_EXTENSIONS:
            data.append(path)
    return images, data

def count_case_diagnosis_assets(case: Dict[str, object]) -> Tuple[int, int]:
    images, data = discover_case_diagnosis_assets(case)
    return len(images), len(data)

def render_metric_cards(section: Dict[str, object], case_configs: List[Dict[str, object]]) -> str:
    metrics: list[tuple[str, str]] = []
    alpha_count = 0
    field_count = 0
    snapshot_count = 0
    diagnosis_count = 0
    for case in case_configs:
        case_alpha_count, case_field_count, case_snapshot_count = count_case_key_snapshots(case)
        diagnosis_images, diagnosis_data = count_case_diagnosis_assets(case)
        alpha_count += case_alpha_count
        field_count = max(field_count, case_field_count)
        snapshot_count += case_snapshot_count
        diagnosis_count += diagnosis_images + diagnosis_data

    error_count = 0
    if section.get("dynamic_gallery") == "tc1":
        error_count = len(tc1_error_contour_items())
    elif section.get("dynamic_gallery") == "tc2":
        error_count = len(tc2_error_norm_items())
    elif section.get("dynamic_gallery") == "tc3":
        error_count = len(tc3_error_norm_items())
    media_count = sum(1 for media in section.get("media", []) if Path(media["source"]).exists())
    if snapshot_count:
        key_day_labels = [
            ", ".join(format_number(day) for day in snapshot_key_days(case))
            for case in case_configs
        ]
        key_day_value = key_day_labels[0] if len(set(key_day_labels)) == 1 else "; ".join(key_day_labels)
        metrics = [
            ("Alphas", str(alpha_count)),
            ("Fields", str(field_count)),
            ("Key Days", key_day_value),
            ("Diagnosis Assets", str(diagnosis_count or error_count)),
        ]
    else:
        metrics = [
            ("Main Panels", str(media_count)),
            ("Diagnosis Assets", str(diagnosis_count)),
            ("Equations", "Shown"),
            ("Expected Output", "Included"),
        ]

    cards = "".join(
        "<div class='metric-card'>"
        f"<span>{html.escape(label)}</span>"
        f"<strong>{html.escape(value)}</strong>"
        "</div>"
        for label, value in metrics
    )
    return f"<div class='metric-grid'>{cards}</div>"

def render_experiment_summary(section: Dict[str, object]) -> str:
    summary = str(section.get("summary", "")).strip()
    if not summary:
        return ""
    return f"<p class='experiment-summary'>{html.escape(summary)}</p>"

def detail_sort_key(indexed_detail: Tuple[int, Dict[str, str]]) -> Tuple[int, int]:
    original_index, detail = indexed_detail
    title = str(detail["title"]).lower()
    for token, rank in DETAIL_ORDER:
        if token in title:
            return rank, original_index
    return len(DETAIL_ORDER), original_index

def detail_kind(title: str) -> str:
    normalized = title.lower()
    if "definition" in normalized:
        return "definition"
    if "what this test measures" in normalized:
        return "measure"
    if "expected output" in normalized:
        return "expected"
    if "initial conditions" in normalized or "mitgcm equations" in normalized:
        return "equations"
    if "parameters" in normalized:
        return "parameters"
    return "notes"

def render_williamson_details(section: Dict[str, object]) -> str:
    details = WILLIAMSON_DETAILS.get(str(section["slug"]))
    if not details:
        return ""

    blocks = []
    for _index, detail in sorted(enumerate(details), key=detail_sort_key):
        kind = detail_kind(str(detail["title"]))
        blocks.append(
            f"<section class='description-block detail-{html.escape(kind)}'>"
            f"<h3>{html.escape(detail['title'])}</h3>"
            f"<div class='description-copy'>{detail['body']}</div>"
            "</section>"
        )
    return f"<article class='experiment-description' aria-label='Experiment details'>{''.join(blocks)}</article>"

def render_media_card(media: Dict[str, object], slug: str) -> str:
    source = Path(media["source"])
    relative_path = ensure_docs_asset(source, slug)
    label = html.escape(str(media["label"]))
    kind = str(media["kind"])

    if kind == "video":
        body = (
            f"<video controls preload='metadata'>"
            f"<source src='{html.escape(relative_path)}' type='video/webm'>"
            "</video>"
        )
    else:
        body = f"<img src='{html.escape(relative_path)}' alt='{label}' />"

    return (
        "<div class='media-card'>"
        f"<h4>{label}</h4>"
        f"{body}"
        "</div>"
    )

def render_subfigure_gallery(title: str, items: List[Dict[str, object]], slug: str) -> str:
    figures = []
    for item in items:
        figure = render_plot_figure(item, slug)
        if figure:
            figures.append(figure)

    if not figures:
        return ""

    return (
        f"<details id='{html.escape(slug)}-errors' class='media-section results-section error-details'>"
        f"<summary>{html.escape(title)}</summary>"
        "<div class='results-grid subfigure-grid'>"
        f"{''.join(figures)}"
        "</div>"
        "</details>"
    )

def render_error_comparison_block(title: str, items: List[Dict[str, object]], slug: str) -> str:
    figures = []
    for item in items:
        figure = render_plot_figure(item, slug, "subfigure comparison-figure")
        if figure:
            figures.append(figure)
    if not figures:
        return ""
    return (
        f"<details id='{html.escape(slug)}-errors' class='media-section results-section error-block error-details'>"
        f"<summary>{html.escape(title)}</summary>"
        "<div class='results-grid comparison-grid'>"
        f"{''.join(figures)}"
        "</div>"
        "</details>"
    )

def render_dynamic_gallery(section: Dict[str, object], slug: str) -> str:
    gallery = section.get("dynamic_gallery")
    if gallery == "tc1":
        return render_error_comparison_block("Results: signed-error contours", tc1_error_contour_items(), slug)
    if gallery == "tc2":
        return "".join(
            [
                render_subfigure_gallery("Results: steady-state error norms", tc2_error_norm_items(), slug),
            ]
        )
    if gallery == "tc3":
        return render_subfigure_gallery("Results: steady-state drift norms", tc3_error_norm_items(), slug)
    if gallery == "tc6":
        return render_subfigure_gallery(
            "Results: vorticity, Q, energy, and mass diagnostics",
            tc6_postprocessing_items(),
            slug,
        )
    return ""

def render_diagnosis_data_links(paths: List[Path], root: Path, slug: str) -> str:
    if not paths:
        return ""

    links: list[str] = []
    for path in paths:
        relative_path = ensure_docs_asset(path, slug)
        label = diagnosis_asset_label(path, root)
        file_type = path.suffix.lower().lstrip(".").upper() or "DATA"
        links.append(
            "<a class='asset-link' "
            f"href='{html.escape(relative_path)}' target='_blank' rel='noopener'>"
            f"<span>{html.escape(label)}</span>"
            f"<small>{html.escape(file_type)}</small>"
            "</a>"
        )

    return (
        "<div class='diagnosis-data-list'>"
        "<h4>Ready data files</h4>"
        f"<div class='asset-link-grid'>{''.join(links)}</div>"
        "</div>"
    )

def render_case_diagnosis_assets(case: Dict[str, object], slug: str) -> str:
    root = case_diagnosis_root(case)
    if root is None or not root.exists():
        return ""

    images, data_paths = discover_case_diagnosis_assets(case)
    if not images and not data_paths:
        return ""

    label = case_display_label(case)
    panel_id = f"{slug}-{slug_token(label)}-diagnosis-assets"
    figures = "".join(
        render_plot_figure(
            {
                "source": path,
                "caption": f"{label}: {diagnosis_asset_label(path, root)}",
            },
            slug,
            "subfigure diagnosis-figure",
        )
        for path in images
    )
    figures_html = (
        "<div class='results-grid subfigure-grid diagnosis-grid'>"
        f"{figures}"
        "</div>"
        if figures
        else ""
    )
    data_html = render_diagnosis_data_links(data_paths, root, slug)
    counts = f"{len(images)} plot assets, {len(data_paths)} data files"

    return (
        f"<details id='{html.escape(panel_id)}' class='media-section results-section diagnosis-assets' open>"
        f"<summary>{html.escape(label)}: diagnosis assets and ready data</summary>"
        f"<p class='section-note'>{html.escape(counts)} mirrored into the published assets.</p>"
        f"{figures_html}"
        f"{data_html}"
        "</details>"
    )

def render_section(section: Dict[str, object]) -> str:
    slug = str(section["slug"])
    title = html.escape(str(section["title"]))
    section_label = "Williamson test" if slug.startswith("testcase") else "Experiment"
    case_configs = section_case_configs(section)
    metrics_html = render_metric_cards(section, case_configs)
    summary_html = render_experiment_summary(section)
    details_html = render_williamson_details(section)
    tables_html = "".join(
        render_case_tables(case, show_label=len(case_configs) > 1)
        for case in case_configs
    )

    media_items = [
        render_media_card(media, slug)
        for media in section.get("media", [])
        if Path(media["source"]).exists()
    ]
    media_blocks = []
    if media_items:
        media_blocks.append(f"<div class='main-panels'>{''.join(media_items)}</div>")
    for case in case_configs:
        snapshot_html = render_case_snapshot_galleries(case, slug)
        if snapshot_html:
            media_blocks.append(snapshot_html)
    dynamic_html = render_dynamic_gallery(section, slug)
    if dynamic_html:
        media_blocks.append(dynamic_html)
    for case in case_configs:
        diagnosis_html = render_case_diagnosis_assets(case, slug)
        if diagnosis_html:
            media_blocks.append(diagnosis_html)
    media_html = "".join(media_blocks) if media_blocks else "<p class='empty'>No media linked yet for this section.</p>"

    return (
        f"<section id='{html.escape(slug)}' class='section-card experiment-section status-{html.escape(status_for_section(section)[1])}'>"
        "<header class='experiment-top'>"
        "<div>"
        f"<span class='section-label'>{html.escape(section_label)}</span>"
        f"<h2>{title} {render_status_chip(section)}</h2>"
        f"{summary_html}"
        "</div>"
        f"{metrics_html}"
        "</header>"
        "<div class='experiment-body'>"
        f"{details_html}"
        f"<div class='settings-stack'>{tables_html}</div>"
        f"<div class='diagnostics-stack'>{media_html}</div>"
        "</div>"
        "</section>"
    )

def section_fragment_path(slug: str) -> Path:
    return FRAGMENT_DIR / f"{slug}.html"

def ordered_sections() -> List[Dict[str, object]]:
    return sorted(
        SECTIONS,
        key=lambda section: (
            0 if str(section["slug"]).startswith("testcase") else 1,
            natural_key(str(section["slug"])),
        ),
    )

def write_section_fragment(section: Dict[str, object]) -> Path:
    slug = str(section["slug"])
    FRAGMENT_DIR.mkdir(parents=True, exist_ok=True)
    path = section_fragment_path(slug)
    path.write_text(render_section(section) + "\n", encoding="utf-8")
    return path

def write_section_fragments(slugs: Optional[Union[List[str], Tuple[str, ...]]] = None) -> List[Path]:
    selected = set(slugs) if slugs is not None else None
    if selected is not None:
        known = {str(section["slug"]) for section in SECTIONS}
        unknown = sorted(selected - known)
        if unknown:
            raise ValueError(f"unknown GitHub Pages section slug(s): {', '.join(unknown)}")

    paths: list[Path] = []
    for section in ordered_sections():
        slug = str(section["slug"])
        if selected is None or slug in selected or not section_fragment_path(slug).exists():
            paths.append(write_section_fragment(section))
    return paths

def render_fragment_inputs() -> str:
    inputs: list[str] = []
    for section in ordered_sections():
        slug = str(section["slug"])
        path = section_fragment_path(slug)
        if path.exists():
            inputs.append(path.read_text(encoding="utf-8").rstrip())
    return "\n      ".join(inputs)


def parse_optional_float(value: object) -> Optional[float]:
    try:
        number = float(str(value))
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def format_metric_value(value: object) -> str:
    number = parse_optional_float(value)
    if number is None:
        return "n/a"
    if number == 0.0:
        return "0"
    abs_number = abs(number)
    if abs_number < 1.0e-3 or abs_number >= 1.0e4:
        return f"{number:.3e}"
    return f"{number:.6g}"


def render_cfl_audit_section() -> str:
    path = FRAGMENT_DIR / "cfl_deltaT_report.html"
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8").rstrip()


def read_conservation_rows() -> List[Dict[str, str]]:
    rows: list[dict[str, str]] = []
    for case_code, output_name in sorted(
        CASE_OUTPUT_NAMES.items(),
        key=lambda item: natural_key(item[0]),
    ):
        if case_code == "TC1_prime":
            continue
        summary = SANDBOX_OUTPUT_ROOT / output_name / "Diagnosis" / "conservation" / "summary.csv"
        if not summary.exists():
            continue
        with summary.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row = dict(row)
                row["case"] = case_code
                rows.append(row)
    return rows


def conservation_plot_items(case_code: str, alpha: str) -> List[Dict[str, object]]:
    output_name = CASE_OUTPUT_NAMES[case_code]
    root = SANDBOX_OUTPUT_ROOT / output_name / "Diagnosis" / "conservation" / f"alpha_{alpha}"
    labels = (
        ("conservation_mass.png", "mass/free-surface drift"),
        ("conservation_tracer.png", "tracer amount drift"),
        ("conservation_energy.png", "energy proxy drift"),
        ("conservation_pv_enstrophy.png", "PV/enstrophy drift"),
    )
    items: list[dict[str, object]] = []
    for filename, label in labels:
        path = root / filename
        if path.exists():
            items.append(
                {
                    "source": path,
                    "caption": f"{case_code} alpha {alpha}, {label}",
                }
            )
    return items


def render_conservation_plot_galleries(rows: List[Dict[str, str]]) -> str:
    case_blocks: list[str] = []
    for case_code in ("TC1", "TC2", "TC3", "TC4", "TC5", "TC6", "TC7"):
        case_rows = [
            row
            for row in rows
            if row.get("case") == case_code and row.get("status") == "available"
        ]
        figures: list[str] = []
        for row in case_rows:
            alpha = row.get("alpha", "")
            for item in conservation_plot_items(case_code, alpha):
                figure = render_plot_figure(item, "conservation-audit")
                if figure:
                    figures.append(figure)
        if not figures:
            continue
        case_blocks.append(
            "<details class='details-panel case-block'>"
            f"<summary>{html.escape(case_code)} conservation plots</summary>"
            "<div class='results-grid subfigure-grid audit-plot-grid'>"
            f"{''.join(figures)}"
            "</div>"
            "</details>"
        )
    return "".join(case_blocks)


def render_conservation_audit_section() -> str:
    fragment = FRAGMENT_DIR / "conservation_report.html"
    if fragment.exists():
        return fragment.read_text(encoding="utf-8").rstrip()

    rows = read_conservation_rows()
    if not rows:
        return ""

    header = (
        "<tr>"
        "<th>Case</th>"
        "<th>Run</th>"
        "<th>Health</th>"
        "<th>Last day</th>"
        "<th>Mass |rel|</th>"
        "<th>Target |rel|</th>"
        "<th>Energy |rel|</th>"
        "<th>PV enstrophy |rel|</th>"
        "<th>Verdict</th>"
        "</tr>"
    )
    body_rows: list[str] = []
    for row in rows:
        status = row.get("status", "")
        health = row.get("health_verdict", "") or ("unavailable" if status != "available" else "finite state fields")
        verdict = "; ".join(
            item
            for item in (
                f"mass {row.get('mass_verdict', 'n/a')}",
                f"target {row.get('quantity_verdict', 'n/a')}",
                f"energy {row.get('energy_verdict', 'n/a')}",
                f"enstrophy {row.get('enstrophy_verdict', 'n/a')}",
            )
            if item
        )
        if status != "available":
            reason = row.get("reason") or "no run output"
            verdict = f"unavailable: {reason}"
        body_rows.append(
            "<tr>"
            f"<td>{html.escape(row.get('case', ''))}</td>"
            f"<td>{html.escape(row.get('alpha', ''))}</td>"
            f"<td>{html.escape(health)}</td>"
            f"<td>{format_metric_value(row.get('last_day'))}</td>"
            f"<td>{format_metric_value(row.get('max_abs_mass_rel_change'))}</td>"
            f"<td>{format_metric_value(row.get('max_abs_tracer_mass_rel_change'))}</td>"
            f"<td>{format_metric_value(row.get('max_abs_mechanical_energy_rel_change'))}</td>"
            f"<td>{format_metric_value(row.get('max_abs_potential_enstrophy_rel_change'))}</td>"
            f"<td>{html.escape(verdict)}</td>"
            "</tr>"
        )

    plot_galleries = render_conservation_plot_galleries(rows)
    plots_html = plot_galleries or "<p class='empty'>No conservation plots generated yet.</p>"
    return (
        "<section id='conservation-audit' class='section-card experiment-section status-pending'>"
        "<header class='experiment-top'>"
        "<div>"
        "<span class='section-label'>Run checks</span>"
        "<h2>Conservation and Output-Health Checks <span class='status-chip status-pending'>Pending validation</span></h2>"
        "<p class='experiment-summary'>This section checks whether each run output stayed finite and whether mass, the tracked testcase quantity, mechanical energy, and derived potential enstrophy behave as expected.</p>"
        "</div>"
        "<div class='metric-grid'>"
        "<div class='metric-card'><span>Cases</span><strong>TC1-TC7</strong></div>"
        "<div class='metric-card'><span>Source</span><strong>MDS output</strong></div>"
        "<div class='metric-card'><span>State health</span><strong>Eta / U / V</strong></div>"
        "<div class='metric-card'><span>Plots</span><strong>Collapsed</strong></div>"
        "</div>"
        "</header>"
        "<div class='experiment-body'>"
        "<article class='experiment-description' aria-label='Conservation audit details'>"
        "<section class='description-block detail-equations'>"
        "<h3>Diagnostics</h3>"
        "<div class='description-copy'>"
        "<p>TC1's target quantity is passive-tracer amount. TC2, TC3, TC5, TC6, and TC7 are checked through mass, mechanical-energy proxy, and derived potential enstrophy. TC4 is forced, so its energy and enstrophy curves are diagnostic forced-response traces rather than unforced conservation invariants.</p>"
        "<p><strong>Derived potential enstrophy</strong> means the quantity is reconstructed in postprocessing from saved MITgcm fields, not read as a native model diagnostic.</p>"
        + math_block(
            "h=H+\\eta,\\qquad M=\\rho_0\\int_\\Omega h\\,dA",
            "E_m=\\rho_0\\int_\\Omega hK\\,dA+\\frac{1}{2}\\rho_0g\\int_\\Omega\\eta^2\\,dA",
            "\\zeta=\\frac{1}{a\\cos\\theta}\\left(\\frac{\\partial v}{\\partial\\lambda}-\\frac{\\partial(u\\cos\\theta)}{\\partial\\theta}\\right)",
            "q=\\frac{\\zeta+f}{h},\\qquad Z=\\frac{1}{2}\\int_\\Omega hq^2\\,dA",
            "\\Delta_r X(t)=\\frac{X(t)-X(0)}{|X(0)|}",
        )
        + "<p>Here <code>K</code> is native <code>momKE</code> when available, otherwise <code>(u^2+v^2)/2</code>. The Coriolis term <code>f</code> comes from <code>fCoriC.bin</code> when a rotated run writes it; otherwise the diagnostic uses <code>2*omega*sin(latitude)</code>.</p>"
        + "</div>"
        "</section>"
        "<section class='description-block detail-expected'>"
        "<h3>Decision Rule</h3>"
        "<div class='description-copy'>"
        "<p>Use these conservation values only when the health column reports finite state fields. TC5 and TC7 outputs that become non-finite are invalid even when files exist; the corrected reruns use smaller steps, viscosity, and stricter postprocessing health checks. TC2 and TC3 rotated runs must be regenerated after the Coriolis-map fix before their conservation rows can be final.</p>"
        "</div>"
        "</section>"
        "</article>"
        "<div class='settings-stack'>"
        "<div class='table-block data-panel'>"
        "<h3>Conservation Summary</h3>"
        f"<div class='table-scroll'><table><thead>{header}</thead><tbody>{''.join(body_rows)}</tbody></table></div>"
        "</div>"
        "</div>"
        f"<div class='diagnostics-stack'>{plots_html}</div>"
        "</div>"
        "</section>"
    )


def render_audit_sections() -> str:
    return "\n      ".join(
        section
        for section in (
            render_cfl_audit_section(),
            render_conservation_audit_section(),
        )
        if section
    )


def nav_label(section: Dict[str, object]) -> str:
    slug = str(section["slug"])
    return NAV_LABELS.get(slug, str(section["title"]))

def status_for_section(section: Dict[str, object]) -> Tuple[str, str]:
    return STATUS_META.get(str(section["slug"]), ("Pending validation", "pending"))

def render_status_chip(section: Dict[str, object]) -> str:
    label, key = status_for_section(section)
    return (
        f"<span class='status-chip status-{html.escape(key)}'>"
        f"{html.escape(label)}</span>"
    )

def render_nav_status(section: Dict[str, object]) -> str:
    label, key = status_for_section(section)
    return (
        f"<span class='nav-status status-{html.escape(key)}'>"
        f"{html.escape(label)}</span>"
    )

def render_homepage() -> str:
    williamson_links = []
    initialization_links = []
    for section in ordered_sections():
        slug = html.escape(str(section["slug"]))
        label = html.escape(nav_label(section))
        group = "Williamson testcase" if str(section["slug"]).startswith("testcase") else "Initialization experiment"
        status = render_nav_status(section)
        link = (
            "<a class='home-experiment-link' "
            f"href='#{slug}'><span>{html.escape(group)}</span><strong>{label}</strong>{status}</a>"
        )
        if str(section["slug"]).startswith("testcase"):
            williamson_links.append(link)
        else:
            initialization_links.append(link)

    return (
        "<section id='home' class='home-section hero section-card'>"
        "<div class='home-copy'>"
        "<span class='section-label'>MITgcm verification website</span>"
        "<h2>Shallow-water experiments on the sphere</h2>"
        "<p>This site presents MITgcm shallow-water verification experiments, including the "
        "Williamson standard test cases and free-surface initialization runs over constant "
        "and real bathymetry. Each section collects the governing equations, expected output, "
        "numerical settings, snapshots, and available result diagnostics.</p>"
        f"<p class='author-credit'>Code and website by <strong>{html.escape(AUTHOR_NAME)}</strong>.</p>"
        "<div class='home-actions'>"
        "<a class='btn-primary' href='#testcase1'>Start with TC1</a>"
        "<a class='btn-secondary' href='#cfl-deltat-audit'>View run checks</a>"
        "<a class='btn-secondary' href='#case1_constant_bathymetry'>View initialization runs</a>"
        f"<a class='btn-secondary' href='{html.escape(REPO_MAIN_URL)}' target='_blank' rel='noopener'>View code on GitHub</a>"
        "</div>"
        "</div>"
        "<div class='home-experiment-list'>"
        "<div class='home-list-group'>"
        "<h3>Williamson tests</h3>"
        f"{''.join(williamson_links)}"
        "</div>"
        "<div class='home-list-group'>"
        "<h3>Initialization runs</h3>"
        f"{''.join(initialization_links)}"
        "</div>"
        "</div>"
        "</section>"
    )

def render_nav_group(title: str, sections: List[Dict[str, object]], open_group: bool = True) -> str:
    links = []
    for section in sections:
        slug = html.escape(str(section["slug"]))
        title_html = html.escape(nav_label(section))
        links.append(
            f"        <a class='nav-link' href='#{slug}'>{title_html}{render_nav_status(section)}</a>"
        )
    open_attr = " open" if open_group else ""
    return "\n".join(
        [
            f"      <details class='nav-group'{open_attr}>",
            f"        <summary>{html.escape(title)}</summary>",
            "        <div class='nav-group-links'>",
            *links,
            "        </div>",
            "      </details>",
        ]
    )


def render_audit_nav_group() -> str:
    links = []
    for slug, title, _group, status in AUDIT_SECTIONS:
        links.append(
            f"        <a class='nav-link' href='#{html.escape(slug)}'>"
            f"{html.escape(title)}"
            f"<span class='nav-status status-pending'>{html.escape(status)}</span></a>"
        )
    return "\n".join(
        [
            "      <details class='nav-group' open>",
            "        <summary>Run checks</summary>",
            "        <div class='nav-group-links'>",
            *links,
            "        </div>",
            "      </details>",
        ]
    )


def render_navigation() -> str:
    sections = ordered_sections()
    williamson_sections = [section for section in sections if str(section["slug"]).startswith("testcase")]
    experiment_sections = [section for section in sections if not str(section["slug"]).startswith("testcase")]
    return "\n".join(
        [
            "<aside id='site-navigation' class='side-nav' aria-label='Experiment navigation' data-nav-drawer>",
            "      <div class='nav-heading'>",
            "        <span>Navigation</span>",
            "        <button class='nav-close' type='button' data-nav-close aria-label='Close navigation'>x</button>",
            "      </div>",
            "      <a class='nav-link nav-home' href='#home'>Overview</a>",
            render_audit_nav_group(),
            render_nav_group("Williamson tests", williamson_sections, True),
            render_nav_group("Initialization runs", experiment_sections, True),
            "    </aside>",
        ]
    )

def build_html() -> str:
    nav_html = render_navigation()
    audit_sections = render_audit_sections()
    sections_html = render_homepage()
    if audit_sections:
        sections_html = f"{sections_html}\n      {audit_sections}"
    sections_html = f"{sections_html}\n      {render_fragment_inputs()}"
    return f"""---
---
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{HTML_TITLE}</title>
  <link rel="icon" href="{FAVICON_HREF}">
  <link rel="stylesheet" href="site.css">
  <script>
    window.MathJax = {{
      tex: {{
        inlineMath: [['\\\\(', '\\\\)']],
        displayMath: [['\\\\[', '\\\\]']],
        processEscapes: true,
        tags: 'ams'
      }},
      options: {{
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
      }},
      svg: {{
        fontCache: 'global',
        scale: 0.96
      }}
    }};
  </script>
  <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
</head>
<body>
  <header class="site-header">
    <div class="site-header-inner">
      <div>
        <span class="site-kicker">MITgcm verification dashboard</span>
        <h1>{SITE_TITLE}</h1>
        <p class="intro">{html.escape(SITE_SUBTITLE)}</p>
        <p class="site-credit">Code and website by {html.escape(AUTHOR_NAME)} &middot; <a href="{html.escape(REPO_MAIN_URL)}" target="_blank" rel="noopener">GitHub main branch</a></p>
      </div>
      <button class="nav-toggle" type="button" data-nav-toggle aria-controls="site-navigation" aria-expanded="false" aria-label="Open navigation">
        Menu
      </button>
    </div>
  </header>
  <button class="nav-backdrop" type="button" data-nav-backdrop aria-label="Close navigation"></button>
  <div class="site-layout">
    {nav_html}
    <main>
      {sections_html}
    </main>
  </div>
  <div class="plot-modal" data-modal hidden>
    <button class="modal-close" type="button" data-modal-close aria-label="Close enlarged plot">x</button>
    <figure>
      <img data-modal-img alt="">
      <figcaption data-modal-caption></figcaption>
    </figure>
  </div>
  <script src="site.js" defer></script>
</body>
</html>
"""

def build_site(section_slugs: Optional[Union[List[str], Tuple[str, ...]]] = None) -> Tuple[Path, List[Path]]:
    COPIED_ASSETS.clear()
    MISSING_SNAPSHOT_ROOTS.clear()
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_ROOT.mkdir(parents=True, exist_ok=True)
    FRAGMENT_DIR.mkdir(parents=True, exist_ok=True)
    sync_ready_output_assets()
    fragment_paths = write_section_fragments(section_slugs)
    index_path = DOCS_DIR / "index.html"
    index_path.write_text(build_html(), encoding="utf-8")
    return index_path, fragment_paths

def main() -> None:
    index_path, fragment_paths = build_site()
    print(f"wrote {index_path}")
    for path in fragment_paths:
        print(f"wrote {path}")
    print(f"copied/updated {len(COPIED_ASSETS)} asset files")
    if MISSING_SNAPSHOT_ROOTS:
        print("missing snapshot roots:")
        for path in sorted(MISSING_SNAPSHOT_ROOTS, key=lambda item: item.as_posix()):
            print(f"  {path.relative_to(REPO_DIR)}")

if __name__ == "__main__":
    main()

# %%
