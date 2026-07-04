#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from github_pages import build_site
from postprocessing_Quantities import PostprocessingSpec, analyze_run
from run_log import write_run_log
from mitgcm_io import lon_lat, read_data_value, read_mds_field
from shared import infer_alpha_label, snapshot_output_dir
from snapshot_plots import SnapshotSpec, plot_cube_net_snapshot, run_snapshots

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC5"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC5"

RUN_DIRS = [
    Path(item).expanduser()
    for item in os.environ.get("TC5_RUN_DIRS", "").split(os.pathsep)
    if item
] or [CASE_DIR / "run_standard_cs32"]
MAKE_POSTPROCESSING = True

POSTPROCESSING_SPEC = PostprocessingSpec(
    case_code=CASE_CODE,
    eta_candidates=("Eta", "ETAN"),
    u_candidates=("U", "UVEL", "UVELMASS"),
    v_candidates=("V", "VVEL", "VVELMASS"),
    kinetic_energy_candidates=(),
    vorticity_candidates=(),
    compute_kinetic_energy_if_missing=True,
    compute_vorticity_if_missing=True,
    compute_potential_vorticity_if_missing=True,
)

SNAPSHOT_FIELDS = (
    SnapshotSpec("Eta", "eta", "eta", "m", center_zero=True),
    SnapshotSpec("ETAN", "etan", "ETAN", "m", center_zero=True),
    SnapshotSpec("PsiVEL", "psi", "PsiVEL", r"m$^3$ s$^{-1}$", center_zero=True),
    SnapshotSpec("PhiVEL", "phi", "PhiVEL", r"m$^2$ s$^{-1}$", center_zero=True),
)


def save_mountain_plot(run_dir: Path) -> None:
    alpha = infer_alpha_label(run_dir, CASE_CODE)
    output_dir = snapshot_output_dir(CASE_CODE, alpha) / "mountain"
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in output_dir.glob("*.png"):
        path.unlink()

    xc, yc = lon_lat(run_dir)
    depth = read_mds_field(run_dir, "Depth")
    h0 = read_data_value(run_dir, "delR", 5960.0)
    mountain = h0 - depth
    title = (
        "TC5 isolated mountain | static bathymetry | "
        f"peak {float(mountain.max()):.1f} m at 270E, 30N"
    )
    out = output_dir / "tc5_mountain_day_00.00_iter_0000000000.pdf"
    plot_cube_net_snapshot(mountain, xc, yc, title, "m", out, vmin=0.0, vmax=2000.0)


def main() -> None:
    for run_dir in RUN_DIRS:
        run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        save_mountain_plot(run_dir)
        if MAKE_POSTPROCESSING:
            analyze_run(run_dir, spec=POSTPROCESSING_SPEC)
    build_site(["testcase5"])
    for run_dir in RUN_DIRS:
        write_run_log(run_dir, case_code=CASE_CODE)


if __name__ == "__main__":
    main()
