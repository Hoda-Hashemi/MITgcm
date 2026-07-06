#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from github_pages import build_site
from postprocessing_Quantities import PostprocessingSpec, analyze_run
from run_log import write_run_log
from snapshot_plots import SnapshotSpec, run_snapshots

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC4"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC4"

RUN_DIRS = [
    Path(item).expanduser()
    for item in os.environ.get("TC4_RUN_DIRS", "").split(os.pathsep)
    if item
] or [CASE_DIR / "run_u0_20", CASE_DIR / "run_u0_40"]
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


def main() -> None:
    for run_dir in RUN_DIRS:
        run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        if MAKE_POSTPROCESSING:
            analyze_run(run_dir, spec=POSTPROCESSING_SPEC)
    build_site(["testcase4"])
    for run_dir in RUN_DIRS:
        write_run_log(run_dir, case_code=CASE_CODE)


if __name__ == "__main__":
    main()
