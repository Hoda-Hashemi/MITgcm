#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from github_pages import build_site
from postprocessing_ErrorAnalysis import ErrorAnalysisSpec, ErrorFieldSpec, analyze_error
from postprocessing_Quantities import PostprocessingSpec, analyze_run
from run_log import write_run_log
from snapshot_plots import SnapshotSpec, run_snapshots

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC3"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC3"

RUN_DIRS = [
    Path(item).expanduser()
    for item in os.environ.get("TC3_RUN_DIRS", "").split(os.pathsep)
    if item
] or [CASE_DIR / "run_alpha_0"]
MAKE_SNAPSHOTS = True
MAKE_ERROR_ANALYSIS = True
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

TC3_ERROR_SPEC = ErrorAnalysisSpec(
    case_code=CASE_CODE,
    fields=(
        ErrorFieldSpec("eta", ("Eta", "ETAN"), None, label="eta", normalize=True, reference_kind="initial"),
        ErrorFieldSpec("u", ("U", "UVEL", "UVELMASS"), None, label="u", normalize=True, reference_kind="initial"),
        ErrorFieldSpec("v", ("V", "VVEL", "VVELMASS"), None, label="v", normalize=False, reference_kind="initial"),
    ),
    title="Williamson TC3 compact-support steady-state drift",
    table_title="Williamson TC3 compact-support steady-state drift",
    plot_stem="TC3_error_norms",
    table_stem="TC3_error_norms",
    log_y=True,
)


def main() -> None:
    for run_dir in RUN_DIRS:
        if MAKE_SNAPSHOTS:
            run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        if MAKE_ERROR_ANALYSIS:
            analyze_error(run_dir, TC3_ERROR_SPEC)
        if MAKE_POSTPROCESSING:
            analyze_run(run_dir, spec=POSTPROCESSING_SPEC)
    build_site(["testcase3"])
    for run_dir in RUN_DIRS:
        write_run_log(run_dir, case_code=CASE_CODE)


if __name__ == "__main__":
    main()
