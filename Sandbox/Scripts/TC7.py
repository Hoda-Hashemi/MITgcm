#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from github_pages import build_site
from run_log import write_run_log
from snapshot_plots import SnapshotSpec, run_snapshots

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC7"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC7"

RUN_DIRS = [
    Path(item).expanduser()
    for item in os.environ.get("TC7_RUN_DIRS", "").split(os.pathsep)
    if item
] or [CASE_DIR / "run_analysis"]

SNAPSHOT_FIELDS = (
    SnapshotSpec("Eta", "eta", "eta", "m", center_zero=True),
    SnapshotSpec("ETAN", "etan", "ETAN", "m", center_zero=True),
    SnapshotSpec("PsiVEL", "psi", "PsiVEL", r"m$^3$ s$^{-1}$", center_zero=True),
    SnapshotSpec("PhiVEL", "phi", "PhiVEL", r"m$^2$ s$^{-1}$", center_zero=True),
)


def main() -> None:
    for run_dir in RUN_DIRS:
        run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
    build_site(["testcase7"])
    for run_dir in RUN_DIRS:
        write_run_log(run_dir, case_code=CASE_CODE)


if __name__ == "__main__":
    main()
