#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from snapshot_plots import SnapshotSpec, run_snapshots
from williamson_errors import analyze_tc2_error

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC2"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC2"

RUN_DIRS = [CASE_DIR / "run_alpha_0"]
MAKE_SNAPSHOTS = True
MAKE_ERROR_ANALYSIS = True

SNAPSHOT_FIELDS = (
    SnapshotSpec("Eta", "eta", "eta", "m", -1905.0, 1905.0, center_zero=True),
    SnapshotSpec("ETAN", "etan", "ETAN", "m", center_zero=True),
    SnapshotSpec("PsiVEL", "psi", "PsiVEL", r"m$^3$ s$^{-1}$", center_zero=True),
    SnapshotSpec("PhiVEL", "phi", "PhiVEL", r"m$^2$ s$^{-1}$", center_zero=True),
)


def main() -> None:
    for run_dir in RUN_DIRS:
        if MAKE_SNAPSHOTS:
            run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        if MAKE_ERROR_ANALYSIS:
            analyze_tc2_error(run_dir)


if __name__ == "__main__":
    main()
