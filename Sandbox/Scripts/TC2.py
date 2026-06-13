#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from postprocessing_ErrorAnalysis import DAY, ErrorAnalysisSpec, ErrorFieldSpec, analyze_error
from run_log import write_run_log
from snapshot_plots import SnapshotSpec, run_snapshots

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC2"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC2"

TC2_R_EARTH = 6.371e6
TC2_OMEGA = 7.292e-5
TC2_G = 9.81
TC2_U0 = 2.0 * math.pi * TC2_R_EARTH / (12.0 * DAY)
TC2_GH0 = 2.94e4
TC2_H0 = TC2_GH0 / TC2_G

RUN_DIRS = [CASE_DIR / "run_alpha_0"]
MAKE_SNAPSHOTS = True
MAKE_ERROR_ANALYSIS = True

SNAPSHOT_FIELDS = (
    SnapshotSpec("Eta", "eta", "eta", "m", -1905.0, 1905.0, center_zero=True),
    SnapshotSpec("ETAN", "etan", "ETAN", "m", center_zero=True),
    SnapshotSpec("PsiVEL", "psi", "PsiVEL", r"m$^3$ s$^{-1}$", center_zero=True),
    SnapshotSpec("PhiVEL", "phi", "PhiVEL", r"m$^2$ s$^{-1}$", center_zero=True),
)


def exact_tc2_eta(xc: np.ndarray, yc: np.ndarray, t_sec: float, alpha: float) -> np.ndarray:
    lon = np.deg2rad(xc)
    lat = np.deg2rad(yc)
    mu = -np.cos(lon) * np.cos(lat) * np.sin(alpha) + np.sin(lat) * np.cos(alpha)
    gh = TC2_GH0 - (TC2_R_EARTH * TC2_OMEGA * TC2_U0 + 0.5 * TC2_U0**2) * mu**2
    return gh / TC2_G - TC2_H0


def exact_tc2_u(xc: np.ndarray, yc: np.ndarray, t_sec: float, alpha: float) -> np.ndarray:
    lon = np.deg2rad(xc)
    lat = np.deg2rad(yc)
    return TC2_U0 * (np.cos(lat) * np.cos(alpha) + np.cos(lon) * np.sin(lat) * np.sin(alpha))


def exact_tc2_v(xc: np.ndarray, yc: np.ndarray, t_sec: float, alpha: float) -> np.ndarray:
    return -TC2_U0 * np.sin(np.deg2rad(xc)) * np.sin(alpha)


TC2_ERROR_SPEC = ErrorAnalysisSpec(
    case_code=CASE_CODE,
    fields=(
        ErrorFieldSpec(
            name="eta",
            label="eta",
            candidates=("Eta", "ETAN"),
            reference=exact_tc2_eta,
            normalize=True,
        ),
        ErrorFieldSpec(
            name="u",
            label="u",
            candidates=("U", "UVEL", "UVELMASS"),
            reference=exact_tc2_u,
            normalize=True,
        ),
        ErrorFieldSpec(
            name="v",
            label="v abs",
            candidates=("V", "VVEL", "VVELMASS"),
            reference=exact_tc2_v,
            normalize=False,
        ),
    ),
    title="Williamson TC2 steady-state error norms",
    table_title="Williamson TC2 steady-state error norms",
    plot_stem="TC2_error_norms",
    table_stem="TC2_error_norms",
    log_y=True,
)


def main() -> None:
    for run_dir in RUN_DIRS:
        if MAKE_SNAPSHOTS:
            run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        if MAKE_ERROR_ANALYSIS:
            analyze_error(run_dir, TC2_ERROR_SPEC)
    for run_dir in RUN_DIRS:
        write_run_log(run_dir, case_code=CASE_CODE)


if __name__ == "__main__":
    main()
