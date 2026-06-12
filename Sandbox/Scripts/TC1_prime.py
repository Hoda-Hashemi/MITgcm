#!/usr/bin/env python3
#%%
from __future__ import annotations

from pathlib import Path

from conservation import ConservationSpec, analyze_conservation
from snapshot_plots import SnapshotSpec, run_snapshots
from williamson_errors import analyze_tc1_error

SCRIPT_DIR = Path(__file__).resolve().parent
CASE_CODE = "TC1_prime"
CASE_DIR = SCRIPT_DIR.parent / "vortexSphere_Williamson_TC1_prime"

RUN_DIRS = [
    CASE_DIR / "run_alpha_0",
    CASE_DIR / "run_alpha_1.57",
]
MAKE_SNAPSHOTS = True
MAKE_ERROR_ANALYSIS = True
MAKE_CONSERVATION = True

CONSERVATION_SPEC = ConservationSpec(
    case_code=CASE_CODE,
    output_kind="conservation",
    report_title="TC1 prime conservation diagnostics",
    mass_stem="tc1_prime_mass_conservation",
    mass_title="TC1 prime mass",
    mass_ylabel=r"$M(t)\sim\int (H+\eta)\,dA$ [m$^3$]",
    energy_stem="tc1_prime_mechanical_energy",
    energy_title="TC1 prime mechanical energy",
    energy_ylabel="Mechanical energy [J]",
    table_stem="tc1_prime_conservation_table",
    table_title="Williamson TC1 prime conservation diagnostics",
    mass_column_label="M[m3]",
    kinetic_column_label="KE[J]",
    potential_column_label="PE[J]",
    energy_column_label="ME[J]",
    eta_candidates=("ETAN", "Eta"),
    u_candidates=("UVEL", "U", "UVELMASS"),
    v_candidates=("VVEL", "V", "VVELMASS"),
)

SNAPSHOT_FIELDS = (
    SnapshotSpec("Eta", "eta", "free-surface height", "m", 0.0, 1000.0, center_zero=True),
    SnapshotSpec("ETAN", "etan", "ETAN", "m", center_zero=True),
    SnapshotSpec("PsiVEL", "psi", "PsiVEL", r"m$^3$ s$^{-1}$", center_zero=True),
    SnapshotSpec("PhiVEL", "phi", "PhiVEL", r"m$^2$ s$^{-1}$", center_zero=True),
)

def main() -> None:
    for run_dir in RUN_DIRS:
        if MAKE_SNAPSHOTS:
            run_snapshots(CASE_CODE, run_dir, SNAPSHOT_FIELDS)
        if MAKE_ERROR_ANALYSIS:
            analyze_tc1_error(run_dir, case_code=CASE_CODE, field="ETAN")
        if MAKE_CONSERVATION:
            analyze_conservation(run_dir, CONSERVATION_SPEC)

if __name__ == "__main__":
    main()

# %%
