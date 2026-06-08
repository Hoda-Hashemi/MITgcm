#!/usr/bin/env python3
"""Save passive-height and velocity snapshots for Williamson TC1."""
#%%
from __future__ import annotations

import importlib.util
from pathlib import Path

# %%

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
CASE_DIR = SCRIPT_DIR.parents[1]
RUN_DIR = CASE_DIR / "run"
OUTPUT_DIR = SCRIPT_DIR.parent / "output" / RUN_DIR.name
REFERENCE_SCRIPT = (
    CASE_DIR.parent
    / "vortexSphere_mitgcm_referenceCase"
    / "Scripts"
    / "HPC"
    / "snapshots.py"
)

spec = importlib.util.spec_from_file_location("reference_snapshots", REFERENCE_SCRIPT)
if spec is None or spec.loader is None:
    raise ImportError(f"Cannot load reference plotter: {REFERENCE_SCRIPT}")

reference_snapshots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference_snapshots)

def run_snapshots(run_dir: Path | None = None) -> None:
    run_dir = RUN_DIR if run_dir is None else Path(run_dir)
    run_dir = run_dir.expanduser().resolve()
    delta_t_sec = reference_snapshots.read_delta_t(run_dir)
    output_root = SCRIPT_DIR.parent / "output" / run_dir.name
    output_root.mkdir(parents=True, exist_ok=True)

    timed = reference_snapshots.ac.discover_timed_variables(run_dir)
    iterations = reference_snapshots.ac.common_iterations(timed, ["S"])
    salt0, _ = reference_snapshots.ac.read_selected_field(run_dir, "S", iterations[0])
    if salt0.ndim == 3:
        salt0 = salt0[0]
    land_mask = reference_snapshots.ac.load_land_mask(run_dir, *salt0.shape)
    height_series = reference_snapshots.ac.read_variable_series(
        run_dir, "S", iterations, land_mask
    )
    height_series = [
        field[0] if field.ndim == 3 else field for field in height_series
    ]
    reference_snapshots.save_field(
        {"TC1 height": height_series},
        iterations,
        delta_t_sec,
        output_root / "Height",
    )
    reference_snapshots.write_velocity_magnitude(
        run_dir, output_root, delta_t_sec
    )

if __name__ == "__main__":
    run_snapshots(RUN_DIR)

# %%
