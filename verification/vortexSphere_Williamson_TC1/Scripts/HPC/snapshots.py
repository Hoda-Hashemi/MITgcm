#!/usr/bin/env python3
"""Save TC1 snapshots under this case while reusing the reference plotter."""

from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
REFERENCE_SCRIPT = (
    SCRIPT_PATH.parents[3]
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

# The reused plotter builds its output path relative to this value.
reference_snapshots.SCRIPT_DIR = SCRIPT_PATH.parents[1]

if __name__ == "__main__":
    reference_snapshots.main()
