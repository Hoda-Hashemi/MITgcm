#!/usr/bin/env python3
from pathlib import Path
import sys

scripts_dir = Path(__file__).resolve().parents[2] / "Scripts"
sys.path.insert(0, str(scripts_dir))

from williamson_setup_checks import main

if __name__ == "__main__":
    raise SystemExit(main("TC6"))
