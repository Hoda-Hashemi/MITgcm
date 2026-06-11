"""Save MITgcm single heatmap snapshots for Eta, |U,V|, PsiVEL, and PhiVEL."""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path
import numpy as np
SCRIPT_DIR = Path(__file__).resolve().parent
ANIMATION_DIR = SCRIPT_DIR.parent / "animation"
if str(ANIMATION_DIR) not in sys.path:
    sys.path.insert(0, str(ANIMATION_DIR))
import animationComponents as ac  # noqa: E402
OUTPUT_FORMAT = "pdf"
SINGLE_HEATMAP_DPI = 360
CMAP = "RdBu_r"

DEFAULT_RUN_DIR = SCRIPT_DIR.parents[2] / "run"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Save MITgcm single heatmap snapshots.")
    parser.add_argument(
        "run_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_RUN_DIR,
        help="MITgcm run directory containing output files.",
    )
    return parser.parse_args(argv)
def read_delta_t(run_dir: Path) -> float:
    for data_file in (run_dir / "data", run_dir.parent / "input" / "data"):
        if data_file.exists():
            match = re.search(r"deltaT\s*=\s*([-+0-9.eEdD]+)", data_file.read_text(errors="ignore"))
            if match:
                return float(match.group(1).replace("D", "E").replace("d", "e"))
    return 60.0
def snapshot_days_from_iterations(iterations: list[int], delta_t_sec: float) -> list[float]:
    return [iteration * delta_t_sec / 86400.0 for iteration in iterations]
def constant_scale(fields: dict[str, list[np.ndarray]]) -> dict[str, str]:
    return {name: "constant" for name in fields}
def save_field(fields: dict[str, list[np.ndarray]], iterations: list[int], delta_t_sec: float, output_dir: Path) -> None:
    snapshots = ac.pick_snapshot_iterations(iterations, delta_t_sec, snapshot_days_from_iterations(iterations, delta_t_sec))
    ac.write_single_heatmap_snapshots(fields, iterations, snapshots, output_dir, CMAP, dpi=SINGLE_HEATMAP_DPI, scale_mode_by_field=constant_scale(fields), file_ext=OUTPUT_FORMAT)
def write_eta(run_dir: Path, output_root: Path, delta_t_sec: float) -> None:
    timed = ac.discover_timed_variables(run_dir)
    iterations = ac.common_iterations(timed, ["Eta"])
    eta0, _ = ac.read_selected_field(run_dir, "Eta", iterations[0])
    land_mask = ac.load_land_mask(run_dir, *eta0.shape)
    eta_series = ac.read_variable_series(run_dir, "Eta", iterations, land_mask)
    save_field({"Eta": eta_series}, iterations, delta_t_sec, output_root / "Eta")
def write_velocity_magnitude(run_dir: Path, output_root: Path, delta_t_sec: float) -> None:
    timed = ac.discover_timed_variables(run_dir)
    iterations = ac.common_iterations(timed, ["U", "V"])
    u0, _ = ac.read_selected_field(run_dir, "U", iterations[0])
    land_mask = ac.load_land_mask(run_dir, *u0.shape)
    u_series = ac.read_variable_series(run_dir, "U", iterations, land_mask)
    v_series = ac.read_variable_series(run_dir, "V", iterations, land_mask)
    velocity_magnitude = [np.sqrt(u * u + v * v) for u, v in zip(u_series, v_series, strict=True)]
    save_field({"|U,V|": velocity_magnitude}, iterations, delta_t_sec, output_root / "VelocityMagnitude")
def write_psi(run_dir: Path, output_root: Path, delta_t_sec: float) -> None:
    iterations, _ = ac.common_iterations_for_records(run_dir, ["PsiVEL"])
    psi0, _ = ac.read_record_any_bundle(run_dir, "PsiVEL", iterations[0])
    land_mask = ac.load_land_mask(run_dir, *psi0.shape)
    psi_series, _ = ac.read_record_series_any_bundle(run_dir, "PsiVEL", iterations, land_mask)
    save_field({"PsiVEL": psi_series}, iterations, delta_t_sec, output_root / "Psi")
def write_phi(run_dir: Path, output_root: Path, delta_t_sec: float) -> None:
    iterations, _ = ac.common_iterations_for_records(run_dir, ["PhiVEL"])
    phi0, _ = ac.read_record_any_bundle(run_dir, "PhiVEL", iterations[0])
    land_mask = ac.load_land_mask(run_dir, *phi0.shape)
    phi_series, _ = ac.read_record_series_any_bundle(run_dir, "PhiVEL", iterations, land_mask)
    save_field({"PhiVEL": phi_series}, iterations, delta_t_sec, output_root / "Phi")
def run_snapshots(run_dir: Path | None = None) -> None:
    run_dir = DEFAULT_RUN_DIR if run_dir is None else Path(run_dir)
    run_dir = run_dir.expanduser().resolve()
    delta_t_sec = read_delta_t(run_dir)
    output_root = SCRIPT_DIR / "output" / run_dir.name
    output_root.mkdir(parents=True, exist_ok=True)
    write_eta(run_dir, output_root, delta_t_sec)
    write_velocity_magnitude(run_dir, output_root, delta_t_sec)
    write_psi(run_dir, output_root, delta_t_sec)
    write_phi(run_dir, output_root, delta_t_sec)


if __name__ == "__main__":
    run_snapshots(parse_args().run_dir)
