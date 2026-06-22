#!/usr/bin/env python3
"""Initial-condition CFL preflight report for Williamson TC2.

This is a read-only helper. It does not compile, submit jobs, or write model
inputs. For TC2 the velocity evolves, so this reports the CFL of the analytic
initial velocity field and the DELTA_T values currently declared in the jobs.
"""

import argparse
import ast
import importlib.util
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

try:
    import numpy as np
except ModuleNotFoundError as exc:
    raise SystemExit(
        "ERROR: numpy is required. Activate the MITgcm virtual environment or "
        "run with /home/hmh85/scratch/MITgcm/.venv/bin/python."
    ) from exc


WARN_THRESHOLD = 0.5
FAIL_THRESHOLD = 1.0


class GridConfig(NamedTuple):
    nx: int
    ny: int
    dlon_deg: float
    dlat_deg: float
    r_sphere: float
    template_delta_t: float


class JobCase(NamedTuple):
    name: str
    alpha: float
    delta_t: float
    path: Path


class CflReport(NamedTuple):
    case: JobCase
    status: str
    max_x: float
    max_y: float
    max_total: float
    location_lon: float
    location_lat: float
    suggested_delta_t: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Report TC2 initial-condition advective CFL before compiling or running MITgcm."
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="cases",
        metavar="NAME",
        help="case name to check, e.g. c1. Can be passed more than once.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="check all job_tc2_*.slurm cases. This is the default.",
    )
    parser.add_argument(
        "--warn-threshold",
        type=float,
        default=WARN_THRESHOLD,
        help="CFL value above which status becomes WARN; default {}.".format(WARN_THRESHOLD),
    )
    parser.add_argument(
        "--fail-threshold",
        type=float,
        default=FAIL_THRESHOLD,
        help="CFL value above which status becomes FAIL; default {}.".format(FAIL_THRESHOLD),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="exit nonzero if any selected case is WARN or FAIL.",
    )
    return parser.parse_args()


def safe_eval_number(expr: str) -> float:
    expr = expr.strip().rstrip(",").replace("D", "E").replace("d", "e")
    tree = ast.parse(expr, mode="eval")

    def eval_node(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return eval_node(node.body)
        if isinstance(node, ast.Num):
            return float(node.n)
        if hasattr(ast, "Constant") and isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = eval_node(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp) and isinstance(
            node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)
        ):
            left = eval_node(node.left)
            right = eval_node(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left / right
        raise ValueError("unsupported numeric expression: {}".format(expr))

    return eval_node(tree)


def parse_repeat_grid(value: str, name: str) -> Tuple[int, float]:
    value = value.strip().rstrip(",")
    match = re.fullmatch(
        r"(\d+)\s*\*\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)",
        value,
    )
    if not match:
        raise ValueError("{} must use MITgcm repeat syntax like 1440*0.25".format(name))
    return int(match.group(1)), safe_eval_number(match.group(2))


def parse_data_file(path: Path) -> GridConfig:
    values: Dict[str, str] = {}
    for line in path.read_text().splitlines():
        clean = line.split("#", 1)[0].strip()
        if not clean or "=" not in clean:
            continue
        key, raw_value = clean.split("=", 1)
        key = key.strip()
        raw_value = raw_value.strip().rstrip(",")
        if key in {"deltaT", "rSphere", "delX", "delY"}:
            values[key] = raw_value

    missing = sorted({"deltaT", "rSphere", "delX", "delY"} - set(values))
    if missing:
        raise ValueError("{} is missing: {}".format(path, ", ".join(missing)))

    nx, dlon_deg = parse_repeat_grid(values["delX"], "delX")
    ny, dlat_deg = parse_repeat_grid(values["delY"], "delY")
    return GridConfig(
        nx=nx,
        ny=ny,
        dlon_deg=dlon_deg,
        dlat_deg=dlat_deg,
        r_sphere=safe_eval_number(values["rSphere"]),
        template_delta_t=safe_eval_number(values["deltaT"]),
    )


def parse_size_file(path: Path) -> Tuple[int, int]:
    params: Dict[str, int] = {}
    pattern = re.compile(r"PARAMETER\s*\(\s*(sNx|sNy|nSx|nSy|nPx|nPy)\s*=\s*(\d+)\s*\)", re.I)
    for line in path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            params[match.group(1)] = int(match.group(2))

    missing = sorted({"sNx", "sNy", "nSx", "nSy", "nPx", "nPy"} - set(params))
    if missing:
        raise ValueError("{} is missing: {}".format(path, ", ".join(missing)))
    return params["sNx"] * params["nSx"] * params["nPx"], params["sNy"] * params["nSy"] * params["nPy"]


def parse_export_value(text: str, name: str) -> Optional[str]:
    pattern = re.compile(r"^\s*export\s+{}=(.+?)\s*$".format(re.escape(name)), re.M)
    match = pattern.search(text)
    if not match:
        return None
    return match.group(1).strip().strip("'\"")


def discover_jobs(jobs_dir: Path, template_delta_t: float) -> List[JobCase]:
    jobs: List[JobCase] = []
    for path in sorted(jobs_dir.glob("job_tc2_*.slurm")):
        text = path.read_text()
        name = parse_export_value(text, "CASE_NAME")
        alpha = parse_export_value(text, "ALPHA_VALUE")
        delta_t = parse_export_value(text, "DELTA_T")
        if name is None or alpha is None:
            continue
        jobs.append(
            JobCase(
                name=name,
                alpha=safe_eval_number(alpha),
                delta_t=safe_eval_number(delta_t) if delta_t is not None else template_delta_t,
                path=path,
            )
        )
    if not jobs:
        raise ValueError("no TC2 job files found in {}".format(jobs_dir))
    return jobs


def load_gendata_module(path: Path):
    spec = importlib.util.spec_from_file_location("tc2_gendata_ref", path)
    if spec is None or spec.loader is None:
        raise ValueError("could not import {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_consistency(grid: GridConfig, size_shape: Tuple[int, int], gendata_module) -> List[str]:
    warnings: List[str] = []
    size_nx, size_ny = size_shape
    if (grid.nx, grid.ny) != (size_nx, size_ny):
        warnings.append(
            "input/data grid {}x{} does not match SIZE.h {}x{}".format(
                grid.nx, grid.ny, size_nx, size_ny
            )
        )
    if (grid.nx, grid.ny) != (getattr(gendata_module, "NX", None), getattr(gendata_module, "NY", None)):
        warnings.append(
            "input/data grid does not match input/gendata_ref.py NX/NY ({}x{})".format(
                getattr(gendata_module, "NX", "?"), getattr(gendata_module, "NY", "?")
            )
        )
    module_radius = float(getattr(gendata_module, "R_EARTH", grid.r_sphere))
    if not math.isclose(grid.r_sphere, module_radius, rel_tol=0.0, abs_tol=1.0e-6):
        warnings.append(
            "input/data rSphere={} differs from gendata_ref.py R_EARTH={}".format(
                grid.r_sphere, module_radius
            )
        )
    return warnings


def select_jobs(jobs: Iterable[JobCase], selected: Optional[List[str]]) -> List[JobCase]:
    jobs_by_name = {job.name: job for job in jobs}
    if not selected:
        return list(jobs_by_name.values())

    missing = sorted(set(selected) - set(jobs_by_name))
    if missing:
        raise ValueError("unknown case(s): {}".format(", ".join(missing)))
    return [jobs_by_name[name] for name in selected]


def make_tc2_uv(grid: GridConfig, gendata_module, alpha: float) -> Tuple[np.ndarray, np.ndarray]:
    r_earth = float(getattr(gendata_module, "R_EARTH"))
    u0 = float(getattr(gendata_module, "U0"))
    dlon_rad = math.radians(grid.dlon_deg)
    dlat_rad = math.radians(grid.dlat_deg)

    lon_u = np.arange(grid.nx, dtype=np.float64) * dlon_rad
    lat_u = -0.5 * np.pi + (np.arange(grid.ny, dtype=np.float64) + 0.5) * dlat_rad
    lon_u_2d, lat_u_2d = np.meshgrid(lon_u, lat_u, indexing="xy")
    u = u0 * (
        np.cos(lat_u_2d) * np.cos(alpha)
        + np.cos(lon_u_2d) * np.sin(lat_u_2d) * np.sin(alpha)
    )

    lon_v = (np.arange(grid.nx, dtype=np.float64) + 0.5) * dlon_rad
    lat_v = -0.5 * np.pi + np.arange(grid.ny, dtype=np.float64) * dlat_rad
    lon_v_2d, _lat_v_2d = np.meshgrid(lon_v, lat_v, indexing="xy")
    v = -u0 * np.sin(lon_v_2d) * np.sin(alpha)
    v[0, :] = 0.0
    v[-1, :] = 0.0

    if not math.isclose(grid.r_sphere, r_earth, rel_tol=0.0, abs_tol=1.0e-6):
        # CFL distances below use input/data rSphere; this import is used only
        # for the analytic U0 and formula constants.
        pass
    return u, v


def compute_report(
    grid: GridConfig,
    job: JobCase,
    gendata_module,
    warn_threshold: float,
    fail_threshold: float,
) -> CflReport:
    u, v = make_tc2_uv(grid, gendata_module, job.alpha)

    lat = -90.0 + (np.arange(grid.ny, dtype=np.float64) + 0.5) * grid.dlat_deg
    lon = (np.arange(grid.nx, dtype=np.float64) + 0.5) * grid.dlon_deg
    dx = grid.r_sphere * np.cos(np.deg2rad(lat)) * np.deg2rad(grid.dlon_deg)
    dy = grid.r_sphere * np.deg2rad(grid.dlat_deg)

    cfl_x = np.abs(u) * job.delta_t / dx[:, None]
    cfl_y = np.abs(v) * job.delta_t / dy
    cfl_total = cfl_x + cfl_y

    max_index = np.unravel_index(int(np.argmax(cfl_total)), cfl_total.shape)
    max_total = float(cfl_total[max_index])
    status = "PASS"
    if max_total > fail_threshold:
        status = "FAIL"
    elif max_total > warn_threshold:
        status = "WARN"

    if max_total > 0.0:
        suggested_delta_t = job.delta_t * warn_threshold / max_total
    else:
        suggested_delta_t = math.inf

    return CflReport(
        case=job,
        status=status,
        max_x=float(np.max(cfl_x)),
        max_y=float(np.max(cfl_y)),
        max_total=max_total,
        location_lon=float(lon[max_index[1]]),
        location_lat=float(lat[max_index[0]]),
        suggested_delta_t=suggested_delta_t,
    )


def print_report(
    grid: GridConfig,
    reports: List[CflReport],
    warnings: List[str],
    warn_threshold: float,
    fail_threshold: float,
) -> None:
    print("TC2 initial-condition CFL preflight")
    print("grid       : {} x {}".format(grid.nx, grid.ny))
    print("spacing    : dlon={:g} deg, dlat={:g} deg".format(grid.dlon_deg, grid.dlat_deg))
    print("template dt: {:g} s".format(grid.template_delta_t))
    print("thresholds : WARN>{:g}, FAIL>{:g}".format(warn_threshold, fail_threshold))
    print("note       : TC2 velocities evolve; this checks the analytic initial field.")
    print()

    if warnings:
        print("Consistency warnings:")
        for warning in warnings:
            print("  - {}".format(warning))
        print()

    header = (
        "case   alpha(rad)  deltaT(s)  status   max_CFL_x   max_CFL_y   "
        "max_total   lon(deg)   lat(deg)  dt_for_warn(s)"
    )
    print(header)
    print("-" * len(header))
    for report in reports:
        dt_text = "inf" if math.isinf(report.suggested_delta_t) else "{:11.6g}".format(report.suggested_delta_t)
        print(
            "{:<5} {:11.6g} {:10.6g}  {:<6} {:11.6g} {:11.6g} {:11.6g} "
            "{:10.4f} {:10.4f} {:>15}".format(
                report.case.name,
                report.case.alpha,
                report.case.delta_t,
                report.status,
                report.max_x,
                report.max_y,
                report.max_total,
                report.location_lon,
                report.location_lat,
                dt_text,
            )
        )


def main() -> int:
    args = parse_args()
    if args.warn_threshold <= 0.0 or args.fail_threshold <= 0.0:
        print("ERROR: thresholds must be positive.", file=sys.stderr)
        return 2
    if args.warn_threshold >= args.fail_threshold:
        print("ERROR: warn-threshold must be smaller than fail-threshold.", file=sys.stderr)
        return 2

    case_dir = Path(__file__).resolve().parents[1]
    grid = parse_data_file(case_dir / "input" / "data")
    size_shape = parse_size_file(case_dir / "code" / "SIZE.h")
    jobs = discover_jobs(case_dir / "jobs" / "large", grid.template_delta_t)
    selected_jobs = select_jobs(jobs, args.cases)
    gendata_module = load_gendata_module(case_dir / "input" / "gendata_ref.py")

    warnings = validate_consistency(grid, size_shape, gendata_module)
    reports = [
        compute_report(grid, job, gendata_module, args.warn_threshold, args.fail_threshold)
        for job in selected_jobs
    ]
    print_report(grid, reports, warnings, args.warn_threshold, args.fail_threshold)

    if args.strict and any(report.status != "PASS" for report in reports):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
