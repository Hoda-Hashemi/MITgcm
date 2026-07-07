from __future__ import annotations

import argparse
import ast
import importlib.util
import math
import re
import sys
from pathlib import Path
from typing import Iterable, NamedTuple

import numpy as np

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
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            left = eval_node(node.left)
            right = eval_node(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left / right
        raise ValueError(f"unsupported numeric expression: {expr}")

    return eval_node(tree)


def parse_repeat_grid(value: str, name: str) -> tuple[int, float]:
    value = value.strip().rstrip(",")
    match = re.fullmatch(
        r"(\d+)\s*\*\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)",
        value,
    )
    if not match:
        raise ValueError(f"{name} must use MITgcm repeat syntax like 1440*0.25")
    return int(match.group(1)), safe_eval_number(match.group(2))


def parse_data_file(path: Path) -> GridConfig:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        clean = line.split("#", 1)[0].strip()
        if not clean or "=" not in clean:
            continue
        key, raw_value = clean.split("=", 1)
        key = key.strip()
        raw_value = raw_value.strip().rstrip(",")
        if key in {"deltaT", "rSphere", "delX", "delY"}:
            values[key] = raw_value

    missing = sorted({"deltaT", "rSphere"} - set(values))
    if missing:
        raise ValueError(f"{path} is missing: {', '.join(missing)}")
    values.setdefault("delX", "1440*0.25")
    values.setdefault("delY", "720*0.25")

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


def parse_size_file(path: Path) -> tuple[int, int]:
    params: dict[str, int] = {}
    pattern = re.compile(r"PARAMETER\s*\(\s*(sNx|sNy|nSx|nSy|nPx|nPy)\s*=\s*(\d+)\s*\)", re.I)
    for line in path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            params[match.group(1)] = int(match.group(2))

    missing = sorted({"sNx", "sNy", "nSx", "nSy", "nPx", "nPy"} - set(params))
    if missing:
        raise ValueError(f"{path} is missing: {', '.join(missing)}")
    return params["sNx"] * params["nSx"] * params["nPx"], params["sNy"] * params["nSy"] * params["nPy"]


def parse_export_value(text: str, name: str) -> str | None:
    pattern = re.compile(rf"^\s*export\s+{re.escape(name)}=(.+?)\s*$", re.M)
    match = pattern.search(text)
    if not match:
        return None
    return match.group(1).strip().strip("'\"")


def discover_jobs(case_dir: Path, case_code: str, template_delta_t: float) -> list[JobCase]:
    jobs: list[JobCase] = []
    for path in sorted((case_dir / "jobs" / "large").glob(f"job_{case_code.lower()}_*.slurm")):
        text = path.read_text()
        name = parse_export_value(text, "CASE_NAME")
        alpha = parse_export_value(text, "ALPHA_VALUE")
        delta_t = parse_export_value(text, "DELTA_T")
        if name is None:
            continue
        jobs.append(
            JobCase(
                name=name,
                alpha=safe_eval_number(alpha) if alpha is not None else 0.0,
                delta_t=safe_eval_number(delta_t) if delta_t is not None else template_delta_t,
                path=path,
            )
        )
    if not jobs:
        raise ValueError(f"no {case_code} job files found in {case_dir / 'jobs' / 'large'}")
    return jobs


def load_gendata_module(path: Path):
    spec = importlib.util.spec_from_file_location("williamson_gendata_ref", path)
    if spec is None or spec.loader is None:
        raise ValueError(f"could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_consistency(grid: GridConfig, size_shape: tuple[int, int], gendata_module) -> list[str]:
    warnings: list[str] = []
    size_nx, size_ny = size_shape
    if (grid.nx, grid.ny) != (size_nx, size_ny):
        warnings.append(f"input/data grid {grid.nx}x{grid.ny} does not match SIZE.h {size_nx}x{size_ny}")
    if (grid.nx, grid.ny) != (getattr(gendata_module, "NX", None), getattr(gendata_module, "NY", None)):
        warnings.append(
            "input/data grid does not match input/gendata_ref.py NX/NY "
            f"({getattr(gendata_module, 'NX', '?')}x{getattr(gendata_module, 'NY', '?')})"
        )
    module_radius = float(getattr(gendata_module, "R_EARTH", grid.r_sphere))
    if not math.isclose(grid.r_sphere, module_radius, rel_tol=0.0, abs_tol=1.0e-6):
        warnings.append(f"input/data rSphere={grid.r_sphere} differs from gendata_ref.py R_EARTH={module_radius}")
    return warnings


def select_jobs(jobs: Iterable[JobCase], selected: list[str] | None) -> list[JobCase]:
    jobs_by_name = {job.name: job for job in jobs}
    if not selected:
        return list(jobs_by_name.values())
    missing = sorted(set(selected) - set(jobs_by_name))
    if missing:
        raise ValueError(f"unknown case(s): {', '.join(missing)}")
    return [jobs_by_name[name] for name in selected]


def compute_report(
    grid: GridConfig,
    job: JobCase,
    gendata_module,
    warn_threshold: float,
    fail_threshold: float,
) -> CflReport:
    if not hasattr(gendata_module, "make_velocity_fields"):
        raise ValueError("gendata_ref.py must define make_velocity_fields(alpha_rad=...)")
    try:
        u, v = gendata_module.make_velocity_fields(alpha_rad=job.alpha)
    except TypeError as exc:
        if "alpha_rad" not in str(exc):
            raise
        u0_value = parse_export_value(job.path.read_text(), "TC4_U0_VALUE")
        if u0_value is not None and hasattr(gendata_module, "U0"):
            gendata_module.U0 = safe_eval_number(u0_value)
        u, v = gendata_module.make_velocity_fields()

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

    suggested_delta_t = job.delta_t * warn_threshold / max_total if max_total > 0.0 else math.inf
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


def parse_args(case_code: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Report {case_code} initial-condition advective CFL.")
    parser.add_argument("--case", action="append", dest="cases", metavar="NAME")
    parser.add_argument("--all", action="store_true", help="check all job files; default behavior")
    parser.add_argument("--warn-threshold", type=float, default=WARN_THRESHOLD)
    parser.add_argument("--fail-threshold", type=float, default=FAIL_THRESHOLD)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def print_report(
    case_code: str,
    grid: GridConfig,
    reports: list[CflReport],
    warnings: list[str],
    warn_threshold: float,
    fail_threshold: float,
) -> None:
    print(f"{case_code} initial-condition CFL preflight")
    print(f"grid       : {grid.nx} x {grid.ny}")
    print(f"spacing    : dlon={grid.dlon_deg:g} deg, dlat={grid.dlat_deg:g} deg")
    print(f"template dt: {grid.template_delta_t:g} s")
    print(f"thresholds : WARN>{warn_threshold:g}, FAIL>{fail_threshold:g}")
    print("note       : this checks the generated initial velocity field only.")
    print()
    if warnings:
        print("Consistency warnings:")
        for warning in warnings:
            print(f"  - {warning}")
        print()

    header = (
        "case   alpha(rad)  deltaT(s)  status   max_CFL_x   max_CFL_y   "
        "max_total   lon(deg)   lat(deg)  dt_for_warn(s)"
    )
    print(header)
    print("-" * len(header))
    for report in reports:
        dt_text = "inf" if math.isinf(report.suggested_delta_t) else f"{report.suggested_delta_t:11.6g}"
        print(
            f"{report.case.name:<5} {report.case.alpha:11.6g} {report.case.delta_t:10.6g}  "
            f"{report.status:<6} {report.max_x:11.6g} {report.max_y:11.6g} {report.max_total:11.6g} "
            f"{report.location_lon:10.4f} {report.location_lat:10.4f} {dt_text:>15}"
        )


def main(case_code: str) -> int:
    args = parse_args(case_code)
    if args.warn_threshold <= 0.0 or args.fail_threshold <= 0.0:
        print("ERROR: thresholds must be positive.", file=sys.stderr)
        return 2
    if args.warn_threshold >= args.fail_threshold:
        print("ERROR: warn-threshold must be smaller than fail-threshold.", file=sys.stderr)
        return 2

    case_dir = Path(__file__).resolve().parent.parent / f"vortexSphere_Williamson_{case_code}"
    grid = parse_data_file(case_dir / "input" / "data")
    size_shape = parse_size_file(case_dir / "code" / "SIZE.h")
    jobs = discover_jobs(case_dir, case_code, grid.template_delta_t)
    selected_jobs = select_jobs(jobs, args.cases)
    gendata_module = load_gendata_module(case_dir / "input" / "gendata_ref.py")
    warnings = validate_consistency(grid, size_shape, gendata_module)
    reports = [
        compute_report(grid, job, gendata_module, args.warn_threshold, args.fail_threshold)
        for job in selected_jobs
    ]
    print_report(case_code, grid, reports, warnings, args.warn_threshold, args.fail_threshold)
    if args.strict and any(report.status != "PASS" for report in reports):
        return 1
    return 0
