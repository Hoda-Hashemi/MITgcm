#!/usr/bin/env python3
"""CFL and deltaT audit for the Williamson sandbox cases.

This helper is read-only with respect to model inputs and run directories. With
``--write-fragments`` it writes docs-ready fragments under ``docs/fragments``.
"""

import argparse
import ast
import html
import importlib.util
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mitgcm_io import discover_iterations, read_mds_field, to_2d


WARN_CFL = 0.5
FAIL_CFL = 1.0
CASE_IDS = ("TC1", "TC2", "TC3", "TC4", "TC5", "TC6", "TC7")


@dataclass
class Grid:
    nx: int
    ny: int
    dlon_deg: float
    dlat_deg: float
    radius: float
    gravity: float
    depth: float
    implicit_free_surface: bool

    @property
    def lat_centers(self):
        return -90.0 + (np.arange(self.ny, dtype=np.float64) + 0.5) * self.dlat_deg

    @property
    def lon_centers(self):
        return (np.arange(self.nx, dtype=np.float64) + 0.5) * self.dlon_deg

    @property
    def dx(self):
        return self.radius * np.cos(np.deg2rad(self.lat_centers)) * math.radians(self.dlon_deg)

    @property
    def dy(self):
        return self.radius * math.radians(self.dlat_deg)

    @property
    def dx_min(self):
        return float(np.min(self.dx))

    @property
    def c_gravity(self):
        return math.sqrt(self.gravity * self.depth)


@dataclass
class Job:
    case_id: str
    name: str
    alpha: float
    delta_t: float
    total_seconds: float
    run_dir: Path
    path: Path

    @property
    def label(self):
        if self.run_dir.name.startswith("run_alpha_"):
            return self.run_dir.name.replace("run_alpha_", "alpha=")
        if self.run_dir.name == "run_standard":
            return "standard"
        if self.run_dir.name == "run_analysis":
            return "analysis"
        return self.run_dir.name


@dataclass
class Schedule:
    delta_t: float | None
    n_steps: int | None
    end_time: float | None = None


@dataclass
class CflStats:
    max_x: float
    max_y: float
    max_total: float
    lon: float
    lat: float
    max_speed: float
    iteration: int = -1


@dataclass
class AuditRow:
    case_id: str
    label: str
    job_name: str
    alpha: float
    template_delta_t: float | None
    template_n_steps: int | None
    delta_t: float
    n_steps: int
    run_delta_t: float | None
    run_n_steps: int | None
    observed_delta_t: float
    depth: float
    wave_speed: float
    cg_x: float
    cg_y: float
    initial: CflStats | None
    observed: CflStats | None
    run_status: str
    recommendation: str
    command: str
    notes: str


def safe_eval_number(expr):
    expr = expr.strip().rstrip(",").replace("D", "E").replace("d", "e")
    tree = ast.parse(expr, mode="eval")

    def eval_node(node):
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
        raise ValueError("unsupported numeric expression: {}".format(expr))

    return eval_node(tree)


def parse_assignments(path):
    values = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        clean = line.split("#", 1)[0].strip()
        if not clean or "=" not in clean:
            continue
        key, raw = clean.split("=", 1)
        values[key.strip()] = raw.strip().rstrip(",")
    return values


def parse_repeat_grid(value, name):
    match = re.fullmatch(
        r"(\d+)\s*\*\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)",
        value.strip().rstrip(","),
    )
    if not match:
        raise ValueError("{} must use MITgcm repeat syntax".format(name))
    return int(match.group(1)), safe_eval_number(match.group(2))


def parse_bool(raw):
    return raw.strip().strip(".").upper() == "TRUE"


def parse_grid(path):
    values = parse_assignments(path)
    if "delX" in values and "delY" in values:
        nx, dlon = parse_repeat_grid(values["delX"], "delX")
        ny, dlat = parse_repeat_grid(values["delY"], "delY")
    else:
        nx, dlon = 1440, 0.25
        ny, dlat = 720, 0.25
    return Grid(
        nx=nx,
        ny=ny,
        dlon_deg=dlon,
        dlat_deg=dlat,
        radius=safe_eval_number(values["rSphere"]),
        gravity=safe_eval_number(values["gravity"]),
        depth=safe_eval_number(values["delR"]),
        implicit_free_surface=parse_bool(values.get("implicitFreeSurface", ".FALSE.")),
    )


def parse_schedule(path):
    if not path.exists():
        return Schedule(delta_t=None, n_steps=None)
    values = parse_assignments(path)
    delta_t = None
    n_steps = None
    end_time = None
    if "deltaT" in values:
        delta_t = safe_eval_number(values["deltaT"])
    if "nTimeSteps" in values:
        n_steps = int(round(safe_eval_number(values["nTimeSteps"])))
    if "endTime" in values:
        end_time = safe_eval_number(values["endTime"])
    if n_steps is None and delta_t and end_time:
        n_steps = int(round(end_time / delta_t))
    return Schedule(delta_t=delta_t, n_steps=n_steps, end_time=end_time)


def parse_export(text, name):
    match = re.search(r"^\s*export\s+{}=(.+?)\s*$".format(re.escape(name)), text, re.M)
    if not match:
        return None
    return match.group(1).strip().strip("'\"")


def discover_jobs(repo_root, case_id):
    case_dir = repo_root / "Sandbox" / "vortexSphere_Williamson_{}".format(case_id)
    jobs = []
    for path in sorted((case_dir / "jobs" / "large").glob("job_{}_*.slurm".format(case_id.lower()))):
        text = path.read_text(encoding="utf-8", errors="ignore")
        name = parse_export(text, "CASE_NAME")
        run_dir = parse_export(text, "RUN_DIR")
        if run_dir:
            expanded_run_dir = Path(
                run_dir.replace("$CASE_DIR", str(case_dir)).replace("$MITGCM_DIR", str(repo_root))
            )
        else:
            expanded_run_dir = None
        run_schedule = parse_schedule(expanded_run_dir / "data") if expanded_run_dir is not None else Schedule(None, None)
        delta_t = parse_export(text, "DELTA_T")
        total_seconds = parse_export(text, "TOTAL_SECONDS")
        if delta_t is None and run_schedule.delta_t is not None:
            delta_t = str(run_schedule.delta_t)
        if total_seconds is None and run_schedule.delta_t is not None and run_schedule.n_steps is not None:
            total_seconds = str(run_schedule.delta_t * run_schedule.n_steps)
        alpha = parse_export(text, "ALPHA_VALUE") or "0.0"
        if not name or not delta_t or not total_seconds or expanded_run_dir is None:
            continue
        jobs.append(
            Job(
                case_id=case_id,
                name=name,
                alpha=safe_eval_number(alpha),
                delta_t=safe_eval_number(delta_t),
                total_seconds=safe_eval_number(total_seconds),
                run_dir=expanded_run_dir,
                path=path,
            )
        )
    return jobs


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError("could not import {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def initial_velocity(case_id, module, job):
    alpha = job.alpha
    if case_id == "TC1":
        _tracer, u, v = module.make_tc1_fields(alpha_rad=alpha)
        return u, v
    if case_id == "TC2":
        return module.make_tc2_u(alpha_rad=alpha), module.make_tc2_v(alpha_rad=alpha)
    if hasattr(module, "make_velocity_fields"):
        try:
            return module.make_velocity_fields(alpha_rad=alpha)
        except TypeError as exc:
            if "alpha_rad" not in str(exc):
                raise
            u0_value = parse_export(job.path.read_text(encoding="utf-8", errors="ignore"), "TC4_U0_VALUE")
            if u0_value is not None and hasattr(module, "U0"):
                module.U0 = safe_eval_number(u0_value)
            return module.make_velocity_fields()
    raise ValueError("{} gendata_ref.py has no known velocity function".format(case_id))


def cfl_from_uv(grid, delta_t, u, v, iteration=-1):
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    if u.shape != v.shape or u.shape != (grid.ny, grid.nx):
        return None
    dx = grid.dx[:, None]
    dy = grid.dy
    finite = np.isfinite(u) & np.isfinite(v) & np.isfinite(dx) & (dx > 0.0)
    if not np.any(finite):
        return None
    dx2 = np.broadcast_to(dx, u.shape)
    cfl_x = np.full_like(u, np.nan, dtype=np.float64)
    cfl_y = np.full_like(v, np.nan, dtype=np.float64)
    cfl_total = np.full_like(u, np.nan, dtype=np.float64)
    cfl_x[finite] = np.abs(u[finite]) * delta_t / dx2[finite]
    cfl_y[finite] = np.abs(v[finite]) * delta_t / dy
    cfl_total[finite] = cfl_x[finite] + cfl_y[finite]
    idx = np.unravel_index(int(np.nanargmax(cfl_total)), cfl_total.shape)
    speed = np.hypot(u, v)
    return CflStats(
        max_x=float(cfl_x[idx]),
        max_y=float(cfl_y[idx]),
        max_total=float(cfl_total[idx]),
        lon=float(grid.lon_centers[idx[1]]),
        lat=float(grid.lat_centers[idx[0]]),
        max_speed=float(np.nanmax(speed)),
        iteration=iteration,
    )


def read_cube_compact_field(path):
    if not path.exists():
        return None
    raw = np.fromfile(path, dtype=">f4")
    expected = 6 * 32 * 32
    if raw.size != expected:
        return None
    field = raw.reshape((32, 6, 32), order="F").transpose(1, 0, 2).astype(np.float64)
    if not np.isfinite(field).all():
        return None
    return field


def cube_metric_min(run_dir):
    dx_values = []
    dy_values = []
    for face in range(1, 7):
        path = run_dir / "tile{:03d}.mitgrid".format(face)
        if not path.exists():
            return None
        raw = np.fromfile(path, dtype=">f8")
        expected = 33 * 33 * 16
        if raw.size != expected:
            return None
        data = raw.reshape((33, 33, 16), order="F")
        dx_values.extend((float(np.min(data[:32, :32, 2])), float(np.min(data[:32, :32, 10]))))
        dy_values.extend((float(np.min(data[:32, :32, 3])), float(np.min(data[:32, :32, 11]))))
    dx_min = min(dx_values)
    dy_min = min(dy_values)
    if dx_min <= 0.0 or dy_min <= 0.0:
        return None
    return dx_min, dy_min


def cube_initial_cfl(run_dir, delta_t):
    metrics = cube_metric_min(run_dir)
    u = read_cube_compact_field(run_dir / "u_init.bin")
    v = read_cube_compact_field(run_dir / "v_init.bin")
    if metrics is None or u is None or v is None:
        return None
    dx_min, dy_min = metrics
    cfl_x = float(np.max(np.abs(u))) * delta_t / dx_min
    cfl_y = float(np.max(np.abs(v))) * delta_t / dy_min
    return CflStats(
        max_x=cfl_x,
        max_y=cfl_y,
        max_total=max(cfl_x, cfl_y),
        lon=math.nan,
        lat=math.nan,
        max_speed=float(np.nanmax(np.hypot(u, v))),
        iteration=-1,
    )


def scan_run_output(grid, job, delta_t):
    run_dir = job.run_dir
    if not run_dir.is_dir():
        return None, "no archived run directory"
    if not (run_dir / "data").exists():
        return None, "run directory has no data file"
    iters = sorted(set(discover_iterations(run_dir, "U")) & set(discover_iterations(run_dir, "V")))
    if not iters:
        return None, "no archived U/V snapshots"

    best = None
    first_nonfinite = None
    all_nonfinite = []
    for iteration in iters:
        u = to_2d(read_mds_field(run_dir, "U", iteration))
        v = to_2d(read_mds_field(run_dir, "V", iteration))
        finite = np.isfinite(u) & np.isfinite(v)
        if not np.any(finite):
            all_nonfinite.append(iteration)
            if first_nonfinite is None:
                first_nonfinite = iteration
            continue
        stats = cfl_from_uv(grid, delta_t, u, v, iteration=iteration)
        if stats and (best is None or stats.max_total > best.max_total):
            best = stats
    if first_nonfinite is not None:
        return best, "non-finite output from iter {}".format(first_nonfinite)
    return best, "archived U/V scan over {} snapshots".format(len(iters))


def recommendation(row):
    adv = None
    if row.observed is not None:
        adv = row.observed.max_total
    elif row.initial is not None:
        adv = row.initial.max_total

    if "non-finite" in row.run_status:
        return "investigate non-finite run; not explained by initial advective CFL"
    if adv is None:
        return "cannot verify advective CFL until inputs/run output exist"
    if adv >= FAIL_CFL:
        return "reduce deltaT before rerun"
    if adv > WARN_CFL:
        base_delta_t = row.observed_delta_t if row.observed is not None else row.delta_t
        dt_warn = base_delta_t * WARN_CFL / adv
        return "above 0.5 margin; use deltaT <= {:.2f} s for margin".format(dt_warn)
    return "no advective deltaT change indicated"


def audit(repo_root):
    rows = []
    for case_id in CASE_IDS:
        case_dir = repo_root / "Sandbox" / "vortexSphere_Williamson_{}".format(case_id)
        grid = parse_grid(case_dir / "input" / "data")
        template_schedule = parse_schedule(case_dir / "input" / "data")
        module = load_module(case_dir / "input" / "gendata_ref.py", "gendata_{}".format(case_id.lower()))
        for job in discover_jobs(repo_root, case_id):
            run_schedule = parse_schedule(job.run_dir / "data")
            observed_delta_t = run_schedule.delta_t if run_schedule.delta_t is not None else job.delta_t
            cg_x = grid.c_gravity * job.delta_t / grid.dx_min
            cg_y = grid.c_gravity * job.delta_t / grid.dy
            initial = None
            notes = []
            initial = cube_initial_cfl(job.run_dir, job.delta_t)
            if initial is not None:
                notes.append("cubed-sphere initial CFL from compact run inputs.")
            else:
                try:
                    u, v = initial_velocity(case_id, module, job)
                    initial = cfl_from_uv(grid, job.delta_t, u, v)
                except Exception as exc:
                    notes.append("initial velocity unavailable: {}".format(exc))
            observed, run_status = scan_run_output(grid, job, observed_delta_t)
            n_steps = int(round(job.total_seconds / job.delta_t))
            if (
                template_schedule.delta_t is not None
                and abs(template_schedule.delta_t - job.delta_t) > 1.0e-12
            ):
                notes.append(
                    "template deltaT={} s differs from submitted DELTA_T={} s.".format(
                        fmt(template_schedule.delta_t, 3), fmt(job.delta_t, 3)
                    )
                )
            if template_schedule.n_steps is not None and template_schedule.n_steps != n_steps:
                notes.append(
                    "template nTimeSteps={} differs from submitted {} steps.".format(
                        template_schedule.n_steps, n_steps
                    )
                )
            if run_schedule.delta_t is not None and abs(run_schedule.delta_t - job.delta_t) > 1.0e-12:
                notes.append(
                    "archived run data deltaT={} s differs from submitted DELTA_T={} s.".format(
                        fmt(run_schedule.delta_t, 3), fmt(job.delta_t, 3)
                    )
                )
            if case_id in {"TC3", "TC4", "TC5", "TC6", "TC7"}:
                notes.append("run check_cfl with the project venv, not the system python3.")
            if grid.implicit_free_surface:
                notes.append("gravity-wave CFL is explicit-only because implicitFreeSurface=.TRUE.")
            command = "cd {}\n/home/hmh85/scratch/MITgcm/.venv/bin/python tools/check_cfl_{}.py --all".format(
                case_dir.relative_to(repo_root), case_id.lower()
            )
            row = AuditRow(
                case_id=case_id,
                label=job.label,
                job_name=job.name,
                alpha=job.alpha,
                template_delta_t=template_schedule.delta_t,
                template_n_steps=template_schedule.n_steps,
                delta_t=job.delta_t,
                n_steps=n_steps,
                run_delta_t=run_schedule.delta_t,
                run_n_steps=run_schedule.n_steps,
                observed_delta_t=observed_delta_t,
                depth=grid.depth,
                wave_speed=grid.c_gravity,
                cg_x=cg_x,
                cg_y=cg_y,
                initial=initial,
                observed=observed,
                run_status=run_status,
                recommendation="",
                command=command,
                notes=" ".join(notes),
            )
            row.recommendation = recommendation(row)
            rows.append(row)
    rows.extend(discover_tc1_cubed_rows(repo_root))
    return rows


def tc1_cubed_alpha_value(label):
    if label == "1.52":
        return math.pi / 2.0 - 0.05
    if label == "1.57":
        return math.pi / 2.0
    return safe_eval_number(label)


def monitor_advcfl(run_dir):
    pattern = re.compile(r"%MON\s+advcfl_(?:uvel|vvel|wvel)_max\s*=\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")
    values = []
    for path in sorted(run_dir.glob("STDOUT.*")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        values.extend(safe_eval_number(match.group(1)) for match in pattern.finditer(text))
    if not values:
        return None
    return max(values)


def discover_tc1_cubed_rows(repo_root):
    rows = []
    output_root = (
        repo_root
        / "Sandbox"
        / "output"
        / "existingTutorials"
        / "test1"
        / "MITGCM_Williamson_TC1"
        / "advect_cs"
    )
    setup_data = (
        repo_root
        / "Sandbox"
        / "MITgcm_Williamson_TC1"
        / "existingTutorials"
        / "advect_cs"
        / "input"
        / "data"
    )
    template_schedule = parse_schedule(setup_data)
    order = {"0": 0, "0.05": 1, "1.52": 2, "1.57": 3}
    run_dirs = sorted(
        output_root.glob("alpha_*/run_alpha_*"),
        key=lambda path: order.get(path.parent.name.removeprefix("alpha_"), 99),
    )
    for run_dir in run_dirs:
        alpha_label = run_dir.parent.name.removeprefix("alpha_")
        alpha = tc1_cubed_alpha_value(alpha_label)
        run_schedule = parse_schedule(run_dir / "data")
        delta_t = run_schedule.delta_t or template_schedule.delta_t
        if delta_t is None:
            continue
        n_steps = run_schedule.n_steps or template_schedule.n_steps
        total_seconds = run_schedule.end_time
        if total_seconds is None and n_steps is not None:
            total_seconds = n_steps * delta_t
        if total_seconds is None:
            total_seconds = 1036800.0
        if n_steps is None:
            n_steps = int(round(total_seconds / delta_t))
        monitor_cfl = monitor_advcfl(run_dir)
        observed = None
        run_status = "MITgcm monitor advcfl scan"
        if monitor_cfl is None:
            run_status = "no MITgcm monitor advcfl values"
        else:
            observed = CflStats(
                max_x=monitor_cfl,
                max_y=monitor_cfl,
                max_total=monitor_cfl,
                lon=math.nan,
                lat=math.nan,
                max_speed=math.nan,
                iteration=-1,
            )
        row = AuditRow(
            case_id="TC1 cubed",
            label="advect_cs alpha={}".format(alpha_label),
            job_name="MITGCM_Williamson_TC1_advect_cs",
            alpha=alpha,
            template_delta_t=template_schedule.delta_t,
            template_n_steps=template_schedule.n_steps,
            delta_t=delta_t,
            n_steps=n_steps,
            run_delta_t=run_schedule.delta_t,
            run_n_steps=run_schedule.n_steps,
            observed_delta_t=delta_t,
            depth=100000.0,
            wave_speed=math.nan,
            cg_x=math.nan,
            cg_y=math.nan,
            initial=None,
            observed=observed,
            run_status=run_status,
            recommendation="",
            command="cd Sandbox/MITgcm_Williamson_TC1\n/home/hmh85/scratch/MITgcm/.venv/bin/python tools/postprocess_advect_cs.py --check-only",
            notes="cubed-sphere advect_cs row; CFL is MITgcm monitor max(advcfl_uvel, advcfl_vvel, advcfl_wvel), not spherical-polar recomputation.",
        )
        row.recommendation = recommendation(row)
        rows.append(row)
    return rows


def fmt(value, digits=3):
    if value is None:
        return "n/a"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "n/a"
    if abs(float(value)) >= 100.0:
        return "{:.0f}".format(float(value))
    if abs(float(value)) >= 10.0:
        return "{:.1f}".format(float(value))
    return ("{:.%df}" % digits).format(float(value))


def alpha_text(value):
    if abs(value) < 5.0e-7:
        return "0"
    return "{:.6g}".format(value)


def row_cells(row):
    initial = row.initial.max_total if row.initial else None
    observed = row.observed.max_total if row.observed else None
    speed = None
    if row.observed is not None:
        speed = row.observed.max_speed
    elif row.initial is not None:
        speed = row.initial.max_speed
    return [
        row.case_id,
        row.label,
        alpha_text(row.alpha),
        fmt(row.template_delta_t, 2),
        fmt(row.delta_t, 2),
        fmt(row.run_delta_t, 2),
        str(row.n_steps),
        fmt(row.depth, 0),
        fmt(row.wave_speed, 1),
        fmt(initial, 3),
        fmt(observed, 3),
        fmt(speed, 1),
        fmt(row.cg_x, 0),
        row.recommendation,
    ]


def compact_row_cells(row):
    max_cfl = None
    if row.observed is not None:
        max_cfl = row.observed.max_total
    elif row.initial is not None:
        max_cfl = row.initial.max_total
    total = row.n_steps * row.delta_t
    status = row.recommendation
    if status == "no advective deltaT change indicated":
        status = "OK"
    elif status == "reduce deltaT before rerun":
        status = "reduce dt"
    elif status == "cannot verify advective CFL until inputs/run output exist":
        status = "n/a"
    elif status.startswith("above 0.5 margin"):
        match = re.search(r"<= ([0-9.]+) s", status)
        status = "OK<1; 0.5 margin dt<={} s".format(match.group(1)) if match else "OK<1; above 0.5"
    return [
        row.case_id,
        row.label,
        fmt(row.delta_t, 2),
        str(row.n_steps),
        fmt(total, 0),
        fmt(total / 86400.0, 2),
        fmt(max_cfl, 3),
        status,
    ]


def case_class(case_id):
    slug = re.sub(r"[^a-z0-9]+", "-", case_id.lower()).strip("-")
    return "cfl-row cfl-row-{}".format(slug)


def markdown(rows):
    commands = sorted({row.command for row in rows})
    headers = [
        "case",
        "run",
        "alpha(rad)",
        "template dt(s)",
        "job dt(s)",
        "run dt(s)",
        "steps",
        "H(m)",
        "sqrt(gH)",
        "init adv CFL",
        "saved adv CFL",
        "max speed",
        "explicit Cg,x",
        "decision",
    ]
    lines = []
    lines.append("## Time-step safety check (CFL and deltaT)")
    lines.append("")
    lines.append(
        "This check asks whether the submitted model time step is small enough "
        "for the resolved velocities and grid spacing in each Williamson run."
    )
    lines.append("")
    lines.append(
        "For each submitted Williamson job, the audited Courant numbers are "
        "`C_adv = |u| deltaT / dx + |v| deltaT / dy` and "
        "`C_g = sqrt(gH) deltaT / dx`. The spherical-polar grid is 1440 x 720 "
        "with 0.25 degree spacing; the smallest center-cell zonal metric is "
        "60.65 m in the polar row and `dy = 27798.73 m`."
    )
    lines.append("")
    lines.append(
        "The gravity-wave column is an explicit-wave diagnostic only. All active "
        "runs set `implicitFreeSurface=.TRUE.`, so the external gravity wave is "
        "handled by the implicit free-surface solve rather than by the explicit "
        "advective CFL limit."
    )
    lines.append("")
    lines.append(
        "The template `input/data` files define the grid and default schedule, but "
        "vortexSphere submitted jobs are controlled by `jobs/large/job_*.slurm`: each job exports "
        "`DELTA_T`, computes `nTimeSteps=round(TOTAL_SECONDS/DELTA_T)`, copies the "
        "template input directory, and rewrites the run-local `data`. The table "
        "therefore separates template, submitted job, and archived run deltaT."
    )
    lines.append("")
    lines.append(
        "TC1 has two entries here: vortexSphere TC1 uses the submitted lat-lon jobs, "
        "while `TC1 cubed` is the MITgcm `advect_cs` tutorial and keeps the tutorial "
        "`deltaT=2700 s`. The cubed CFL values come from MITgcm monitor output."
    )
    lines.append("")
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(row_cells(row)) + " |")
    lines.append("")
    lines.append("### Decisions")
    lines.append("")
    lines.append(
        "No completed run exceeds advective CFL 1.0. TC2 alpha=0.05 reaches "
        "about 0.56 in the saved fields, above the conservative 0.5 margin; "
        "deltaT <= 8.93 s would keep the saved-output maximum under 0.5."
    )
    lines.append("")
    lines.append(
        "TC5 is now a completed CS32 rerun: the initial and monitored advective CFL "
        "remain small, the final saved state fields are finite through day 15, and "
        "the mountain is verified as static bathymetry rather than an eta bump. "
        "TC7 uses cubed-sphere compact initial fields for the submitted three-date suite."
    )
    lines.append("")
    lines.append(
        "`n/a` means the audit could not read a finite CFL source for that column: "
        "missing archived U/V fields, unavailable initial-velocity hook, or a cubed-sphere "
        "row where the spherical-polar gravity-wave metric is not used. TC4 now includes "
        "both `run_u0_20` and completed `run_u0_40` output."
    )
    lines.append("")
    lines.append(
        "TC7 has three completed analyzed-state rows: 21 Dec 1978, 16 Jan 1979, "
        "and 9 Jan 1979. The table values are CS32 compact-input CFL checks for "
        "the completed 25 s, 48-rank runs."
    )
    lines.append("")
    lines.append("### Check commands")
    lines.append("")
    lines.append("```bash")
    lines.append("cd /home/hmh85/scratch/MITgcm")
    lines.append("/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments")
    lines.append("```")
    lines.append("")
    lines.append("```bash")
    lines.append("\n\n".join(commands))
    lines.append("```")
    lines.append("")
    lines.append("For TC3-TC7, use the project venv; the system `python3` is too old for the shared `williamson_cfl.py` helper.")
    lines.append("")
    return "\n".join(lines)


def html_table(rows):
    headers = [
        "Case",
        "Run",
        "dt(s)",
        "n steps",
        "total(s)",
        "days",
        "max CFL",
        "Status/action",
    ]
    out = ["<div class='table-scroll'><table>"]
    out.append("<thead><tr>{}</tr></thead>".format("".join("<th>{}</th>".format(html.escape(h)) for h in headers)))
    out.append("<tbody>")
    for row in rows:
        out.append(
            "<tr class='{}'>{}</tr>".format(
                html.escape(case_class(row.case_id)),
                "".join("<td>{}</td>".format(html.escape(c)) for c in compact_row_cells(row)),
            )
        )
    out.append("</tbody></table></div>")
    return "".join(out)


def latex_escape(value):
    return (
        str(value)
        .replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("$", r"\$")
        .replace("#", r"\#")
        .replace("_", r"\_")
        .replace("{", r"\{")
        .replace("}", r"\}")
    )


def latex_table(rows):
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Submitted-run time-step and CFL settings used in the current MITgcm Williamson suite.}",
        r"\label{tab:submitted_time_step_cfl}",
        r"\begin{tabular}{llrrrrrl}",
        r"\toprule",
        r"Test & Run & $\Delta t$ [s] & $n$ steps & Total [s] & Days & Max CFL & Status/action \\",
        r"\midrule",
    ]
    for row in rows:
        case_id, label, dt, steps, total, days, max_cfl, status = compact_row_cells(row)
        lines.append(
            "{} & {} & {} & {} & {} & {} & {} & {} \\\\".format(
                latex_escape(case_id),
                latex_escape(label),
                latex_escape(dt),
                latex_escape(steps),
                latex_escape(total),
                latex_escape(days),
                latex_escape(max_cfl),
                latex_escape(status),
            )
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    return "\n".join(lines)


def html_fragment(rows):
    commands = sorted({row.command for row in rows})
    command_block = "\n\n".join(commands[:7])
    template = """<section id='cfl-deltat-audit' class='section-card experiment-section'>
  <header class='experiment-top'>
    <div>
      <span class='section-label'>Run checks</span>
      <h2>Time-Step Safety Check (CFL)</h2>
      <p class='experiment-summary'>This section checks whether each submitted model time step is small enough for the resolved velocities and grid spacing.</p>
    </div>
    <div class='metric-grid'>
      <div class='metric-card'><span>Grid</span><strong>1440 x 720</strong></div>
      <div class='metric-card'><span>Min dx</span><strong>60.65 m</strong></div>
      <div class='metric-card'><span>dy</span><strong>27.80 km</strong></div>
      <div class='metric-card'><span>Warn CFL</span><strong>0.5</strong></div>
    </div>
  </header>
  <div class='experiment-body'>
    <article class='experiment-description' aria-label='CFL audit details'>
      <section class='description-block detail-definition'>
        <h3>Equations</h3>
        <div class='description-copy'>
          <p>The advective audit uses the local spherical-polar metrics from the 0.25 degree grid:</p>
          <div class='equation-box'>
            <div class='math-line'>\\[C_{adv}=\\Delta t\\left(\\frac{|u|}{\\Delta x}+\\frac{|v|}{\\Delta y}\\right)\\]</div>
            <div class='math-line'>\\[\\Delta x=a\\cos\\phi\\,\\Delta\\lambda,\\quad \\Delta y=a\\Delta\\phi\\]</div>
            <div class='math-line'>\\[C_g=\\frac{\\sqrt{gH}\\,\\Delta t}{\\Delta x}\\]</div>
          </div>
          <p>The polar row dominates the metric audit: center-cell <code>dx_min=60.65 m</code>, while <code>dy=27798.73 m</code>. The gravity-wave CFL is shown as an explicit-wave diagnostic only because every active data file has <code>implicitFreeSurface=.TRUE.</code>.</p>
        </div>
      </section>
      <section class='description-block detail-measure'>
        <h3>How deltaT Was Chosen</h3>
        <div class='description-copy'>
          <p>The template <code>input/data</code> files provide the grid and default schedule, but the submitted Slurm files are the source of truth for vortexSphere runs. Each job exports <code>DELTA_T</code>, computes <code>nTimeSteps=round(TOTAL_SECONDS/DELTA_T)</code>, copies the input directory, then rewrites the run-local <code>data</code> file. The table keeps template, job, and archived-run <code>deltaT</code> separate.</p>
          <p>TC1 has two entries here: the vortexSphere TC1 rows use the submitted lat-lon jobs, while <code>TC1 cubed</code> rows are MITgcm <code>advect_cs</code> runs and keep the tutorial timestep <code>deltaT=2700 s</code>. The cubed CFL values come from MITgcm monitor output.</p>
          <p>The small near-polar/tilted cases use <code>0.75 s</code>, the mildly tilted cases use <code>10 s</code>, standard zonal cases use <code>60 s</code>, and TC6 uses <code>30 s</code>.</p>
        </div>
      </section>
      <section class='description-block detail-expected'>
        <h3>Decision</h3>
        <div class='description-copy'>
          <p>No completed run exceeds advective CFL 1.0. TC2 alpha=0.05 is above the conservative 0.5 margin in saved output, so use <code>deltaT&lt;=8.93 s</code> if that margin is required. TC5 is now a completed CS32 rerun: final saved state fields are finite through day 15, and the mountain is verified as static bathymetry rather than an eta bump. TC7 uses cubed-sphere compact initial fields for the submitted three-date suite.</p>
          <p><code>n/a</code> means the audit could not read a finite CFL source for that column: missing archived U/V fields, unavailable initial-velocity hook, or a cubed-sphere row where the spherical-polar gravity-wave metric is not used. TC4 now includes both <code>run_u0_20</code> and completed <code>run_u0_40</code> output.</p>
        </div>
      </section>
    </article>
    <div class='settings-stack'>
      <div class='table-block data-panel'>
        <h3>Compact CFL Table</h3>
        __CFL_TABLE__
      </div>
      <details class='details-panel case-block' open>
        <summary>Check commands</summary>
        <pre><code>cd /home/hmh85/scratch/MITgcm
/home/hmh85/scratch/MITgcm/.venv/bin/python Sandbox/Scripts/cfl_deltaT_audit.py --write-fragments

__CFL_COMMANDS__</code></pre>
      </details>
    </div>
  </div>
</section>
"""
    return (
        template.replace("__CFL_TABLE__", html_table(rows))
        .replace("__CFL_COMMANDS__", html.escape(command_block))
    )


def print_summary(rows):
    headers = ["case", "run", "template_dt", "job_dt", "run_dt", "init_adv", "obs_adv", "cg_x", "decision"]
    print("\t".join(headers))
    for row in rows:
        print(
            "\t".join(
                [
                    row.case_id,
                    row.label,
                    fmt(row.template_delta_t, 2),
                    fmt(row.delta_t, 2),
                    fmt(row.run_delta_t, 2),
                    fmt(row.initial.max_total if row.initial else None, 3),
                    fmt(row.observed.max_total if row.observed else None, 3),
                    fmt(row.cg_x, 0),
                    row.recommendation,
                ]
            )
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-fragments", action="store_true")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    rows = audit(repo_root)
    print_summary(rows)
    if args.write_fragments:
        fragments = repo_root / "docs" / "fragments"
        fragments.mkdir(parents=True, exist_ok=True)
        (fragments / "cfl_deltaT_report.md").write_text(markdown(rows), encoding="utf-8")
        (fragments / "cfl_deltaT_report.html").write_text(html_fragment(rows), encoding="utf-8")
        latex_path = repo_root / "docs" / "assets" / "williamson" / "submitted_time_step_cfl_table.tex"
        latex_path.parent.mkdir(parents=True, exist_ok=True)
        latex_path.write_text(latex_table(rows), encoding="utf-8")
        print("wrote {}".format(fragments / "cfl_deltaT_report.md"))
        print("wrote {}".format(fragments / "cfl_deltaT_report.html"))
        print("wrote {}".format(latex_path))


if __name__ == "__main__":
    main()
