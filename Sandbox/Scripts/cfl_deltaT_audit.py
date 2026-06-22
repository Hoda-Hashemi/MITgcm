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
    nx, dlon = parse_repeat_grid(values["delX"], "delX")
    ny, dlat = parse_repeat_grid(values["delY"], "delY")
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
    if "deltaT" in values:
        delta_t = safe_eval_number(values["deltaT"])
    if "nTimeSteps" in values:
        n_steps = int(round(safe_eval_number(values["nTimeSteps"])))
    return Schedule(delta_t=delta_t, n_steps=n_steps)


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
        delta_t = parse_export(text, "DELTA_T")
        total_seconds = parse_export(text, "TOTAL_SECONDS")
        run_dir = parse_export(text, "RUN_DIR")
        alpha = parse_export(text, "ALPHA_VALUE") or "0.0"
        if not name or not delta_t or not total_seconds or not run_dir:
            continue
        jobs.append(
            Job(
                case_id=case_id,
                name=name,
                alpha=safe_eval_number(alpha),
                delta_t=safe_eval_number(delta_t),
                total_seconds=safe_eval_number(total_seconds),
                run_dir=Path(run_dir.replace("$CASE_DIR", str(case_dir))),
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


def initial_velocity(case_id, module, alpha):
    if case_id == "TC1":
        _tracer, u, v = module.make_tc1_fields(alpha_rad=alpha)
        return u, v
    if case_id == "TC2":
        return module.make_tc2_u(alpha_rad=alpha), module.make_tc2_v(alpha_rad=alpha)
    if hasattr(module, "make_velocity_fields"):
        return module.make_velocity_fields(alpha_rad=alpha)
    raise ValueError("{} gendata_ref.py has no known velocity function".format(case_id))


def cfl_from_uv(grid, delta_t, u, v, iteration=-1):
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
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
            try:
                u, v = initial_velocity(case_id, module, job.alpha)
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
        "submitted jobs are controlled by `jobs/large/job_*.slurm`: each job exports "
        "`DELTA_T`, computes `nTimeSteps=round(TOTAL_SECONDS/DELTA_T)`, copies the "
        "template input directory, and rewrites the run-local `data`. The table "
        "therefore separates template, submitted job, and archived run deltaT."
    )
    lines.append("")
    lines.append(
        "TC1 has the most visible discrepancy: its legacy `tools/check_cfl_tc1.py` "
        "prints the template `deltaT=1 s`, while the submitted jobs use 60 s, 10 s, "
        "or 0.75 s and the archived run-local `data` files confirm those values. "
        "Use this audit table for submitted-run CFL decisions."
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
        "TC5 does not look like a simple CFL failure: the initial advective CFL "
        "is 0.043, but archived fields become non-finite after iteration 1440 "
        "and the CG residuals later print NaN. Treat TC5 as needing a run-health "
        "fix or a targeted shorter-deltaT rerun before using later-day plots."
    )
    lines.append("")
    lines.append(
        "TC4 and TC7 should be read from the job schedule, not from their template "
        "`nTimeSteps`: both submitted jobs target 5 days at 60 s, i.e. 7200 steps. "
        "TC4 now has completed `run_u0_20` output and its saved advective CFL can "
        "be read from archived U/V fields, while TC7 has input data staged but needs "
        "a smaller deltaT before a final rerun."
    )
    lines.append("")
    lines.append(
        "TC7 cannot be fully audited yet because completed `run_analysis` output is not "
        "archived. The staged analyzed input is present, and the preflight advective CFL "
        "indicates the current 60 s job schedule is too large for final validation."
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
        "Decision",
    ]
    out = ["<div class='table-scroll'><table>"]
    out.append("<thead><tr>{}</tr></thead>".format("".join("<th>{}</th>".format(html.escape(h)) for h in headers)))
    out.append("<tbody>")
    for row in rows:
        out.append("<tr>{}</tr>".format("".join("<td>{}</td>".format(html.escape(c)) for c in row_cells(row))))
    out.append("</tbody></table></div>")
    return "".join(out)


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
          <p>The template <code>input/data</code> files provide the grid and default schedule, but the submitted Slurm files are the source of truth for runs. Each job exports <code>DELTA_T</code>, computes <code>nTimeSteps=round(TOTAL_SECONDS/DELTA_T)</code>, copies the input directory, then rewrites the run-local <code>data</code> file. The table keeps template, job, and archived-run <code>deltaT</code> separate.</p>
          <p>TC1 has the most visible discrepancy: its legacy <code>tools/check_cfl_tc1.py</code> prints the template <code>deltaT=1 s</code>, while the submitted jobs use <code>60 s</code>, <code>10 s</code>, or <code>0.75 s</code>. Use this audit table for submitted-run CFL decisions.</p>
          <p>The small near-polar/tilted cases use <code>0.75 s</code>, the mildly tilted cases use <code>10 s</code>, standard zonal cases use <code>60 s</code>, and TC6 uses <code>30 s</code>.</p>
        </div>
      </section>
      <section class='description-block detail-expected'>
        <h3>Decision</h3>
        <div class='description-copy'>
          <p>No completed run exceeds advective CFL 1.0. TC2 alpha=0.05 is above the conservative 0.5 margin in saved output, so use <code>deltaT&lt;=8.93 s</code> if that margin is required. TC5 becomes non-finite after iteration 1440 despite an initial advective CFL of 0.043, so that case needs a run-health investigation rather than a CFL-only explanation.</p>
          <p>TC4 has completed <code>run_u0_20</code> output, so its saved advective CFL is read from archived U/V fields. TC7 has staged analyzed input, but the preflight audit flags <code>deltaT=60 s</code> as too large before a final rerun.</p>
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
        print("wrote {}".format(fragments / "cfl_deltaT_report.md"))
        print("wrote {}".format(fragments / "cfl_deltaT_report.html"))


if __name__ == "__main__":
    main()
