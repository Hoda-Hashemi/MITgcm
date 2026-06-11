from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np

SANDBOX_DIR = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = SANDBOX_DIR / "output"

CASE_OUTPUT_NAMES = {
    "TC1": "TestCase1",
    "TC2": "TestCase2",
    "TC3": "TestCase3",
}

RUN_LABELS = {
    "TC1": {"c1": "0", "c2": "0.05", "c3": "1.52", "c4": "1.57"},
    "TC2": {"c1": "0", "c2": "0.05", "c3": "1.52", "c4": "1.57"},
    "TC3": {"c1": "0", "c2": "1.0472"},
}


def read_text_value(path: Path, pattern: str) -> float | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
    if match is None:
        return None
    return float(match.group(1).replace("D", "e").replace("d", "e"))


def format_alpha_label(value: float) -> str:
    if not math.isfinite(value):
        return "unknown"
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text or "0"


def detect_case_code(run_dir: Path) -> str:
    text = str(run_dir)
    for case_code in CASE_OUTPUT_NAMES:
        if case_code in text:
            return case_code
    raise ValueError(f"Could not infer test case from {run_dir}")


def case_output_name(case_code: str) -> str:
    return CASE_OUTPUT_NAMES[case_code]


def alpha_from_run_name(run_dir: Path, case_code: str) -> str | None:
    run_name = run_dir.name
    direct_match = re.search(r"run_alpha_([0-9.]+)", run_name)
    if direct_match:
        return direct_match.group(1)

    case_match = re.search(r"alpha_(c\d+)", run_name)
    if case_match:
        return RUN_LABELS.get(case_code, {}).get(case_match.group(1))

    return None


def alpha_from_gendata(run_dir: Path) -> float | None:
    candidates = list(run_dir.glob("gendata*.py"))
    parent_input = run_dir.parent / "input"
    if parent_input.is_dir():
        candidates.extend(sorted(parent_input.glob("gendata*.py")))

    for path in candidates:
        value = read_text_value(path, r"^\s*ALPHA_RAD\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return value
        value = read_text_value(path, r"^\s*ALPHA\s*=\s*np\.deg2rad\(([+\-0-9.eEdD]+)\)")
        if value is not None:
            return math.radians(value)
        value = read_text_value(path, r"^\s*ALPHA_DEG\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return math.radians(value)
    return None


def alpha_from_data_mypackage(run_dir: Path) -> float | None:
    for path in (run_dir / "data.mypackage", run_dir.parent / "input" / "data.mypackage"):
        value = read_text_value(path, r"myPa_param1\s*=\s*([+\-0-9.eEdD]+)")
        if value is not None:
            return value
    return None


def infer_alpha_label(run_dir: Path, case_code: str, override: str | None = None) -> str:
    if override:
        return override

    label = alpha_from_run_name(run_dir, case_code)
    if label is not None:
        return label

    value = alpha_from_data_mypackage(run_dir)
    if value is not None:
        return format_alpha_label(value)

    value = alpha_from_gendata(run_dir)
    if value is not None:
        return format_alpha_label(value)

    return run_dir.name


def snapshot_output_dir(case_code: str, alpha_label: str) -> Path:
    return OUTPUT_ROOT / case_output_name(case_code) / "Snapshots" / f"alpha_{alpha_label}"


def diagnosis_output_dir(case_code: str, kind: str, alpha_label: str) -> Path:
    return OUTPUT_ROOT / case_output_name(case_code) / "Diagnosis" / kind / f"alpha_{alpha_label}"


def comparison_output_dir(case_code: str, alpha_1: str, alpha_2: str) -> Path:
    return (
        OUTPUT_ROOT
        / case_output_name(case_code)
        / "Diagnosis"
        / "comparison"
        / f"alpha_{alpha_1}"
        / f"vs_alpha_{alpha_2}"
    )


def write_manifest(output_dir: Path, payload: dict[str, object]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def save_figure_variants(
    figure: object,
    output_path: Path,
    dpi: int | None = None,
    formats: tuple[str, ...] = ("pdf", "png"),
) -> list[Path]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    base_path = output_path.with_suffix("")
    saved_paths: list[Path] = []

    for fmt in formats:
        target = base_path.parent / f"{base_path.name}.{fmt}"
        save_kwargs: dict[str, object] = {"bbox_inches": "tight"}
        if dpi is not None:
            save_kwargs["dpi"] = dpi
        figure.savefig(target, **save_kwargs)
        saved_paths.append(target)

    return saved_paths


def resolve_color_limits(
    field: object,
    vmin: float | None = None,
    vmax: float | None = None,
    symmetric: bool = False,
) -> tuple[float, float]:
    values = np.asarray(field, dtype=np.float64)
    finite = np.isfinite(values)

    if np.any(finite):
        finite_values = values[finite]
        data_min = float(np.min(finite_values))
        data_max = float(np.max(finite_values))
    else:
        data_min = 0.0
        data_max = 1.0

    if symmetric:
        candidates: list[float] = []
        if vmin is not None:
            candidates.append(abs(float(vmin)))
        if vmax is not None:
            candidates.append(abs(float(vmax)))
        if vmin is None or vmax is None:
            candidates.extend((abs(data_min), abs(data_max)))

        limit = max(candidates) if candidates else 1.0
        if not math.isfinite(limit) or limit == 0.0:
            limit = 1.0
        return -limit, limit

    resolved_min = float(vmin) if vmin is not None else data_min
    resolved_max = float(vmax) if vmax is not None else data_max
    if resolved_min == resolved_max:
        resolved_max = resolved_min + 1.0
    return resolved_min, resolved_max
