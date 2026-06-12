#!/usr/bin/env python3
"""Generate the GitHub Pages homepage under docs/index.html."""
# %%
from __future__ import annotations

import html
import csv
import re
import shutil
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parents[2]
SANDBOX_DIR = REPO_DIR / "Sandbox"
DOCS_DIR = REPO_DIR / "docs"
ASSET_ROOT = DOCS_DIR / "assets" / "williamson"
FRAGMENT_DIR = DOCS_DIR / "fragments"
SANDBOX_OUTPUT_ROOT = SANDBOX_DIR / "output"
SITE_TITLE = "SHALLOW WATER MITGCM VERIFICATION CASES"

SECTIONS = [
    {
        "slug": "case1_constant_bathymetry",
        "title": "Gaussian Patch Free Surface Initialization with Constant Bathymetry",
        "media": [
            {
                "label": "Prognostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Constant_Bathymetry_d_-4000_m" / "prognostic.webm",
            },
            {
                "label": "Diagnostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Constant_Bathymetry_d_-4000_m" / "diagnostic.webm",
            },
        ],
    },
    {
        "slug": "case2_real_bathymetry",
        "title": "Gaussian Patch Free Surface Initialization with Real Bathymetry",
        "media": [
            {
                "label": "Prognostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Real_Bathymetry" / "prognostic.webm",
            },
            {
                "label": "Diagnostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Gaussian_Free-Surface_Patch_with_Real_Bathymetry" / "diagnostic.webm",
            },
        ],
    },
    {
        "slug": "case3_geostrophic_adjustment",
        "title": "Geostrophic Free Surface Initialization with Real Bathymetry",
        "media": [
            {
                "label": "Prognostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Geostrophic_Adjustment_with_Real_Bathymetry" / "prognostic.webm",
            },
            {
                "label": "Diagnostic",
                "kind": "video",
                "source": DOCS_DIR / "videos" / "Geostrophic_Adjustment_with_Real_Bathymetry" / "diagnostic.webm",
            },
        ],
    },
    {
        "slug": "testcase1",
        "title": "Advection of Cosine Bell",
        "cases": [
            {
                "label": "TC1 passive tracer",
                "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC1",
                "snapshot_fields": [
                    ("tracer", "Passive Tracer"),
                    ("etan", "ETAN"),
                    ("psi", "PsiVEL"),
                    ("phi", "PhiVEL"),
                ],
            },
            # {
            #     "label": "TC1 prime free surface",
            #     "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC1_prime",
            #     "snapshot_fields": [
            #         ("eta", "Eta"),
            #         ("etan", "ETAN"),
            #         ("psi", "PsiVEL"),
            #         ("phi", "PhiVEL"),
            #         ("velocity_magnitude", "Velocity Magnitude"),
            #     ],
            # },
        ],
        "media": [],
        "dynamic_gallery": "tc1",
    },
    # {
    #     "slug": "testcase2",
    #     "title": "Global Steady-State Nonlinear Zonal Geostrophic Flow",
    #     "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC2",
    #     "media": [],
    #     "dynamic_gallery": "tc2",
    # },
    # {
    #     "slug": "testcase3",
    #     "title": "Steady Geostrophic Flow with Compact Support",
    #     "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC3",
    #     "media": [],
    # },
]

DATA_COLUMNS = [
    "rhoConst",
    "gravity",
    "deltaT",
    "nTimeSteps",
    "dumpFreq",
    "monitorFreq",
    "delR",
    "delX",
    "delY",
]
SIZE_COLUMNS = ["sNx", "sNy", "nPx", "nPy", "Nx", "Ny", "mpi_ranks"]
DATA_COLUMN_META = {
    "rhoConst": ("rhoConst", "kg m<sup>-3</sup>"),
    "gravity": ("gravity", "m s<sup>-2</sup>"),
    "deltaT": ("deltaT", "s"),
    "nTimeSteps": ("nTimeSteps", "steps"),
    "dumpFreq": ("dumpFreq", "s"),
    "monitorFreq": ("monitorFreq", "s"),
    "delR": ("delR", "m"),
    "delX": ("delX", "deg"),
    "delY": ("delY", "deg"),
}
SIZE_COLUMN_META = {
    "sNx": ("sNx", "cells"),
    "sNy": ("sNy", "cells"),
    "nPx": ("nPx", "MPI tiles"),
    "nPy": ("nPy", "MPI tiles"),
    "Nx": ("Nx", "cells"),
    "Ny": ("Ny", "cells"),
    "mpi_ranks": ("mpi_ranks", "ranks"),
}
DAY_PATTERN = re.compile(r"_day_([0-9]+(?:\.[0-9]+)?)")
KEY_SNAPSHOT_DAYS = (0.0, 3.0, 6.0, 9.0, 12.0)
CASE_OUTPUT_NAMES = {
    "TC1_prime": "TestCase1Prime",
    "TC1": "TestCase1",
    "TC2": "TestCase2",
    "TC3": "TestCase3",
}
SNAPSHOT_FIELD_LABELS = {
    "tracer": "Passive Tracer",
    "eta": "Eta",
    "etan": "ETAN",
    "psi": "PsiVEL",
    "phi": "PhiVEL",
    "velocity_magnitude": "Velocity Magnitude",
}
SNAPSHOT_FIELD_UNITS = {
    "tracer": "m",
    "eta": "m",
    "etan": "m",
    "psi": "m<sup>3</sup> s<sup>-1</sup>",
    "phi": "m<sup>2</sup> s<sup>-1</sup>",
    "velocity_magnitude": "m s<sup>-1</sup>",
}
SNAPSHOT_FIELD_ORDER = {
    field: index
    for index, field in enumerate(("tracer", "eta", "etan", "psi", "phi", "velocity_magnitude"))
}
FIELD_TAB_ORDER = (
    ("etan", "ETAN"),
    ("tracer", "Tracer"),
    ("psi", "PsiVEL"),
    ("phi", "PhiVEL"),
)
COPIED_ASSETS: set[Path] = set()
MISSING_SNAPSHOT_ROOTS: set[Path] = set()

def parse_data_file(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values: dict[str, str] = {}
    for key in DATA_COLUMNS:
        match = re.search(rf"\b{re.escape(key)}\s*=\s*([^,\n/]+)", text, re.IGNORECASE)
        values[key] = match.group(1).strip() if match else "—"
    return values

def parse_size_file(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values: dict[str, int] = {}
    for key in ("sNx", "sNy", "nSx", "nSy", "nPx", "nPy"):
        match = re.search(rf"\b{re.escape(key)}\s*=\s*([0-9]+)", text)
        values[key] = int(match.group(1)) if match else 0

    nx = values["sNx"] * values["nSx"] * values["nPx"]
    ny = values["sNy"] * values["nSy"] * values["nPy"]
    mpi_ranks = values["nPx"] * values["nPy"]
    return {
        "sNx": str(values["sNx"]),
        "sNy": str(values["sNy"]),
        "nPx": str(values["nPx"]),
        "nPy": str(values["nPy"]),
        "Nx": str(nx),
        "Ny": str(ny),
        "mpi_ranks": str(mpi_ranks),
    }

def relative_to_root(path: Path, root: Path) -> Path | None:
    try:
        return path.relative_to(root)
    except ValueError:
        return None

def ensure_docs_asset(source: Path, slug: str) -> str:
    source = source.resolve()
    docs_relative = relative_to_root(source, DOCS_DIR.resolve())
    if docs_relative is not None:
        return docs_relative.as_posix()

    sandbox_output_relative = relative_to_root(source, SANDBOX_OUTPUT_ROOT.resolve())
    if sandbox_output_relative is not None:
        target = ASSET_ROOT / sandbox_output_relative
    else:
        sandbox_relative = relative_to_root(source, SANDBOX_DIR.resolve())
        if sandbox_relative is not None:
            target = ASSET_ROOT / sandbox_relative
        else:
            target = ASSET_ROOT / slug / source.name

    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    COPIED_ASSETS.add(target)
    return target.relative_to(DOCS_DIR).as_posix()

def natural_key(text: str) -> list[str]:
    return [part.zfill(12) if part.isdigit() else part for part in re.split(r"(\d+)", text)]

def path_sort_key(path: Path) -> tuple[object, ...]:
    day = extract_day(path)
    day_key = f"{day:012.6f}" if day is not None else ""
    return (*natural_key(path.as_posix()), day_key)

def extract_day(path: Path) -> float | None:
    match = DAY_PATTERN.search(path.name)
    return float(match.group(1)) if match else None

def format_number(value: float) -> str:
    return f"{value:.2f}".rstrip("0").rstrip(".")

def day_caption(path: Path) -> str:
    day = extract_day(path)
    return "Day unknown" if day is None else f"Day {format_number(day)}"

def path_preference(path: Path) -> tuple[int, int]:
    return (1 if "_iter" in path.name else 0, len(path.name))

def alpha_from_name(name: str) -> str:
    if name.startswith("alpha_"):
        return name.removeprefix("alpha_")
    if name.startswith("run_alpha_"):
        return name.removeprefix("run_alpha_")
    return name

def label_from_token(token: str) -> str:
    return SNAPSHOT_FIELD_LABELS.get(token, token.replace("_", " ").title())

def unit_from_token(token: str) -> str:
    return SNAPSHOT_FIELD_UNITS.get(token, "")

def label_with_unit(label: str, unit: str) -> str:
    label_html = html.escape(label)
    if not unit:
        return label_html
    return f"{label_html} <span class='unit'>[{unit}]</span>"

def infer_case_output_name(case_dir: Path) -> str | None:
    case_name = case_dir.name.lower()
    for case_code, output_name in sorted(CASE_OUTPUT_NAMES.items(), key=lambda item: len(item[0]), reverse=True):
        if case_code.lower() in case_name:
            return output_name
    return None

def section_case_configs(section: dict[str, object]) -> list[dict[str, object]]:
    configured_cases = section.get("cases")
    if configured_cases:
        return [dict(case) for case in configured_cases]  # type: ignore[arg-type]

    case_dir = section.get("case_dir")
    if case_dir is None:
        return []

    case: dict[str, object] = {"case_dir": case_dir}
    for key in ("label", "output_name", "snapshot_fields"):
        if key in section:
            case[key] = section[key]
    return [case]

def case_display_label(case: dict[str, object]) -> str:
    label = case.get("label")
    if label:
        return str(label)
    case_dir = case.get("case_dir")
    return Path(case_dir).name if case_dir is not None else ""

def case_snapshot_root(case: dict[str, object]) -> Path | None:
    case_dir = case.get("case_dir")
    if case_dir is None:
        return None
    case_path = Path(case_dir)
    output_name = case.get("output_name") or infer_case_output_name(case_path)
    if output_name is None:
        return None
    return SANDBOX_OUTPUT_ROOT / str(output_name) / "Snapshots"

def normalize_snapshot_field(field: object) -> tuple[str, str]:
    if isinstance(field, dict):
        folder = str(field["folder"])
        return folder, str(field.get("label", label_from_token(folder)))
    if isinstance(field, (list, tuple)):
        folder = str(field[0])
        label = str(field[1]) if len(field) > 1 else label_from_token(folder)
        return folder, label
    folder = str(field)
    return folder, label_from_token(folder)

def snapshot_field_sort_key(field: tuple[str, str]) -> tuple[int, list[str]]:
    folder, label = field
    return (SNAPSHOT_FIELD_ORDER.get(folder, len(SNAPSHOT_FIELD_ORDER)), natural_key(label))

def discover_snapshot_fields(snapshot_root: Path) -> list[tuple[str, str]]:
    fields: set[str] = set()
    for alpha_dir in snapshot_root.glob("alpha_*"):
        if not alpha_dir.is_dir():
            continue
        for field_dir in alpha_dir.iterdir():
            if field_dir.is_dir() and any(field_dir.glob("*.png")):
                fields.add(field_dir.name)
    return sorted(((field, label_from_token(field)) for field in fields), key=snapshot_field_sort_key)

def case_snapshot_fields(case: dict[str, object], snapshot_root: Path) -> list[tuple[str, str]]:
    configured_fields = case.get("snapshot_fields")
    if configured_fields:
        return [normalize_snapshot_field(field) for field in configured_fields]  # type: ignore[union-attr]
    return discover_snapshot_fields(snapshot_root)

def reduce_to_one_image_per_day(paths: list[Path]) -> list[Path]:
    by_day: dict[float | str, Path] = {}
    for path in paths:
        key: float | str = extract_day(path)
        if key is None:
            key = path.name
        old_path = by_day.get(key)
        if old_path is None or path_preference(path) > path_preference(old_path):
            by_day[key] = path
    return sorted(by_day.values(), key=path_sort_key)

def is_key_snapshot_day(path: Path) -> bool:
    day = extract_day(path)
    return day is not None and any(abs(day - key_day) < 0.01 for key_day in KEY_SNAPSHOT_DAYS)

def case_snapshot_items(snapshot_root: Path, folder: str, field_label: str) -> list[dict[str, object]]:
    items: list[dict[str, object]] = []
    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        field_dir = alpha_dir / folder
        if not field_dir.is_dir():
            continue
        alpha = alpha_from_name(alpha_dir.name)
        paths = sorted(field_dir.glob("*.png"), key=path_sort_key)
        items.extend(
            {
                "source": path,
                "caption": f"alpha {alpha}, {field_label}, {day_caption(path)}",
            }
            for path in reduce_to_one_image_per_day(paths)
        )
    return items

def case_alpha_snapshot_items(
    snapshot_root: Path,
    alpha_dir: Path,
    folder: str,
    field_label: str,
) -> list[dict[str, object]]:
    field_dir = alpha_dir / folder
    if not field_dir.is_dir():
        return []

    alpha = alpha_from_name(alpha_dir.name)
    paths = [
        path
        for path in reduce_to_one_image_per_day(sorted(field_dir.glob("*.png"), key=path_sort_key))
        if is_key_snapshot_day(path)
    ]
    return [
        {
            "source": path,
            "caption": f"alpha {alpha}, {field_label}, {day_caption(path)}",
        }
        for path in paths
    ]

def render_plot_figure(item: dict[str, object], slug: str, figure_class: str = "subfigure") -> str:
    source = Path(item["source"])
    if not source.exists():
        return ""
    relative_path = ensure_docs_asset(source, slug)
    caption = html.escape(str(item["caption"]))
    return (
        f"<figure class='{html.escape(figure_class)}'>"
        "<button class='plot-open' type='button' "
        f"data-modal-src='{html.escape(relative_path)}' "
        f"data-modal-caption='{caption}'>"
        f"<img src='{html.escape(relative_path)}' alt='{caption}' loading='lazy' />"
        "</button>"
        f"<figcaption>{caption}</figcaption>"
        "</figure>"
    )

def render_snapshot_grid(items: list[dict[str, object]], slug: str) -> str:
    figures = "".join(render_plot_figure(item, slug) for item in items)
    if not figures:
        return "<p class='empty'>No key-day snapshots available.</p>"
    return f"<div class='snapshot-grid'>{figures}</div>"

def render_field_units_note(fields: list[tuple[str, str]]) -> str:
    entries = [
        label_with_unit(label, unit_from_token(folder))
        for folder, label in fields
        if unit_from_token(folder)
    ]
    if not entries:
        return ""
    return f"<p class='field-units'>Units: {', '.join(entries)}.</p>"

def render_case_snapshot_browser(case: dict[str, object], slug: str) -> str:
    snapshot_root = case_snapshot_root(case)
    if snapshot_root is None:
        return ""
    if not snapshot_root.exists():
        MISSING_SNAPSHOT_ROOTS.add(snapshot_root)
        return ""

    configured_fields = dict(case_snapshot_fields(case, snapshot_root))
    alpha_dirs = [
        alpha_dir
        for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name))
        if alpha_dir.is_dir()
    ]
    if not alpha_dirs:
        return ""

    alpha_blocks: list[str] = []
    for alpha_index, alpha_dir in enumerate(alpha_dirs):
        alpha = alpha_from_name(alpha_dir.name)
        tab_base = f"{slug}-{alpha_dir.name.replace('.', '_')}"
        buttons: list[str] = []
        panels: list[str] = []
        first_panel = True

        for folder, button_label in FIELD_TAB_ORDER:
            field_label = configured_fields.get(folder, label_from_token(folder))
            items = case_alpha_snapshot_items(snapshot_root, alpha_dir, folder, field_label)
            if not items:
                continue

            tab_id = f"{tab_base}-{folder}"
            active_class = " is-active" if first_panel else ""
            selected = "true" if first_panel else "false"
            hidden = "" if first_panel else " hidden"
            button_label_html = label_with_unit(button_label, unit_from_token(folder))
            buttons.append(
                "<button class='field-tab"
                f"{active_class}' type='button' data-tab-target='{html.escape(tab_id)}' "
                f"aria-selected='{selected}'>{button_label_html}</button>"
            )
            panels.append(
                f"<div class='tab-panel{active_class}' data-tab-panel='{html.escape(tab_id)}'{hidden}>"
                f"{render_snapshot_grid(items, slug)}"
                "</div>"
            )
            first_panel = False

        if not panels:
            continue

        buttons.append(
            f"<a class='field-tab field-tab-link' href='#{html.escape(slug)}-errors'>Error</a>"
        )
        alpha_blocks.append(
            f"<details class='alpha-panel' {'open' if alpha_index == 0 else ''}>"
            f"<summary>&alpha;={html.escape(alpha)}</summary>"
            "<div class='field-tabs'>"
            f"<div class='tab-list'>{''.join(buttons)}</div>"
            f"{''.join(panels)}"
            "</div>"
            "</details>"
        )

    if not alpha_blocks:
        return ""

    days = ", ".join(format_number(day) for day in KEY_SNAPSHOT_DAYS)
    label = html.escape(case_display_label(case))
    field_units_note = render_field_units_note(case_snapshot_fields(case, snapshot_root))
    return (
        "<div class='media-section snapshot-browser'>"
        f"<h3>{label}: key-day snapshots</h3>"
        f"<p class='section-note'>Rendered days: {html.escape(days)}.</p>"
        f"{field_units_note}"
        f"{''.join(alpha_blocks)}"
        "</div>"
    )

def render_case_snapshot_galleries(case: dict[str, object], slug: str) -> str:
    return render_case_snapshot_browser(case, slug)

def comparison_caption(path: Path) -> str:
    parts = list(path.parts)
    for part in parts:
        match = re.match(r"compare_run_alpha_([^_]+)_vs_run_alpha_([^/]+)$", part)
        if match:
            return f"alpha {match.group(1)} vs alpha {match.group(2)}"

    if "comparison" in parts:
        index = parts.index("comparison")
        if index + 2 < len(parts):
            alpha_1 = alpha_from_name(parts[index + 1])
            alpha_2 = alpha_from_name(parts[index + 2].removeprefix("vs_"))
            return f"alpha {alpha_1} vs alpha {alpha_2}"

    return "TC1 tracer comparison"

def final_day_from_error_table(error_dir: Path) -> float | None:
    table = error_dir / "tc1_error_table.csv"
    if not table.exists():
        return None
    with table.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    value = rows[-1].get("day")
    return float(value) if value not in (None, "") else None

def tc1_snapshot_items() -> list[dict[str, object]]:
    snapshot_root = SANDBOX_OUTPUT_ROOT / "TestCase1" / "Snapshots"
    items: list[dict[str, object]] = []
    if not snapshot_root.exists():
        return items

    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        alpha = alpha_from_name(alpha_dir.name)
        for folder in ("tracer", "eta"):
            paths = sorted((alpha_dir / folder).glob("*.png"), key=path_sort_key)
            if paths:
                by_day: dict[float | str, Path] = {}
                for path in paths:
                    key: float | str = extract_day(path)
                    if key is None:
                        key = path.name
                    old_path = by_day.get(key)
                    if old_path is None or path_preference(path) > path_preference(old_path):
                        by_day[key] = path
                items.extend(
                    {
                        "source": path,
                        "caption": f"alpha {alpha}, {day_caption(path)}",
                    }
                    for path in sorted(by_day.values(), key=path_sort_key)
                )
                break
    return items

def tc1_tracer_overlay_items() -> list[dict[str, object]]:
    central = (
        SANDBOX_OUTPUT_ROOT
        / "TestCase1"
        / "Diagnosis"
        / "comparison"
    ).glob("alpha_*/vs_alpha_*/passive_tracer_overlay/*.png")
    legacy = (
        SANDBOX_DIR
        / "vortexSphere_Williamson_TC1"
        / "Scripts"
        / "output"
    ).glob("compare_*/PassiveTracerOverlay/tc1_tracer_overlay_day_*.png")
    paths = sorted([*central, *legacy], key=path_sort_key)
    by_case_day: dict[tuple[str, float | str], Path] = {}
    for path in paths:
        day_key: float | str = extract_day(path)
        if day_key is None:
            day_key = path.name
        key = (comparison_caption(path), day_key)
        old_path = by_case_day.get(key)
        if old_path is None or path_preference(path) > path_preference(old_path):
            by_case_day[key] = path
    return [
        {
            "source": path,
            "caption": f"{comparison_caption(path)}, {day_caption(path)}",
        }
        for path in sorted(by_case_day.values(), key=path_sort_key)
    ]

def tc1_error_contour_items() -> list[dict[str, object]]:
    items: list[dict[str, object]] = []
    seen_alpha: set[str] = set()

    central_paths = sorted(
        (
            SANDBOX_OUTPUT_ROOT
            / "TestCase1"
            / "Diagnosis"
            / "error"
        ).glob("alpha_*/tc1_signed_error_contours.png"),
        key=path_sort_key,
    )
    legacy_paths = sorted(
        (
            SANDBOX_DIR
            / "vortexSphere_Williamson_TC1"
            / "Scripts"
            / "output"
        ).glob("run_alpha_*/TC1Error/tc1_signed_error_contours.png"),
        key=path_sort_key,
    )

    for path in [*central_paths, *legacy_paths]:
        alpha = alpha_from_name(path.parent.name)
        if path.parent.name == "TC1Error":
            alpha = alpha_from_name(path.parent.parent.name)
        if alpha in seen_alpha:
            continue
        seen_alpha.add(alpha)
        final_day = final_day_from_error_table(path.parent)
        day_text = f", Day {format_number(final_day)}" if final_day is not None else ""
        items.append(
            {
                "source": path,
                "caption": f"alpha {alpha}, signed error contours{day_text}",
            }
        )
    return items

def final_day_from_tc2_table(error_dir: Path) -> float | None:
    table = error_dir / "TC2_error_norms.csv"
    if not table.exists():
        return None
    with table.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    value = rows[-1].get("day")
    return float(value) if value not in (None, "") else None

def tc2_snapshot_items() -> list[dict[str, object]]:
    snapshot_root = SANDBOX_OUTPUT_ROOT / "TestCase2" / "Snapshots"
    items: list[dict[str, object]] = []
    if not snapshot_root.exists():
        return items

    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        alpha = alpha_from_name(alpha_dir.name)
        for folder, label in (("eta", "eta"), ("etan", "ETAN")):
            paths = sorted((alpha_dir / folder).glob("*.png"), key=path_sort_key)
            if not paths:
                continue
            by_day: dict[float | str, Path] = {}
            for path in paths:
                key: float | str = extract_day(path)
                if key is None:
                    key = path.name
                old_path = by_day.get(key)
                if old_path is None or path_preference(path) > path_preference(old_path):
                    by_day[key] = path
            items.extend(
                {
                    "source": path,
                    "caption": f"alpha {alpha}, {label}, {day_caption(path)}",
                }
                for path in sorted(by_day.values(), key=path_sort_key)
            )
            break
    return items

def tc2_error_norm_items() -> list[dict[str, object]]:
    error_root = SANDBOX_OUTPUT_ROOT / "TestCase2" / "Diagnosis" / "error"
    items: list[dict[str, object]] = []
    if not error_root.exists():
        return items

    for path in sorted(error_root.glob("alpha_*/TC2_error_norms.png"), key=path_sort_key):
        alpha = alpha_from_name(path.parent.name)
        final_day = final_day_from_tc2_table(path.parent)
        day_text = f", final day {format_number(final_day)}" if final_day is not None else ""
        items.append(
            {
                "source": path,
                "caption": f"alpha {alpha}, error norms{day_text}",
            }
        )
    return items

def render_column_header(column: str, column_meta: dict[str, tuple[str, str]]) -> str:
    label, unit = column_meta.get(column, (column, ""))
    unit_html = f"<span class='table-unit'>[{unit}]</span>" if unit else ""
    return f"<th><span class='table-symbol'>{html.escape(label)}</span>{unit_html}</th>"

def render_table(
    title: str,
    columns: list[str],
    values: dict[str, str],
    column_meta: dict[str, tuple[str, str]],
) -> str:
    header = "".join(render_column_header(column, column_meta) for column in columns)
    row = "".join(f"<td>{html.escape(values.get(column, '—'))}</td>" for column in columns)
    return (
        f"<div class='table-block'><h3>{html.escape(title)}</h3>"
        f"<table><thead><tr>{header}</tr></thead><tbody><tr>{row}</tr></tbody></table></div>"
    )

def render_case_tables(case: dict[str, object], show_label: bool) -> str:
    case_dir = case.get("case_dir")
    if case_dir is None:
        return ""

    case_path = Path(case_dir)
    data_path = case_path / "input" / "data"
    size_path = case_path / "code" / "SIZE.h"
    if not data_path.exists() or not size_path.exists():
        return ""

    data_values = parse_data_file(data_path)
    size_values = parse_size_file(size_path)
    heading = (
        f"<h3 class='case-heading'>{html.escape(case_display_label(case))}</h3>"
        if show_label
        else ""
    )
    return (
        "<details class='details-panel case-block'>"
        "<summary>Numerical settings</summary>"
        f"{heading}"
        "<div class='tables'>"
        f"{render_table('Numerical Settings', DATA_COLUMNS, data_values, DATA_COLUMN_META)}"
        f"{render_table('Grid And MPI Layout', SIZE_COLUMNS, size_values, SIZE_COLUMN_META)}"
        "</div>"
        "</details>"
    )

def count_case_key_snapshots(case: dict[str, object]) -> tuple[int, int, int]:
    snapshot_root = case_snapshot_root(case)
    if snapshot_root is None or not snapshot_root.exists():
        return 0, 0, 0

    alphas: set[str] = set()
    fields: set[str] = set()
    count = 0
    configured_fields = dict(case_snapshot_fields(case, snapshot_root))
    for alpha_dir in sorted(snapshot_root.glob("alpha_*"), key=lambda path: natural_key(path.name)):
        if not alpha_dir.is_dir():
            continue
        alpha_has_items = False
        for folder in configured_fields:
            paths = [
                path
                for path in reduce_to_one_image_per_day(sorted((alpha_dir / folder).glob("*.png"), key=path_sort_key))
                if is_key_snapshot_day(path)
            ]
            if paths:
                fields.add(folder)
                count += len(paths)
                alpha_has_items = True
        if alpha_has_items:
            alphas.add(alpha_from_name(alpha_dir.name))
    return len(alphas), len(fields), count

def render_metric_cards(section: dict[str, object], case_configs: list[dict[str, object]]) -> str:
    metrics: list[tuple[str, str]] = []
    alpha_count = 0
    field_count = 0
    snapshot_count = 0
    for case in case_configs:
        case_alpha_count, case_field_count, case_snapshot_count = count_case_key_snapshots(case)
        alpha_count += case_alpha_count
        field_count = max(field_count, case_field_count)
        snapshot_count += case_snapshot_count

    error_count = len(tc1_error_contour_items()) if section.get("dynamic_gallery") == "tc1" else 0
    media_count = sum(1 for media in section.get("media", []) if Path(media["source"]).exists())
    if snapshot_count:
        metrics = [
            ("Alphas", str(alpha_count)),
            ("Fields", str(field_count)),
            ("Key Days", ", ".join(format_number(day) for day in KEY_SNAPSHOT_DAYS)),
            ("Error Plots", str(error_count)),
        ]
    else:
        metrics = [
            ("Main Panels", str(media_count)),
            ("Primary", "Prognostic / Diagnostic" if media_count >= 2 else "Linked media"),
            ("Format", "WEBM" if media_count else "None"),
            ("Details", "Collapsible"),
        ]

    cards = "".join(
        "<div class='metric-card'>"
        f"<span>{html.escape(label)}</span>"
        f"<strong>{html.escape(value)}</strong>"
        "</div>"
        for label, value in metrics
    )
    return f"<div class='metric-grid'>{cards}</div>"

def render_experiment_summary(section: dict[str, object]) -> str:
    summary = str(section.get("summary", "")).strip()
    if not summary:
        return ""
    return f"<p class='experiment-summary'>{html.escape(summary)}</p>"

def render_media_card(media: dict[str, object], slug: str) -> str:
    source = Path(media["source"])
    relative_path = ensure_docs_asset(source, slug)
    label = html.escape(str(media["label"]))
    kind = str(media["kind"])

    if kind == "video":
        body = (
            f"<video controls preload='metadata'>"
            f"<source src='{html.escape(relative_path)}' type='video/webm'>"
            "</video>"
        )
    else:
        body = f"<img src='{html.escape(relative_path)}' alt='{label}' />"

    return (
        "<div class='media-card'>"
        f"<h4>{label}</h4>"
        f"{body}"
        "</div>"
    )

def render_subfigure_gallery(title: str, items: list[dict[str, object]], slug: str) -> str:
    figures = []
    for item in items:
        figure = render_plot_figure(item, slug)
        if figure:
            figures.append(figure)

    if not figures:
        return ""

    return (
        "<div class='media-section'>"
        f"<h3>{html.escape(title)}</h3>"
        "<div class='subfigure-grid'>"
        f"{''.join(figures)}"
        "</div>"
        "</div>"
    )

def render_error_comparison_block(title: str, items: list[dict[str, object]], slug: str) -> str:
    figures = []
    for item in items:
        figure = render_plot_figure(item, slug, "subfigure comparison-figure")
        if figure:
            figures.append(figure)
    if not figures:
        return ""
    return (
        f"<div id='{html.escape(slug)}-errors' class='media-section error-block'>"
        f"<h3>{html.escape(title)}</h3>"
        "<div class='comparison-grid'>"
        f"{''.join(figures)}"
        "</div>"
        "</div>"
    )

def render_dynamic_gallery(section: dict[str, object], slug: str) -> str:
    gallery = section.get("dynamic_gallery")
    if gallery == "tc1":
        return render_error_comparison_block("Error Contours", tc1_error_contour_items(), slug)
    if gallery == "tc2":
        return "".join(
            [
                render_subfigure_gallery("Steady-State Error Norms", tc2_error_norm_items(), slug),
            ]
        )
    return ""

def render_section(section: dict[str, object]) -> str:
    slug = str(section["slug"])
    title = html.escape(str(section["title"]))
    case_configs = section_case_configs(section)
    metrics_html = render_metric_cards(section, case_configs)
    summary_html = render_experiment_summary(section)
    tables_html = "".join(
        render_case_tables(case, show_label=len(case_configs) > 1)
        for case in case_configs
    )

    media_items = [
        render_media_card(media, slug)
        for media in section.get("media", [])
        if Path(media["source"]).exists()
    ]
    media_blocks = []
    if media_items:
        media_blocks.append(f"<div class='main-panels'>{''.join(media_items)}</div>")
    for case in case_configs:
        snapshot_html = render_case_snapshot_galleries(case, slug)
        if snapshot_html:
            media_blocks.append(snapshot_html)
    dynamic_html = render_dynamic_gallery(section, slug)
    if dynamic_html:
        media_blocks.append(dynamic_html)
    media_html = "".join(media_blocks) if media_blocks else "<p class='empty'>No media linked yet for this section.</p>"

    return (
        f"<section id='{html.escape(slug)}' class='section-card'>"
        "<header class='experiment-top'>"
        "<div>"
        "<span class='section-label'>Experiment</span>"
        f"<h2>{title}</h2>"
        f"{summary_html}"
        "</div>"
        f"{metrics_html}"
        "</header>"
        f"{tables_html}"
        f"{media_html}"
        "</section>"
    )

def section_fragment_path(slug: str) -> Path:
    return FRAGMENT_DIR / f"{slug}.html"

def write_section_fragment(section: dict[str, object]) -> Path:
    slug = str(section["slug"])
    FRAGMENT_DIR.mkdir(parents=True, exist_ok=True)
    path = section_fragment_path(slug)
    path.write_text(render_section(section) + "\n", encoding="utf-8")
    return path

def write_section_fragments(slugs: list[str] | tuple[str, ...] | None = None) -> list[Path]:
    selected = set(slugs) if slugs is not None else None
    if selected is not None:
        known = {str(section["slug"]) for section in SECTIONS}
        unknown = sorted(selected - known)
        if unknown:
            raise ValueError(f"unknown GitHub Pages section slug(s): {', '.join(unknown)}")

    paths: list[Path] = []
    for section in SECTIONS:
        slug = str(section["slug"])
        if selected is None or slug in selected or not section_fragment_path(slug).exists():
            paths.append(write_section_fragment(section))
    return paths

def render_fragment_inputs() -> str:
    inputs: list[str] = []
    for section in SECTIONS:
        slug = str(section["slug"])
        path = section_fragment_path(slug).relative_to(DOCS_DIR).as_posix()
        inputs.append(f"{{% include_relative {path} %}}")
    return "\n      ".join(inputs)

def render_navigation() -> str:
    links = []
    for section in SECTIONS:
        slug = html.escape(str(section["slug"]))
        title = html.escape(str(section["title"]))
        links.append(f"      <a href='#{slug}'>{title}</a>")
    return "\n".join(
        [
            "<aside class='side-nav'>",
            "      <strong>Experiments</strong>",
            *links,
            "    </aside>",
        ]
    )

def build_html() -> str:
    nav_html = render_navigation()
    sections_html = render_fragment_inputs()
    return f"""---
---
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{SITE_TITLE}</title>
  <link rel="stylesheet" href="site.css">
</head>
<body>
  <header class="site-header">
    <h1>{SITE_TITLE}</h1>
    <p class="intro">Verification homepage for the shallow-water MITgcm experiments. Use the section navigation to inspect current results.</p>
  </header>
  <div class="site-layout">
    {nav_html}
    <main>
      {sections_html}
    </main>
  </div>
  <div class="plot-modal" data-modal hidden>
    <button class="modal-close" type="button" data-modal-close aria-label="Close enlarged plot">x</button>
    <figure>
      <img data-modal-img alt="">
      <figcaption data-modal-caption></figcaption>
    </figure>
  </div>
  <script src="site.js" defer></script>
</body>
</html>
"""

def build_site(section_slugs: list[str] | tuple[str, ...] | None = None) -> tuple[Path, list[Path]]:
    COPIED_ASSETS.clear()
    MISSING_SNAPSHOT_ROOTS.clear()
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_ROOT.mkdir(parents=True, exist_ok=True)
    FRAGMENT_DIR.mkdir(parents=True, exist_ok=True)
    fragment_paths = write_section_fragments(section_slugs)
    index_path = DOCS_DIR / "index.html"
    index_path.write_text(build_html(), encoding="utf-8")
    return index_path, fragment_paths

def main() -> None:
    index_path, fragment_paths = build_site()
    print(f"wrote {index_path}")
    for path in fragment_paths:
        print(f"wrote {path}")
    print(f"copied/updated {len(COPIED_ASSETS)} asset files")
    if MISSING_SNAPSHOT_ROOTS:
        print("missing snapshot roots:")
        for path in sorted(MISSING_SNAPSHOT_ROOTS, key=lambda item: item.as_posix()):
            print(f"  {path.relative_to(REPO_DIR)}")

if __name__ == "__main__":
    main()

# %%
