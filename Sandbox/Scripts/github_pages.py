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
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC1",
        "media": [],
        "dynamic_gallery": "tc1",
    },
    {
        "slug": "testcase2",
        "title": "Global Steady-State Nonlinear Zonal Geostrophic Flow",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC2",
        "media": [],
        "dynamic_gallery": "tc2",
    },
    {
        "slug": "testcase3",
        "title": "Steady Geostrophic Flow with Compact Support",
        "case_dir": SANDBOX_DIR / "vortexSphere_Williamson_TC3",
        "media": [],
    },
]

DATA_COLUMNS = ["rhoConst", "gravity", "deltaT", "nTimeSteps", "dumpFreq", "monitorFreq", "delR"]
SIZE_COLUMNS = ["sNx", "sNy", "nPx", "nPy", "Nx", "Ny", "mpi_ranks"]
DAY_PATTERN = re.compile(r"_day_([0-9]+(?:\.[0-9]+)?)")


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


def render_table(title: str, columns: list[str], values: dict[str, str]) -> str:
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    row = "".join(f"<td>{html.escape(values.get(column, '—'))}</td>" for column in columns)
    return (
        f"<div class='table-block'><h3>{html.escape(title)}</h3>"
        f"<table><thead><tr>{header}</tr></thead><tbody><tr>{row}</tr></tbody></table></div>"
    )


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
        source = Path(item["source"])
        if not source.exists():
            continue
        relative_path = ensure_docs_asset(source, slug)
        caption = html.escape(str(item["caption"]))
        figures.append(
            "<figure class='subfigure'>"
            f"<img src='{html.escape(relative_path)}' alt='{caption}' loading='lazy' />"
            f"<figcaption>{caption}</figcaption>"
            "</figure>"
        )

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


def render_dynamic_gallery(section: dict[str, object], slug: str) -> str:
    gallery = section.get("dynamic_gallery")
    if gallery == "tc1":
        return "".join(
            [
                render_subfigure_gallery("Passive Tracer Snapshots", tc1_snapshot_items(), slug),
                render_subfigure_gallery("Passive Tracer Overlay Comparisons", tc1_tracer_overlay_items(), slug),
                render_subfigure_gallery("Error Contours", tc1_error_contour_items(), slug),
            ]
        )
    if gallery == "tc2":
        return "".join(
            [
                render_subfigure_gallery("Eta Snapshots", tc2_snapshot_items(), slug),
                render_subfigure_gallery("Steady-State Error Norms", tc2_error_norm_items(), slug),
            ]
        )
    return ""


def render_section(section: dict[str, object]) -> str:
    slug = str(section["slug"])
    title = html.escape(str(section["title"]))
    tables_html = ""

    case_dir = section.get("case_dir")
    if case_dir is not None:
        case_path = Path(case_dir)
        data_values = parse_data_file(case_path / "input" / "data")
        size_values = parse_size_file(case_path / "code" / "SIZE.h")
        tables_html = (
            "<div class='tables'>"
            f"{render_table('Numerical Settings', DATA_COLUMNS, data_values)}"
            f"{render_table('Grid And MPI Layout', SIZE_COLUMNS, size_values)}"
            "</div>"
        )

    media_items = [
        render_media_card(media, slug)
        for media in section.get("media", [])
        if Path(media["source"]).exists()
    ]
    media_blocks = []
    if media_items:
        media_blocks.append(f"<div class='media-grid'>{''.join(media_items)}</div>")
    dynamic_html = render_dynamic_gallery(section, slug)
    if dynamic_html:
        media_blocks.append(dynamic_html)
    media_html = "".join(media_blocks) if media_blocks else "<p class='empty'>No media linked yet for this section.</p>"

    return (
        f"<section id='{html.escape(slug)}' class='section-card'>"
        f"<h2>{title}</h2>"
        f"{tables_html}"
        f"{media_html}"
        "</section>"
    )


def render_navigation() -> str:
    cards = []
    for section in SECTIONS:
        slug = html.escape(str(section["slug"]))
        title = html.escape(str(section["title"]))
        cards.append(
            "<a class='nav-card' "
            f"href='#{slug}'>"
            f"<span class='nav-kicker'>Experiment</span>"
            f"<strong>{title}</strong>"
            "</a>"
        )
    return "<section class='nav-grid'>" + "".join(cards) + "</section>"


def build_html() -> str:
    nav_html = render_navigation()
    sections_html = "".join(render_section(section) for section in SECTIONS)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{SITE_TITLE}</title>
  <style>
    :root {{
      --ink: #122033;
      --paper: #f6f7fb;
      --card: #ffffff;
      --accent: #0a6c74;
      --accent-soft: #dff1f3;
      --line: #d6dbe6;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(10,108,116,0.16), transparent 30%),
        linear-gradient(180deg, #eef4f7 0%, var(--paper) 100%);
    }}
    main {{
      max-width: 1200px;
      margin: 0 auto;
      padding: 48px 20px 64px;
    }}
    h1 {{
      margin: 0 0 12px;
      font-size: clamp(2rem, 4vw, 3.2rem);
      letter-spacing: 0.04em;
    }}
    .intro {{
      margin: 0 0 28px;
      font-size: 1.05rem;
      max-width: 70ch;
    }}
    .nav-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 14px;
      margin: 0 0 28px;
    }}
    .nav-card {{
      display: block;
      text-decoration: none;
      color: var(--ink);
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 16px 18px;
      box-shadow: 0 12px 34px rgba(18, 32, 51, 0.06);
      transition: transform 120ms ease, border-color 120ms ease;
    }}
    .nav-card:hover {{
      transform: translateY(-2px);
      border-color: var(--accent);
    }}
    .nav-kicker {{
      display: block;
      margin-bottom: 6px;
      color: var(--accent);
      font-size: 0.78rem;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }}
    .section-card {{
      scroll-margin-top: 24px;
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 24px;
      margin-bottom: 24px;
      box-shadow: 0 14px 40px rgba(18, 32, 51, 0.08);
    }}
    .tables {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
      gap: 16px;
      margin: 18px 0 20px;
    }}
    .table-block h3 {{
      margin: 0 0 10px;
      font-size: 1rem;
      color: var(--accent);
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.92rem;
      table-layout: fixed;
    }}
    th, td {{
      border: 1px solid var(--line);
      padding: 8px 10px;
      text-align: center;
      word-wrap: break-word;
    }}
    th {{
      background: var(--accent-soft);
    }}
    .media-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 16px;
    }}
    .media-card {{
      border: 1px solid var(--line);
      border-radius: 14px;
      padding: 14px;
      background: #fbfcfe;
    }}
    .media-card h4 {{
      margin: 0 0 10px;
      font-size: 1rem;
    }}
    .media-section {{
      margin-top: 18px;
    }}
    .media-section h3 {{
      margin: 0 0 12px;
      font-size: 1rem;
      color: var(--accent);
    }}
    .subfigure-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 14px;
    }}
    .subfigure {{
      margin: 0;
      border: 1px solid var(--line);
      border-radius: 10px;
      padding: 10px;
      background: #fbfcfe;
    }}
    .subfigure figcaption {{
      margin-top: 8px;
      font-size: 0.9rem;
      text-align: center;
    }}
    video, img {{
      width: 100%;
      display: block;
      border-radius: 10px;
      background: #000;
    }}
    .empty {{
      margin: 0;
      color: #556273;
      font-style: italic;
    }}
  </style>
</head>
<body>
  <main>
    <h1>{SITE_TITLE}</h1>
    <p class="intro">Verification homepage for the shallow-water MITgcm experiments. Use the clickable section cards below to jump to each experiment and inspect its current results.</p>
    {nav_html}
    {sections_html}
  </main>
</body>
</html>
"""


def main() -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_ROOT.mkdir(parents=True, exist_ok=True)
    index_path = DOCS_DIR / "index.html"
    index_path.write_text(build_html(), encoding="utf-8")
    print(f"wrote {index_path}")


if __name__ == "__main__":
    main()

# %%
