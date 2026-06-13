#!/usr/bin/env python3
"""Write persistent run inventories for MITgcm sandbox jobs."""
from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_DIR = Path(__file__).resolve().parents[2]
SANDBOX_DIR = REPO_DIR / "Sandbox"
DOCS_DIR = REPO_DIR / "docs"
OUTPUT_ROOT = SANDBOX_DIR / "output"
DOCS_ASSET_ROOT = DOCS_DIR / "assets" / "williamson"
DEFAULT_LOG_NAME = "run_inventory.log"

CASE_OUTPUT_NAMES = {
    "TC1_prime": "TestCase1Prime",
    "TC1": "TestCase1",
    "TC2": "TestCase2",
    "TC3": "TestCase3",
}
RUN_LABELS = {
    "TC1_prime": {"c1": "0", "c2": "0.05", "c3": "1.52", "c4": "1.57"},
    "TC1": {"c1": "0", "c2": "0.05", "c3": "1.52", "c4": "1.57"},
    "TC2": {"c1": "0", "c2": "0.05", "c3": "1.52", "c4": "1.57"},
    "TC3": {"c1": "0", "c2": "1.0472"},
}

ITERATED_STEM_RE = re.compile(r"^(?P<stem>.+)\.(?P<iteration>[0-9]{10})$")
REFERENCE_RE = re.compile(r"""(?:src|href|data-modal-src)=['"]([^'"]+)['"]""")
PLOT_EXTENSIONS = {".png", ".pdf", ".svg", ".webm", ".mp4", ".mov"}
COMPUTED_EXTENSIONS = {".csv", ".txt", ".tex", ".json", ".nc", ".npy", ".npz"}


@dataclass(frozen=True)
class DiagnosticMeta:
    name: str
    levels: str
    mate: str
    code: str
    units: str
    title: str


@dataclass(frozen=True)
class DiagnosticStream:
    kind: str
    stream: str
    file_name: str
    frequency: str
    phase: str
    fields: tuple[str, ...]


@dataclass(frozen=True)
class RunOutputGroup:
    stem: str
    data_files: int
    meta_files: int
    iterations: tuple[int, ...]
    bytes_total: int


@dataclass(frozen=True)
class RunProduct:
    output_dir: Path
    manifest_path: Path
    manifest: dict[str, Any]
    files: tuple[Path, ...]


def natural_key(text: str) -> tuple[object, ...]:
    return tuple(int(part) if part.isdigit() else part for part in re.split(r"([0-9]+)", text))


def relpath(path: Path, root: Path = REPO_DIR) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def same_path(left: Path, right: Path) -> bool:
    return left.expanduser().resolve(strict=False) == right.expanduser().resolve(strict=False)


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


def alpha_from_run_name(run_dir: Path, case_code: str) -> str | None:
    direct_match = re.search(r"run_alpha_([0-9.]+)", run_dir.name)
    if direct_match:
        return direct_match.group(1)

    case_match = re.search(r"alpha_(c\d+)", run_dir.name)
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


def human_bytes(value: int) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    amount = float(value)
    for unit in units:
        if abs(amount) < 1024.0 or unit == units[-1]:
            return f"{amount:.1f} {unit}" if unit != "B" else f"{int(amount)} B"
        amount /= 1024.0
    return f"{value} B"


def utc_time_from_timestamp(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, tz=timezone.utc).isoformat(timespec="seconds")


def strip_namelist_comments(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        quoted = False
        out: list[str] = []
        for char in line:
            if char == "'":
                quoted = not quoted
                out.append(char)
                continue
            if not quoted and char in ("#", "!"):
                break
            out.append(char)
        lines.append("".join(out))
    return "\n".join(lines)


def parse_namelist_assignments(text: str) -> list[tuple[str, str | None, str]]:
    clean = strip_namelist_comments(text)
    pattern = re.compile(
        r"(?ms)^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(\([^=\n]*\))?\s*=\s*"
        r"(.*?)(?=^\s*[A-Za-z_][A-Za-z0-9_]*\s*(?:\([^=\n]*\))?\s*=|^\s*&|\Z)"
    )
    return [
        (match.group(1), match.group(2), match.group(3).strip())
        for match in pattern.finditer(clean)
    ]


def stream_index(index_text: str | None) -> str | None:
    if not index_text:
        return None
    parts = [part.strip() for part in index_text.strip("()").split(",") if part.strip()]
    return parts[-1] if parts else None


def quoted_values(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in re.findall(r"'([^']*)'", value) if item.strip())


def scalar_value(value: str) -> str:
    text = value.strip().rstrip(",").strip()
    strings = quoted_values(text)
    if strings:
        return strings[0]
    parts = [part.strip() for part in text.split(",") if part.strip()]
    return parts[0] if parts else ""


def diagnostics_file_for_run(run_dir: Path) -> Path | None:
    for path in (run_dir / "data.diagnostics", run_dir.parent / "input" / "data.diagnostics"):
        if path.exists():
            return path
    return None


def available_diagnostics_file_for_run(run_dir: Path) -> Path | None:
    path = run_dir / "available_diagnostics.log"
    return path if path.exists() else None


def parse_available_diagnostics(path: Path | None) -> dict[str, DiagnosticMeta]:
    if path is None or not path.exists():
        return {}

    diagnostics: dict[str, DiagnosticMeta] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) < 7 or not parts[0].isdigit():
            continue
        meta = DiagnosticMeta(
            name=parts[1],
            levels=parts[2],
            mate=parts[3],
            code=parts[4],
            units=parts[5],
            title=parts[6],
        )
        diagnostics[meta.name] = meta
    return diagnostics


def parse_requested_diagnostics(data_path: Path | None) -> tuple[DiagnosticStream, ...]:
    if data_path is None or not data_path.exists():
        return ()

    streams: dict[str, dict[str, Any]] = defaultdict(dict)
    stat_streams: dict[str, dict[str, Any]] = defaultdict(dict)

    for key, index_text, value in parse_namelist_assignments(data_path.read_text(encoding="utf-8", errors="ignore")):
        index = stream_index(index_text)
        if index is None:
            continue

        key_lower = key.lower()
        if key_lower == "fields":
            streams[index]["fields"] = quoted_values(value)
        elif key_lower in {"filename", "frequency", "timephase", "averagingfreq", "averagingphase", "repeatcycle"}:
            streams[index][key_lower] = scalar_value(value)
        elif key_lower == "stat_fields":
            stat_streams[index]["fields"] = quoted_values(value)
        elif key_lower in {"stat_fname", "stat_freq", "stat_phase", "stat_region"}:
            stat_streams[index][key_lower] = scalar_value(value)

    parsed: list[DiagnosticStream] = []
    for index, values in sorted(streams.items(), key=lambda item: natural_key(item[0])):
        parsed.append(
            DiagnosticStream(
                kind="diagnostics",
                stream=index,
                file_name=str(values.get("filename", "")),
                frequency=str(values.get("frequency", "")),
                phase=str(values.get("timephase", "")),
                fields=tuple(values.get("fields", ())),
            )
        )
    for index, values in sorted(stat_streams.items(), key=lambda item: natural_key(item[0])):
        parsed.append(
            DiagnosticStream(
                kind="statistics",
                stream=index,
                file_name=str(values.get("stat_fname", "")),
                frequency=str(values.get("stat_freq", "")),
                phase=str(values.get("stat_phase", "")),
                fields=tuple(values.get("fields", ())),
            )
        )
    return tuple(parsed)


def inspect_run_dir_outputs(
    run_dir: Path,
    *,
    exclude_names: tuple[str, ...] = (DEFAULT_LOG_NAME,),
) -> tuple[tuple[RunOutputGroup, ...], tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    groups: dict[str, dict[str, Any]] = defaultdict(lambda: {"data": 0, "meta": 0, "iterations": set(), "bytes": 0})
    file_manifest: list[dict[str, Any]] = []
    directories: list[dict[str, Any]] = []

    for path in sorted(run_dir.iterdir(), key=lambda item: natural_key(item.name)):
        if path.name in exclude_names:
            continue

        stat = path.stat()
        if path.is_dir():
            directories.append({"name": path.name, "mtime_utc": utc_time_from_timestamp(stat.st_mtime)})
            continue

        size = stat.st_size
        file_manifest.append(
            {
                "name": path.name,
                "bytes": size,
                "mtime_utc": utc_time_from_timestamp(stat.st_mtime),
            }
        )

        if path.suffix not in {".data", ".meta"}:
            continue
        base = path.with_suffix("").name
        match = ITERATED_STEM_RE.match(base)
        stem = match.group("stem") if match else base
        if path.suffix == ".data":
            groups[stem]["data"] += 1
        else:
            groups[stem]["meta"] += 1
        if match:
            groups[stem]["iterations"].add(int(match.group("iteration")))
        groups[stem]["bytes"] += size

    output_groups = tuple(
        RunOutputGroup(
            stem=stem,
            data_files=int(values["data"]),
            meta_files=int(values["meta"]),
            iterations=tuple(sorted(values["iterations"])),
            bytes_total=int(values["bytes"]),
        )
        for stem, values in sorted(groups.items(), key=lambda item: natural_key(item[0]))
    )
    return output_groups, tuple(file_manifest), tuple(directories)


def find_generated_products(run_dir: Path, output_root: Path = OUTPUT_ROOT) -> tuple[RunProduct, ...]:
    if not output_root.exists():
        return ()

    products: list[RunProduct] = []
    for manifest_path in sorted(output_root.rglob("manifest.json"), key=lambda item: natural_key(item.as_posix())):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue

        source_run = manifest.get("source_run")
        if not source_run or not same_path(Path(str(source_run)), run_dir):
            continue

        output_dir = manifest_path.parent
        files = tuple(
            path
            for path in sorted(output_dir.rglob("*"), key=lambda item: natural_key(item.as_posix()))
            if path.is_file()
        )
        products.append(RunProduct(output_dir=output_dir, manifest_path=manifest_path, manifest=manifest, files=files))
    return tuple(products)


def read_docs_references(docs_dir: Path = DOCS_DIR) -> dict[str, set[Path]]:
    html_paths: list[Path] = []
    index_path = docs_dir / "index.html"
    if index_path.exists():
        html_paths.append(index_path)
    fragment_dir = docs_dir / "fragments"
    if fragment_dir.exists():
        html_paths.extend(sorted(fragment_dir.glob("*.html"), key=lambda item: natural_key(item.name)))

    references: dict[str, set[Path]] = defaultdict(set)
    for html_path in html_paths:
        text = html_path.read_text(encoding="utf-8", errors="ignore")
        for reference in REFERENCE_RE.findall(text):
            references[reference].add(html_path)
    return references


def docs_asset_for_output(path: Path, output_root: Path = OUTPUT_ROOT, docs_asset_root: Path = DOCS_ASSET_ROOT) -> Path | None:
    try:
        relative = path.resolve().relative_to(output_root.resolve())
    except ValueError:
        return None
    return docs_asset_root / relative


def find_indexed_assets(
    products: tuple[RunProduct, ...],
    *,
    docs_dir: Path = DOCS_DIR,
    output_root: Path = OUTPUT_ROOT,
) -> tuple[dict[str, Any], ...]:
    references = read_docs_references(docs_dir)
    indexed: list[dict[str, Any]] = []

    for product in products:
        for output_file in product.files:
            docs_asset = docs_asset_for_output(output_file, output_root=output_root)
            if docs_asset is None:
                continue
            docs_reference = docs_asset.relative_to(docs_dir).as_posix()
            referring_html = references.get(docs_reference)
            if not referring_html:
                continue
            indexed.append(
                {
                    "output_file": output_file,
                    "docs_asset": docs_asset,
                    "referenced_by": tuple(sorted(referring_html, key=lambda item: natural_key(item.as_posix()))),
                }
            )

    return tuple(indexed)


def classify_product_files(product: RunProduct) -> tuple[tuple[Path, ...], tuple[Path, ...], tuple[Path, ...]]:
    plotted: list[Path] = []
    computed: list[Path] = []
    other: list[Path] = []
    for path in product.files:
        if path.name == "manifest.json":
            other.append(path)
        elif path.suffix.lower() in PLOT_EXTENSIONS:
            plotted.append(path)
        elif path.suffix.lower() in COMPUTED_EXTENSIONS:
            computed.append(path)
        else:
            other.append(path)
    return tuple(plotted), tuple(computed), tuple(other)


def format_iteration_span(iterations: tuple[int, ...]) -> str:
    if not iterations:
        return "none"
    if len(iterations) == 1:
        return str(iterations[0])
    return f"{iterations[0]}..{iterations[-1]} ({len(iterations)})"


def format_run_log(
    run_dir: Path,
    *,
    case_code: str | None = None,
    include_file_manifest: bool = True,
    output_root: Path = OUTPUT_ROOT,
    docs_dir: Path = DOCS_DIR,
    log_name: str = DEFAULT_LOG_NAME,
) -> str:
    run_dir = run_dir.expanduser().resolve()
    if case_code is None:
        try:
            case_code = detect_case_code(run_dir)
        except ValueError:
            case_code = "unknown"
    try:
        alpha_label = infer_alpha_label(run_dir, case_code) if case_code != "unknown" else "unknown"
    except Exception:
        alpha_label = "unknown"

    diagnostic_path = diagnostics_file_for_run(run_dir)
    available_path = available_diagnostics_file_for_run(run_dir)
    diagnostic_meta = parse_available_diagnostics(available_path)
    requested = parse_requested_diagnostics(diagnostic_path)
    groups, file_manifest, directories = inspect_run_dir_outputs(run_dir, exclude_names=(log_name,))
    products = find_generated_products(run_dir, output_root=output_root)
    indexed_assets = find_indexed_assets(products, docs_dir=docs_dir, output_root=output_root)

    run_file_count = len(file_manifest)
    run_dir_bytes = sum(int(item["bytes"]) for item in file_manifest)
    product_files = tuple(path for product in products for path in product.files)
    plotted_files = tuple(path for product in products for path in classify_product_files(product)[0])
    computed_files = tuple(path for product in products for path in classify_product_files(product)[1])

    lines: list[str] = [
        "MITgcm run inventory",
        f"generated_utc: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"run_dir: {relpath(run_dir)}",
        f"case_code: {case_code}",
        f"alpha_label: {alpha_label}",
        f"log_file: {log_name}",
        "",
        "Summary",
        f"  run_dir_files: {run_file_count}",
        f"  run_dir_bytes: {run_dir_bytes} ({human_bytes(run_dir_bytes)})",
        f"  run_output_groups: {len(groups)}",
        f"  requested_diagnostic_streams: {len(requested)}",
        f"  generated_product_dirs: {len(products)}",
        f"  generated_product_files: {len(product_files)}",
        f"  plotted_files: {len(plotted_files)}",
        f"  computed_files: {len(computed_files)}",
        f"  docs_indexed_assets: {len(indexed_assets)}",
        "",
        "Diagnostics request",
        f"  data_diagnostics: {relpath(diagnostic_path) if diagnostic_path else 'not found'}",
        f"  available_diagnostics: {relpath(available_path) if available_path else 'not found'}",
    ]

    if requested:
        for stream in requested:
            lines.append(
                f"  [{stream.kind} stream {stream.stream}] "
                f"fileName={stream.file_name or 'not set'} frequency={stream.frequency or 'not set'} "
                f"phase={stream.phase or 'not set'}"
            )
            for field in stream.fields:
                meta = diagnostic_meta.get(field.strip())
                if meta is None:
                    lines.append(f"    {field.strip()}")
                else:
                    lines.append(f"    {meta.name}: units={meta.units}; levels={meta.levels}; title={meta.title}")
    else:
        lines.append("  no requested diagnostic streams parsed")

    lines.extend(["", "Run directory outputs"])
    if groups:
        for group in groups:
            lines.append(
                f"  {group.stem}: data={group.data_files} meta={group.meta_files} "
                f"iterations={format_iteration_span(group.iterations)} bytes={group.bytes_total} "
                f"({human_bytes(group.bytes_total)})"
            )
    else:
        lines.append("  no .data/.meta outputs found")

    if directories:
        lines.extend(["", "Run directory subdirectories"])
        for directory in directories:
            lines.append(f"  {directory['name']} mtime_utc={directory['mtime_utc']}")

    lines.extend(["", "Generated products"])
    if products:
        for product in products:
            plotted, computed, other = classify_product_files(product)
            manifest = product.manifest
            lines.append(f"  output_dir: {relpath(product.output_dir)}")
            lines.append(
                "    "
                f"case={manifest.get('case', 'unknown')} "
                f"product={manifest.get('product', 'unknown')} "
                f"alpha={manifest.get('alpha', 'unknown')}"
            )
            saved = manifest.get("saved", ())
            if saved:
                lines.append(f"    saved: {', '.join(str(item) for item in saved)}")
            fields = manifest.get("fields")
            if fields:
                lines.append(f"    fields: {json.dumps(fields, sort_keys=True)}")
            lines.append(f"    plotted_files: {len(plotted)}")
            for path in plotted:
                lines.append(f"      {relpath(path)}")
            lines.append(f"    computed_files: {len(computed)}")
            for path in computed:
                lines.append(f"      {relpath(path)}")
            if other:
                lines.append(f"    other_files: {len(other)}")
                for path in other:
                    lines.append(f"      {relpath(path)}")
    else:
        lines.append("  no product manifest matched this run")

    lines.extend(["", "Sent to docs index"])
    if indexed_assets:
        for item in indexed_assets:
            refs = ", ".join(relpath(path) for path in item["referenced_by"])
            lines.append(f"  output: {relpath(item['output_file'])}")
            lines.append(f"    docs_asset: {relpath(item['docs_asset'])}")
            lines.append(f"    referenced_by: {refs}")
    else:
        lines.append("  no generated assets from this run are referenced by docs fragments")

    if include_file_manifest:
        lines.extend(["", "Run directory file manifest"])
        if file_manifest:
            for item in file_manifest:
                lines.append(
                    f"  {item['name']} bytes={item['bytes']} ({human_bytes(int(item['bytes']))}) "
                    f"mtime_utc={item['mtime_utc']}"
                )
        else:
            lines.append("  no files found")

    return "\n".join(lines) + "\n"


def write_run_log(
    run_dir: Path,
    *,
    case_code: str | None = None,
    log_name: str = DEFAULT_LOG_NAME,
    include_file_manifest: bool = True,
    output_root: Path = OUTPUT_ROOT,
    docs_dir: Path = DOCS_DIR,
) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")
    log_path = run_dir / log_name
    log_text = format_run_log(
        run_dir,
        case_code=case_code,
        include_file_manifest=include_file_manifest,
        output_root=output_root,
        docs_dir=docs_dir,
        log_name=log_name,
    )
    log_path.write_text(log_text, encoding="utf-8")
    return log_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write persistent MITgcm run inventory logs.")
    parser.add_argument("run_dirs", nargs="+", type=Path, help="Run directory or directories to inspect.")
    parser.add_argument("--case-code", default=None, help="Case code override, for example TC1 or TC2.")
    parser.add_argument("--log-name", default=DEFAULT_LOG_NAME, help=f"Output log name inside each run dir.")
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT, help="Sandbox product output root.")
    parser.add_argument("--docs-dir", type=Path, default=DOCS_DIR, help="GitHub Pages docs directory.")
    parser.add_argument("--no-file-manifest", action="store_true", help="Omit full per-file run_dir manifest.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for run_dir in args.run_dirs:
        path = write_run_log(
            run_dir,
            case_code=args.case_code,
            log_name=args.log_name,
            include_file_manifest=not args.no_file_manifest,
            output_root=args.output_root,
            docs_dir=args.docs_dir,
        )
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
