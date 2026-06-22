#!/usr/bin/env python3
"""Safe cleaner for the TC1 build and run directories."""

import argparse
import shutil
import sys
from pathlib import Path
from typing import List


RUN_OUTPUT_SUFFIXES = {".data", ".meta", ".txt", ".png", ".pdf", ".jpg", ".jpeg", ".npz", ".csv"}
RUN_OUTPUT_PREFIXES = (
    "dynDiag",
    "dyn_Aux",
    "dynStDiag",
    "Eta.",
    "PH.",
    "PHL.",
    "S.",
    "T.",
    "U.",
    "V.",
    "W.",
    "pickup",
    "core",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Clean TC1 build and/or run directories. Dry-run by default."
    )
    target = parser.add_argument_group("targets")
    target.add_argument("--all", action="store_true", help="clean both build and run directories")
    target.add_argument("--build", action="store_true", help="clean build directories")
    target.add_argument("--runs", action="store_true", help="clean run_alpha_* directories")
    target.add_argument("--logs", action="store_true", help="also clean logs/slurm-* files")
    target.add_argument(
        "--for-git",
        action="store_true",
        help="clean generated build, run, log, and Python cache artifacts before committing or pushing",
    )

    select = parser.add_argument_group("selection")
    select.add_argument(
        "--builds",
        nargs="+",
        default=["all"],
        metavar="NAME",
        help="build subdirectories to clean, e.g. c1 c2; default: all",
    )
    select.add_argument(
        "--alphas",
        nargs="+",
        default=["all"],
        metavar="ALPHA",
        help="run_alpha values to clean, e.g. 0 0.05 1.57; default: all",
    )
    select.add_argument(
        "--run-mode",
        choices=("outputs", "all"),
        default="outputs",
        help="outputs removes model products but keeps inputs/logs/executable; all empties run dirs",
    )

    action = parser.add_argument_group("action")
    action.add_argument("--yes", action="store_true", help="actually delete files")
    action.add_argument("--dry-run", action="store_true", help="preview only; this is the default")
    return parser.parse_args()


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def require_safe_case_path(path: Path, case_dir: Path) -> None:
    path = path.resolve()
    build_dir = (case_dir / "build").resolve()
    logs_dir = (case_dir / "logs").resolve()
    jobs_dir = (case_dir / "jobs").resolve()
    tools_dir = (case_dir / "tools").resolve()
    if is_relative_to(path, build_dir):
        return
    if is_relative_to(path, logs_dir):
        return
    if path.name == "__pycache__" and is_relative_to(path, jobs_dir):
        return
    if path.name == "__pycache__" and is_relative_to(path, tools_dir):
        return
    if path.parent == case_dir and path.name.startswith("run_alpha_"):
        return
    raise SystemExit(f"Refusing unsafe cleanup target: {path}")


def delete_path(path: Path, dry_run: bool) -> None:
    print(f"{'would remove' if dry_run else 'removing'} {path}")
    if dry_run:
        return
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def clean_dir_contents(path: Path, dry_run: bool) -> None:
    for child in sorted(path.iterdir()):
        delete_path(child, dry_run)


def is_run_output(path: Path) -> bool:
    if path.is_dir():
        return path.name == "__pycache__"
    if path.suffix in RUN_OUTPUT_SUFFIXES:
        return True
    return path.name.startswith(RUN_OUTPUT_PREFIXES)


def clean_run_outputs(path: Path, dry_run: bool) -> None:
    for child in sorted(path.iterdir()):
        if is_run_output(child):
            delete_path(child, dry_run)


def select_build_dirs(case_dir: Path, names: List[str]) -> List[Path]:
    build_root = case_dir / "build"
    if not build_root.exists():
        return []
    if names == ["all"]:
        return sorted(p for p in build_root.iterdir() if p.is_dir())
    return [build_root / name for name in names]


def select_run_dirs(case_dir: Path, alphas: List[str]) -> List[Path]:
    if alphas == ["all"]:
        return sorted(p for p in case_dir.glob("run_alpha_*") if p.is_dir())
    return [case_dir / f"run_alpha_{alpha}" for alpha in alphas]


def clean_logs(case_dir: Path, dry_run: bool) -> None:
    logs_dir = case_dir / "logs"
    if not logs_dir.exists():
        return
    for path in sorted(logs_dir.glob("slurm-*")):
        delete_path(path, dry_run)


def clean_python_caches(case_dir: Path, dry_run: bool) -> None:
    for root in (case_dir / "jobs", case_dir / "tools"):
        if not root.exists():
            continue
        for path in sorted(p for p in root.rglob("__pycache__") if p.is_dir()):
            require_safe_case_path(path, case_dir)
            delete_path(path, dry_run)


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    case_dir = script_dir.parent.resolve()
    dry_run = not args.yes

    clean_build = args.for_git or args.all or args.build
    clean_runs = args.for_git or args.all or args.runs
    clean_logs_target = args.for_git or args.logs
    clean_caches = args.for_git
    run_mode = "all" if args.for_git else args.run_mode
    if not (clean_build or clean_runs or clean_logs_target or clean_caches):
        print(
            "Choose at least one target: --build, --runs, --logs, --for-git, or --all.",
            file=sys.stderr,
        )
        return 2

    if dry_run:
        print("DRY RUN: pass --yes to actually delete files.")

    if clean_build:
        for build_dir in select_build_dirs(case_dir, args.builds):
            require_safe_case_path(build_dir, case_dir)
            if build_dir.exists():
                clean_dir_contents(build_dir, dry_run)
            else:
                print(f"missing build directory: {build_dir}")

    if clean_runs:
        for run_dir in select_run_dirs(case_dir, args.alphas):
            require_safe_case_path(run_dir, case_dir)
            if not run_dir.exists():
                print(f"missing run directory: {run_dir}")
                continue
            if run_mode == "all":
                clean_dir_contents(run_dir, dry_run)
            else:
                clean_run_outputs(run_dir, dry_run)

    if clean_logs_target:
        require_safe_case_path(case_dir / "logs", case_dir)
        clean_logs(case_dir, dry_run)

    if clean_caches:
        clean_python_caches(case_dir, dry_run)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
