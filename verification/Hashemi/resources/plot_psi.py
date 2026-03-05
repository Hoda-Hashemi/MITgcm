#!/usr/bin/env python3
#%%
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Hashemi PSI output from PsiOcean.csv."
    )
    parser.add_argument(
        "--input",
        default="verification/Hashemi/output/Fortran/solution/PsiOcean.csv",
        help="Path to PsiOcean.csv",
    )
    parser.add_argument(
        "--outdir",
        default="verification/Hashemi/results",
        help="Directory to write figures and stats",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap for heatmap",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively in addition to saving",
    )
    return parser.parse_args()

def load_psi_csv(path: Path) -> tuple[np.ndarray, int, int]:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")

    rows: list[tuple[int, int, float]] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header row.")

        required = {"iLat", "iLon", "psi"}
        if not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"CSV header must contain {sorted(required)}; got {reader.fieldnames}"
            )

        for index, row in enumerate(reader, start=2):
            try:
                i_lat = int(row["iLat"])
                i_lon = int(row["iLon"])
                psi = float(row["psi"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid value at CSV line {index}: {row}"
                ) from exc
            if i_lat < 1 or i_lon < 1:
                raise ValueError(
                    f"indices must be 1-based positive at line {index}: {row}"
                )
            rows.append((i_lat, i_lon, psi))

    if not rows:
        raise ValueError("Input CSV is empty (no data rows).")

    n_lat = max(item[0] for item in rows)
    n_lon = max(item[1] for item in rows)

    psi_grid = np.full((n_lat, n_lon), np.nan, dtype=float)
    for i_lat, i_lon, psi in rows:
        psi_grid[i_lat - 1, i_lon - 1] = psi

    return psi_grid, n_lat, n_lon

def write_stats(path: Path, psi_grid: np.ndarray, n_lat: int, n_lon: int) -> None:
    finite = psi_grid[np.isfinite(psi_grid)]
    missing = np.count_nonzero(~np.isfinite(psi_grid))
    text = "\n".join(
        [
            "PsiOcean statistics",
            f"nLat={n_lat}",
            f"nLon={n_lon}",
            f"samples={finite.size}",
            f"missing_cells={missing}",
            f"min={np.min(finite):.16e}",
            f"max={np.max(finite):.16e}",
            f"mean={np.mean(finite):.16e}",
            f"std={np.std(finite):.16e}",
        ]
    )
    path.write_text(text + "\n")

def main() -> int:
    args = parse_args()

    if not args.show:
        import matplotlib

        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    psi_grid, n_lat, n_lon = load_psi_csv(input_path)
    finite = psi_grid[np.isfinite(psi_grid)]

    heatmap_path = outdir / "psi_ocean_heatmap.png"
    hist_path = outdir / "psi_ocean_hist.png"
    stats_path = outdir / "psi_ocean_stats.txt"

    fig, ax = plt.subplots(figsize=(10, 5))
    image = ax.imshow(
        psi_grid,
        origin="lower",
        aspect="auto",
        cmap=args.cmap,
        interpolation="nearest",
    )
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("psi")
    ax.set_title("PsiOcean Field")
    ax.set_xlabel("iLon")
    ax.set_ylabel("iLat")
    fig.tight_layout()
    fig.savefig(heatmap_path, dpi=180)

    fig_hist, ax_hist = plt.subplots(figsize=(8, 4.5))
    ax_hist.hist(finite, bins=80)
    ax_hist.set_title("PsiOcean Distribution")
    ax_hist.set_xlabel("psi")
    ax_hist.set_ylabel("count")
    fig_hist.tight_layout()
    fig_hist.savefig(hist_path, dpi=180)

    write_stats(stats_path, psi_grid, n_lat, n_lon)

    print(f"Wrote {heatmap_path}")
    print(f"Wrote {hist_path}")
    print(f"Wrote {stats_path}")

    if args.show:
        plt.show()

    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

# %%
