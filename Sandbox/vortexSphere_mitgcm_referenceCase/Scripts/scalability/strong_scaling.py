#%%
"""Strong-scaling plot for the Octopus MITgcm case."""

from __future__ import annotations

import sys
from datetime import timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from scripts import set_style  # noqa: E402

RANKS = np.array([8, 12, 16, 24, 32, 48], dtype=float)

WALL_TIMES = [
    "03:17:18",
    "03:10:31",
    "02:44:48",
    "02:19:05",
    "01:43:18",
    "01:02:44",
]

CASE_TITLE = "Strong Scaling of the 10-Day Octopus MITgcm Run"
CASE_NOTE = "Fixed workload: 10 days, dt = 60 s, 14,400 iterations"
REFERENCE_NOTE = "Partition used: large, nodes onode[13-16]"

OUTPUT_DIR = Path(__file__).resolve().parent / "results"
OUTPUT_STEM = "octopus_strong_scaling"

SUPTITLE_SIZE = 18
TITLE_SIZE = 16
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 16
LEGEND_SIZE = 16
ANNOTATION_SIZE = 16
FOOTNOTE_TEXT_SIZE = 16

def configure_plotting() -> None:
    """Apply shared publication-style plot settings."""

    set_style()

    plt.rcParams.update(
        {
            "font.size": TICK_LABEL_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "figure.titlesize": SUPTITLE_SIZE,
            "axes.titlepad": 10,
            "axes.labelpad": 6,
            "axes.grid": True,
            "grid.alpha": 0.30,
            "grid.linestyle": "--",
            "grid.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "legend.handlelength": 2.0,
            "legend.handletextpad": 0.6,
            "legend.borderaxespad": 0.4,
            "legend.labelspacing": 0.3,
            "lines.linewidth": 3.0,
            "lines.markersize": 9,
            "lines.markeredgewidth": 1.1,
        }
    )

def parse_hhmmss(text: str) -> float:
    """Convert HH:MM:SS string to seconds."""

    hours_str, minutes_str, seconds_str = text.split(":")

    duration = timedelta(
        hours=int(hours_str),
        minutes=int(minutes_str),
        seconds=int(seconds_str),
    )

    return duration.total_seconds()

def format_hhmmss(seconds: float) -> str:
    """Format seconds as HH:MM:SS."""

    total_seconds = int(round(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)

    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def annotate_runtime_points(
    ax: plt.Axes,
    ranks: np.ndarray,
    wall_hours: np.ndarray,
    labels: list[str],
    color: str,
) -> None:
    """Add compact runtime labels without overlapping the plotted data."""

    label_offsets = [
        (0, 12, "center", "bottom"),
        (0, -16, "center", "top"),
        (0, 12, "center", "bottom"),
        (0, -16, "center", "top"),
        (0, 12, "center", "bottom"),
        (0, 12, "center", "bottom"),
    ]

    for rank, hours, label, offset in zip(
        ranks,
        wall_hours,
        labels,
        label_offsets,
        strict=True,
    ):
        dx, dy, ha, va = offset

        ax.annotate(
            label,
            xy=(rank, hours),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha,
            va=va,
            fontsize=ANNOTATION_SIZE,
            color=color,
            bbox=dict(
                boxstyle="round,pad=0.18",
                fc="white",
                ec="0.85",
                lw=0.4,
                alpha=0.90,
            ),
            clip_on=False,
        )

def make_scaling_figure() -> tuple[plt.Figure, np.ndarray]:
    """Build the runtime, speedup, and cost panels."""

    configure_plotting()

    wall_seconds = np.array(
        [parse_hhmmss(value) for value in WALL_TIMES],
        dtype=float,
    )
    wall_hours = wall_seconds / 3600.0

    reference_rank = RANKS[0]
    reference_time = wall_seconds[0]

    ideal_runtime_hours = (reference_time * reference_rank / RANKS) / 3600.0

    speedup = reference_time / wall_seconds
    ideal_speedup = RANKS / reference_rank

    cost_core_hours = RANKS * wall_hours
    ideal_cost_core_hours = np.full_like(
        cost_core_hours,
        reference_rank * reference_time / 3600.0,
    )

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(13.2, 14.4),
        sharex=True,
    )

    fig.suptitle(
        CASE_TITLE,
        fontsize=SUPTITLE_SIZE,
        fontweight="bold",
        y=0.985,
    )

    fig.text(
        0.5,
        0.015,
        f"{CASE_NOTE} | {REFERENCE_NOTE}",
        ha="center",
        va="bottom",
        fontsize=FOOTNOTE_TEXT_SIZE,
        color="dimgray",
    )

    runtime_color = "#1f77b4"
    speedup_color = "#9467bd"
    cost_color = "#2ca02c"
    ideal_color = "#4D4D4D"

    # ------------------------------------------------------------
    # Panel 1: Wall-clock runtime
    # ------------------------------------------------------------
    ax = axes[0]

    ax.plot(
        RANKS,
        wall_hours,
        marker="o",
        color=runtime_color,
        label="Measured runtime",
    )

    ax.plot(
        RANKS,
        ideal_runtime_hours,
        linestyle="--",
        color=ideal_color,
        label="Perfect 1/p",
    )

    annotate_runtime_points(
        ax=ax,
        ranks=RANKS,
        wall_hours=wall_hours,
        labels=WALL_TIMES,
        color=runtime_color,
    )

    ax.set_ylabel("Wall time [h]")
    ax.set_title("Wall-clock time")
    ax.set_ylim(0.0, wall_hours.max() * 1.35)
    ax.legend(loc="upper right")

    # ------------------------------------------------------------
    # Panel 2: Speedup
    # ------------------------------------------------------------
    ax = axes[1]

    ax.plot(
        RANKS,
        speedup,
        marker="o",
        color=speedup_color,
        label="Measured speedup",
    )

    ax.plot(
        RANKS,
        ideal_speedup,
        linestyle="--",
        color=ideal_color,
        label="Perfect speedup",
    )

    ax.set_ylabel("Speedup")
    ax.set_title("Strong-scaling speedup")
    ax.set_ylim(0.75, ideal_speedup.max() * 1.05)
    ax.legend(loc="upper left")

    # ------------------------------------------------------------
    # Panel 3: Computational cost
    # ------------------------------------------------------------
    ax = axes[2]

    ax.plot(
        RANKS,
        cost_core_hours,
        marker="o",
        color=cost_color,
        label="Measured cost",
    )

    ax.plot(
        RANKS,
        ideal_cost_core_hours,
        linestyle="--",
        color=ideal_color,
        label="Constant core-hours",
    )

    ax.set_ylabel("Core-hours")
    ax.set_title("Computational cost")
    ax.set_ylim(
        cost_core_hours.min() * 0.95,
        cost_core_hours.max() * 1.05,
    )
    ax.legend(loc="upper left")

    # ------------------------------------------------------------
    # Shared axis formatting
    # ------------------------------------------------------------
    for ax in axes:
        ax.grid(True, which="major", alpha=0.30)
        ax.grid(True, which="minor", alpha=0.18, linestyle=":")
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_SIZE)

    axes[-1].set_xlabel("MPI ranks")
    axes[-1].set_xticks(RANKS)
    axes[-1].set_xlim(RANKS.min() - 2, RANKS.max() + 2)

    fig.tight_layout(rect=[0, 0.06, 1, 0.965], h_pad=2.0)

    print("\nStrong-scaling summary")
    print(f"{'ranks':>6} {'wall time':>12} {'speedup':>10} {'core-hours':>12}")

    for rank, wall_time, spd, cost in zip(
        RANKS,
        wall_seconds,
        speedup,
        cost_core_hours,
        strict=True,
    ):
        print(
            f"{int(rank):>6d} "
            f"{format_hhmmss(wall_time):>12} "
            f"{spd:>10.2f} "
            f"{cost:>12.2f}"
        )

    return fig, axes

def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    """Save the figure as PNG and PDF."""

    output_dir.mkdir(parents=True, exist_ok=True)

    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"

    fig.savefig(png_path, dpi=400, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=400, bbox_inches="tight")

    print(f"\nSaved figure to:\n  {png_path}\n  {pdf_path}")

def main() -> None:
    fig, _ = make_scaling_figure()
    save_figure(fig, OUTPUT_DIR, OUTPUT_STEM)
    plt.close(fig)

if __name__ == "__main__":
    main()

# %%

