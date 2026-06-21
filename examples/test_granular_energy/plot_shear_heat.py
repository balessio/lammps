#!/usr/bin/env python3.11
"""Plot generated heat for pressure-loaded granular shear-cell tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_table(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--avewindow", type=int, required=True)
    parser.add_argument("--compact-steps", type=int, required=True)
    parser.add_argument("--title", default="granular shear heat")
    args = parser.parse_args()

    heat = load_table(args.outdir / "heat.dat")
    step = heat[:, 0]
    time = step * args.dt
    window_time = args.dt * args.avewindow
    compact_time = args.compact_steps * args.dt

    heat_rate = heat[:, 1] + heat[:, 2] + heat[:, 3]
    phase = heat[:, 7]

    total_heat = np.cumsum(heat_rate * window_time)
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    ax.plot(time, total_heat, label="heat")
    if compact_time > 0.0:
        ax.axvline(compact_time, color="0.3", lw=1.0, ls="--", label="shear start")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("generated heat (J)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(args.outdir / "heat_vs_time.png", dpi=180)
    plt.close(fig)

    shear_mask = phase >= 0.5
    if compact_time == 0.0:
        heat_at_shear = 0.0
    else:
        heat_at_shear = float(total_heat[shear_mask][0]) if np.any(shear_mask) else float(total_heat[-1])
    summary = {
        "time_unit": "s",
        "heat_unit": "J",
        "heat_rate_unit": "W",
        "compact_steps": args.compact_steps,
        "final_step": int(step[-1]),
        "compact_time_s": float(compact_time),
        "final_time_s": float(time[-1]),
        "shear_duration_s": float(time[-1] - compact_time),
        "final_heat": float(total_heat[-1]),
        "heat_at_shear_start": heat_at_shear,
        "shear_heat": float(total_heat[-1] - heat_at_shear),
    }
    (args.outdir / "heat_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
