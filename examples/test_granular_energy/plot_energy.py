#!/usr/bin/env python3.11
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--avewindow", type=float, required=True)
    parser.add_argument("--title", default="")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    energy_path = outdir / "energy.dat"
    energy = np.loadtxt(energy_path, comments="#")
    if energy.ndim == 1:
        energy = energy.reshape(1, -1)

    step = energy[:, 0]
    time = step * args.dt
    kin = energy[:, 1]
    rot = energy[:, 2]
    contact_strain = energy[:, 3]
    gravity = energy[:, 4]
    wall_norm = energy[:, 5]
    wall_tang = energy[:, 6]

    diss_path = outdir / "dissipation.dat"
    if diss_path.exists():
        diss = np.loadtxt(diss_path, comments="#")
        if diss.ndim == 1:
            diss = diss.reshape(1, -1)
        if energy[0, 0] == 0 and diss[0, 0] != 0:
            zero_diss = np.zeros((1, diss.shape[1]))
            zero_diss[0, 0] = energy[0, 0]
            diss = np.vstack((zero_diss, diss))
        if len(diss) != len(energy) or np.any(diss[:, 0] != step):
            raise ValueError(f"{diss_path} is not aligned with {energy_path}")
        if diss.shape[1] >= 8:
            pair_diss_parts_inc = diss[:, 1:4]
            pair_diss_inc = np.sum(pair_diss_parts_inc, axis=1)
            wall_diss_inc = diss[:, 4]
            pair_atom_parts_rate = diss[:, 5:8]
            pair_atom_rate = np.sum(pair_atom_parts_rate, axis=1)
            if diss.shape[1] >= 11:
                wall_diss_parts_inc = diss[:, 8:11]
            else:
                wall_diss_parts_inc = np.column_stack(
                    (wall_diss_inc, np.zeros_like(wall_diss_inc), np.zeros_like(wall_diss_inc))
                )
        else:
            raise ValueError(f"{diss_path} has {diss.shape[1]} columns; expected at least 8")
    else:
        pair_diss_parts_inc = energy[:, 7:10]
        pair_diss_inc = np.sum(pair_diss_parts_inc, axis=1)
        wall_diss_inc = energy[:, 10]
        pair_atom_parts_rate = energy[:, 11:14]
        pair_atom_rate = np.sum(pair_atom_parts_rate, axis=1)
        wall_diss_parts_inc = np.column_stack(
            (wall_diss_inc, np.zeros_like(wall_diss_inc), np.zeros_like(wall_diss_inc))
        )

    wall_diss = np.cumsum(wall_diss_inc * args.avewindow)
    atom_pair_diss = np.cumsum(pair_atom_rate * args.dt * args.avewindow)
    pair_local_diss = np.cumsum(pair_diss_inc * args.avewindow)
    pair_diss = atom_pair_diss
    mechanical = kin + rot + contact_strain + gravity
    heat_total = pair_diss + wall_diss
    corrected = mechanical + heat_total
    e0 = corrected[0]
    observed_energy_scale = max(
        abs(e0),
        float(np.max(np.abs(corrected))) if corrected.size else 0.0,
        float(np.max(np.abs(mechanical))) if mechanical.size else 0.0,
    )
    energy_scale = observed_energy_scale if observed_energy_scale > 0.0 else 1.0

    summary = {
        "time_unit": "s",
        "energy_unit": "J",
        "rows": int(len(energy)),
        "first_step": int(step[0]),
        "last_step": int(step[-1]),
        "first_time": float(time[0]),
        "last_time": float(time[-1]),
        "initial_energy": float(e0),
        "final_total_energy": float(corrected[-1]),
        "final_relative_drift": float((corrected[-1] - e0) / e0) if e0 else None,
        "max_relative_excursion": float(np.max(np.abs(corrected - e0)) / abs(e0)) if e0 else None,
        "final_mechanical_energy": float(mechanical[-1]),
        "final_pair_dissipation": float(pair_diss[-1]),
        "final_pair_local_dissipation": float(pair_local_diss[-1]),
        "final_wall_dissipation": float(wall_diss[-1]),
        "energy_scale_for_drift": float(energy_scale),
        "pair_local_vs_per_atom_pair_dissipation_mismatch": float(pair_local_diss[-1] - atom_pair_diss[-1]),
        "max_contact_strain_energy": float(np.max(contact_strain)),
        "max_wall_strain_energy": float(np.max(wall_norm + wall_tang)),
        "max_pair_dissipation_increment": float(np.max(pair_diss_inc)),
        "max_wall_dissipation_increment": float(np.max(wall_diss_inc)),
        "min_pair_dissipation_increment": float(np.min(pair_diss_inc)),
        "min_wall_dissipation_increment": float(np.min(wall_diss_inc)),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"matplotlib unavailable, wrote {outdir / 'summary.json'} only: {exc}")
        return

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(time, kin, label="kinetic")
    axes[0].plot(time, rot, label="rotational kinetic")
    axes[0].plot(time, contact_strain, label="contact strain")
    axes[0].plot(time, gravity, label="gravitational")
    axes[0].plot(time, pair_diss, label="pair accumulated dissipation")
    axes[0].plot(time, wall_diss, label="wall accumulated dissipation")
    axes[0].set_ylabel("energy (J)")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(time, np.full_like(time, e0), label="initial")
    axes[1].plot(time, corrected, label="total")
    axes[1].set_ylabel("energy (J)")
    axes[1].legend(loc="best", fontsize=8)
    axes[1].grid(True, alpha=0.3)

    axes[1].set_xlabel("time (s)")

    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(outdir / "energy.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
