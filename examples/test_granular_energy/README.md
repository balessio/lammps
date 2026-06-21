# Granular Energy Tests

This folder contains four parameterized LAMMPS tests plus runners:

- `in.two_grains_periodic`: two grains in a small periodic box undergoing an oblique collision.
- `in.settling_wall`: 64 slightly polydisperse grains with random initial velocities settling under gravity between bottom and top granular elastic walls.
- `in.shear_cell_2d`: about 60 grains in a pressure-loaded 2D rough-clump shear cell.
- `in.shear_cell_3d`: about 300 grains in a pressure-loaded 3D rough-clump shear cell.
- `run_energy_tests.sh`: runs the two-grain and settling-wall energy tests across several contact laws and plots energy histories.
- `run_shear_heat_tests.sh`: runs the 2D and 3D shear heat tests for `hertz_mindlin_rescale/force`.
- `plot_energy.py`: writes `summary.json` and `energy.png` for each case.
- `plot_shear_heat.py`: writes `heat_summary.json` and `heat_vs_time.png`.

Run the default sweep:

```sh
DUMP=1 ./run_energy_tests.sh
```

Useful options:

```sh
QUICK=1 DUMP=1 ./run_energy_tests.sh
CASES="elastic_hooke elastic_history hooke_history hertz_history hertz_mindlin hertz_mindlin_force hertz_mindlin_rescale_force" ./run_energy_tests.sh
CASES="hooke_normal_damping hooke_history_damping_only hooke_history_friction_only" ./run_energy_tests.sh
SCENARIOS="two_grains" ./run_energy_tests.sh
PYTHON=python3.11 ./run_energy_tests.sh
```

Run the pressure-loaded shear heat cases:

```sh
DUMP=1 PYTHON=python3.11 ./run_shear_heat_tests.sh
```

The default full sweep runs `300000` two-grain steps at `dt = 1e-7` and `2000000` settling-wall steps at `dt = 1e-7`. Override with `STEPS_TWO`, `STEPS_SETTLE`, `DT_TWO`, `DT_SETTLE`, `DUMP_EVERY`, or `PAIR_CUTOFF` when needed. The default `PAIR_CUTOFF=0.006` is slightly larger than the 0.005 m contact diameter so deleted history contacts remain visible for one-step energy accounting. The `hooke_*_only` cases are diagnostics for separating normal damping, tangential damping, and Coulomb slip.

Each result folder contains `energy.dat`, `dissipation.dat`, `energy.png`,
`summary.json`, the generated `pair.inc`/`wall.inc`, and LAMMPS logs.
`energy.dat` stores instantaneous endpoint mechanical-energy components.
`dissipation.dat` stores window-averaged per-step heat increments/rates plus a
signed residual-work diagnostic. For wall cases, `compute pe` includes wall
strain energy via `fix_modify wall energy yes`; the wall normal/tangential
columns remain diagnostics from `fix wall/gran contacts`. `energy.png` reports
energy in J, time in s, and includes residual work normalized by the initial or
largest observed total energy scale.

The shear heat runs first prepare a compacted packing from a no-contact top
platen gap and write `packed.restart` plus inspectable `packed.data`. The shear
phase then reads the restart, resets the heat counters, and plots shear heat
only. The nonrotating kinematic rough-clump platens are pressure controlled in
the vertical direction, while the horizontal velocity imposes total shear 1.
`heat_vs_time.png` is generated heat in J versus time in s. Dump files include
per-grain cumulative heat in J, per-grain dissipation rate in W, the four
dissipation-rate components in W, and particle temperature in K.
