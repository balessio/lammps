
Relevant variables in `src/GRANULAR/granular_model.h:74-95`:

```cpp
int beyond_contact, limit_damping, history_update, synchronized_verlet, dissipative_heat;
double heat_norm_damp, heat_tang_damp, heat_tang_fric;
double dq_conduct, dq_dissipate;
double dq_damp_hold, dq_friction_hold;
double StrainEnergyNorm, StrainEnergyTang;
```

- `dissipative_heat` turns heat accounting on.
- `heat_norm_damp`, `heat_tang_damp`, and `heat_tang_fric` say what fraction of
  each loss becomes heat.
- `StrainEnergyNorm` and `StrainEnergyTang` are contact strain energy.
- `dq_dissipate` is the total heat made by one contact during one step.
- `dq_damp_hold` and `dq_friction_hold` are temporary tangential heat pieces.




`src/GRANULAR/granular_model.cpp:292-299` rejects heat scale factors outside `[0,1]`:

```cpp
if (heat_norm_damp < 0.0 || heat_norm_damp > 1.0)
  error->all(...);
if (heat_tang_damp < 0.0 || heat_tang_damp > 1.0)
  error->all(...);
if (heat_tang_fric < 0.0 || heat_tang_fric > 1.0)
  error->all(...);
```

For restarts they're written in `src/GRANULAR/granular_model.cpp:379-383` and read in `src/GRANULAR/granular_model.cpp:419-428`.





Every force calculation clears old values first
(`src/GRANULAR/granular_model.cpp:464-471`):

```cpp
StrainEnergyNorm = 0.0;
StrainEnergyTang = 0.0;
dq_conduct = 0.0;
dq_dissipate = 0.0;
dq_damp_hold = 0.0;
dq_friction_hold = 0.0;
```

then heat is added in `src/GRANULAR/granular_model.cpp:530-532`:

```cpp
dq_dissipate = damping_model->calculate_heat();
dq_dissipate += tangential_model->calculate_heat();
```





Pair contacts store three heat values in neighbor history. The size is set in `src/GRANULAR/pair_granular.cpp:46`:

```cpp
static constexpr int ENERGY_HISTORY_SIZE = 3;
```

The extra history starts at `energy_history_offset` (`src/GRANULAR/pair_granular.cpp:519-522`):

```cpp
energy_history_offset = size_history;
size_history += ENERGY_HISTORY_SIZE;
```

After `model->calculate_forces()`, the three heat values are saved in `src/GRANULAR/pair_granular.cpp:284-290`:

```cpp
history[energy_history_offset] = model->svector[0];
history[energy_history_offset + 1] = model->svector[1];
history[energy_history_offset + 2] = model->svector[2];
```

and are put into `pair/local p13 p14 p15` in `src/GRANULAR/pair_granular.cpp:892-895`:

```cpp
svector[12] = history[energy_history_offset];
svector[13] = history[energy_history_offset + 1];
svector[14] = history[energy_history_offset + 2];
```

- `p13`: normal damping heat
- `p14`: tangential damping heat
- `p15`: tangential friction/deleted-contact heat


Deleted tangential history is added to friction heat in
`src/GRANULAR/pair_granular.cpp:220-231`:

```cpp
const double lost_energy =
    model->heat_tang_fric * model->tangential_model->elastic_potential();
heatflow[i] += 0.5 * lost_power;
heatflow[j] += 0.5 * lost_power;
```


Pair strain energy is added to `compute pe` in
`src/GRANULAR/pair_granular.cpp:317-320`:

```cpp
const double strain_energy = model->StrainEnergyNorm + model->StrainEnergyTang;
ev_tally_xyz(..., strain_energy, ...);
```




`compute dissipationrate/atom` now has three columns (`src/GRANULAR/compute_dissipationrate_atom.cpp:37-39`):

```cpp
size_peratom_cols = 3;
comm_reverse = 3;
```

It reads `svector[12]`, `svector[13]`, and `svector[14]`, then gives half to each
grain (`src/GRANULAR/compute_dissipationrate_atom.cpp:155-164`):

```cpp
dissipationrate[i][0] += 0.5 * force->pair->svector[12] / update->dt;
dissipationrate[i][1] += 0.5 * force->pair->svector[13] / update->dt;
dissipationrate[i][2] += 0.5 * force->pair->svector[14] / update->dt;
```




Wall contacts use the same heat scale factors as pair contacts (`src/GRANULAR/fix_wall_gran.cpp:123-129`). Deleted wall-contact history is assigned fully to the particle (there is no pair for `src/GRANULAR/fix_wall_gran.cpp:527-536`):

```cpp
const double lost_energy =
    model->heat_tang_fric * model->tangential_model->elastic_potential();
heatflow[i] += lost_energy / update->dt;
array_atom[i][10] = lost_energy;
array_atom[i][13] = lost_energy;
```

Active wall-contact heat is also assigned fully to the particle (`src/GRANULAR/fix_wall_gran.cpp:559-570`):

```cpp
wallstrain += model->StrainEnergyNorm + model->StrainEnergyTang;
heatflow[i] += model->dq_conduct + model->dq_dissipate / update->dt;
```

The per-contact wall output stores strain and heat in `src/GRANULAR/fix_wall_gran.cpp:573-587`:

```cpp
array_atom[i][8] = model->StrainEnergyNorm;
array_atom[i][9] = model->StrainEnergyTang;
array_atom[i][10] = model->dq_dissipate;
array_atom[i][11 + n] = model->svector[n];
```

Region walls use the same full wall-heat rule in `src/GRANULAR/fix_wall_gran_region.cpp:246-309`.




You can `bash examples/test_granular_energy/run_energy_tests.sh` and ` bash examples/test_granular_energy/run_shear_heat_tests.sh` and check the output for validation.