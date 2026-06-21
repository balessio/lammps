# Granular Energy Accounting Summary

## Code Changes

The implementation keeps heat production nonnegative and stores energy
increments during the real contact force evaluation instead of reconstructing
them later from a diagnostic `pair->single()` call.

- `PairGranular` stores per-contact energy increments in neighbor history:
  normal damping heat, tangential damping heat, tangential friction/deleted
  contact heat, and signed residual work.
- `compute dissipationrate/atom` reports the same four components as per-grain
  rates, split half to each grain for grain-grain contacts.
- `fix wall/gran` and `fix wall/gran/region` expose matching wall-contact
  energy components, including wall normal/tangential strain energy.
- Pair strain is tallied into `compute pe` as `U_n + U_t`; wall strain is exposed
  through LAMMPS fix energy and included when the input uses
  `fix_modify <wallfix> energy yes`.
- `linear_history` and Mindlin tangential models now clip Coulomb-limited
  contacts by scaling both the stored tangential elastic state and tangential
  damping work, matching the force split used by `DEM_Ben`.
- `examples/test_granular_energy` contains the regression harness, plotting,
  and a documented plan for future cleanup.

## Edited Code Locations

- `src/GRANULAR/pair_granular.cpp` and `.h`: added the contact-history offset
  for four energy channels, store them immediately after
  `model->calculate_forces()`, carry them through history transfer, and expose
  them from `single()` as `p13` through `p16`.
- `src/GRANULAR/compute_dissipationrate_atom.cpp`: expanded the per-atom
  compute to four columns and tallies the stored pair heat/residual channels
  onto both grains.
- `src/GRANULAR/fix_wall_gran.cpp` and
  `src/GRANULAR/fix_wall_gran_region.cpp`: added wall normal/tangential strain
  energy, heat, heat subcomponents, and residual work to the wall contact
  output; deleted wall contacts release their stored tangential strain energy
  into heat. `fix_wall_gran` also reports wall strain through `compute_scalar()`
  for `fix_modify energy yes`.
- `src/GRANULAR/pair_granular.cpp`: tallies `StrainEnergyNorm + StrainEnergyTang`
  as the pair potential-energy contribution, leaving the Coulomb energy slot at
  zero.
- `src/GRANULAR/gran_sub_mod_tangential.cpp` and `.h`: split tangential heat
  into damping, friction/deleted-history heat, and residual work; `linear_history`
  and Mindlin Coulomb clipping now scale the stored tangential state and
  tangential damping work by the same force scale.
- `src/GRANULAR/gran_sub_mod_normal.cpp`: corrected Hertz normal strain energy
  to the integrated potential `U_n = 2/5 F_n delta_n`.
- `src/GRANULAR/granular_model.cpp`: writes and reads the `dissipative_heat`
  flag plus its normal damping, tangential damping, and tangential friction heat
  scale factors in binary restarts, so prepared granular packings can restart
  with contact history and heat accounting intact.
- `examples/test_granular_energy/*`: added standalone LAMMPS inputs, runners,
  plotting, and this summary for reproducible energy-accounting checks.

## Physical Bookkeeping

For an active contact, the plotted total is

```text
E_total = K_trans + K_rot + U_g + U_n + U_t + Q_heat + W_res
```

where `Q_heat` is monotone and `W_res` is a signed diagnostic, not heat.

Normal strain energy is

```text
U_n = 1/2 k_n delta_n^2       Hooke
U_n = 2/5 F_n delta_n         Hertz
```

Tangential strain energy is

```text
U_t = 1/2 k_t |delta_t|^2
```

The tangential work balance is decomposed as

```text
dW_t = -F_t . ddelta_t
dW_t = dU_t + dQ_t + dW_res
```

with

```text
dQ_t = dQ_t,damp + dQ_t,fric >= 0
```

and the nonmonotone remainder stored separately as `dW_res`.

## Coulomb Clipping Fix

The important discrepancy with `DEM_Ben` was not the frame projection. The
elastic no-damping history tests close well, so preserving the tangential
history magnitude after projection is acceptable.

The problem was Coulomb clipping. The old LAMMPS history path effectively
absorbed an unscaled tangential dashpot contribution into the stored history.
The updated rule uses the contact force scale factor

```text
s = min(1, mu F_n / |F_t,elastic + F_t,damp|)
```

and applies it consistently:

```text
delta_t <- s delta_t
F_t <- s F_t
dQ_t,damp <- s c_t |v_t|^2 dt
```

This makes the heat term monotone while preventing the clipped dashpot force
from being double-counted as recoverable tangential strain.

## Current Demonstration

The full post-fix sweep `runs_final_demclip` gave settling-wall relative drift
below about `1%` for the main dissipative cases:

```text
Hooke linear_history          7.1e-3
Hertz linear_history          3.0e-3
Hertz-Mindlin                -4.3e-3
Hertz-Mindlin force          -7.7e-3
Hertz-Mindlin rescale force  -9.5e-3
```

The remaining drift is now consistent with discrete integration/model-splitting
error rather than a missing heat channel.
