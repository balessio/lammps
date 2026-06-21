# Granular Dissipation Accounting Plan

## Goal

Make granular dissipation and strain-energy accounting close the mechanical
energy budget for grain-grain and grain-wall contacts without requiring a
negative heat-production term.

The most important change is architectural: record dissipated increments during
the real force evaluation, when the previous contact history, updated contact
history, forces, and wall state are all simultaneously available. Avoid
reconstructing dissipated energy later from `pair->single()`.

## Findings From `DEM_Ben`

`DEM_Ben` differs from the current LAMMPS branch in ways that matter for energy
closure:

- It computes `disrateDash`, `disrateFric`, `peN`, and `peT` inside the contact
  resolution loop, at the same time contact history is advanced.
- Its Hertz-Mindlin model is closer to LAMMPS `mindlin_rescale/force` than to
  plain LAMMPS `mindlin`: it rescales the previous tangential elastic force on
  unloading and then derives an effective tangential stiffness.
- Its frictional work term is signed:
  `-slip . 0.5 * (f_te_prev + f_te) / dt`. That closes energy well, but it is
  not a monotone heat-production rate.
- Its Hertz normal strain energy uses `0.5 * k_n * delta_n^2`; LAMMPS should keep
  the physically integrated Hertz potential, `0.4 * F_n * delta`, for the
  granular Hertz normal law.

## Current LAMMPS Problem

The test scripts currently get pair dissipation from:

- `compute pair/local p13 p14 p15`
- `compute dissipationrate/atom`

Both routes depend on `PairGranular::single()`. That diagnostic call recomputes
the force model after the real force calculation, with `history_update = 0`.
For tangential-history models, this is not equivalent to the force-loop state:
the previous tangential history and updated tangential history are no longer both
available in the same way.

This is especially problematic for Mindlin-style models, where the tangential
elastic state depends on normal overlap/contact radius.

## Source-Code Plan

### 1. Add Persistent Per-Contact Energy Outputs For Pair Contacts

Add a small per-neighbor-history payload, or an adjacent fix-owned array, to
store the energy increments produced during `PairGranular::compute()`.

Required per-contact fields:

- normal damping dissipated increment
- tangential damping dissipated increment
- tangential friction/slip dissipated increment
- optional nonheat residual/work diagnostic
- normal strain energy
- tangential strain energy

These should be written exactly once per real contact evaluation, immediately
after `model->calculate_forces()`.

### 2. Stop Using `pair->single()` For Dissipation Accounting

Keep `pair->single()` useful for force/geometry diagnostics, but do not treat it
as the authoritative source for dissipated energy.

Replace or augment `compute dissipationrate/atom` with a compute that reads the
stored pair-contact energy increments from the real force loop. Likewise, prefer
a pair-local output path backed by the stored contact-loop data rather than a
recomputed `single()` call.

### 3. Mirror The Wall Contact Strategy For Pair Contacts

The wall path already has the right general shape: the wall fix stores contact
forces, strain energies, and dissipated increments after the real wall contact
calculation. Pair contacts should expose the same style of data:

- current contact flag
- contact force/contact vector
- normal strain energy
- tangential strain energy
- dissipated increment
- dissipated subparts

This will also make the test scripts simpler because pair and wall accounting
will have parallel semantics.

### 4. Implement Mindlin Heat As Energy-Balance Diagnostics, Not Recomputed Work

For `linear_history`, the nonnegative per-step tangential heat can be computed
from the force-loop work and strain-energy change:

```text
dq_tangential = max(0, -(F_t . dx_t + U_t,new - U_t,old))
```

For Mindlin variants, do not try to retrofit this inside `pair->single()`.
Compute it only in the real force loop, where these are available:

- previous history
- previous normal overlap/contact radius, if needed
- previous tangential elastic force
- updated history
- updated tangential elastic force
- actual tangential displacement increment used by the force update

For `mindlin` displacement history, explicitly decide how to report the energy
change caused by changing contact radius at fixed tangential displacement. That
term is not monotone heat. It should be a separate diagnostic residual unless a
physically justified irreversible interpretation is chosen.

For `mindlin_rescale/force`, use that model as the closest LAMMPS analogue to
`DEM_Ben`.

### 5. Treat Deleted Contacts In The Same Stored-Data System

When a history contact disappears, any remaining recoverable tangential strain
energy should be handled in the same force-loop/contact-cleanup path that owns
the history deletion.

Do not depend on `pair->single()` finding a just-deleted contact later. Instead:

- detect contact deletion during the neighbor/history update
- read the stored tangential elastic state before clearing it
- record the deleted-contact energy increment/residual in the stored output
- then clear the history

### 6. Preserve Per-Grain Attribution

For grain-grain contacts, split stored increments half to each grain unless a
future model has a reason to do otherwise. For grain-wall contacts, assign the
whole wall-contact increment to the grain.

The per-grain compute should report rates by dividing stored increments by
`update->dt`, but the stored canonical quantity should be an energy increment.

## Test Plan

Use `examples/test_granular_energy` as the main regression suite.

Required cases:

- `elastic_hooke`: no tangential force and no damping; total mechanical energy
  should remain nearly constant.
- `hooke_history`: checks linear tangential history accounting.
- `hertz_history`: checks Hertz normal energy plus linear tangential history.
- `hertz_mindlin`: documents the plain displacement-history Mindlin behavior.
- `hertz_mindlin_rescale_force`: closest LAMMPS analogue to `DEM_Ben`.

Required scenarios:

- two grains in periodic domain, started close enough to collide quickly
- 64 grains settling between bottom and top elastic walls

For each case, compare:

- mechanical energy
- pair normal/tangential strain energy
- wall normal/tangential strain energy
- accumulated pair dissipation
- accumulated wall dissipation
- total = mechanical + strain + accumulated dissipation

Acceptance targets should be model-specific:

- elastic cases: total relative drift near timestep integration error
- Hooke/Hertz with linear history: total drift decreases with timestep
- Mindlin variants: total drift decreases after contact-loop storage replaces
  `pair->single()` reconstruction; any nonmonotone residual is reported
  separately from heat

## Implementation Order

1. Add stored contact-loop pair energy increments without changing existing
   public computes.
2. Add a temporary debug compute or local output to compare stored force-loop
   increments against current `pair->single()`-based increments.
3. Switch `compute dissipationrate/atom` to the stored increments.
4. Add a pair-local output path for stored increments and update the tests.
5. Rework Mindlin tangential heat using force-loop previous/current state.
6. Add deleted-contact accounting to the pair history cleanup path.
7. Remove temporary debug output once the tests close.

## Risks

- Pair neighbor/history storage is more delicate than wall storage because of
  Newton on/off, ghost atoms, and neighbor-list rebuilds.
- `compute pair/local` currently expects data from `pair->single()`. Replacing
  that path may require a new compute style or new `pair/local` keywords rather
  than overloading existing `p13`-style outputs.
- Plain `mindlin` has an inherent overlap-coupled tangential energy change. It
  may need a separate residual diagnostic, even after the heat accounting is
  fixed.

## Implemented In This Branch

- Pair contacts now store per-contact energy increments during the real
  `PairGranular::compute()` force loop and expose those stored values through
  `pair/local p13 p14 p15 p16`.
- `compute dissipationrate/atom` now reports four columns: normal damping heat,
  tangential damping heat, tangential friction/deleted-contact heat, and signed
  residual work.
- Wall contact arrays now include the same residual-work column through
  `f_wall[15]`, with wall normal strain, tangential strain, and heat included in
  the settling tests.
- Linear-history tangential heat is nonnegative. The nonmonotone part is
  reported separately as residual work and is not added to `dq_dissipate` or
  heatflow.
- The `examples/test_granular_energy` harness now runs the two-grain and
  settling-wall tests, writes component time series, optionally writes dumps,
  and plots `total = mechanical + accumulated heat + residual work`.

## Current Test Outcome

The contact-loop storage fixed the pair-local versus per-atom mismatch, and
the two-grain tests are now useful regressions:

- elastic Hooke: near machine-level drift for the sampled run
- Mindlin variants: about `3e-4` relative drift in the two-grain test
- Hertz linear-history: about `5e-3` relative drift in the two-grain test
- Hooke linear-history: about `1e-2` relative drift in the two-grain test

The 64-grain settling-wall test still exposes unresolved model/accounting
issues until the tangential-history clipping rule was changed to match the
`DEM_Ben` force split more closely. The important source difference was not the
frame projection itself. Both codes project tangential history into the current
tangent plane while preserving magnitude for grain-grain contacts. The large
drift came from Coulomb clipping: `DEM_Ben` scales both the tangential spring
state and the tangential dashpot force, whereas the original LAMMPS path folded
an unscaled dashpot contribution back into the history variable.

After switching `linear_history` and Mindlin clipping to scale the stored
tangential state and the tangential damping work consistently, the full
2,000,000-step 64-grain settling-wall sweep produced these relative drifts:

- elastic Hooke: `3.1e-6`
- elastic tangential history: `-5.0e-4`
- Hooke `linear_history`: `7.1e-3`
- Hertz `linear_history`: `3.0e-3`
- plain Hertz-Mindlin: `-4.3e-3`
- Hertz-Mindlin force: `-7.7e-3`
- Hertz-Mindlin rescale force: `-9.5e-3`

The two-grain tests remain tighter, with the Mindlin variants below `5e-4` and
the linear-history cases at roughly `5e-3` to `1e-2`. These residuals now look
like ordinary discrete-integration/model-splitting error rather than a missing
energy channel.
