# Patches on this fork

This is `michsos/thermochimica-patched`, a fork of
[ORNL-CEES/thermochimica](https://github.com/ORNL-CEES/thermochimica).

`master` tracks upstream. All local work lives on separate branches so rebasing
onto new upstream releases stays cheap.

## Branches

### `bench-patches`

Branched from upstream `c3f5958` (Merge PR #185). Contains Fortran solver
improvements discovered while building the FactSage/Thermochimica benchmark
harness in `michsos/thermochimica-bench` (now archived).

Commits (oldest first):

1. `6cef638` Fortran solver fixes: SIGSEGV guard in
   `RetryCalculationFirstPhase`, `sizeof`→`size` in `SaveReinitData`,
   best-driving-force seeding in `InitGemCheckSolnPhase`, debug symbols,
   `-O2 → -O3`.
2. `246b30e` Additional changes — new files `BestAssemblageSnapshot.f90`,
   `SeedSolutionPhases.f90`, `PostConvergenceCheck.f90`,
   `DumpDatabaseJSON.f90` and related integration edits. **Marked as
   "unknown patches need verification" — individual commits should be
   audited and potentially split before merging to master.**

### `thermochimica-py-patches`

Branched from upstream `c3f5958`. Contains solver and API additions needed
by the `thermochimica-py` Python wrapper and discovered during the aventurine
glaze investigation.

Commits (oldest first):

1. `4d5b731` Add API helpers and QKTOM support for thermochimica-py
2. `088e338` Fix pure condensed driving force indexing
3. `0eb3669` feat: add dormant phase diagnostics API
4. `61b5a95` extend Gibbs parser support for four extra terms
5. `98e4ae9` fix nine solver and API bugs discovered during aventurine
   investigation

## Strategy

- Keep `master` clean and fast-forward onto upstream.
- Rebase `bench-patches` and `thermochimica-py-patches` onto new `master`
  after each upstream sync. Both branches are expected to need occasional
  manual conflict resolution.
- When a patch stabilizes and is upstreamable, open an upstream PR and
  drop it from the branch once merged.
- Eventually collapse the two branches into one integrated `patched`
  branch once both streams are audited and no longer overlap.

## Provenance

These patches were reconstructed from the previous `michsos/thermochimica`
fork (deleted 2026-04-21) so that the new fork could be created fresh
against current upstream. A local bare-mirror backup of the old fork is
kept at `_legacy/thermochimica-backup.git` in the workspace that drove the
migration.
