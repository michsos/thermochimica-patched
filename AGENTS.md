# thermochimica-patched agent guide

## Scope

Fork of [ORNL-CEES/thermochimica](https://github.com/ORNL-CEES/thermochimica)
carrying two patch branches. `master` tracks upstream. All local work lives
on separate branches so rebasing onto new upstream releases stays cheap.

See `PATCHES.md` for the branch strategy, commit list, and patch provenance.

## Branch overview

| Branch | Purpose |
|--------|---------|
| `master` | Tracks upstream ORNL-CEES/thermochimica cleanly |
| `bench-patches` | Fortran solver fixes from the bench investigation (needs audit) |
| `thermochimica-py-patches` | QKTOM support, dormant phase API, aventurine bug fixes |

## Building

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Or use the existing `Makefile_MacOS` for the macOS shared-library build
(produces `lib/libthermochimica.dylib` used by the Python benchmark harness
in `michsos/thermo-julia`).

## Patch discipline

- Do **not** commit directly to `master`. Keep it fast-forwardable to upstream.
- Rebase `bench-patches` and `thermochimica-py-patches` onto `master` after
  each upstream sync.
- When a patch is ready to upstream, open a PR to ORNL-CEES/thermochimica
  and drop it from the branch once merged.
- `bench-patches` commit `246b30e` is marked "needs verification" — audit
  and split before treating as stable.

## Knowledge base

Cross-project domain knowledge lives in
[thermo-kb](https://github.com/michsos/thermo-kb):

- FToxid database reference (phase names, system coverage) useful when
  tracing solver behaviour against specific oxide phases
- Scheil/viscosity design notes for the deferred kinetic simulation work

## Relationship to thermo-julia

`thermochimica-patched` is the Fortran reference solver. The Julia rewrite
(`michsos/thermo-julia`) runs alongside it for regression comparison. Neither
replaces the other until the Julia solver clears all four staircase stages.

## Commits

Commit when a change leaves the repo in a clearer state. Tag commits that
touch patched Fortran with the branch they belong to (`bench-patches` or
`thermochimica-py-patches`). Do not batch unrelated changes.
