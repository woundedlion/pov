# HyperLattice experimental presets (2026-09-27)

The simulator's HyperLattice preset list adds three previews after Cubic Flight
and Hypercube Flight. Configuration also selects them directly. All three are
marked Experimental and remain outside standard Teensy admission.

| Preset ID | Configuration | Work bound |
| --- | --- | --- |
| `experimental-triangular-flight` | Experimental / Triangular | 64 plane candidates, 32 layers |
| `experimental-cosine-surface` | Experimental / Cosine | 256 total field queries |
| `experimental-gyroid-surface` | Experimental / Gyroid | 512 total field queries |

The triangular-prism framework uses the existing nonorthogonal analytic event
adapter and approximate footprint coverage. Cosine and gyroid reuse the shared
certified first-boundary marcher and exact trigonometric nodal definitions.
They use isovalue zero, a 1.5-unit period, a six-unit far distance, and
period / 10,000 position tolerance. Refinement and ordinary probes both consume
the total query budget. Step and refinement caps equal that budget.

All three begin at radial offset zero with depth coloring, near fade 0.08,
flight phase speed 0.008 and 3D spin 0.0024 per frame. The framework's far
distance is 4.5, with wire radius 0.055 cells. The camera follows a bounded,
continuous orbit around (0.17, 0.31, 0.43) cells, without cell-wrap jumps.
Pause stops preset choreography; spatial motion and palette cycling continue.

Unresolved, invalid, unsupported, or exhausted searches increment the read-only
Unfinished Rays count for the rendered frame. A failed surface search contributes
no surface; it is not relabeled as empty geometry or accepted on proximity.
Analytic traversal preserves already-consumed layers if a work cap is reached.
The previews do not claim error-free coverage or the Teensy's 62.5 ms deadline.
The [earlier device sweep](spherical_periodic_device_2026-09-26.md) already
exceeded that deadline at much smaller periodic-field budgets.

Changing Configuration adopts its complete default geometry, preserving common
color and near-fade controls. Automatic transitions switch incompatible geometry
at their midpoint. Schema 11 rejects older snapshots without mutation. Lattice
Planes and Softness are read-only for every experiment; Wire Radius and AA
Strength are additionally read-only for periodic surfaces. 4D Spin is read-only
outside the 4D cubic configuration. Cell Size changes the framework spacing or
surface period; camera distances retain their world-unit convention after the
effect converts its cell-relative radial start.

## Build availability

`HS_ENABLE_HYPERLATTICE_EXPERIMENTS` defaults to one for Emscripten and native
test-oracle builds, and zero for firmware. Explicit zero or one overrides it.
With zero, experimental renderer calls, selectors, telemetry, and preset IDs
are absent. No new effect is registered.

An opt-in device capture uses the normal locked wrapper, for example:

```powershell
$env:HS_PROFILE_TREE = 'C:/work/Holosphere'
$env:HS_TEENSY_PORT = 'COM3'
bash tools/profile_one.sh HyperLattice profile 70 4 '-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=4'
```

Preset indices 2, 3, and 4 select triangular, cosine, and gyroid respectively.
Fixed selection avoids mistaking a partial slow experimental cycle for a full
preset sweep. Device size/layout gates still apply to opt-in builds.

## Validation

Native tests cover experimental selection, visible output and continued motion,
schema restore, invalid pattern values, inactive controls, and atomic geometry
transitions. The shared backend suites cover candidate ordering, verified hit
acceptance, and failure/query-budget accounting. Shipping firmware and the
simulator are built separately to exercise both feature-switch branches.
