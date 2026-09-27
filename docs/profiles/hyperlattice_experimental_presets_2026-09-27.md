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
continuous world-space orbit around (0.255, 0.465, 0.645), without cell-wrap jumps.
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
surface period without moving the camera. Experimental Sphere Radius, near
fade, and far distance use world units; the legacy cubic radial start remains
cell-relative for compatibility.

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

The full native suite passes (100 tests, one intentional replay skip), followed
by six passing focused tests after the camera/world-ray changes. The release
WASM smoke selects and renders every HyperLattice preset at every supported
resolution. [Native full run](evidence/hyperlattice_experiments_2026-09-27/native-full.txt),
[final focused run](evidence/hyperlattice_experiments_2026-09-27/native-final-focused.txt),
and [WASM smoke](evidence/hyperlattice_experiments_2026-09-27/wasm-smoke.txt)
preserve the results.

The standard Phantasm build passes its size/layout gates with 195,544 bytes of
RAM1 code (unchanged), 314,784 bytes of RAM1 variables (unchanged), and 725,452
bytes of FLASH data (+4). The zero-reserve ITCM ceiling remains 196,608 bytes,
leaving 1,064 bytes. [Build evidence](evidence/hyperlattice_experiments_2026-09-27/shipping-build.txt)
uses the integrated preset source before the subsequent experimental-only
camera and query improvements.

The following previews use the final renderer source `c012be4a8`, at 288 by
144 after 24 frames from each preset's fresh initialization. From top to bottom:
triangular framework, cosine, gyroid. RGB16 linear output is converted to sRGB
for this flat equirectangular preview; it is not a photograph of the device.

![Triangular, cosine and gyroid previews](evidence/hyperlattice_experiments_2026-09-27/preview.png)

The last rendered frames report zero, 104, and 69 unfinished rays respectively;
the maximum observed WASM stack watermark is 944 bytes. The
[preview summary](evidence/hyperlattice_experiments_2026-09-27/wasm-preview-summary.json)
also records host timings, including the first draw. Those host timings are not
Teensy measurements or admission evidence.

## Opt-in device validation

Final cold builds at clean source `596cdc1fd` pass the default and opt-in
Phantasm size/layout gates with zero first-party warnings. Their runtime source
is identical to landed `920a9f622`. The
[layout/symbol audit](evidence/hyperlattice_experiments_2026-09-27/device-layout.txt),
[default build](evidence/hyperlattice_experiments_2026-09-27/default-final-build.txt),
[opt-in build](evidence/hyperlattice_experiments_2026-09-27/optin-final-build.txt),
[warning audit](evidence/hyperlattice_experiments_2026-09-27/warnings.txt), and
[7,550 passing native assertions](evidence/hyperlattice_experiments_2026-09-27/native-placement.txt)
retain the final evidence.

| Full Phantasm | Standard | Experiments enabled | Change |
| --- | ---: | ---: | ---: |
| RAM1 code | 195,544 | 196,584 | +1,040 |
| RAM1 variables | 314,784 | 314,784 | 0 |
| FLASH code | 502,112 | 515,120 | +13,008 |
| FLASH data | 725,452 | 725,728 | +276 |
| ITCM headroom | 1,064 | 24 | -1,040 |

Both retain 12,896 bytes of stack space and 4,224 free bytes in RAM2. The
standard ELF contains no experimental renderer symbols. Four measured
`HS_HOT_FLASH_MEMBER` placements move the surface search, surface shading,
outlined evaluation lambda, and framework validity check to cached flash.
This reduces opt-in ITCM from 199,560 to 196,584 bytes without changing the
zero-reserve ceiling. Small emit and step helpers retain the compiler's inline
choices. These annotations are no-ops for Clang; the installed WASM binary
matches the tested preview binary byte for byte.

A 30-second locked capture on COM4 selects gyroid preset 4 with four-frame
counter windows, shipping optimization, and live segmented-driver interrupts.
It validates execution, not realtime admission: runtime frames 2 through 38
average **731.814 ms**, with **629.698 / 820.923 ms** minimum/maximum and
**37/37 deadline spills**. Setup frame 1 takes 1,439.063 ms and is excluded.
Nine complete counter windows cover frames 1 through 36; the runtime figures
also retain the two trailing individual frame records.

[Raw capture](evidence/hyperlattice_experiments_2026-09-27/gyro.txt),
[provenance](evidence/hyperlattice_experiments_2026-09-27/gyro.provenance),
[profile build](evidence/hyperlattice_experiments_2026-09-27/gyro-build.txt),
[exact flags](evidence/hyperlattice_experiments_2026-09-27/gyro-envdump.json),
[summary](evidence/hyperlattice_experiments_2026-09-27/gyro-summary.txt), and
[parser validation](evidence/hyperlattice_experiments_2026-09-27/gyro-validate.txt)
identify the image and frame ranges. This timing capture does not log per-ray
failure counts; the native/WASM diagnostic results above are separate evidence.
The wrapper's paired Phantasm attestation is the default image; the full opt-in
image is separately identified by the layout audit and build log above.
