# HyperLattice experimental presets (2026-09-27)

The simulator's HyperLattice preset list adds Octet Flight and Octet 4D Slice
after Cubic Flight and Hypercube Flight. **Pattern** selects Cubic or
Experimental / Octet Truss; **View** independently selects 3D perspective or
4D slice. Both patterns define both views. Octet remains outside standard
Teensy admission. Triangular, Cosine, and Gyroid have been removed from the
effect's preset list; their earlier measurements are retained below.

| Preset ID | Geometry | Work bound |
| --- | --- | --- |
| `experimental-octet-flight` | 3D D3/FCC edges | 64 plane candidates, 32 layers |
| `experimental-octet-4d-slice` | 4D D4 edges | 64 plane candidates, 32 layers |

Octet's vertices are integer lattice coordinates whose sum is even. Edges join
nearest neighbors, differing by one unit in two coordinates. In 3D this is the
D3/FCC lattice, with 12 neighbors per vertex and regular tetrahedral and
octahedral cells. In 4D the same rule gives D4, with 24 neighbors per vertex.
Cell Size is the nearest-neighbor strut length in world units in either view.

The 3D adapter visits four equally spaced plane families with tetrahedral
normals. Their pairwise acute angles are approximately 70.53 degrees. Its
coverage calculation uses the three intersecting line families in the crossed
plane. The 4D adapter visits eight diagonal hyperplane families and measures
distance to D4's edge graph. Both reuse monotone plane cursors and the shared
layer compositor. Coverage is approximate; neither adapter claims certified
surface intersections.

Both presets begin at radial offset zero with depth coloring and continuous
bounded world-space camera motion. Pause stops preset choreography; spatial
motion and palette cycling continue. A 4D view rotates the three-dimensional
ray domain inside the four-dimensional lattice. It does not project a 3D
octet image or extrude that image along another axis.

Unresolved, invalid, unsupported, or exhausted searches increment the read-only
Unfinished Rays count for the rendered frame. Analytic traversal preserves
already-consumed layers if a work cap is reached.
The previews do not claim error-free coverage or the Teensy's 62.5 ms deadline.
The [earlier device sweep](spherical_periodic_device_2026-09-26.md) already
exceeded that deadline at much smaller periodic-field budgets.

Changing Pattern or View adopts the selected tuple's complete default geometry,
preserving common color and near-fade controls. Automatic transitions switch
incompatible geometry at their midpoint. Schema 12 rejects older snapshots
without mutation. Lattice Planes and Softness are read-only for Octet;
Wire Radius and AA Strength remain editable. 4D Spin is read-only in 3D.
Cell Size changes framework spacing without moving the camera.
Experimental Sphere Radius, near
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
bash tools/profile_one.sh HyperLattice profile 70 16 '-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2'
```

Preset indices 2 and 3 select Octet in 3D and 4D respectively.
Fixed selection avoids mistaking a partial slow experimental cycle for a full
preset sweep. Device size/layout gates still apply to opt-in builds.

## Historical validation of the replaced previews

The following validation and images describe the earlier Triangular, Cosine,
and Gyroid presets, before their replacement by Octet. They do not measure
the current Octet implementation.

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

## Historical opt-in device validation

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
