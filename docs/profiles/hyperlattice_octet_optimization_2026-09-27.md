# HyperLattice Octet optimization ledger

Measured on Teensy 4.0 at 600 MHz, COM4, with the real segmented driver and live interrupts. All changes retain the authored geometry, camera motion, ray budget and coverage model. The experimental presets remain opt-in.

[Shipping profile](shipping/profile_hyperlattice_octet_teensy_2026-09-27.md) and [global-O3 reference](O3/profile_hyperlattice_octet_teensy_2026-09-27.md) contain scope trees, ISR costs, cadence, memory and capture validation.

## Matched comparisons

Startup frame 1 is excluded. Comparisons use the same animation frame numbers; the faster 4D captures otherwise cover more camera motion during their fixed recording duration.

| Preset | Frames | Baseline mean / peak ms | Final mean / peak ms | Speedup | Render reduction |
|---|---:|---:|---:|---:|---:|
| 3D | 2–549 | 89.562 / 94.478 | 81.398 / 87.013 | 1.100× | 9.12% |
| 4D | 2–77 | 524.586 / 549.729 | 178.304 / 183.605 | 2.942× | 66.01% |

Independent rebuilt 3D repeat, frames 2–549: baseline 89.564/94.469 ms mean/peak; final 81.401/87.032 ms, 1.100× speedup. Both runs used the same board and capture parameters.

Both final presets still exceed the 62.5 ms display interval. The shipping captures deliver about 8 fps in 3D and 5.3 fps in 4D; the original 4D capture delivered about 1.8 fps. Render-time reduction does not translate continuously into fps because display flips quantize cadence.

## Experiment ledger

All rows are shipping selective-O3 images. Means and peaks below include every recorded live frame after startup, so use the matched table above for the headline 4D comparison. Every trial has a raw capture, provenance, parser validation and numeric summary in the evidence directory.

| Trial | Result | Frames | Mean ms | Peak ms | Raw evidence |
|---|---|---:|---:|---:|---|
| 3D baseline | baseline | 2–549 | 89.562 | 94.478 | [octet_base_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_base_ship.txt) |
| 3D grouped owner lists | rejected: slower | 2–549 | 96.693 | 101.307 | [octet_grouped_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_grouped_ship.txt) |
| 3D squared support test | accepted | 2–549 | 87.278 | 93.000 | [octet_squared_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_squared_ship.txt) |
| 3D unrolled incidence switch | rejected: slower and larger | 2–549 | 98.494 | 107.519 | [octet_incidence_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_incidence_ship.txt) |
| 3D normalized plane residuals | accepted | 2–549 | 86.554 | 92.170 | [octet_normalized_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_normalized_ship.txt) |
| 3D reciprocal cursor setup | rejected: slower | 2–549 | 90.237 | 95.650 | [octet_reciprocal_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_reciprocal_ship.txt) |
| 3D prepared camera projection | accepted | 2–549 | 85.942 | 91.856 | [octet_camera_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_camera_ship.txt) |
| 4D baseline | baseline | 2–77 | 524.586 | 549.729 | [octet4_base_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet4_base_ship.txt) |
| 4D parity / dominant-coordinate selection | accepted | 2–185 | 189.733 | 198.877 | [octet4_parity_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet4_parity_ship.txt) |
| 4D scalar distance specialization | accepted | 2–232 | 180.867 | 192.335 | [octet4_scalar_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet4_scalar_ship.txt) |
| 4D cached normalized ray | accepted | 2–232 | 177.299 | 188.377 | [octet4_cached_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet4_cached_ship.txt) |
| 3D combined final | retained | 2–550 | 81.407 | 87.013 | [octet_final_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_final_ship.txt) |
| 4D combined final | retained | 2–232 | 177.308 | 188.387 | [octet4_final_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet4_final_ship.txt) |
| 3D rebuilt baseline | repeat | 2–549 | 89.564 | 94.469 | [octet_repeat_base_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_repeat_base_ship.txt) |
| 3D rebuilt final | repeat | 2–550 | 81.410 | 87.032 | [octet_repeat_final_ship](evidence/hyperlattice_octet_optimization_2026-09-27/octet_repeat_final_ship.txt) |

These are sequential whole-image experiments, not isolated instruction microbenchmarks. In particular, the final 3D improvement over the camera-only capture followed the subsequent 4D changes and a code-layout shift. Its repeat validates the combined image; it does not establish a separate arithmetic win in the unchanged 3D loop.

## Geometry and arithmetic

### D4: choose the strut without searching twelve directions

Normalize the query into lattice coordinates and write q = n + r, where n is componentwise round(q) and each residual lies in [-1/2, 1/2]. Let R² be the squared residual length and let a ≥ b be its two largest absolute components. The parity of the integer-coordinate sum determines the nearest D4 edge class.

For even parity, the best squared distance for a coordinate pair is R² − (a+b)²/2. For odd parity it is R² + 1/2 − a − b − (a−b)²/2. Both decrease as either selected magnitude increases on this residual domain. Thus the two largest magnitudes select the minimizing pair; candidates obtained by changing an unselected integer coordinate cannot improve the minimum. Residual signs and parity then determine the strut orientation.

The implementation uses stable comparisons and explicitly resolves the multiple-half-coordinate boundary case. This preserves the original feature identity and normal at ties, not only the minimum distance. The ordinary query drops the twelve-direction search and reduces round operations from 28 to 5.

Distance-only callers evaluate the two transverse residual squares plus half the squared along-pair residual. They avoid constructing a four-component offset and then taking its norm. Offset and normal callers retain their vector specialization. The 4D event stream caches the normalized ray origin and direction once, avoiding repeated world-to-lattice conversions at each crossing.

### 3D: exploit the tetrahedral Gram matrix

The four unit plane normals have pairwise dot product −1/3. At a crossing of plane i, the squared distance to the ray/strut closest approach for pair (i,j) is the squared residual of plane j times

```text
weight = a_i² / (a_i² + a_j² + (2/3) a_i a_j)
a_k = dot(plane_normal_k, ray_direction)
```

Store this squared weight rather than its square root. For each candidate, compare the minimum squared distance against (wire_radius + footprint_width/2)². Crossings beyond that support contribute zero, so they require no square root. A contributing crossing takes one square root before the existing coverage calculation. Plane ownership, event ordering and generic cursor initialization remain unchanged.

Normalize plane positions and speeds once per ray. Candidate evaluation becomes U − round(U), with spacing² absorbed into the stored weight. The camera embedding and tetrahedral projections are composed once per frame, leaving four direction dot products per ray instead of rebuilding world-space projections at every pixel.

The prepared path relies on the effect’s existing validated camera/settings contract; the public general constructor retains its validation. These transformations introduce small floating-point reassociation differences, tested against independent geometry and rendering oracles.

## ARM code generation

Toolchain: pinned Teensy GCC 15.2.1, Cortex-M7 hard-float. [Symbol sizes](evidence/hyperlattice_octet_optimization_2026-09-27/codegen_symbols.json), [baseline disassembly](evidence/hyperlattice_octet_optimization_2026-09-27/octet_base_ship_codegen.txt), and [final disassembly](evidence/hyperlattice_octet_optimization_2026-09-27/octet_final_ship_codegen.txt) are preserved.

| Symbol / hot path | Baseline bytes | Final bytes |
|---|---:|---:|
| 3D shade | 3456 | 3236 |
| 3D ray constructor (general → prepared) | 1376 | 1746 |
| 4D nearest edge (vector → normalized scalar hot specialization) | 2108 | 568 |
| 4D shade | 4068 | 4122 |

The edge symbols have different return contracts; the size comparison describes the actual hot path, not equivalent standalone APIs. All these symbols execute from cached flash. roundf already compiles to VRINTA.F32 and square root to VSQRT.F32: the wins come from eliminating operations and candidate work, not replacing library calls. Color clamps already lower to six VMINNM/VMAXNM instructions per visible layer.

The unrolled incidence switch grew the 3D shader to 4366 bytes and regressed runtime. Grouped loops reduced code size but still slowed the image. Replacing cursor divisions with reciprocal setup also regressed. These were discarded. Residual recurrences were not adopted because accumulated drift would complicate crossing/tie equivalence for a hardware round instruction that is already cheap.

## Memory

Instrumented shipping single-effect images, baseline → final:

| Region | Baseline bytes | Final bytes | Delta |
|---|---:|---:|---:|
| RAM1 code | 20376 | 20376 | +0 |
| RAM1 variables | 315008 | 315008 | +0 |
| FLASH data | 148708 | 148708 | +0 |
| FLASH code | 76136 | 75144 | -992 |

The default full-roster Phantasm image passes the size/layout gates with RAM1 code 195544 bytes. Only 1064 bytes remain before the next 32 KiB ITCM allocation boundary, so moving the whole shader into ITCM is unsuitable. A separate full-roster build with HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 also passes: RAM1 code 196200, RAM1 variables 314784, FLASH data 725644, stack space 12896, and RAM2 free 4224 bytes. These are full-roster figures, distinct from the instrumented single-effect table.

## Validation and reproducibility

- Native CTest: 100 passed, zero failures, one unrelated replay fixture skipped because it requires Clang 22 and the host has Clang 23.
- Strict and optimized (-Os, -ffast-math, -fno-finite-math-only) builds each pass all four targeted geometry/rendering modules. The SDF module has 769002 passing assertions. The exact-grid fixture runs in normalized coordinates to avoid an irrational scale round-trip; zero offset tolerance and exact feature checks remain enabled.
- Geometry tests compare the D4 selector with the original exhaustive search across an 83521-point exact grid and 10000 random points, plus boundary/scaled/normal and cached-event checks.
- 3D coverage uses an independent ray/line cross-product oracle across 600 rays and multiple scales; crossing times are checked against the generic plane streams.
- Prepared-camera rendering is compared over 3456 rays, including rotations, scale, radial offset, clip interval and depth changes; 144 invalid direction cases retain validation. Alpha tolerance is 3e-4 and RGB tolerance is two code values.
- The WASM release build passes. The simulator snapshot is installed from the landed source.
- Firmware builds and size/layout gates pass for default and experimental full-roster images; every capture passes the profile parser validator.

Baseline source is a3b4c5265 (continuous experimental camera motion). The final measured implementation is e1b2bb9d1, rebased onto the peer max-speed default change 251f5a951. The authored preset used here has its own speed; frame-matched camera comparisons remain applicable. A later test-only fixture correction and this documentation do not change the measured firmware.

Each raw filename has matching .provenance, _validate.txt and _summary.json files. The evidence directory also preserves source patches from the baseline, dirty-tree snapshot patches, SHA-256-wrapped build logs/environment dumps, and selected disassembly. Original ELF/map artifacts and complete working logs are retained locally under C:/work/Holosphere/build/prof/octet_opt_20260927. This makes discarded and pre-rebase trials inspectable without depending on their temporary branches.

Reproduce a shipping 3D capture from the repository root (use preset 3 and 45 seconds for 4D; profile_o3 selects the global-O3 reference):

```sh
export HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM4
bash tools/profile_one.sh HyperLattice profile 70 16 \
  "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2"
```

No whole-preset cycle, later camera trajectory, or unchanged output bit pattern is claimed. The report measures held presets with active camera/palette motion and bounded floating-point differences.
