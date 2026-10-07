# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot. Raw capture: `build/prof/hyperlattice_ship.log`, captured 2026-10-06 20:43 America/Los_Angeles on COM4.
Replaces the earlier October 6 snapshot at `e2f5b0a3d` and the September 29 canonical report. The October 1 architecture supplements retain their historical three-preset measurements.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel and DMA ISRs live, COM4 |
| Image | `profile`: `-Os` base, selective `HS_O3`; cached shader `HS_HOT_FLASH_MEMBER`, no ITCM relocation |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, clean committed worktree `a13c149d4788836e6f0071db16a78243d28a4bc2` |
| Method | 345 s, 16-frame windows, all nine shipping presets, `-D HS_PROFILE_EPOCH_REVS=2800`, deep counters disabled |
| Compiler | Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 |
| Reproduce | `HS_PROFILE_TREE=/c/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 345 16 "-D HS_PROFILE_EPOCH_REVS=2800"` |

Profile image: `FLASH: code:130140, data:153344, headers:8352` / `RAM1: variables:315008, code:15352, padding:17416, free:176512` / `RAM2: variables:520064, free:4224`.

The full-roster, non-profile Phantasm build passes the region budget and layout gate:

| Region | Baseline `e2f5b0a3d` bytes | Optimized bytes | Delta bytes |
|---|--:|--:|--:|
| FLASH code | 637628 | 636636 | -992 |
| FLASH data | 846452 | 846452 | 0 |
| RAM1 code | 177656 | 177656 | 0 |
| RAM1 variables | 314784 | 314784 | 0 |
| RAM1 free | 12896 | 12896 | 0 |
| RAM2 variables | 520064 | 520064 | 0 |

Capture exactness and all-nine-preset wrap are checked by `tools/parse_profile.py ... validate`: VALID, root cyc/600 versus wall sum within **0.1 ppm**, frames 3089–3104. The verbatim output is retained at `build/prof/hyperlattice_ship.validate.txt`.

## Frame cadence

**Full-cycle peak render: 51.774 ms, frame 3209, down from 60.525 ms (14.5% lower).** This clears the requested 54 ms by 2.226 ms and the 62.5 ms display window by 10.726 ms. **0/5487 spilled frames**, across all nine presets and their transitions.

`hl_shader_draw` averages 26.35 ms/frame; its worst window is 41.97 ms/frame, frames 3089–3104. Baseline values were 27.52 and 52.54 ms. Setup frame 1 (42.811 ms) is reported separately and excluded from cadence buckets. Even including setup, this capture stays below 54 ms.

The shader retains 146×73 = 10,658 rays/frame, the same presets, support, trace budgets and appearance. `canvas_buffer_wait` waits for the next display flip. Render timing excludes this idle time; wall timing includes it. Cadence remains 16 fps.

## Phase-by-phase readout

Nine presets cycle through Cubic, Octet Truss and Shell families. Each has a 320-frame hold followed by a 240-frame segue. Within a pattern and view, parameters lerp; family changes fade through black. Markers visit all nine and wrap to the first before this capture ends.

### Peak window (frames 3201–3216)

```
frame                   62.06 ms  37.24 Mcyc  100.0%
  hl_shader_draw          41.26 ms  24.76 Mcyc   66.5%
  pov_preserve_half        0.14 ms   0.08 Mcyc    0.2%
  hl_timeline_step         0.03 ms   0.02 Mcyc    0.0%
  canvas_clear             0.09 ms   0.05 Mcyc    0.1%
  canvas_buffer_wait      18.73 ms  11.24 Mcyc   30.2%
```

Wall min/avg/max = 51.65/62.06/74.63 ms. Tree figures are per-frame window averages, distinct from the peak frame. The limiting preset remains Octet 4D Flight.

### Per-preset table

Cadence buckets use per-frame ownership, including transitions. Setup frame 1 is excluded. Clean shader cost is the costliest modal-call-count window that does not straddle a preset advance. `#` is zero-based; serial markers are one-based. The comparison is the same nine-preset COM4 capture at `e2f5b0a3d`, earlier October 6.

| # | Preset | Peak render ms | Spilled/frames | Clean shader ms/f | Clean windows | Baseline peak ms |
|---|---|--:|--:|--:|--:|--:|
| 5 | Octet 4D Flight | 🟢 51.77 | 0/559 | 41.97 | 34/35 | 60.52 |
| 8 | Shell 4D Flight | 🟢 37.19 | 0/559 | 32.07 | 34/35 | 37.22 |
| 2 | Hypercube Flight | 🟢 36.83 | 0/559 | 28.33 | 34/35 | 36.95 |
| 7 | Shell Close Flight | 🟢 31.95 | 0/679 | 27.59 | 42/42 | 31.76 |
| 3 | Octet Flight | 🟢 31.90 | 0/439 | 27.20 | 26/27 | 31.84 |
| 4 | Octet Wide Flight | 🟢 31.25 | 0/679 | 28.53 | 42/43 | 31.23 |
| 1 | Cubic Wide Flight | 🟢 29.50 | 0/817 | 25.11 | 50/51 | 29.40 |
| 6 | Shell Flight | 🟢 28.44 | 0/439 | 23.81 | 26/27 | 28.49 |
| 0 | Cubic Flight | 🟢 26.22 | 0/757 | 23.81 | 46/48 | 26.29 |

### Per-pixel figures

The mean shader cost is approximately 1483 cycles/ray. Pixels are written premultiplied; there is no `filter_blend` scope.

## Deep profile and code generation

The regression is concentrated in Octet 4D Flight: the earlier nine-preset peak rose from 55.87 ms on September 29 to 60.525 ms on October 6. The historical captures differ in source and board; this investigation does not identify a specific regressing commit. Optimizations below were measured against the current October 6 source on COM4.

A separate baseline diagnostic capture on COM3 instruments whole 4D rays, each family walk, plane crossings, exact evaluation, insertion and compositing. In window 417–432 it observes **4 family walks, 16.039 plane crossings, 2.671 exact evaluations and 1.108 covered inserts per ray**. Only **16.65%** of candidates reach exact evaluation. Walk setup and rejected candidates therefore consume substantial time before useful shading.

The instrumented window reports 24.67 ms/frame of walk setup (walk scopes minus crossing scopes), 32.69 ms/frame in crossings, 12.46 ms/frame in nested exact evaluations, 0.59 ms/frame in insertion, and 3.24 ms/frame in compositing. These scopes overlap; do not add them. Deep counter overhead and layout changes make these diagnostic times unsuitable for comparison to shipping peaks. Duplicate family-walk scope names are summed explicitly.

The final canonical path:

- Delays class denominator and transverse calculations until a conservative bound accepts a class, removing persistent per-class arrays.
- Rejects a crossing using the minimum squared free-coordinate residual before computing every class's across-coordinate bound.
- Uses positive canonical owner speeds and known non-parallel ownership for all but the largest-coordinate sum.
- Reuses inverse plane speeds in the threshold calculation, eliminating four divide instructions.
- Truncates fixed-point starts and advances with a wider conservative rounding margin; pools the one-crossing advance guard.

Exact residuals, parity handling, closest-line selection and compositing retain their existing arithmetic. The independent generic event tracer checks the optimized result, including axis, near-axis, equal-magnitude and parallel boundary directions.

`shade_octet<true>` disassembly agrees between the effect-only profile and full shipping Phantasm ELF:

| Static codegen metric | Baseline | Final |
|---|--:|--:|
| Function bytes | 11824 | 11072 |
| Instructions | 3350 | 3124 |
| Instructions addressing `[sp]` | 356 | 325 |
| Local stack reservation, excluding register saves | 1172 B | 1084 B |
| `vdiv.f32` instructions | 20 | 16 |
| `vmrs` instructions | 113 | 109 |

These are static counts over the entire function, including cold fallback paths, not executed instruction counts. The shipping symbol remains in cached flash (`0x60056828`); baseline was `0x60056918`.

### Controlled trials

Each trial uses COM4, preset index 5, 70 s and 16-frame windows; startup is excluded. Mean render here covers all post-startup raw frame records, including the final partial window. The table is for optimization analysis, not the ranked archive.

| Trial | Peak ms | Mean render ms | Outcome |
|---|--:|--:|---|
| Current-source baseline | 60.900 | 52.301 | Reference |
| Separate non-inlined family walks | 79.523 | 70.705 | Rejected |
| Separate walks with compile-time pair tables | 62.921 | 53.683 | Rejected |
| Lazy transverse data and free-coordinate rejection | 59.611 | 49.215 | Kept |
| Lazy denominators and pooled advance guard | 55.908 | 45.977 | Kept |
| O2 shader specialization | 92.855 | 80.390 | Rejected |
| Reciprocal threshold and conservative truncation | 55.906 | 45.225 | Kept |
| Canonical class-ownership shortcut | 54.383 | 43.202 | Kept |
| Positive-speed setup and minimum free-pair bound | 51.777 | 40.957 | Final |

The final fixed-preset trial reduces mean render by 21.7%. Lower static instruction count alone did not predict performance: the O2 trial was much slower, and splitting the walk lost the benefit of inlining and specialization.

## Column-ISR / DMA marshaling cost

Peak-window ISR measurements; times are min/avg/max:

```
isr_wake          1145.9/frame  0.6/1.6/11.2 us  cpu 3.04%
isr_pack           143.2/frame  6.2/6.7/9.4 us  cpu 1.53%
isr_dma_submit     143.2/frame  0.8/0.9/1.0 us  cpu 0.21%
```

CYCCNT free-runs, so shader scopes absorb ISR time. DMA wire transfer is asynchronous; CPU submission cost is measured separately.

## Summary ranking

1. `hl_shader_draw`: 41.26 ms/frame in the peak window; the optimized 4D traversal remains dominant.
2. `pov_preserve_half`: 0.14 ms/frame.
3. `canvas_clear`: 0.09 ms/frame.
4. `hl_timeline_step`: 0.03 ms/frame.

README cells: peak 🟢 51.77 (9), spilled 🟢 0/5487 (0.00%).

Historical O3 columns were not re-measured and retain their own source pairing.

## Validation and caveats

- Native Debug: 105 passed, one existing skipped replay-form test, zero failures; `HS_SMOKE_FRAMES=120`, 106 CTest nodes.
- Optimized assertions-off IEEE and fast-math HyperLattice harnesses: 2/2 passed, including the 16,000-ray canonical oracle and new boundary-direction test.
- WASM Release builds; smoke validation passes all 84 effect/resolution combinations at 120 frames per effect, with no all-black result or dynamic code generation.
- Mutation check: replacing `free_nearest > THRESHOLD` with `free_nearest > 0u` causes the new boundary test to fail 60 assertions. Production source was restored and both optimized harnesses passed again.
- Profile and non-profile Phantasm builds pass; shipping memory/layout gate passes. Commit hooks pass formatting, whitespace, documentation, build pins and license checks.
- Shader scopes are inclusive; duplicate labels are not exclusive subtree totals.
- Partial trailing telemetry outside complete windows is retained in the raw capture but excluded from archive bucket totals.
- This validates the captured deterministic cycle and the independent oracle tolerances; timing remains dependent on source, compiler and device configuration.

## Harness

`targets/Profile/Profile.ino` plus `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`, supported lock-owning `tools/profile_one.sh`. No permanent deep instrumentation was added.

Source/build hashes, compiler and archived ELF paths are in `build/prof/hyperlattice_ship.provenance`. Final profile ELF SHA-256: `ecccd66d081a84591334e8a4ae91d1a82f7a2046c380d9e9d95553c8e9a7cf58`; full shipping ELF SHA-256: `c1ec21db559356349545fe986b463cf24da1e10883c887267cdbd3ee2599b052`.

Baseline raw capture: `build/prof/hyperlattice_baseline_e2f5b0a3d_ship.log`. Investigation evidence, trial ELFs, disassembly, diagnostic header/diff, native and mutation logs remain in `C:/work/temp/hl-peak-1006/`; final ELF artifacts are archived under `build/prof/artifacts/hyperlattice_ship_a13c149d4788_ecccd66d081a/`.
