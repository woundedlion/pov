# HyperLattice optimization plan — sub-59 ms/frame

## Outcome

The quality-neutral implementation pass cuts the full-cycle shipping peak from
**172.06 ms to 81.21 ms**, a **52.8% reduction**. It preserves the render
algorithm and deterministic output signature, but does not reach the strict
59 ms target. The remaining full-cycle gap is **22.21 ms**, or **27.3% of the
current render time**.

Work stopped at that boundary. No sample reduction, shell reduction,
approximate distance function, altered coverage curve, perceptual pruning, or
mixed-resolution rendering has been applied. Continuing toward 59 ms now
requires explicit approval for a visual-error budget.

### Exact-only implementation result

All measurements below are shipping builds on the same Teensy 4.0 at 600 MHz,
using 16-frame windows and the real segmented POV driver.

| Capture | Before | After | Reduction |
|---|---:|---:|---:|
| Full four-preset cycle, worst-window mean | 172.06 ms | 81.21 ms | 52.8% |
| Fixed preset 3, worst-window mean | 149.24 ms | 71.84 ms | 51.9% |

The final fixed-preset sweep separates steady-state cost from transition cost:

| Preset | Worst-window mean | Maximum frame |
|---:|---:|---:|
| 0 `cubic-flight` | 41.63 ms | 41.92 ms |
| 1 `deep-grid` | 54.84 ms | 55.83 ms |
| 2 `dimensional-rift` | 56.28 ms | 60.11 ms |
| 3 `hypercube-flight` | 71.84 ms | 72.35 ms |
| Full cycle, including transitions | **81.21 ms** | **81.87 ms** |

Raw final captures:

- `build/prof/hyperlattice_exact_preset0_ship.log`
- `build/prof/hyperlattice_exact_preset1_ship.log`
- `build/prof/hyperlattice_exact_preset2_ship.log`
- `build/prof/hyperlattice_exact_final_preset3_ship.log`
- `build/prof/hyperlattice_exact_final_full_cycle_ship.log`

The accepted changes remove redundant construction, inline exact coverage and
coordinate transforms, use the opaque palette path, fuse transitional metrics,
reuse direction reciprocals, specialize periodic coordinate evaluation, unroll
the fixed-size edge metrics and event minimum, precompute frame invariants, and
fuse layer compositing into the cached-flash shader. The deterministic native
render signature covers every preset with multiple origins, rotations, and ray
directions; the dedicated HyperLattice suite passes. The final integrated
Phantasm release image also passes the Teensy memory gate at 191,080 bytes of
RAM1 code, leaving 2,456
bytes below the 193,536-byte derived ceiling.

Code generation, not source-level operation count alone, constrained the final
pass. Inlining layer compositing moved the complete hot shader into
`Scan::Shader::draw_cached` in cached flash and produced the largest late-stage
gain, from 80.74 ms to 72.74 ms on fixed preset 3. Conversely, compile-time
dispatch of the 4D kernel removed runtime branches but expanded the generated
path and regressed the same capture to 88.86 ms. Two alternative ordered-event
schedulers, a sentinel cursor layout, partial metric bounds, forced outlining,
and broader optimization attributes also regressed or were neutral and were
discarded.

### Permission-gated continuation plan

If a visual-error budget is approved, evaluate changes one at a time in this
order:

1. Replace the square-root wire-distance ramp with a fitted squared-distance
   ramp. This removes roughly one square root per evaluated plane while keeping
   geometry and sampling unchanged.
2. Add a conservative contribution cutoff based on remaining transmittance,
   fog, and the minimum encodable output alpha. This may omit contributions
   that cannot normally survive quantization, but it is not mathematically
   identical to the current floating-point composite.
3. Render only distant shells at 2x2 resolution and bilinearly reconstruct
   their premultiplied contribution. Keep the nearest shell at full resolution
   and widen its distant AA footprint consistently.

For each candidate, compare all presets and transitions against the pinned
reference with per-channel error histograms, changed-pixel counts, temporal
difference captures, and side-by-side on-device video. Accept a candidate only
after the visual threshold is agreed and two full 380-second captures both
remain below 59.0 ms. Uniform whole-frame resolution reduction is a fallback,
not the first recommendation.

## Original analysis and implementation rationale

The original profile showed that HyperLattice needed an algorithmic
optimization, not a compiler-flag change. The full-cycle shipping capture
peaked at 172.06 ms render time and spilled all 2,672 frames. Global `-O3`
lowered the peak to 148.24 ms and still spilled all 2,816 frames.

The deep capture attributes the work to two multiplicative terms:

1. HyperLattice shades 10,658 sphere samples per segmented frame.
2. The worst preset expands each sample into 10.58 ordered lattice events,
   9.74 plane evaluations, and 3.82 visible-layer composites.

The first implementation pass should remove compiler-visible redundancy and
hot flash calls. The second must replace the event scheduler and reduce the
cost or count of plane evaluations. If those exact-rendering changes do not
reach 65 ms, use mixed-resolution rendering for distant shells; a uniform
2x2 coarse fallback has enough theoretical margin but carries greater visual
risk.

## Evidence

### Acceptance baseline: full four-preset cycle

The existing 2026-08-24 full-cycle reports remain the acceptance baseline.
HyperLattice's render source has not changed since those captures.

| Configuration | Peak render | Slowest clean shader window | Spilled |
|---|---:|---:|---:|
| Shipping selective-`O3` | 172.06 ms | 167.43 ms | 2,672/2,672 |
| Global `-O3` | 148.24 ms | 143.80 ms | 2,816/2,816 |

Global `-O3` buys only 13.8% at the full-cycle peak. The display threshold is
62.5 ms, but this project target is stricter: **peak render below 59.0 ms**.

### Matched fixed-preset codegen pair

Commit `4bf6551041ed` was captured on COM3 with preset 3
(`hypercube-flight`) held, window 16, for 70 seconds per configuration.
Both logs validate to 1.2 ppm or better and carry attested GCC 15.2.1 ELFs.

| Configuration | Peak render | Worst shader mean | Spilled |
|---|---:|---:|---:|
| Shipping | 149.24 ms | 145.47 ms | 384/384 |
| Global `-O3` | 125.36 ms | 123.97 ms | 544/544 |

The matched peak speedup is 1.19x. These shorter fixed-preset runs do not cover
enough flight phase to replace the 172.06 ms full-cycle acceptance peak; their
purpose is a controlled codegen comparison.

Raw captures:

- `build/prof/hyperlattice_p3_ship.log`
- `build/prof/hyperlattice_p3_o3.log`

### Deep worst-preset workload

The corrected deep capture uses the same source commit and fixed preset. Its
worst mean window is frames 273–288:

| Quantity | Per frame | Per pixel |
|---|---:|---:|
| Sphere samples | 10,658 | 1.00 |
| Ordered events | 112,802 | 10.58 |
| Plane evaluations | 103,791 | 9.74 |
| Visible-layer composites | 40,671 | 3.82 |

Instrumentation raises the same-window shader mean from 145.47 ms to
182.21 ms, so deep absolute times are not release timings. Their ratios are
stable enough for prioritization:

| Deep bucket | Inclusive/exclusive share of deep shader |
|---|---:|
| Plane evaluation | 41.3% |
| Event search/cursor work, exclusive of plane and composite scopes | 33.8% |
| Ray setup and trace-loop residual | 13.8% |
| Palette/layer composite | 8.1% |
| Remaining shade/scan work | 3.0% |

The four-preset count sweep shows that event multiplication, not pixel count,
drives preset cost:

| Preset | Events/pixel | Planes/pixel | Layers/pixel |
|---:|---:|---:|---:|
| 0 `cubic-flight` | 6.20 | 5.20 | 2.27 |
| 1 `deep-grid` | 8.68 | 7.68 | 3.51 |
| 2 `dimensional-rift` | 7.42 | 6.42 | 2.13 |
| 3 `hypercube-flight` | 10.58 | 9.74 | 3.82 |

The preset 0–2 runs used an earlier setup sub-scope that perturbed setup
codegen. Use their call counts and relative scaling, not their absolute times.
The corrected preset-3 log is
`build/prof/hyperlattice_p3_ship_deep.log`.

## ELF and generated-code findings

### Size and placement

| Item | Shipping | Global `-O3` | Delta |
|---|---:|---:|---:|
| Profile FLASH code | 53,328 B | 64,296 B | +10,968 B |
| Profile RAM1 code | 22,056 B | 30,728 B | +8,672 B |
| `trace_layers` symbol | 1,152 B | 3,108 B | +1,956 B |

`Scan::Shader::draw_cached` stays in flash; the dominant `trace_layers`
function runs from ITCM. Global `-O3` expands that function 2.70x and inlines
`Mat4::apply`, `fast_rsqrt`, the layer callback, and `BakedPalette::get`.
Shipping leaves them out of line. The larger trace is faster but nowhere near
the required speedup.

The full Phantasm image currently uses 192,488 B of RAM1 code against the
193,536 B derived ceiling: only **1,048 B** of gate headroom. Applying global
`-O3` behavior to the whole trace would add about 1,956 B before secondary
effects, so it is not shippable without reclaiming code or narrowing the
optimized region.

The landed profiling scopes and fixed-preset hook are release-inert: the full
Phantasm ELF and HEX are byte-identical to their parent build.

### Dynamic pathologies at the worst preset

The shipping disassembly exposes work that `-Os` does not eliminate:

- `trace_plane` zero-initializes a 16-byte `Vec4` with `memset`, then overwrites
  all four values. That is about **103,791 small memsets/frame**.
- `Mat4::apply` builds a 52-byte stack frame, calls `memset`, loops over four
  rows, and copies the result back through the ABI. It runs twice per pixel:
  about **21,316 calls and memsets/frame**.
- The cursor array adds one 48-byte `memset` per pixel: another 10,658/frame.
- The generic `smooth_ramp` is a flash-resident function reached through an
  ITCM veneer. It contains a hardware divide. Wire and near coverage call it
  for every plane, and shell horizons call it again: approximately
  **0.21–0.25 million flash calls/frame**.
- Source-count inference puts the frame near **0.4–0.45 million scalar
  divides**, including cursor setup, projected width, smooth ramps, horizon
  normalization, and final alpha normalization.
- The 4D worst preset performs about **311,000 periodic floor operations** and
  **104,000 square roots/frame**.
- `BakedPalette::get` runs once per visible layer, about 40,700/frame, even
  though both HyperLattice palettes are opaque and the shader already produces
  a unit-range coordinate.

The divide and ramp counts are inferred from the measured call counts and the
source path; validate them with dedicated count counters before relying on an
individual expected saving.

## Optimization plan

### Phase 0 — guardrails and repeatable measurement

1. Keep the landed `HS_PROFILE_DEEP` scopes and fixed-preset hook.
2. Add deterministic native snapshots for representative rays from all four
   presets, including coincident-plane events, horizon fades, and long-flight
   origins/rotations.
3. Add a low-overhead count-only mode for ramp, divide-class, horizon, early
   exit, and coincident-event counts. Do not use RAII deep timings as absolute
   release costs.
4. After every optimization batch, run a fixed-preset A/B first, then the full
   380-second four-preset cycle. Keep compiler, board, and capture window fixed.

### Phase 1 — exact-semantics codegen cleanup

Do these as separately measurable commits.

1. **Remove redundant zero initialization.** Return `Mat4::apply` through a
   direct four-value aggregate and construct the plane point from four direct
   expressions. Initialize only the cursor fields that are read. Require
   native output equality.
2. **Replace hot generic ramps.** Add HyperLattice-local always-inline kernels
   with precomputed inverse spans. Near coverage has frame-constant edges;
   horizon ramp width is exactly one. This removes flash veneers and avoids
   generic divisions without changing the cubic curve.
3. **Reuse cursor reciprocals.** Store both `magnitude` and `step`.
   `distance / abs(direction)` becomes `distance * step`, and
   `distance / step` becomes `distance * magnitude`.
4. **Use the opaque palette path.** Resolve the active palette once per frame,
   call `get_color_unit`, and remove alpha-LUT interpolation and the multiply
   by a known-one palette alpha. Preserve hit coverage as the composite alpha.
5. **Fuse transitional metrics.** Preset 2 currently computes 4D periodic
   distances and then recomputes the 3D subset. Calculate the component
   distances once and derive both metrics.
6. **Apply narrow inlining/selective `-O3` last.** Measure forced inlining of
   the two matrix transforms, palette lookup, and layer callback independently.
   Do not apply `HS_O3_FN` to the whole trace unless the full Phantasm memory
   gate passes; the current image cannot absorb the observed 1,956-byte trace
   expansion.

Phase-1 go/no-go target: full-cycle shipping peak at or below **135 ms** with
no material image difference and no RAM1 gate regression. These changes are
necessary but are not expected to reach 59 ms alone.

### Phase 2 — event and plane-kernel redesign

1. **Replace the two-scan event loop with a fixed four-stream merge.** Select
   the nearest active cursor with an unrolled 4-way minimum, advance only the
   selected stream, and handle the rare tolerance-grouped ties explicitly.
   Preserve front-to-back order and the max-coverage rule at coincident planes.
2. **Specialize the plane kernel.** Compute only coordinates other than the
   plane axis. Test `abs(x - nearbyint(x))` as the periodic-distance spelling
   and require bitwise/native equivalence at half-cell boundaries before use.
3. **Cache per-stream increments.** For the second and third shell on one axis,
   advance the other lattice coordinates by precomputed direction ratios
   instead of rebuilding the full 4D point and refolding every component.
4. **Add conservative far-event pruning.** Bound the maximum remaining
   contribution from fog, remaining alpha, and remaining event count. Stop only
   when the bound is below the minimum encodable alpha. Record pruned events and
   compare images over long-flight phases.
5. **Evaluate squared-distance wire coverage.** Replacing `sqrt(metric_sq)`
   with a ramp in squared-distance space removes one square root per plane but
   is not exactly the same curve. Land it only if image-error and motion tests
   pass a separately agreed visual threshold.

Phase-2 target: full-cycle shipping peak at or below **85 ms**, with event
scheduling at most 12 ms/frame and plane evaluation at most 25 ms/frame in a
low-overhead attribution build.

### Phase 3 — quality-scalable path if exact rendering remains above 59 ms

Prefer mixed resolution by depth, not a uniform blur:

1. Trace and shade the nearest shell at full 288x144 segmented resolution.
2. Evaluate distant shells on a 2x2 grid and bilinearly reconstruct their
   premultiplied contribution before front-to-back composition.
3. Increase projected AA width for the coarse distant pass and dither the
   full/coarse boundary over time.
4. If that architecture is too invasive, prototype uniform 2x2 shading as the
   feasibility bound. Quartering the dominant sample work gives ample timing
   margin from a 172 ms baseline, but it is the highest visual-risk option.

Trigger Phase 3 if the exact Phase-2 path remains above **65 ms** after two
independent full-cycle captures.

## Target budget

Allocate the 59 ms frame as follows:

| Bucket | Peak budget |
|---|---:|
| Ray setup and scan traversal | 7 ms |
| Event scheduling/cursor maintenance | 12 ms |
| Plane metric and coverage | 25 ms |
| Palette and front-to-back composite | 8 ms |
| ISR, canvas, state, and safety margin | 7 ms |
| **Total** | **59 ms** |

These are engineering gates, not measurements from the distortion-heavy deep
capture. Update them when low-overhead counters are available.

## Acceptance criteria

An optimization series is complete only when all of the following hold:

1. Shipping selective-`O3` completes a full four-preset wrap with
   `HS_PROFILE_WINDOW=16`, an epoch of at least 3,600 revolutions, and no epoch
   boundary in the capture.
2. Peak per-frame **render**, including transitions and live ISR time, is below
   **59.0 ms**; spilled frames are **0/N**.
3. A second full-cycle run reproduces the result on the same board. A global
   `-O3` twin records the remaining compiler ceiling.
4. HyperLattice native tests, deterministic image snapshots, the complete
   native suite, `pio run -e phantasm`, and all Teensy memory gates pass.
5. The full Phantasm image remains under the 193,536-byte RAM1-code ceiling;
   performance code is not accepted by silently consuming the 3,072-byte
   boundary safety ratchet.
6. Any approximate coverage, pruning, or mixed-resolution change includes
   long-flight visual comparisons for all presets and their transitions.
