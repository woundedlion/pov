# MindSplatter sub-59-ms optimization plan

**Status: PARTIALLY LANDED PLAN (2026-07-25).** The live path now uses
`draw_fused_vertex`, covering part of §1.1/§3.2, and the hole-kernel multiply
barrier from §2.5 is landed. Remaining phases are proposals and must be
revalidated against the current profile.

Date: 2026-07-25  
Target: Teensy 4.0, i.MX RT1062 Cortex-M7 at 600 MHz  
Primary objective: every measured MindSplatter render frame under 59.0 ms  
Secondary objective: increase resident particle count without giving back the
59-ms guarantee

## Decision

The next optimization must replace MindSplatter's generic particle
draw/raster/AA stack with a MindSplatter-specific streaming renderer. Small
compiler annotations are exhausted, global `-O3` is both insufficient and too
large, and the full-roster image has no room for an additive specialization.

The intended end state is:

1. A fused decode, transform, gate, classify, shade, and splat pipeline.
2. A fixed-point one-dot path for the roughly 91% of classified edges that do
   not need the adaptive geodesic sampler.
3. A compact, out-of-line long-edge fallback.
4. Paired signed-axis physics.
5. A work-budgeted particle representation supporting 2,048 visible heads with
   either short trails or a stable full-tail subset.

The 59-ms requirement is a hard maximum, not the display driver's 62.5-ms spill
boundary. "Zero spills" alone is not a pass.

## Current baseline

### Device profile (superseded baseline)

The capture below is the 2026-07-25 planning baseline. The later 2026-07-27
capture at `ece0955b` measured a 60.93 ms peak and 0/1,728 spilled frames, so it
supersedes the cadence and peak figures below while leaving the stricter 59 ms
objective open.

The latest matched shipping capture is:

- Source report: 2026-07-25 shipping capture (no longer tracked)
- Raw log: `build/prof/mindsplatter_ship.log`
- Profile revision: `32478115`
- Resolution and driver: 288 x 144, four-segment POV path
- Capacity: 1,088 particles
- Trail: 23 retained snorm16 positions
- Worst individual render: **73.09 ms**
- Frames at the worst preset: 5/449 spills

The worst clean-hold window averaged:

| Scope | Time/frame | Share of particle scan |
|---|---:|---:|
| `plot_ps_raster` | 34.16 ms | 66% |
| `plot_ps_gate` | 7.78 ms | 15% |
| `plot_ps_vertex` | 4.69 ms | 9% |
| `plot_ps_tween` | 3.38 ms | 6% |
| `plot_ps_deferred` | 1.20 ms | 2% |
| `msp_particle_scan` | 51.74 ms | 100% |
| `msp_particle_step` | 2.43 ms | — |
| Timeline and clear | 0.16 ms | — |

The window's average render was 54.22 ms, while its individual peak was 73.09
ms. The workload has high view-dependent variance; optimizing an average frame
to 58 ms is not enough.

The global-`-O3` twin reduces raster cost from approximately 137 to 126
microseconds per rasterized trail, only about 8%, and still reaches an 81.53-ms
peak on a different RandomWalk trajectory. Its ITCM cost is +28,864 bytes.
Global `-O3` is not the solution.

### Current full-roster ELF

The disassembly audit used the current Phantasm ELF:

```text
revision                         810db7c6f0b0
ELF                              .pio/build/phantasm/firmware.elf
.text.itcm                       196,048 B
next ITCM bank boundary          196,608 B
literal ITCM headroom                560 B
.data + .bss                     314,176 B
allocated DTCM                   327,680 B
DTCM remaining for stack/locals   13,504 B
.bss.dma                         519,552 B
OCRAM remaining                    4,736 B
```

There is no seventh ITCM bank available. Crossing 196,608 bytes takes a 32-KiB
bank from DTCM and makes variables plus the required stack floor fail to fit.

### Disassembly findings

| Hot symbol | Code | Local/total stack | `[sp]` references | Calls |
|---|---:|---:|---:|---:|
| MindSplatter particle `draw_impl` | 2,736 B | 220/328 B | 122 | 35 |
| Direct rasterizer | 3,264 B | 332/440 B | 241 | 52 |
| Direct AA sink | 2,388 B | 24/72 B | 17 | 2 |
| Signed-axis particle physics | 1,070 B | 52/144 B | 8 | 8 |

The direct rasterizer has 95 `vldr`, 87 `vstr`, 93 integer `ldr`, and 241
explicit stack references. Its one-dot path still copies complete 48-byte
`Fragment` values, invokes the fragment shader indirectly, and enters a
440-byte general-raster stack frame.

The direct AA sink is 744 instructions. It recomputes clip state and executes a
large four-tap branch tree for every plotted dot. `ShapeShifter` also uses
`DirectAntiAliasSink`, so replacing MindSplatter's renderer cannot remove that
shared specialization.

The fragment shader calls the 220-byte `BakedPalette::get` for every sample.
That lookup interpolates two 12-byte `Color4` entries and uses a 20-byte local
temporary. The hole shader calls `hole_quintic_kernel`; its explicit multiply
barrier keeps GCC `-Os` from introducing a `__powisf2` veneer.

Physics still contains eight floating divides and three square roots in its
1,070-byte body. The six signed-axis attractors are traversed individually even
though each `+axis/-axis` pair shares perpendicular magnitude and most
intermediates.

### Required gain

The hard gap from the latest observed peak is:

```text
60.93 - 59.00 = 1.93 ms, or 3.2% of the peak
```

Because live RandomWalk captures diverge between builds, the engineering target
is stricter:

- deterministic replay maximum: <= 55.0 ms;
- repeated live-capture maximum: <= 58.5 ms;
- zero frames at or above 59.0 ms;
- zero 62.5-ms display spills.

That margin is needed before increasing particle count.

## Performance budgets

### Stage A: current 1,088-particle look

| Component | Current worst clean hold | Target |
|---|---:|---:|
| Particle scan | 51.74 ms | <= 40.0 ms |
| Raster/splat portion | 34.16 ms | <= 23.5 ms |
| Decode + transform + gate + deferred | 17.05 ms | <= 14.0 ms |
| Physics | 2.43 ms | <= 2.0 ms |
| Individual render maximum | 73.09 ms | <= 55.0 ms replay, <59.0 ms live |

These are diagnostic budgets, not independently additive promises. The fused
kernel moves work between the old scopes.

### Stage B: high-density mode

The preferred stretch target is 2,048 resident particles:

- all 2,048 heads visible;
- no more than about 20,000-22,000 retained/rendered trail edges per frame;
- <= 55-ms deterministic maximum;
- <59-ms live maximum;
- zero spills;
- full-roster ITCM below 196,608 bytes;
- at least the existing 12-KiB DTCM stack floor and 4-KiB OCRAM floor.

If 2,048 cannot meet the visual threshold, land the highest verified count
above 1,088, but keep the representation capable of 2,048.

## Phase 0 — make performance deterministic

Do this before changing the renderer. Ordinary ship/O3 comparisons follow
different feedback-driven walks and cannot isolate codegen.

### 0.1 Capture an expensive replay corpus

Add a host/device test seam that snapshots:

- particle position, velocity, life, color seed, and ring state;
- all retained snorm16 history;
- orientation quaternion;
- `MobiusParams`;
- active preset and live parameters;
- canvas clip;
- baked palette.

Capture at least:

1. The worst known preset-2 view.
2. The worst view from each of the other three presets.
3. A north/south-pole stress view.
4. A seam-crossing view.
5. A saturated-population view at each proposed particle count.

Replay exactly the same bytes after one shared physics step. Never compare
independently evolved builds for a go/no-go decision.

### 0.2 Add path counters

Build a count-only instrumentation image recording:

- resident and live particles;
- full, partial, and draining histories;
- Cartesian latitude rejects;
- Cartesian meridian rejects;
- exact-gate fallbacks;
- visible trails;
- one-dot and long edges;
- plotted samples;
- accepted AA taps;
- tap writes split by one, two, three, and four taps;
- interior versus clip-boundary splats;
- palette endpoint versus interpolated lookups;
- hole-kernel early-outs;
- framebuffer cache-line addresses or, at minimum, unique 32-byte lines per
  bounded batch.

The existing 91% one-dot figure is strong evidence, but it must be refreshed on
the current four-preset, 1,088-particle workload.

### 0.3 Add cycle and stall scopes

Use short batches so the auxiliary DWT counters do not overflow. Record:

- `CYCCNT`;
- `CPICNT`;
- `LSUCNT`;
- `EXCCNT`;

around:

- history decode;
- Mobius plus rotation;
- trail reject;
- edge classification;
- palette and alpha;
- one-dot AA;
- pixel blending;
- long-edge fallback;
- signed-axis physics.

Keep instrumentation builds separate from timing builds.

### 0.4 Automate codegen inspection

For every candidate, archive:

- `.text.itcm`, `.data`, `.bss`, and `.bss.dma`;
- symbol size deltas from `arm-none-eabi-nm -C -S --size-sort`;
- disassembly of the streaming kernel, one-dot splat, fallback, and physics;
- local and total stack;
- instruction mix;
- indirect calls in inner loops;
- stack references and obvious spill/reload chains;
- unexpected `memset`, `memcpy`, `__powisf2`, divide, and math veneers.

Reject a candidate when a hot one-dot iteration contains an indirect call or a
full `Fragment` copy.

## Phase 1 — replace the particle renderer

This is the main sub-59-ms lever.

### 1.1 Build one fixed-resolution streaming kernel

Implement a fixed 288 x 144 MindSplatter renderer instead of another broad
`Plot::ParticleSystem` template variant. Keep the public effect template, but
route the device specialization to one non-cloned kernel.

For each particle:

1. Read ring metadata directly.
2. Decode one snorm16 sample.
3. Apply the fused homogeneous Mobius map.
4. Apply the precomputed rotation matrix.
5. Store only transformed xyz and compact projected coordinates.
6. Run the whole-trail quadrant reject.
7. Classify visible edges.
8. Shade and splat one-dot edges directly.
9. Dispatch only long edges to the out-of-line fallback.

Do not construct an `ArenaVector<Fragment>`. Do not construct the parallel
original-position vector. Re-decode an original snorm16 point for the hole
shader only after the trail survives; that second six-byte load is cheaper than
materializing 12 extra bytes for every point of every rejected trail.

A candidate compact point is:

```text
transformed xyz       12 B
column Q8.8            2 B
row Q8.8               2 B
edge flags             1 B per outgoing edge, separate array
```

At a 23-point maximum this needs about 391 bytes of scratch, but it is arena
scratch rather than a nested C stack frame. Benchmark a two-point rolling
version as well; prefer it if it can preserve the whole-trail gate without
repeated transforms.

### 1.2 Split the common and uncommon paths

The common function should contain only:

- one-dot visible edge;
- non-polar finite coordinates;
- direct palette/alpha evaluation;
- direct framebuffer splat.

Outline:

- exact geodesic fallback;
- pole and antipodal cases;
- clip-boundary splats;
- assertions and diagnostic traps;
- short-history setup.

Use measured branch probabilities to place `__builtin_expect` only after the
counter run. Verify the result in disassembly; source-level hints are not
success criteria.

### 1.3 Reuse the already-linked generic rasterizer for long edges

Evaluate two fallbacks:

1. Construct only the two endpoint fragments and route the edge through the
   already-linked `PipelineRef` geodesic rasterizer.
2. Implement a compact MindSplatter-only long-edge walker accepting endpoint
   xyz, age, seed, and hole alpha.

The first option has indirect-call overhead but near-zero marginal rasterizer
code because the generic `PipelineRef` implementation already exists for other
effects. The second may run faster but must justify every added ITCM byte.

The fallback is not allowed to retain the current 3,264-byte direct rasterizer
just for roughly 9% of edges. If it does, the replacement has failed the
code-size objective.

### 1.4 Make the replacement net-negative in ITCM

Routing MindSplatter away from the current path should remove approximately:

```text
MindSplatter draw_impl                 2,736 B
Direct-AA rasterizer                   3,264 B
                                      -------
gross effect-specific code              6,000 B
```

The shared direct-AA sink remains for ShapeShifter. Other shared helpers may
also remain. Require at least **2,048 bytes net ITCM
reclaimed** after the new kernel and fallback link. Preferred reclaim is 4 KiB
or more. This creates room for later physics work and reduces the risk of
layout noise crossing the 196,608-byte boundary.

## Phase 2 — replace the one-dot shade and AA path

Run these experiments in order and retain the smallest visual change that meets
the time budget.

### 2.1 Pre-bind frame and clip state

Build a small frame context once:

```text
framebuffer base
row start/end
column start/end and wrap
row stride
palette base
global opacity Q16
Mobius coefficients
rotation rows
```

The one-dot loop must not call `Canvas::clip()`,
`ClipRegion::x_clip()`, or resolve the framebuffer base for every sample.

Split splats into:

- interior 2 x 2 footprint: no clip branches;
- boundary footprint: explicit tap checks.

The trail gate proves possible visibility, not that all four taps are inside, so
the boundary variant remains required.

### 2.2 Add a Q8.8-coordinate, Q16-alpha splat

Convert projected row and column once per transformed point:

```text
integer pixel = coord >> 8
fraction      = coord & 0xff
```

Use a 256- or 257-entry smootherstep table to map the fractional byte to a Q16
axis weight. Place the constant table in flash, then test whether copying its
small hot subset to DTCM is beneficial; DTCM is too constrained to assume that
it is.

Compute:

```text
wx0 = 65535 - wx1
wy0 = 65535 - wy1
tap weight = ((wx * wy + 32768) >> 16)
tap alpha  = ((base_alpha * tap weight + 32768) >> 16)
```

This removes four float-to-integer alpha conversions, repeated clamps, two
quintic float polynomials, and `floorf` from the hot splat. Q8.8 introduces
less than 1/256-pixel coordinate error.

Benchmark three variants:

1. Existing float coordinates with base alpha converted to Q16 once.
2. Q8.8 coordinates with a smootherstep LUT.
3. Q12.4 or a small interpolated LUT if Q8.8 code/table traffic is not a win.

Do not assume fixed point wins on Cortex-M7; select by device cycles and
`LSUCNT`.

### 2.3 Remove `BakedPalette::get` from the one-dot path

The particle already owns a uint16 color seed and trail progress is known from
the ring index. Test:

1. Exact two-entry interpolation, expressed directly in the kernel without
   `Color4`.
2. Nearest 8-bit LUT entry:

   ```text
   palette_index = high_byte(color_seed + progress_q16)
   ```

3. A 9-bit phase with one-bit interpolation if nearest-entry banding is visible.

The nearest-entry error is at most one palette bin and removes the 220-byte
lookup call, its local `Color4`, float wrap, and per-sample color interpolation.
Require image-difference evidence before accepting it.

Precompute `progress_q16` for all possible history lengths 2..23 in flash. This
removes the per-trail float reciprocal and repeated float index conversion.

### 2.4 Make alpha integer end-to-end

For a one-dot point:

```text
fade = min(progress, hole_alpha)
base_alpha = palette_alpha * fade * fade * opacity
```

Evaluate a Q16 implementation using widening 32-bit multiplies. Keep hole alpha
as float initially and convert once. If the palette is opaque in the generated
palette configuration, prove that fact and remove its alpha multiply.

The framebuffer blend remains the existing q16 source-over operation so channel
semantics and draw order do not change.

### 2.5 Replace the hole kernel only if it still matters

The deferred hole pass is currently about 1.2 ms/frame, so this is not a
first-line lever.

The explicit `t2=t*t; t3=t2*t` barrier is already present in
`hole_quintic_kernel`. If the pass remains material after fusion:

- test a 128- or 256-entry LUT keyed by `max(abs(x),abs(y),abs(z))`;
- test a polynomial directly in `1-m`, where
  `m=max(abs(x),abs(y),abs(z))`, avoiding `fast_acos`;
- preserve the `m < cos(EVENT_HORIZON)` exact early-out.

Measure framebuffer error specifically around the six event horizons.

## Phase 3 — reduce front-end work

### 3.1 Pull the quadrant reject through the Mobius transform

A Möbius map sends circles to circles. Each output clip halfspace boundary is a
circle on the sphere; after composing inverse orientation and inverse Möbius,
it is another source-space circle.

Once per frame:

1. Express the active latitude and meridian clip boundaries as sphere planes.
2. Pull those planes back through the rotation and Möbius transform.
3. Represent each pulled-back circle as a source-space plane.
4. Test decoded source points against those planes before the full Mobius
   transform.

The renderer connects transformed endpoints with geodesic edges, so endpoint
classification alone is not sufficient. Derive a conservative slack bound from
the maximum local Möbius scale over the trail, or fall back whenever that bound
is uncertain.

This experiment is successful only if it rejects far trails before their
per-point division. The already-rejected compact-materialization experiment
performed the transform first and saved only 0.08-0.20 ms; do not repeat it.

### 3.2 Fuse projection with the gate

For non-rejected trails, compute row and column once and retain Q8.8 values for
both gate and splat. No later path may call `vector_to_pixel` for a classified
one-dot edge.

The present projection uses `fast_acos` for row and `fast_atan2` for column.
Evaluate gate-only approximations only if they are conservative. Rendered
coordinates still need the current visual accuracy or the validated Q8.8
equivalent.

### 3.3 Consider a rolling transform window

The gate currently wants whole-trail extrema while the draw path wants adjacent
pairs. Compare:

- compact full-trail transformed scratch;
- two passes over snorm16 history with a rolling two-point window;
- block processing four points at a time to expose independent divides and
  improve Cortex-M7 scheduling.

The block-of-four variant is especially worth disassembling: four independent
Mobius denominators may give GCC enough work to hide VFP latency without
creating the 440-byte spill frame seen in the generic rasterizer.

## Phase 4 — pair signed-axis physics

All six attractors are fixed signed unit axes with shared live strength and
horizon parameters. Replace the generic six-attractor loop in this
specialization with three axis pairs.

For each axis coordinate `q`:

```text
cross_sq = |p|^2 - q^2
dist_plus_sq  = |p|^2 + 1 - 2q
dist_minus_sq = |p|^2 + 1 + 2q
```

The common far-from-horizon path should:

1. Handle `+axis` and `-axis` kill/horizon exceptional cases first.
2. Share `cross_sq` and its reciprocal square root.
3. Compute both distance reciprocals with one pair-product divide:

   ```text
   inv_pair  = 1 / (dist_plus_sq * dist_minus_sq)
   inv_plus  = dist_minus_sq * inv_pair
   inv_minus = dist_plus_sq  * inv_pair
   ```

4. Accumulate both tangent impulses before writing velocity.

This can reduce six distance divides to three. Guard very small distances before
forming the product and retain the current exact horizon behavior.

Then evaluate `fast_rsqrt(cross_sq)` versus `1/sqrtf(cross_sq)`. Accept the
approximation only if long-run trajectory and framebuffer metrics remain within
bounds; a small per-step force error can accumulate over 160 frames.

Also compute `speed_sq` and `axis_sq` before square roots so stationary particles
skip normalization without two `Vector::magnitude()` calls.

Target physics reduction:

- 1,088 particles: 2.43 -> <= 1.8 ms;
- 2,048 particles: projected <= 3.5 ms.

## Phase 5 — increase particle count without increasing unbounded work

Do not simply change `NUM_PARTICLES`. Capacity, resident population, trail
memory, and rendered work are separate controls.

### 5.1 First land a safe intermediate count

With the current emitters and `max_life=160`, the natural resident ceiling is
about:

```text
8 emitters * (159 live frames + 22 draining frames) = 1,448 particles
```

After the 1,088-particle path is below 55 ms, test 1,280 and then 1,448
residents. Hold the population saturated in the replay; do not infer support
from pool capacity.

### 5.2 Implement a total trail-work budget

For 2,048 particles, a ten-sample trail occupies 100 bytes per particle and the
pool is about 200 KiB. Its maximum edge count is:

```text
2,048 * 9 = 18,432 edges
```

The current 1,088 x 23 configuration can present:

```text
1,088 * 22 = 23,936 edges
```

Thus 2,048 x 10 has 23% fewer maximum edges, but almost twice the head and
front-end count. It is the first full-population configuration to test after
the streaming kernel.

Drive and sustain 2,048 residents by increasing emission rate or lifetime in
the benchmark. A capacity of 2,048 with the present spawn/lifetime schedule
does not exercise 2,048 particles.

### 5.3 Test "all heads, stable half tails"

If 2,048 x 10 shortens the visual trails too much, separate particle state from
trail ownership:

- state for all 2,048 particles;
- a head rendered for every particle;
- full 23-sample trails for a deterministic 768-1,024-particle subset;
- stable ownership keyed by particle identity/seed, not frame number.

This preserves current tail density while doubling moving heads. Do not rotate
the tail subset every frame; that creates temporal flicker.

The representation should use:

- a dense 28-byte head-state pool;
- a separate fixed trail slab;
- a compact owner/free list;
- no per-particle heap or arena allocation.

### 5.4 Optional multi-resolution history

If ten consecutive samples look too short, test ten samples spanning a longer
time interval, with denser recent history and sparse old history. This may
increase long-edge subdivision, so it is a visual experiment, not a presumed
CPU win.

Reject it if plotted sample count rises back toward the 23-sample baseline.

## Phase 6 — codegen experiments after the structural win

These are bounded A/Bs, not substitutes for the new renderer.

### 6.1 Function-specific compiler flags

For the streaming translation unit, compare:

- `-Os`;
- selective `-O2`;
- selective `-O3 -fno-unswitch-loops`;
- `-O3 -fno-unswitch-loops -fno-tree-loop-distribute-patterns`;
- explicit `noinline` on fallback and diagnostics;
- explicit `always_inline` only on the tiny q16 multiply and common splat leaf.

The `-fno-tree-loop-distribute-patterns` candidate specifically guards against
the previously observed per-sample 48-byte `memset` regression.

Select on device cycles, stack, and net ITCM. Source size or host timing is not
evidence.

### 6.2 Block scheduling

Try processing two or four independent particles/points together where it
exposes independent VFP divide, square-root, or framebuffer address chains.
The Cortex-M7 is in-order; independent chains can matter more than fewer scalar
operations.

Reject any unroll that materially raises spills. The disassembly threshold is
no more than 64 total bytes of local stack for the hot one-dot leaf and no
growth in stack references per dot.

### 6.3 Toolchain bake-off

Only after the kernel is deterministic and isolated, compile that translation
unit with a newer compatible Arm GCC and compare assembly. Do not migrate the
firmware toolchain merely because a synthetic benchmark improves.

A toolchain candidate must pass:

- full native tests;
- full Phantasm link and memory gate;
- deterministic framebuffer replay;
- interrupt/DMA device run;
- complete four-preset captures.

### 6.4 No PGO without device-representative branches

Host PGO is likely to learn the wrong cache and branch costs. If PGO is tried,
feed it the deterministic device workload and use it primarily for branch
layout/inlining decisions. A hand-split common path is expected to be safer.

## Correctness and visual gates

### Exact stages

The first streaming implementation must be pixel-identical for:

- snorm16 decode;
- Mobius and rotation;
- gate decisions;
- one-dot versus long-edge classification;
- hole alpha;
- palette interpolation;
- AA weights and tap order;
- q16 source-over blending.

This isolates structural/codegen gain.

### Approximate stages

For Q8.8 coordinates, palette quantization, hole LUTs, or reciprocal
approximations, report:

- changed pixels/frame;
- maximum and mean per-channel q16 delta;
- total absolute channel delta;
- total luminance delta;
- coverage flips;
- maximum coordinate displacement;
- edge-classification changes;
- trail rejection false positives, which must remain zero;
- temporal difference over at least one complete preset wrap.

Suggested initial limits:

- no source-space false reject;
- coordinate error <= 1/256 pixel;
- maximum channel delta <= 256 q16 units for AA-only differences;
- total luminance delta <= 0.1%;
- no visible banding, popping, seam, pole, or event-horizon regression in the
  captured video loop.

The numerical limits may be tightened after observing the exact renderer's
natural snorm16 sensitivity.

## Landing sequence and kill criteria

Land only measured wins, one independent step at a time.

1. Deterministic corpus, counters, and codegen report.
2. Exact streaming renderer with current float shade and AA.
3. Compact long-edge fallback and removal of the direct raster specialization.
4. Pre-bound frame/clip state.
5. Q16 alpha path.
6. Q8.8 smootherstep AA path.
7. Direct palette phase lookup.
8. Source-space reject, if it rejects enough trails.
9. Paired signed-axis physics.
10. 1,280 and 1,448 resident qualification.
11. 2,048 x 10 and all-heads/half-tails bake-off.
12. Optional compiler/toolchain experiments.

For each step:

- reject if the deterministic worst case improves by less than 1% unless it
  also reclaims at least 512 ITCM bytes;
- reject if ITCM crosses 196,608 bytes;
- reject if DTCM stack headroom falls below 12 KiB;
- reject if OCRAM free space falls below 4 KiB;
- reject if a hot-loop indirect call, full `Fragment` copy, or per-sample
  `memset` reappears;
- revert immediately on a device regression even if host timing improves.

## Final qualification

A renderer is sub-59 qualified only after all of the following pass:

1. Native unit and framebuffer tests.
2. Deterministic corpus at 1,088 particles:
   - maximum <= 55.0 ms;
   - all correctness/visual gates pass.
3. At least five complete four-preset live wraps across multiple starting
   orientations:
   - every individual render < 59.0 ms;
   - aggregate maximum <= 58.5 ms;
   - zero spills.
4. ISR and DMA enabled exactly as shipping.
5. Full Phantasm roster build:
   - `.text.itcm < 196,608`;
   - DTCM stack floor >= 12 KiB;
   - OCRAM floor >= 4 KiB.
6. The same qualification repeated at the selected higher-particle mode.
7. A 30-minute soak confirming no slow trajectory, pool churn, or tail-owner
   fragmentation escapes the shorter captures.

## Expected outcome

The credible gain stack for 1,088 particles is:

| Lever | Expected gain | Confidence |
|---|---:|:---:|
| Streaming one-dot renderer and removal of generic `Fragment`/callbacks | 6-10 ms | High |
| Pre-bound, integer-oriented direct AA | 4-8 ms | Medium-high |
| Direct palette phase and integer alpha | 1-3 ms | Medium |
| Front-end fusion/source reject | 1-4 ms | Medium-low |
| Paired signed-axis physics | 0.5-0.9 ms | Medium-high |

The ranges overlap and must not be summed naively. The first two levers are
expected to cover the 14.09-ms hard gap together; the remaining work creates
variance margin and particle-count headroom.

For density, the preferred order is:

```text
1,088 x 23 exact look
    -> sub-55-ms engineering margin
    -> 1,448 natural residents
    -> 2,048 x 10
    -> 2,048 heads plus a stable full-tail subset
```

Keeping 2,048 particles with 23 independent four-tap-AA trails is neither a
memory-feasible nor a credible CPU target on this device. Supporting 2,048
visible particles is credible when history work is explicitly budgeted.

## Known dead ends

Do not spend another device cycle on these without new evidence:

- global `-O3`;
- broad additive template specialization;
- the tested compact-materialization split after full transform;
- moving the particle pool to OCRAM;
- a DTCM particle SoA justified as a cache optimization;
- further framebuffer clear work;
- physics-only optimization as the sub-59 solution;
- host-only timing;
- accepting zero spills as proof of <59 ms.

## Evidence

- `effects/MindSplatter.h`
- `core/animation/sprites.h`
- `core/animation/trails.h`
- `core/render/plot.h`
- `core/render/filter.h`
- `core/color/composition.h`
- `docs/mindsplatter_particle_scaling_analysis.md`
- 2026-07-25 shipping and O3 profile captures (no longer tracked)
- `build/prof/mindsplatter_ship.log`
- `build/prof/mindsplatter_o3.log`
- `.pio/build/phantasm/firmware.elf`
