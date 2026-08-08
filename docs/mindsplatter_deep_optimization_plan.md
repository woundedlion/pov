# MindSplatter deep optimization plan

**Status: PLAN (2026-08-01).** This supersedes the execution strategy in
`docs/mindsplatter_optimization.md` without erasing that document's experiment
history. The next optimization effort targets a durable 8-12 ms reduction, not
another marginal improvement near the display deadline.

Target: Teensy 4.0, i.MX RT1062 Cortex-M7 at 600 MHz  
Effect configuration: 288 x 144, segmented POV driver, 1,088 particles,
23-sample quantized trails, four presets

## Decision

The next pass must optimize MindSplatter's long-edge trail renderer. It must
not begin with physics, clearing, compiler annotations, or a universal one-dot
fast path.

The evidence now shows that the old "roughly 91% one-dot" premise was not true
for every preset. In the captured heavy regime, only about 3% of classified
edges were one-dot and approximately 5,000 long edges were considered per
frame. The previous streaming renderer regressed because it dispatched those
long edges individually through fresh generic raster and scratch setup.

The intended end state is:

1. A deterministic, production-resolution device replay harness.
2. A single-pass or otherwise compact long-edge traversal.
3. One trail-wide typed MindSplatter kernel with no full `Fragment` copies or
   inner-loop indirect calls.
4. A conservative source-space reject before Mobius transformation.
5. A measured fixed-point shade and anti-alias path where it is beneficial.
6. A full-roster image that remains inside the current FlexRAM allocation.

## Current evidence

### Latest standard shipping capture

The latest published shipping report is the 2026-07-27 capture at `ece0955b`:

- peak frame render: **60.93 ms**;
- spilled frames: **0/1,728**;
- clean-hold particle draw: **32.65-42.34 ms/frame** across presets;
- particle physics: approximately **1.95 ms/frame**;
- timeline and clear: approximately **0.15 ms/frame** combined;
- ISR share: approximately **4.98%**.

Zero spills is not sufficient. The peak left only about 1.6 ms before the
62.5 ms display-window boundary and missed the existing 59 ms engineering
target.

This report also predates the single-pass rasterizer and later Plot changes.
Current `HEAD` must be recaptured before any new A/B result is accepted.

### Heavy-path classification

The deep count capture recorded the following workload spread. These figures
come from an older live trajectory and must be refreshed, but they are enough
to reject a one-dot-only strategy.

| Preset | One-dot edges | Long edges/frame | Raster ms/frame |
|---:|---:|---:|---:|
| 0 | 53.8% | 1,722 | 9.25 |
| 1 | 62.0% | 2,118 | 11.81 |
| 2 | 36.6% | 3,628 | 16.56 |
| 3 | **3.2%** | **5,008** | **20.48** |

The gate and materialization front end is also material:

| Preset | Tween/transform ms/frame | Gate ms/frame | Raster ms/frame |
|---:|---:|---:|---:|
| 0 | 6.68 | 6.24 | 11.71 |
| 1 | 7.73 | 8.44 | 12.06 |
| 2 | 7.08 | 8.29 | 19.64 |
| 3 | 7.49 | 7.29 | 17.97 |

Within the gate, the per-edge tier costs about 2.8-3.9 ms/frame, while row and
column projection each cost about 1.0-1.4 ms/frame. Exact gate fallbacks are
small. This makes early whole-trail rejection more valuable than tuning the
exact fallback.

### Why the first streaming renderer failed

The reverted `528e0337` implementation correctly avoided the generic path for
one-dot edges, but passed every long edge separately to a helper that:

- opened a new scratch scope;
- allocated a two-`Fragment` vector;
- recomputed endpoint hole alpha;
- rebuilt a fragment-shader closure;
- called the generic rasterizer for one edge.

It regressed the shipping peak from 73.09 to 101.83 ms and the global-O3 peak
to 120.71 ms. A successor must retain one persistent traversal and one reused
scratch/cache context across all edges of a trail.

### Current constraints

- The full-roster Phantasm image has only about 1.5-2.4 KiB of recently
  reported ITCM headroom before the 196,608-byte bank boundary. Rebuild at the
  execution tip and use that build's exact number.
- Crossing the boundary consumes another 32 KiB FlexRAM bank and violates the
  required DTCM stack floor.
- OCRAM has only about 4.7 KiB free, so a large command or tile buffer is not
  available.
- The existing host replay seam uses a small test resolution. It proves state
  restoration and determinism but is not a production-resolution device
  benchmark.

## Performance and correctness targets

For the current 1,088-particle, 23-sample effect:

- deterministic device replay maximum: **at most 55.0 ms**;
- repeated live-capture maximum: **below 58.5 ms**;
- no frame at or above 59.0 ms;
- zero 62.5 ms display spills;
- stretch goal: **at most 50 ms peak**;
- full-roster ITCM below 196,608 bytes;
- DTCM stack reserve at least 12 KiB;
- OCRAM reserve at least 4 KiB.

Preserve particle count, trail length, blend order, palette behavior, hole
behavior, and anti-alias coverage unless an approximation passes a separately
defined image-error gate.

Do not increase particle count until the current look clears the deterministic
55 ms gate.

## Phase 0 - establish a current deterministic baseline

### 0.1 Recapture current `HEAD`

Run the standard full-cycle captures through the locked profiling workflow:

```text
bash tools/profile_one.sh MindSplatter profile    110 16
bash tools/profile_one.sh MindSplatter profile_o3 110 16
```

Requirements:

- all four presets visited;
- wrap back to preset 0 confirmed;
- no epoch reset;
- per-frame render telemetry present;
- cycle/wall exactness validated;
- both configurations captured on the same board where possible;
- two shipping passes to quantify live-trajectory variance.

The old shipping and global-O3 logs follow different RandomWalk trajectories.
Their peak difference is not a valid code-generation comparison.

### 0.2 Add a production-resolution replay target

Create a single-effect benchmark mode that embeds one frozen corpus at a time
in flash and repeatedly renders it without advancing physics. Reuse the
existing `MindSplatterWhiteBox::ReplaySnapshot` representation and serialization
logic, but do not duplicate an entire live pool in DTCM.

Maintain separate corpora for:

1. the worst observed state from each preset;
2. a saturated-pool state;
3. a seam-crossing state;
4. a north/south-pole stress state;
5. the worst long-edge state;
6. the worst one-dot-heavy state.

Build one selected corpus into a profiling image through a build flag. This
keeps flash and RAM bounded while making every A/B render the exact same state.

Record:

- framebuffer hash;
- changed-pixel and channel-error statistics against the reference;
- individual-frame render cycles;
- every named stage counter;
- corpus identifier and source revision.

### 0.3 Archive code-generation evidence

For every baseline and candidate, archive:

- `.text.itcm`, `.data`, `.bss`, and `.bss.dma`;
- full-roster and single-effect image sizes;
- size-sorted symbols;
- disassembly of particle draw, gate, rasterizer, direct AA, and physics;
- local and total stack;
- indirect calls and veneers in hot loops;
- load/store and explicit stack-reference counts;
- unexpected `memcpy`, `memset`, divide, sqrt, and math-library calls.

## Phase 1 - refresh path and stall measurements

Add count-only instrumentation for:

- live and resident particles;
- full, partial, and draining histories;
- Cartesian whole-trail rejects;
- row/column prologue rejects;
- per-edge rejects;
- visible trails;
- one-dot and long edges;
- adaptive samples per long edge;
- exact-gate fallbacks;
- fragment-shader calls;
- AA tap masks with zero through four accepted taps;
- interior and clip-boundary splats;
- palette endpoint and interpolated samples;
- hole-kernel early-outs.

Add short-batch `CYCCNT`, `CPICNT`, `LSUCNT`, and `EXCCNT` measurements around:

- history decode plus vertex transform;
- trail gate;
- long-edge strategy setup;
- adaptive simulation;
- normalized replay;
- shade and palette sample;
- projection;
- AA weight generation;
- framebuffer blending.

Keep instrumentation images separate from timing images. Per-pixel RAII scopes
are prohibited because they perturb the workload.

## Phase 2 - test single-pass geodesic rasterization

The current rasterizer can emit adaptive samples while stepping instead of
first storing a step cache and replaying it. This facility landed after the
last MindSplatter capture and is the smallest credible long-edge experiment.

### Implementation

- Route only MindSplatter's direct-AA geodesic raster instantiation through
  `rasterize<..., true>`.
- Preserve one raster call per whole trail.
- Keep particle, edge, sample, and tap order unchanged where the single-pass
  stepping permits it.
- Add geodesic single-pass tests for endpoints, omit-end behavior, gaps, poles,
  seams, long arcs, and active quadrant clips.

### Acceptance

- at least 1.5 ms deterministic improvement in the heavy corpus;
- no missing lit pixels or AA coverage;
- bounded, documented channel error if normalized sample placement changes;
- no full-roster ITCM-bank crossing;
- no larger stack frame.

Reject it if tangent evaluation costs more than the eliminated replay or if
sample-placement changes are visually material.

## Phase 3 - build a trail-wide typed MindSplatter kernel

This is the primary architectural lever. It must retain long edges inside one
trail traversal instead of reproducing the failed per-edge fallback.

### 3.1 Compact control points

Replace the generic `Fragment` materialization with only the state the renderer
uses:

```text
transformed xyz                  12 B
progress                         Q16 or float
particle color phase             Q16
hole/life factor                 Q16 or float
projected endpoint row/column    only when required by the gate
```

Compare a compact full-trail array against a rolling or block-of-four version.
The full array is acceptable when it lives in arena scratch and produces
better scheduling; a large nested stack frame is not.

### 3.2 One ordered trail pipeline

For each particle:

1. Decode the quantized history once.
2. Apply the Mobius transform and cached rotation.
3. Compute the whole-trail gate and every edge classification once.
4. Return immediately when no edge survives.
5. Mark points incident to surviving edges.
6. Evaluate the hole kernel only for marked points.
7. Render surviving edges in original order.
8. Send one-dot edges directly to a preprojected splat.
9. Process long edges with one compact adaptive sampler whose scratch and
   stepping state are reused across the trail.
10. Shade directly to `{Pixel, alpha}` and write through direct AA.

The hot path must not:

- copy or zero a complete `Fragment`;
- call a type-erased fragment shader;
- allocate scratch per edge;
- reconstruct a two-point `Fragments` vector per edge;
- re-enter the generic rasterizer per edge;
- repeat clip preparation or framebuffer-base lookup per sample.

### 3.3 Code-size contract

The kernel must replace the present MindSplatter-specific `draw_impl` and
rasterizer instantiations in the final ELF. It may reuse out-of-line shared
helpers for uncommon pole, antipodal, and exact-gate cases.

Accept only if:

- deterministic peak improves by at least 4 ms;
- the inner long-edge loop has no indirect call;
- stack usage materially decreases;
- full-roster ITCM remains below its derived ceiling;
- preferably at least 2 KiB net ITCM is reclaimed.

## Phase 4 - reject trails before the Mobius transform

The current Cartesian gate rejects many trails only after every history point
has been decoded, Mobius-transformed, rotated, and materialized. Move a
conservative tier ahead of that work.

### 4.1 Pull clip circles back into source space

Each output quadrant boundary is a circle on the sphere. The inverse
orientation and inverse Mobius map take it to another sphere circle, which can
be represented as:

```text
n dot p + d = 0
```

Once per frame:

1. Build the active output clip boundaries including AA margin.
2. Apply inverse rotation.
3. Apply the inverse Mobius matrix `(d,-b;-c,a)`.
4. Fit or derive each pulled-back circle plane.
5. Compute conservative numeric slack for projection, geodesic bow, and float
   error.

For each source trail, scan cheap plane evaluations before doing its Mobius
transforms. Reject only when every point lies outside one expanded boundary.
Singular, near-pole, or uncertain cases fall through to the exact current
path.

### 4.2 Proof and kill criterion

Before device timing, test millions of randomized combinations covering:

- all four clips;
- identity and animated MindSplatter Mobius parameters;
- arbitrary orientations;
- seams and both poles;
- maximum-length and short histories;
- edges at the conservative margin.

False rejection is forbidden. A false accept is safe and only reduces the
performance benefit.

Do not land the tier unless it removes enough transform plus gate work to save
at least 2 ms in a deterministic heavy corpus.

## Phase 5 - specialize shade and anti-alias arithmetic

Run these experiments only after Phase 3 scopes show what remains hot.

### 5.1 Q16 color phase

The particle already stores a `uint16_t` color seed. Store trail progress as a
Q16 phase and add the seed with natural wrap. Add a `BakedPalette` sampler that
derives LUT index and interpolation weight directly from the phase, avoiding:

- float seed normalization;
- float wrap;
- multiply by 255;
- float-to-integer index conversion.

Keep interpolated RGB initially. Test nearest-entry or reduced interpolation
only as a separate visual approximation.

### 5.2 Integer alpha propagation

The palette is proven opaque for MindSplatter. Compute the remaining life,
hole, progress, opacity, and squared fade as Q16 where measured beneficial.
Convert the hole result once per control point rather than once per raster
sample.

### 5.3 Fixed-point splat coordinates

Compare:

1. current float coordinates with one Q16 base-alpha conversion;
2. Q8.8 coordinates plus a 257-entry smootherstep table;
3. a smaller interpolated table if flash traffic dominates.

Retain the existing full-interior and masked-boundary branches. Preserve tap
order and `Pixel::lerp16` source-over semantics.

For every approximation report:

- changed pixels;
- added or dropped taps;
- maximum RGB16 channel delta;
- total channel delta;
- seam and pole subsets;
- device cycles and LSUCNT;
- code-size change.

The combined shade/splat work should save at least 2 ms to justify a
non-bit-exact result.

## Phase 6 - residual front-end and physics work

Only after the renderer clears 55 ms:

### 6.1 Front-end cleanup

- Precompute Q16 progress tables for history lengths 2-23.
- Hoist the framebuffer base, clip tables, Mobius coefficients, rotation rows,
  opacity, and palette base into one frame context.
- Compare scalar, pair, and block-of-four transforms. Retain blocking only when
  disassembly shows improved scheduling without added spills.
- Outline pole, antipodal, exact-gate, diagnostic, and short-history paths.

### 6.2 Signed-axis physics

Paired signed-axis physics is already landed and the whole step is about
1.95 ms/frame. Remaining experiments are secondary:

- derive perpendicular magnitude from `1-q^2` where unit-vector error permits;
- test a bounded reciprocal-square-root approximation;
- interleave independent axis pairs;
- reduce repeated horizon/kill constants.

Require trajectory and framebuffer-error bounds over long deterministic runs.
A realistic remaining gain is 0.3-0.8 ms.

## Phase 7 - optional workload scaling

Do not begin this phase until the current look has a deterministic peak at or
below 55 ms.

The present architecture cannot store 2,048 particles with 23 history samples:
the pool alone would require about 360 KiB against a roughly 266 KiB persistent
budget. Raising capacity also does not guarantee 2,048 residents because the
current eight emitters and lifetime naturally top out around 1,448 including
draining trails.

Evaluate, in order:

1. the highest safe resident count with the existing 23-sample look;
2. 2,048 particles with approximately 10 retained samples;
3. 2,048 visible heads with a deterministic 768-1,024-particle full-tail
   subset;
4. a total per-frame trail-edge budget rather than equal history for every
   particle.

Every density mode must force and sustain its intended resident population in
the replay corpus. Pool capacity alone is not evidence.

## Qualification and landing protocol

For every candidate:

1. Run the production replay correctness comparison.
2. Run a shipping deterministic device capture.
3. Run the matched global-O3 capture as a code-generation ceiling.
4. Run two complete live four-preset cycles.
5. Run the native test suite.
6. Build the complete `phantasm` image.
7. Run the memory and stack gates.
8. Audit the ELF and hot disassembly.
9. Update both configuration reports and all ranked profile tables.
10. Land one measured lever at a time.

Use `profile_one.sh` for every flash and capture so the shared device lock,
port pinning, stale-image checks, and log verification remain in force.

Small, exact changes may land for a repeatable improvement of at least 0.5 ms
when they are code-size neutral. Structural changes must deliver at least
1.5-2 ms, and the trail-wide kernel must deliver at least 4 ms.

## Known dead ends

Do not repeat these without materially new evidence:

- the reverted per-edge streaming renderer;
- optimizing only the one-dot path;
- the compact materialization split after the Mobius transform;
- global O3 as a shipping solution;
- broad additive template specialization;
- effect-local O3 annotations already measured nearly dead;
- physics-first optimization;
- framebuffer-clear work;
- moving the particle pool to OCRAM;
- a DTCM particle SoA justified as a cache optimization;
- a large tile-command buffer;
- host timing as a landing criterion;
- accepting zero spills as proof of sufficient margin;
- increasing particle count before the current renderer is stable;
- 2,048 particles with 23 complete retained trails.

## Expected outcome

The most credible winning stack is:

```text
single-pass long-edge traversal
  + trail-wide compact typed kernel
  + conservative pre-transform trail rejection
  + measured Q16 shade/AA cleanup
```

The gains overlap and must not be added mechanically. Together they have a
credible 8-15 ms ceiling, which would move the last observed 60.93 ms peak into
approximately the 46-53 ms range if the hypotheses hold. The first hard stop
is a deterministic maximum at or below 55 ms; the effort should continue toward
50 ms while safe, measured levers remain.

## Evidence

- `effects/MindSplatter.h`
- `core/render/plot.h`
- `core/render/filter.h`
- `core/animation/sprites.h`
- `core/color/composition.h`
- `tests/test_effects.h`
- `docs/profiles/shipping/profile_mindsplatter_teensy_2026-07-27.md`
- `docs/mindsplatter_optimization.md`
- `docs/MINDSPLATTER_PARTICLE_SCALING_ANALYSIS.md`
- `build/prof/mindsplatter_ship.log`
- `build/prof/mindsplatter_ship_gate_breakdown_deep.log`
- `build/prof/mindsplatter_ship_raster_peak_counts.log`
- reverted streaming implementation `528e0337` / `805f2d0a`
