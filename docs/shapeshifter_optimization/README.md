# ShapeShifter deep optimization plan

Status: executed and validated on 2026-07-31. The selected implementation is
commit `02ccaca8`; this document retains the rejected-candidate ledger and the
remaining optional research directions.

## Mission

Turn the current ShapeShifter from a cadence-flapping 93-98 ms renderer into a
reliably sub-window effect on the Teensy 4.0, then continue toward enough compute
headroom for denser shapes or a future 32 fps mode. The work is allowed to change
the renderer's architecture, add an effect-specialized Plot kernel, replace
generic type-erased paths, and use bounded mathematical approximations. It is not
allowed to silently reduce the preset workload, drop shapes, weaken segmented
rendering, or consume another FlexRAM bank.

The primary target is the shipping Phantasm image, not a single-effect global
`-O3` build. Global `-O3` is a diagnostic ceiling and code-generation oracle.

## Outcome

The shipping renderer now peaks at **43.52 ms**, down from the audited
**98.79 ms** baseline: **55.27 ms / 55.95% faster**. It spills **0/1,728**
frames, leaving **15.48 ms** against the 59 ms owner target. The global-`-O3`
twin peaks at 42.39 ms with 0/1,728 spills, so selective code generation is now
within 1.13 ms of the compiler ceiling.

| Gate | Final result |
|---|---:|
| Shipping peak render | 43.52 ms |
| Global-`-O3` peak render | 42.39 ms |
| Shipping spills | 0/1,728 |
| Full-roster RAM1 code | 189,416 B / 196,608 B |
| ITCM headroom | 7,192 B |
| Native tests | 59/59 passed |
| Firmware environments | `phantasm`, `holosphere`, `holosphere_dma` passed |
| WASM smoke | 24 effects at 96x20 and 288x144 passed |

The final Phase C preserves the shader architecture. ShapeShifter passes a
compile-time typed callable with full `Fragment` semantics through the raster;
constant color is merely the current shader body, not the renderer API. The
specialization removes runtime shader/pipeline erasure, replay, and redundant
projection work without preventing future position-, varying-, texture-, or
noise-driven shading.

## Current ground truth

The reference implementation was `effects/ShapeShifter.h` at `2c0a7a27`. The
validated implementation is `02ccaca8`. Neither is the older Plot+Scan
implementation discussed in the 2026-07-29 analysis.

The current effect:

- draws only through Plot with terminal four-tap anti-aliasing;
- draws 43-75 concentric primitives per frame in the shipping preset range;
- has four interactive shape types, but its three shipping endpoints select
  SphericalPolygon, Star, and Flower;
- continuously lerps three presets for 240 frames per leg, so a full preset loop
  is 720 frames, not the obsolete 192-frame/four-hard-cut cycle;
- uses SINE at all three preset endpoints, while all four waveforms remain part
  of the interactive control contract;
- renders a 72x144 segmented quadrant, or 10,368 pixels, per board.

### Measured baseline

Source and raw data:

- `docs/profiles/shipping/profile_shapeshifter_teensy_2026-07-31.md`
- `docs/profiles/O3/profile_shapeshifter_teensy_2026-07-31.md`
- `build/prof/shapeshifter_ship.log`
- `build/prof/shapeshifter_o3.log`
- attested source state: `5e57cf0b` plus the diff subsequently landed as
  `2c0a7a27`

| Metric | Shipping selective `-O3` | Global `-O3` reference |
|---|---:|---:|
| Peak render | 97.71 ms | 92.65 ms |
| Average `ss_draw_all` | 54.57 ms | 51.08 ms |
| Worst-window `ss_draw_all` | 93.52 ms | 88.24 ms |
| Spilled frames | 421/1296 (32.5%) | 399/1328 (30.0%) |
| Worst-window `ss_plot_dispatch` | 93.46 ms | 88.19 ms |
| Worst-window `filter_blend` | 3.14 ms | 2.89 ms |
| Profile-image RAM1 code | 43,320 B | 66,248 B |

An immediate same-source baseline recapture used for candidate comparisons
measured 98.79 ms peak and 423/1,296 spills. The small difference from the
97.71 ms report is normal bench variance; optimization deltas below use the
98.79 ms audited recapture.

Endpoint costs in the shipping capture:

| Endpoint | Workload | Peak render | Cadence |
|---|---|---:|---:|
| SphericalPolygon | count 74, sides 3 | 97.71 ms | 8 fps |
| Flower | count 70, sides 3 | 92.13 ms | 8 fps |
| Star | count 43, sides 6 | 82.43 ms | 16/8 fps |

The final blend is only about 3% of the worst frame. RandomWalk, preset
interpolation, palette lookup, clear, and bookkeeping are collectively below
0.2 ms. The optimization target is therefore curve construction, edge setup,
adaptive sampling, projection, and Plot dispatch.

### Memory baseline

The attested full-roster Phantasm image from the same source state reports:

- RAM1 code: 185,016 B;
- derived ceiling: 196,608 B (six 32 KiB ITCM banks);
- headroom: 11,592 B;
- RAM1 local-variable headroom: 13,504 B;
- RAM2 free: 4,736 B.

No accepted optimization may cross the 196,608 B RAM1-code ceiling. Prefer to
finish at or below 192,512 B, preserving at least 4 KiB inside the current bank.
Do not buy speed with a large persistent RAM2 allocation; that region is already
effectively full.

### Disassembly baseline

The attested single-effect ELFs are:

- `build/prof/artifacts/shapeshifter_ship_5e57cf0b53aa_cd9d0ed8d19d/profile.elf`
- `build/prof/artifacts/shapeshifter_o3_5e57cf0b53aa_29c338f4d87e/profile.elf`

The relevant symbol sizes are:

| Symbol | Shipping | Global `-O3` |
|---|---:|---:|
| shared `Plot::rasterize` | 7,572 B | 7,712 B |
| O3 planar process lambda | folded into rasterizer | 5,320 B |
| ShapeShifter `dispatch_plot` | 728 B | 4,200 B |
| ShapeShifter `draw_all` | 428 B | 1,764 B |

Static instruction census of the shipping rasterizer finds 24 indirect branch
sites, 35 `vdiv` sites, 13 `vsqrt` sites, 268 `vldr` sites, and 177 `vstr`
sites. These are static code counts, not executed-event counts, but they expose
the code-generation pressure. The generic rasterizer retains indirect shader
and pipeline calls in its hottest replay loops.

The sampler audit also finds `sinf`/`cosf` calls through flash veneers inside
the per-vertex loops for Ring, Star, and Flower. At high count, the same regular
angular topology is rebuilt with libm hundreds to thousands of times per frame.

Shipping already applies the selective-`O3` Plot-raster and blend regions. The
global-`O3` image spends another 22,928 B of RAM1 code for only a 5.2% peak win.
That is strong evidence that the next large gains require less work and a more
specialized data path, not a broader optimization flag.

## Success criteria

### Performance gates

Use render time, never wall time, for all cadence decisions.

1. Cadence gate: peak render below 59.0 ms and zero frames above 62.5 ms over a
   complete three-preset loop in the shipping build.
2. Deep-optimization target: peak render at or below 45.0 ms, with every
   endpoint's worst clean window below 42.0 ms.
3. Stretch target: peak render at or below 31.25 ms without reducing Count,
   Sides, anti-alias quality, or the interactive control surface.
4. No endpoint may regress more than 2% while another endpoint improves unless
   the aggregate peak and cadence target both improve.
5. ISR/DMA measurements remain reported separately. They are included in DWT
   foreground scopes but are not a ShapeShifter optimization target unless the
   renderer falls below 45 ms and ISR cost becomes material.

### Correctness and visual gates

Classify each candidate before implementation:

- Exact: must be framebuffer byte-identical to the reference over the full
  matrix below.
- Reassociated: same algorithm with changed float evaluation order; must have
  identical lit-pixel topology and a bounded RGB16 error.
- Approximate: sampling or trig is changed; must pass an explicit perceptual
  error budget and device review. Never merge it under an “exact” label.

Every accepted candidate must preserve:

- outer-to-inner ring order and source-over blend order;
- all four Shape selections and all four Function selections;
- Count 1..144 and Sides 3..16 at 288x144;
- pole, antipode, longitude-seam, and segmented clip behavior;
- pause behavior, preset timing, strobe behavior, and parameter names/ranges;
- deterministic native rendering for a fixed seed;
- non-profile WASM and Teensy builds.

For approximate candidates, define the budget relative to a natural sub-pixel
motion oracle: compare the candidate/reference delta against the reference
renderer moved by 1/16 pixel. Candidate total RGB16 error must be no greater
than 25% of that oracle, it may introduce no isolated holes, and seam/pole masks
must show no concentrated error. Also report differing pixels, coverage XOR,
maximum channel delta, total absolute RGB delta, and SSIM or an equivalent image
metric. The numeric budget must be frozen before the candidate is timed.

### Size and code-generation gates

- Full-roster RAM1 code at or below 196,608 B; preferred at or below 192,512 B.
- No additional FlexRAM bank and no reduction of the 12,288 B stack floor.
- No per-frame heap allocation.
- Scratch allocations must remain within the existing arena contract and be
  checked at the maximum Count/Sides combination.
- The final hot loop must contain no indirect shader or pipeline branch.
- The final regular-shape vertex loop should contain no flash `sinf`/`cosf`
  call, no general `fmodf`, and no redundant normalize of an analytically unit
  vector.
- Record symbol size, stack frame, call graph, `vdiv`/`vsqrt`, and spill counts
  for every accepted code-generation change.

## Measurement program

### 1. Preserve an oracle before changing the renderer

Add a test-only/reference entry point for the current ShapeShifter path. A
candidate kernel and the reference kernel must render identical initial canvas
contents with identical parameters into separate buffers in one process. Do not
use a stored golden hash: dual rendering keeps the oracle useful as unrelated
palette and engine code evolves.

The production-resolution matrix must cover:

- Shape: all four values;
- Function: all four values;
- Sides: 3, 4, 5, 6, 8, 12, 16;
- Count: 1, 4, 16, 43, 70, 75, 144;
- phase: 0, values just below/above wrap, and eight evenly spaced turns;
- orientation: identity, both display poles, longitude seam, and at least eight
  fixed RandomWalk seeds;
- segmented clips: all four 72x144 segment-band placements with the production
  margin;
- initial canvases: black, deterministic noise, and a non-black gradient so
  alpha/order errors are visible.

Keep the existing slider smoke test, but treat it as coverage rather than
parity evidence.

### 2. Recapture the unmodified baseline

Use the device lock and the current 720-frame preset loop. A valid pass must see
all three preset markers, wrap back to preset 0, stay inside one effect epoch,
and pass `tools/parse_profile.py ... validate`.

```bash
bash tools/profile_one.sh ShapeShifter profile 110 16
bash tools/profile_one.sh ShapeShifter profile_o3 110 16
python tools/parse_profile.py build/prof/shapeshifter_ship.log validate
python tools/parse_profile.py build/prof/shapeshifter_o3.log validate
```

Run the two configurations on the same pinned board when comparing small
deltas. For a landing candidate, use an ABBA sequence on one board and report
the median endpoint/window delta. Do not force a live device lock.

### 3. Deep phase capture

The existing `HS_PROFILE_DEEP` scopes split segment cull, adaptive simulation,
and draw replay. Capture both configurations before adding more probes:

```bash
HS_PROFILE_DEEP=1 bash tools/profile_one.sh ShapeShifter profile 110 16
HS_PROFILE_DEEP=1 bash tools/profile_one.sh ShapeShifter profile_o3 110 16
```

Deep instrumentation perturbs the absolute time. Use it to rank sub-phases and
derive ratios, then use normal captures for acceptance.

Add outer-granularity scopes or accumulators for:

- waveform/palette preparation;
- control-point generation;
- planar arc pre-pass;
- edge clip classification;
- adaptive-step simulation;
- draw replay/projection;
- fragment setup/shader;
- terminal AA/sink.

Never put a `CycleScope` around every sample or AA tap. For inner operations,
time a whole ring/edge batch and divide by count.

### 4. Count-only workload pass

Add a separately flagged, count-only build. Its timings are invalid by design;
its job is to establish the exact work amplification:

- rings submitted and whole-ring culled;
- edges submitted, culled, one-dot, planar, geodesic, and antipodal fallback;
- planar unprojects and arc-table samples;
- simulated and replayed substeps;
- vector normalizations;
- shader calls, projected samples, and AA taps;
- `steps_cache` capacity and backstop hits;
- source-over blends per output pixel.

Report each count per frame and per ring for the three shipping endpoints and
the exhaustive control cases. This pass supplies the denominator for
cycles/edge, cycles/substep, and cycles/plotted sample.

### 5. Cortex-M7 event sampling

For one frozen worst-case frame per endpoint, snapshot the M7 DWT auxiliary
counters around the large phases: `CPICNT`, `LSUCNT`, `FOLDCNT`, `EXCCNT`, and
`SLEEPCNT`, together with `CYCCNT`. Reset and sample one phase at a time because
the auxiliary counters are narrow and saturate. This distinguishes dependency
and load/store pressure from pure instruction count.

Run diagnostic-only cache experiments after the phase split:

- warm versus deliberately cold flash/libm calls;
- current RAM2 framebuffer versus a small synthetic DTCM sink;
- normal ISR load versus a foreground-only host/device micro-harness.

These experiments are attribution tools, not shippable configurations.

### 6. Disassembly ledger

For each candidate, archive the profile and full-roster ELF/map and run:

```text
arm-none-eabi-size -A firmware.elf
arm-none-eabi-nm -C -S --size-sort firmware.elf
arm-none-eabi-objdump -drwC -S firmware.elf
arm-none-eabi-readelf -S -s firmware.elf
```

Track, by named hot symbol:

- section and byte size;
- direct and indirect calls in the inner loops;
- stack-frame size and GPR/VFP saves;
- `vdiv`, `vsqrt`, `vldr`, `vstr`, and flash veneer calls;
- duplicated template instantiations;
- ITCM-to-DTCM granule position in the full roster.

Fail the candidate if the source-level hypothesis is absent in assembly. An
`inline`, `restrict`, `HS_O3`, or `constexpr` annotation is not a result until
the targeted call, spill, or instruction disappears.

## Architecture direction

The recommended end state is an effect-specialized regular-stroke batch, not a
second general Plot framework.

Conceptually:

```text
ShapeShifter::draw_all
  prepare frame: basis, clip, typed AA sink, shape/function/sides once
  build ring descriptors: radius, phase, shader parameters
  switch shape once
    draw_planar_regular_stack(shader, ...)    // polygon/star/flower
    draw_geodesic_regular_stack(shader, ...)  // spherical polygon
```

The batch owns only the semantics ShapeShifter actually uses:

- regular closed topology;
- one shape, side count, basis, and clip per frame;
- a typed compact stroke-sample shader interface carrying position, ring
  coordinate, edge coordinate, radius, and explicitly requested varyings;
- a compile-time constant-color specialization for the current presets,
  without making constant color the architectural boundary;
- no type-erased shader call or interpolation of unused generic Fragment
  registers;
- outer-to-inner ordered source-over rendering;
- a typed terminal anti-alias sink.

Keep the generic Plot path intact for other effects. This bounds code-size risk
and lets the ShapeShifter oracle compare old and new paths directly.

### Planar regular stack

PlanarPolygon, Star, and Flower ultimately draw straight segments in an
azimuthal-equidistant chart. The generic path currently:

1. creates sphere-space Fragment control points;
2. projects the points back into the chart;
3. performs a four-chord arc pre-pass per edge;
4. reconstructs an edge sampler;
5. simulates adaptive steps;
6. replays and re-unprojects the samples;
7. creates/interpolates full Fragments for a flat-color shader.

The specialized path should generate chart-space endpoints directly and stream
them into a typed procedural raster kernel. All edges of a regular polygon are
congruent in the rotationally symmetric chart; Star's alternating edges are
also congruent under reversal. Compute the arc table once per ring, not once per
edge, and prove the reuse against the oracle.

The first implementation must retain the current adaptive sample positions.
After exact/reassociated reuse is measured, a second candidate may replace the
arc-normalized replay with an error-bounded screen-space stepper.

### Geodesic regular stack

For SphericalPolygon, every edge in one ring has the same length:

```text
dot_edge = cos(radius)^2 + sin(radius)^2 * cos(2*pi/sides)
edge_angle = acos(clamp(dot_edge, -1, 1))
```

Compute it once per ring. Generate vertices by rotating one tangent direction
with a complex recurrence, then stream endpoints without `ArenaVector<Fragment>`
or the second `angle_between` pass. Reuse per-ring geodesic coefficients in the
edge sampler.

### Typed shader and sink

Evaluate `Filter::Screen::DirectAntiAliasSink<W,H>` first. It is already proven
pixel-identical to the generic terminal AntiAlias pipeline and is used by
MindSplatter. Prepare it once per frame. The regular-stroke kernel takes a typed
shader callable and a compact `RegularStrokeSample`; the shader returns Pixel
and alpha. A `ConstantStrokeShader` specialization may bypass sample fields used
only by dynamic shading while remaining one model of the same interface.

The sample interface must leave room for interesting ShapeShifter shading:

- normalized sphere position;
- normalized ring index and phase;
- normalized edge position and side index;
- radius and shape-local coordinates;
- optional explicitly declared varyings for gradients, noise, texture lookup,
  edge effects, and animated modulation.

Do not replace this with a flat-color-only API. The optimization target is
unused generic state and runtime erasure, not shader expressiveness.

This removes:

- `PipelineRef` indirection;
- `FragmentShaderFn` indirection;
- per-sample zeroing and interpolation of varyings the selected typed shader
  does not request;
- repeated clip decoding in the terminal chain.

The experiment must include the full-roster size delta. A typed instantiation
that wins only a few milliseconds while consuming most of the 11,592 B ITCM
headroom is not acceptable; the compact typed kernel plus its constant-shader
specialization should be smaller than instantiating the entire generic typed
planar rasterizer.

## Optimization sequence

Each numbered candidate is a separate measurable change. Do not stack an
unmeasured candidate under a later win.

### Phase A - invariant hoisting and topology

1. Move the shape switch outside the ring loop. Dispatch one shape-specific
   stack loop per frame rather than calling `dispatch_plot` 43-75 times.
2. Move the waveform switch outside the ring loop. Keep a compact generic
   fallback if duplicating four loops grows code excessively.
3. Precompute `1/count`, radius step, effective alpha, side step, and palette
   indexing increments once per frame.
4. Replace per-vertex `sinf`/`cosf` pairs with one initial sin/cos plus a complex
   rotation recurrence. Periodically renormalize only if the visual/error audit
   proves drift requires it.
5. Compare newlib `sincosf`, the existing Bhaskara fast trig, and recurrence.
   Exact libm is the reference; fast trig is an approximate candidate and must
   use the visual gate.
6. For SINE waveform, evaluate a bounded recurrence or lookup only after the
   geometry recurrence is measured. Triangle, saw, and square are already cheap.

Expected result: eliminate hundreds or thousands of flash libm calls per dense
frame and remove the repeated shape branch. Kill any sub-candidate that saves
less than 2% peak time while adding more than 1 KiB ITCM.

### Phase B - congruent-edge reuse

1. Compute SphericalPolygon edge length once per ring.
2. Compute one planar arc table per ring and reuse it across congruent edges.
3. Cache chart/basis setup once per frame and antipode setup once per ring.
4. Add a conservative whole-ring clip reject before allocating scratch or
   generating vertices.
5. Measure whether exact per-edge clip culling costs more than it saves for
   rings that broadly intersect a quadrant. Select between whole-ring-only and
   whole-ring-plus-edge culling from a cheap workload threshold.

Expected result: reduce edge setup and cull work by a factor close to Sides for
the dominant three- and six-sided frames. Preserve segment-level counters so a
reduction in tested work cannot hide undersampling.

### Phase C - procedural typed-shader raster

1. Add a small `RegularStroke`/`RegularStrokeBatch` kernel in `core/render/plot.h`
   or a focused adjacent header. It accepts procedural endpoints, edge kind, a
   typed shader callable, and a typed sink.
2. Stream control points; remove per-ring scratch allocation and
   `ArenaVector<Fragment>` construction.
3. Define a compact shader sample with semantic coordinates. Instantiate and
   interpolate only fields requested by the selected shader traits; remove the
   general `Fragment` register bank and runtime shader erasure from this path.
4. Preserve current sample order and AA tap order first, targeting exact or
   reassociated parity.
5. Add a constant-color shader specialization for the current presets, and a
   test shader using position/edge/ring coordinates so the fast path cannot
   accidentally become the only supported architecture.
6. Specialize only the two edge families: planar regular and geodesic regular.
   Avoid a template axis for Shape, Function, Sides, and Count simultaneously.

This is the main architectural lever. Target at least a 20% peak reduction and
no more than 6 KiB net full-roster ITCM growth. If a generic typed rasterizer
exceeds that budget, keep the compact typed kernel, one constant-shader
specialization, and route cold/error branches to shared noinline helpers.

Executed choice: keep `Plot::rasterize` as the shared semantic boundary and add
a compile-time `SinglePass` specialization instead of introducing a second
regular-stroke framework. ShapeShifter supplies a typed shader and a prepared
`DirectAntiAliasSink`; the generic replay path remains available to every other
effect. This delivered the architectural win while retaining the established
Fragment model for interesting shading.

### Phase D - eliminate simulation/replay duplication

The generic rasterizer simulates adaptive step lengths, scales them to land on
the endpoint, then recomputes positions during draw replay. Prototype these in
order:

1. Exact replay cache: cache any expensive frame/edge intermediates that are
   identical between simulation and replay; keep final sample positions exact.
2. Position cache: store final projected positions rather than generic step
   metadata where they can be derived without a second unproject/normalize.
3. One-pass truncated stepper: march in screen-error-bounded steps and emit a
   shorter final step instead of globally rescaling all steps.
4. Forward-difference/recurrence sampler: advance geodesic sin/cos and planar
   chart coordinates incrementally, re-anchoring at each control point.
5. Adaptive subdivision alternative: recursively or iteratively split until
   projected chord length and midpoint deviation meet the AA error bound.

Candidates 3-5 are approximate even if visually indistinguishable. Judge them
against the frozen sub-pixel oracle. Require a further 15% reduction from the
best Phase-C candidate or reject the added complexity.

### Phase E - stack-level batching and locality

Only after the single-ring kernel is efficient:

1. Build compact ring descriptors for the current frame in scratch memory.
2. Bin work into small screen tiles or clip bands while retaining outer-to-inner
   compositing order within every pixel.
3. Reuse projected basis/clip state across rings.
4. Evaluate per-pixel coverage accumulation for multiple hits from the same
   ring, then blend once. Do not reorder different rings.
5. Measure RAM2 framebuffer traffic with DWT `LSUCNT`; terminal blending is only
   3 ms in the baseline, so reject batching that merely moves blend work around.

This is a stretch lever, not the first architecture. It must show a measurable
projection/sampling locality win, not just fewer blend calls.

### Phase F - code-generation closeout

Once the algorithm and data path are stable:

1. Put only the final regular-stroke hot loop under `HS_O3_FN`/the existing
   selective-O3 mechanism.
2. Use `always_inline` only for tiny arithmetic leaves whose missed inline is
   visible in disassembly; use noinline cold helpers to reduce hot register
   pressure and code size.
3. Add `restrict`/alignment facts only where ownership proves them.
4. Test `-O3` with and without loop unrolling on the isolated kernel. Cortex-M7
   has no NEON; favor scheduling and fewer dependencies over code-size-heavy
   pseudo-vectorization.
5. Test a whole-program LTO build as an experiment, not as an assumed release
   change. It must pass every firmware environment and size gate.
6. Inspect branch layout and mark only measured dominant cases likely. The
   shipping presets exercise SphericalPolygon/Flower/Star, but PlanarPolygon and
   non-SINE functions remain required controls and may not become broken cold
   paths.
7. Re-run the full disassembly ledger and chase any remaining inner-loop flash
   veneer, indirect `blx`, redundant divide/square-root, or large stack spill.

The codegen phase succeeds only if the shipping build moves materially toward
the global-`O3` kernel while staying within the current FlexRAM bank. Do not land
global `-O3` for the full roster.

### Phase G - optional quality/performance envelope

Keep these behind explicit evidence and visual gates:

- raise `SCREEN_STEP_PX` from 0.9 in small increments;
- use Bhaskara fast trig for control-point generation;
- quantize phase/radius enough to share more geometry across adjacent rings;
- apply a pose-aware sampling density based on projected curvature;
- use a lower-cost AA kernel only where projected coverage proves equivalence.

Do not reduce Count or Sides in the shipping presets and call it optimization.
If a dynamic quality controller is ever considered, report it as a separate
product/visual decision after the full-quality renderer meets the 59 ms gate.

## Candidate ledger

Fill this table as work proceeds. “Peak” is render peak from normal captures;
deep-instrumented absolute times do not go here.

| Candidate | Class | Ship peak | Spherical | Flower | Star | Spilled | ITCM code | Visual gate | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Baseline `2c0a7a27` recapture | reference | 98.79 | — | — | — | 423/1296 | 185,016 | reference | keep as oracle |
| C1 direct typed AA sink | exact | 80.08 | — | — | — | — | 190,888 | pixel-identical | combine |
| D one-pass replay removal | exact sampling order | 86.17 | — | — | — | 391/1344 | — | oracle pass | combine |
| C2 + D + analytic planar sampler | exact/perceptual-gated | **43.52** | 43.52 | 41.30 | 37.54 | **0/1728** | **189,416** | exhaustive pass | **land `02ccaca8`** |
| F move all specialized hot code to flash | exact | 60.10 | — | — | — | — | 184,488 | unchanged | reject: misses target |
| F flash sampler + runtime mode | exact | 45.71 | — | — | — | — | 193,480 | unchanged | superseded |
| F compile-time mode + flash sampler | exact | 50.29 | — | — | — | — | 188,568 | unchanged | reject: slower |

For every row, link the source commit, raw logs, provenance artifacts, profile
reports, map diff, and framebuffer-diff output.

## Landing and rejection rules

- One hypothesis per candidate commit; no mystery bundles.
- Native tests and the ShapeShifter oracle run before device time.
- Build `phantasm` after every template/codegen candidate, not only at the end.
- A speedup smaller than 1 ms or 2% must repeat in an ABBA same-board run.
- Reject a candidate that adds more than 1 KiB ITCM for less than 2% peak gain.
- Reject a candidate that crosses the six-bank ceiling even if the
  single-effect profile fits.
- Reject a candidate that improves average time while leaving the 97 ms phase
  unchanged; cadence is set by the peak.
- Reject an approximate candidate whose quality budget was invented after the
  diff was seen.
- Remove count-only instrumentation after attribution. Permanent
  `HS_PROFILE` scopes may remain only when they compile to nothing outside the
  profile build.
- Preserve rejected experiments and measurements in this ledger so failed
  fused/batched designs are not repeated without a new hypothesis.

## Final verification

All completion gates passed at `02ccaca8`:

1. Shipping and global-`O3` 110-second captures each cover and wrap the full
   three-preset cycle, validate cleanly, and have current provenance.
2. Shipping peak is below the selected target, with zero spilled frames.
3. The exhaustive four-shape/four-function native matrix passes its declared
   exact or perceptual gate.
4. Production-resolution deterministic effect tests, the complete native suite,
   WASM build/smoke tests, `pio run -e phantasm`, and all Teensy size/layout gates
   pass.
5. The full-roster image stays in six ITCM banks with at least 4 KiB preferred
   headroom and no new arena/stack high-water failure.
6. Disassembly confirms no indirect shader/pipeline call or flash libm call in
   the final inner sample loop, and the documented instruction/spill reductions
   match the implemented hypothesis.
7. The two ShapeShifter profile reports and the three ranked profile READMEs are
   regenerated from the final raw captures.

The final disassembly of the 3,836-byte specialized raster body contains zero
indirect `blx`, zero `PipelineRef`/`FunctionRef` references, zero calls to the
slow planar sampler, and zero calls to `azimuthal_unproject`. It has two direct
static calls to the 568-byte analytic `one_pass` sampler and direct calls to the
typed AA sink. Count-only profiling confirms replay work fell from 48.0 to 0
steps/edge and unprojects fell from 159.8 to 15.9/edge while plotted samples
changed by only 0.2% and the backstop remained zero.

## First execution slice

Start with a bounded, high-information slice:

1. Add the dual-render framebuffer oracle and frozen profile-case controls.
2. Recapture normal and deep shipping/global-`O3` baselines on one board.
3. Add count-only Plot metrics and identify cycles per edge/substep/sample.
4. Prototype switch hoisting plus trig recurrence; keep or reject independently.
5. Prototype `DirectAntiAliasSink` with a typed full-Fragment shader call path
   and measure both speed and full-roster ITCM.
6. Use those results to size the procedural planar/geodesic batch before writing
   it.

This slice resolves the largest uncertainties—where the 90 ms Plot cost splits,
how much flash trig matters, and how expensive a typed sink is in the full
roster—before committing to the deeper raster rewrite.
