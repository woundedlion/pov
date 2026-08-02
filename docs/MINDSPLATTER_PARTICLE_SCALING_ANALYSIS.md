# MindSplatter 2,048-Particle Scaling Analysis

**Status: ANALYSIS (2026-07-22).** 2,048 full particles cannot be reached by
raising `NUM_PARTICLES` — it fails on both the arena and the CPU.

Date: 2026-07-22  
Target: Teensy 4.0, i.MX RT1062 Cortex-M7 at 600 MHz  
Render target: ideally no more than 59 ms peak, with zero spilled frames

## Executive conclusion

The current architecture cannot support 2,048 full 23-sample particles under
59 ms by merely increasing `NUM_PARTICLES`. It fails on both memory and CPU:

1. The pool would require 368,640 bytes, exceeding the 272,384-byte persistent
   arena.
2. A workload-matched doubling projects to roughly 94 ms average and 110-120 ms
   peak render.

The most credible route is:

> 2,048 particles, roughly 10 retained trail samples each, plus a
> MindSplatter-specific streaming one-dot renderer.

That combination fits RAM and has a credible 53-58 ms peak target. Keeping 23
samples on every particle probably requires a more radical rendering-quality
tradeoff.

## Current baseline

The current four-preset raw capture, `build/prof/mindsplatter_ship.log`, is
newer than the five-preset July 20 shipping report. It validates a complete
four-preset cycle:

- 1,408 frames
- 58.27 ms peak
- 0 spills
- Only 0.73 ms below the requested 59 ms limit

The worst window, frames 1249-1264, averaged 47.25 ms render but contained the
58.27 ms peak:

| Scope | Average/frame | Calls/frame |
|---|---:|---:|
| Particle draw | 44.91 ms | 1 |
| Raster | 26.05 ms | 326.2 trails |
| Gate | 9.24 ms | 827.6 trails |
| Vertex transform | 4.36 ms | 827.6 trails |
| Trail decode/materialization | 3.15 ms | 827.6 trails |
| Physics | 2.28 ms | 1 |
| Deferred hole shading | 1.59 ms | 326.2 trails |
| Timeline | 0.06 ms | 1 |

Only 39.4% of front-end trails reach rasterization. About 60.6% pay decoding,
Mobius/orientation transform, and gating before being rejected.

Normalized costs from the cycle counters are approximately:

- Trail decode/materialization: 2,286 cycles/trail
- Vertex transform: 3,158 cycles/trail
- Gate: 6,699 cycles/trail
- Deferred shading: 2,918 cycles/surviving trail
- Raster: 47,921 cycles/surviving trail, roughly 2,180 cycles per edge at 23
  samples

The ISR consumed approximately 5% CPU in this window. For robust sub-59-ms
operation, the effect should target 55-56 ms peak rather than 58.9 ms.

### Baseline provenance caveat

The detailed July 20 report records commit `5afc18b7`, five presets, 1,728
frames, and a 59.88 ms peak. The newer raw log reflects the current four-preset
table and validates a full wrap, but its serial header does not embed a commit
hash. Use the newer log for workload estimates and capture a fresh explicitly
versioned baseline before landing an optimization.

## Capacity is presently RAM-blocked

The current definition is 1,088 particles with 23 samples. Each particle is
180 bytes:

```text
position + velocity + seed/life       28 B
23 x snorm16 xyz                     138 B
alignment                              2 B
ring head/tail/count                  12 B
                                      ----
                                      180 B
```

At 2,048 particles:

```text
2048 x 180 B = 368,640 B
persistent arena   = 272,384 B
```

That is impossible before considering the palette and other persistent
allocations.

Useful configurations are:

| Representation | Particle | 2,048-particle pool | Fits? |
|---|---:|---:|:---:|
| Current, 23 samples | 180 B | 360 KiB | No |
| Current snorm16, 12 samples | 112 B | 224 KiB | Yes |
| Current snorm16, 10 samples | 100 B | 200 KiB | Yes |
| Packed 32-bit sphere point, 23 samples | 132 B | 264 KiB | No; only 2 KiB remains for all other persistent state |

A packed 32-bit trail solves memory but not CPU. It may also add 0.5-1.5 ms of
decoding work unless reconstruction is particularly cheap.

### Capacity is not the same as resident population

The current eight emitters spawn eight particles per frame and `max_life` is
160. Including 22 draining-history frames, the theoretical resident ceiling is
approximately:

```text
8 x (159 live + 22 draining) = 1,448 particles
```

Raising capacity to 2,048 will therefore not produce 2,048 residents without
increasing emission rate or lifetime. The final benchmark must force and
sustain the intended population.

## Disassembly findings

The disassembly was taken from
`build/prof/phantasm_setup3_dead_normals_final.elf`. The MindSplatter hot code
is unchanged by the later preset-only commits.

Static opcode counts are useful as register-pressure indicators. Stack
references include explicit local objects as well as genuine spills.

| Symbol | Code | Local/total stack | `[sp]` references | FP/int divides | sqrt | Indirect calls |
|---|---:|---:|---:|---:|---:|---:|
| Particle draw front end | 2,736 B | 220 / 328 B | 122 | 6 / 2 | 4 | 1 |
| Trail tween | 260 B | 28 / 64 B | 7 | 1 / 0 | 0 | 2 sites |
| Particle physics | 1,070 B | 52 / 144 B | 8 | 8 / 3 | 3 | 0 |
| Direct rasterizer | 3,264 B | 332 / 440 B | 241 | 10 / 1 | 4 | 11 |
| Direct AA sink | 2,388 B | 24 / 72 B | 17 | 0 | 0 | 0 |

### Trail traversal

Current trail traversal no longer divides or performs `% 23` per sample. It
computes one reciprocal per trail and uses conditional wrapping:

```asm
25370: vdiv.f32 s17, s14, s15       ; 1 / (length - 1), once per trail
...
25386: ldrsh    r3, [r4, r3]        ; decode x
25388: ldrsh.w  r1, [r2, #2]        ; decode y
25390: ldrsh.w  r2, [r2, #4]        ; decode z
25398: vcvt.f32.s32 s13, s13
253aa: vmul.f32 s13, s13, s16
253b6: vstr     s13, [sp, #12]      ; materialize callback Vector
253d2: blx      r6                   ; indirect callback
253d4: cmp      r5, #23
253d8: moveq    r5, #0
```

The remaining cost is the decoded temporary, stack stores, callback ABI, and
construction of the full downstream representation.

### One-dot raster path

The one-dot path still copies complete 48-byte `Fragment` objects onto its
large stack frame before invoking the shader:

```asm
153e4: ldmia.w ip!, {r0-r3}
153e8: stmia.w r5!, {r0-r3}
153ea: ldmia.w ip!, {r0-r3}
153ee: stmia.w r5!, {r0-r3}
153f0: ldmia.w ip,  {r0-r3}
153f4: stmia.w r5,  {r0-r3}         ; 48-byte Fragment copy
...
15408: blx      r9                   ; fragment shader callback
...
153d4: bl       DirectAntiAliasSink::plot
```

This is the main remaining compiler-visible structural problem: the dominant
one-dot path still travels through a general polyline rasterizer, general
`Fragment`, callback ABI, and a 440-byte stack frame.

### Particle physics

The physics disassembly contains dependent `VSQRT` to `VDIV` chains:

```asm
23568: vsqrt.f32 s14, s15
2357a: vdiv.f32  s15, s25, s14
23582: vdiv.f32  s15, s23, s14
2358a: vdiv.f32  s15, s22, s14
```

These chains are difficult to hide on an in-order processor. The Cortex-M7
can dual-issue many independent instructions and has dual TCM load channels,
but dependencies, callback boundaries, and load/store pressure restrict that
advantage.

## Microarchitecture and cache analysis

The existing capture contains `CYCCNT` timing only. It does not contain a
hardware D-cache-miss count, so an exact miss number cannot be derived from it.

What can be established:

- The hot functions are in `.text.itcm`; they do not suffer ordinary I-cache
  misses.
- Particle state, trails, scratch allocations, palette, and stack are in DTCM.
  A particle SoA conversion is therefore not a D-cache optimization.
- The framebuffer is in OCRAM and is the meaningful D-cache-pressure source.
- The RT1062 has 32 KiB instruction and data caches.

A framebuffer is 288 x 144 x 6 = 248,832 bytes. The active
segment-plus-margin is approximately 64 KiB, about twice the D-cache.

With six-byte pixels:

- 12.5% of individual pixels straddle a 32-byte line.
- A horizontal two-pixel AA pair crosses a line for 31.25% of alignments.
- A 2x2 AA splat touches 2.625 cache lines on average.
- Randomly interleaved trails reduce reuse before dirty lines are evicted.

The exact follow-up measurement should record `DWT_CPICNT`, `DWT_LSUCNT`, and
`DWT_EXCCNT` alongside `CYCCNT` around the compact front end and AA sink. These
are useful stall proxies even though none is a dedicated cache-miss event.

### ITCM constraint

The analyzed full-roster ELF has 189,696 bytes of ITCM code, only 6,912 bytes
below the 196,608-byte bank boundary. Crossing that boundary allocates another
32 KiB FlexRAM bank to ITCM and removes it from DTCM. With 314,176 bytes of
RAM1 variables, that would make the image fail to fit.

Any specialized kernel must replace or outline existing hot code rather than
simply adding another broad specialization.

## Prioritized recommendations

Estimates below are for the doubled target. Overlapping recommendations should
not be summed blindly.

| Rank | Recommendation | Estimated gain | Confidence |
|---:|---|---:|:---:|
| 1 | Cap effective trails at about 10 samples at 2,048 particles | 40-50 ms versus naive 23-sample doubling | High |
| 2 | Replace the generic particle draw/raster path with a streaming one-dot kernel | 3-6 ms with the 10-sample edge budget; 6-12 ms with full trails | Medium-high |
| 3 | Add a true source-space pre-transform reject; do not use the tested compact-materialization split | Unmeasured; the compact split saved only 0.08-0.20 ms at 1,088 | Low |
| 4 | Pair `+axis/-axis` physics and reuse perpendicular magnitude | 1.1-1.8 ms at 2,048 | Medium-high |
| 5 | Improve framebuffer locality with order-preserving tile command bins | 2-5 ms with short trails; 3-8 ms full workload | Medium-low |
| 6 | Use packed 32-bit trail points if 23-sample storage is mandatory | Enables RAM fit; likely costs 0.5-1.5 ms | High for memory |
| 7 | Dense-mode reduced or stochastic AA if full trails remain mandatory | 12-30 ms | Low visual confidence |

### 1. Use a total trail-work budget

At 2,048 particles and 10 samples:

```text
current edges:  1088 x 22 = 23,936
target edges:   2048 x  9 = 18,432
```

The target has 23% fewer maximum edges than today while nearly doubling particle
heads. The current representation becomes 100 bytes per particle, so the pool
occupies only 200 KiB.

Using the measured per-trail costs, the predicted saturated target is:

- 48-52 ms average before additional optimization
- 59-64 ms peak given current peak/average variation
- 53-58 ms peak after the specialized renderer below

This is the highest-confidence route to doubling within the existing visual
algorithm.

A perceptual alternative is "all heads, half tails":

- Store state for all 2,048 particles.
- Give only a deterministic 1,024-particle subset full 23-sample trails.
- Render every particle's head.

That occupies approximately 208 KiB and preserves the present tail density
while doubling visible moving points.

### 2. Build a true MindSplatter streaming kernel

Do not repeat the rejected fully inlined `Fragment` implementation. That
generated a 48-byte `memset` per sample and raised the prior peak to 70.18 ms.

Instead, use compact SoA scratch containing only:

```text
decoded/transformed position
screen row/column
hole alpha
edge classification byte
```

Then:

1. Decode directly into the compact arrays.
2. Apply Mobius and the precomputed rotation matrix.
3. Gate the trail.
4. For one-dot edges, compute `v0`, seed phase, palette, alpha, and AA directly.
5. Send only the roughly 9% long-edge cases to the existing general
   rasterizer.

This removes the 48-byte fragment copies, indirect shader call, repeated color
zeroing, general raster state, and much of the 440-byte spill-heavy frame from
the dominant path.

### 3. Move rejection before expensive representation

The theoretical ceiling is the 60.6% of trails currently discarded after
front-end processing.

A Mobius transform maps circles to circles, so segment clip boundaries can in
principle be transformed back into source-space circle/plane tests once per
frame. A conservative source-space scan can then reject trails before fragment
construction and full vertex transformation.

Even a less ambitious fused decode/gate pass should avoid building the 48-byte
fragment and 12-byte original-position arrays for rejected trails.

#### Compact-materialization experiment (2026-07-22): rejected

An exact prototype stored only `{original position, transformed position}`
before the existing gate, then constructed full fragments only for surviving
trails. Native pixel-parity tests passed, and both Teensy captures completed a
four-preset wrap on the same COM3 board.

| Preset | Shipping baseline | Shipping compact | Delta |
|---:|---:|---:|---:|
| 2 | 49.76 ms | 49.56 ms | -0.20 ms |
| 3 | 44.90 ms | 44.92 ms | +0.02 ms |
| 0 | 41.40 ms | 41.22 ms | -0.18 ms |
| 1 | 34.45 ms | 34.37 ms | -0.08 ms |

These are the worst clean-hold `msp_draw_particles` windows per preset. Peak
shipping render changed from 71.70 to 71.52 ms, while spilled frames remained
4/1,728. In the worst window, compact decode reduced `plot_ps_tween` by 1.02
ms/frame, but survivor materialization cost 0.43 ms/frame and the gate became
0.33 ms/frame slower. Raster and deferred-shader variation consumed most of
the remainder.

Global `-O3` improved the four clean-hold costs by 1.78-2.19 ms, but that does
not help the selective-O3 shipping image and still leaves the worst preset far
above the required headroom.

The full-roster ELF grew by 2,272 bytes of ITCM, from 189,784 to 192,056 bytes,
leaving only 4,552 bytes before the next FlexRAM bank boundary. Disassembly
shows why: the 2,736-byte integrated baseline `draw_impl` became a 1,288-byte
driver plus two 1,524-byte gate specializations. Nested stack demand also rose
from 328 bytes to 496 bytes. The measured shipping gain does not justify either
cost, so this implementation was not landed.

The remaining version of recommendation 3 must reject in source space before
the Mobius transform itself. Merely changing the intermediate representation
cannot recover the original 4-7 ms estimate.

### 4. Process signed-axis attractors in pairs

For each axis pair, reuse:

```text
cross_sq(+X) == cross_sq(-X)
dist^2(+X) = |p|^2 + 1 - 2px
dist^2(-X) = |p|^2 + 1 + 2px
```

This permits three perpendicular-magnitude square roots rather than six
dynamically, removes runtime axis-selection branching, and shares several
products. Physics becomes about 4.6 ms after doubling, making a 25-40%
reduction worthwhile.

### 5. Address framebuffer locality after fusion

Order-preserving tile bins are valid: append commands in original draw order,
then process one tile at a time. Operations targeting the same pixel retain
their original source-over order; operations to unrelated pixels may be
reordered safely.

The tradeoff is command-buffer RAM. This becomes practical only if
segment-local framebuffer storage or packed history frees OCRAM.

## Low-value or unsuitable levers

- **Global `-O3`:** the measured peak improved only 1.20 ms while RAM1 code
  grew by 28,864 bytes. It is not roster-safe.
- **Broad new hot specializations:** only 6,912 bytes remain before the ITCM
  bank boundary.
- **Further clear optimization:** clear is already approximately 0.10 ms and
  is outside the reported render peak.
- **Moving the existing pool to OCRAM:** RAM2 has only 4,736 bytes free and
  trail reads would compete with framebuffer lines in D-cache.
- **Physics-only optimization:** even eliminating physics would not cover the
  doubling deficit.
- **A simple particle SoA rewrite for cache:** current particle data is already
  in DTCM.

## Recommended landing sequence

1. Add a deterministic saturated-population benchmark that actually holds
   2,048 residents.
2. Try `TRAIL_LEN=10`, preserving snorm16 and all current math.
3. Capture all four presets; reject unless peak is already near 60-62 ms.
4. Land the compact streaming one-dot renderer.
5. Pair the signed-axis physics.
6. Add pre-transform/source-space rejection if more headroom is required.
7. Require 55-56 ms measured peak, zero spills, and full-roster ITCM below
   196,608 bytes.

## Measurement and correctness protocol

For every candidate:

1. Replay identical particle state after one shared simulation step. Ordinary
   independent captures can diverge through floating-point trajectories.
2. Record rendered trails, rejected trails, one-dot edges, long edges, exact
   gate fallbacks, raster samples, and accepted AA taps.
3. Record `CYCCNT`, `CPICNT`, `LSUCNT`, and `EXCCNT` for front end, raster, and
   AA sink. Use sufficiently short batches to avoid auxiliary-counter
   overflow.
4. Compare deterministic framebuffer hashes. If arithmetic reordering prevents
   bit identity, report changed pixels, coverage changes, maximum channel
   delta, and total channel delta.
5. Run a full four-preset wrap with the population held at the target, not just
   a pool whose capacity happens to be 2,048.
6. Build the complete Phantasm roster and reject any landing that crosses the
   196,608-byte ITCM boundary or leaves unsafe DTCM headroom.

## Final assessment

"2,048 particles" is realistic. "2,048 particles each with the present
23-sample, four-tap-AA trail" probably is not realistic under a 59 ms peak on
the current 600 MHz renderer.

The 10-sample or all-heads/half-tails design spends the fixed frame budget on
particle count instead of duplicating historical raster work. A compact,
MindSplatter-specific one-dot renderer is then the most credible way to regain
enough margin to remain solidly green.

## Evidence

- `build/prof/mindsplatter_ship.log`
- `docs/profiles/shipping/profile_mindsplatter_teensy_2026-07-20.md`
- `docs/profiles/O3/profile_mindsplatter_teensy_2026-07-20.md`
- `effects/MindSplatter.h`
- `core/animation/sprites.h`
- `core/animation/trails.h`
- `core/render/plot.h`
- `core/render/filter.h`
- `core/engine/memory.h`
- `build/prof/phantasm_setup3_dead_normals_final.elf`
- Arm Cortex-M7 Processor Technical Reference Manual, ARM DDI 0489D
- NXP i.MX RT1060 Crossover Processors for Consumer Products datasheet
