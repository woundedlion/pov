# ShaderBall red-preset optimization plan

Status: **complete on the current 23-preset bank as of 2026-08-12.** The final
shipping cycle has 23 green buckets, zero spilled frames, and a 49.52 ms worst
render. Earlier 29- and 22-preset investigations remain below as historical
evidence; the final closure ledger supersedes their open holdout status.

The target is a **peak render time below 59.00 ms** for every authored preset
in the shipping selective-O3 image, followed by a full-cycle capture with no
spilled frames. The 59 ms threshold reserves the measured column-ISR/DMA share
inside the 62.5 ms display window; it is stricter than merely avoiding the
next 125 ms cadence tier.

## Measurement basis

The existing shipping cycle at source `0f50fe97` marked presets 0, 1, 14, 19,
20, 26, 27, and 28 red. A fixed-preset capture was then taken for each of those
indices at source `080aa9b9`; the commits differ only in profile documents, not
runtime source.

Each fixed run used the shipping `profile` environment on COM3:

```sh
HS_PROFILE_OUT=build/prof/shaderball_ship_pNN.log \
  bash tools/profile_one.sh ShaderBall profile 70 32 \
  "-D HS_PROFILE_PRESET=NN"
```

All eight captures pass `tools/parse_profile.py ... validate`: monotonically
increasing frames, the requested `Profile preset:` marker, real render
telemetry, a complete counter tree, and root-cycle/wall agreement. The
attested full Phantasm image also passes the shipping size gate.

| Preset | Fixed regime | Shader window range, ms/f | Peak render, ms | Spilled | Speedup to 59 ms |
|---:|---|--:|--:|--:|--:|
| 26 | Peirce primitive lattice | 344.01-370.63 | **391.55** | 160/160 | **6.64x** |
| 19 | diagnostic Peirce | 160.98-161.22 | **166.78** | 352/352 | **2.83x** |
| 27 | kaleidoscope edge-fade liquid | 57.76-67.51 | **76.22** | 502/544 | **1.29x** |
| 28 | dodecahedral grid + mirror/noise | 59.96-68.28 | **73.92** | 544/544 | **1.25x** |
| 0 | kaleidoscope liquid stereo | 51.69-60.18 | **65.45** | 148/928 | **1.11x** |
| 14 | liquid stereo 14 | 45.61-53.89 | **59.56** | 0/1088 | **1.01x** |
| 20 | diagnostic Airocean | 41.35-41.40 | 47.55 | 0/1088 | already below |
| 1 | liquid stereo 1 | 37.47-40.38 | 45.02 | 0/1088 | already below |

Raw captures are `build/prof/shaderball_ship_p00.log`, `p01`, `p14`, `p19`,
`p20`, `p26`, `p27`, and `p28`, with the matching `.provenance`, attested ELF,
map, build log, and environment dump beside each log.

### What “red preset” means here

The cycle parser assigns a transition frame to the preset whose marker owns
that frame. ShaderBall's through-clear transition renders only one endpoint on
each non-clear frame, so a cheap incoming preset can inherit the outgoing
endpoint's cost.

- Presets 1 and 20 are transition-attribution false positives. Fixed in place,
  both have more than 11 ms of margin.
- Preset 14 is green by cadence but misses the requested engineering target by
  0.56 ms.
- Presets 0, 19, 26, 27, and 28 are intrinsically over budget and form the
  implementation set.

The fixed runs are deliberately longer than a cycle hold. They expose
orientation, palette, and wrapped-noise phases that a two-window compressed
hold can miss; optimization acceptance must use the fixed peak, not only the
clean-hold average from the roster sweep.

## Generated ARM code

Analysis used the GCC 15.2.1 tools that built the attested profile ELF:

```text
C:/Users/gabe/.platformio/packages/
  toolchain-gccarmnoneeabi-teensy/bin/arm-none-eabi-{nm,objdump,size}.exe
```

The representative shipping ELF is the preset-26 artifact named in
`build/prof/shaderball_ship_p26.provenance`. The rendered call chain is:

```text
Scan::Shader::draw
  FrameShader::operator()
    apply_lens
    project_branch
      project_nonstereographic
        project_peirce
          projections::peirce_projection
    shade_projected
      planar_warp_lookup
      source-coordinate conditioning
      source sampling
      material shaping
      colorize
```

Every arrow below `Scan::Shader::draw` is taken per pixel. A quadrant contains
10,368 pixels, so even a 100-cycle leaf costs about 1.73 ms/frame.

### Peirce kernel

`projections::peirce_projection` is flash-resident, 0x478 bytes, and 340 ARM
instructions in the shipping ELF. A square-layout, edge-distance-enabled call
executes these out-of-line leaves per pixel:

```text
atan2f, fmodf, sinf, cosf,
asinf x2, peirce_elliptic_integral x2,
acosf, and a conditional third asinf
```

That is 10368 executions of the projection and roughly 93k-104k
transcendental/helper calls per frame before source and color work. The two
Peirce presets both require edge distance, so the existing
`calculate_edge_distance=false` shortcut does not apply to them.

The diagnostic controls isolate this cost unusually well: presets 19 and 20
share the diagnostic pipeline and differ principally in Peirce versus Airocean
projection. Their fixed peaks are 166.78 and 47.55 ms, a 119.23 ms delta.

### Shipping versus global-O3 topology

The same runtime source produces materially different hot-code topology in the
shipping and global-O3 profile ELFs:

| Symbol | Shipping placement/bytes | Global-O3 placement/bytes |
|---|---:|---:|
| `shade_projected` | flash / 980 | flash / 212 |
| `sample_source` | inlined into caller | ITCM / 888 |
| `FrameShader::operator()` | flash / 310 | flash / 314 |
| `peirce_projection` | flash / 1,144 | flash / 1,152 |
| `colorize` | ITCM / 776 | ITCM / 1,728 |
| `Scan::Shader::draw` | ITCM / 516 | ITCM / 672 |

The important difference is not Peirce's own byte count. Global O3 outlines
the large source switch into an ITCM leaf and leaves a much smaller
`shade_projected` dispatcher. Shipping inlines the switch into a branch-heavy
flash block. The global-O3 cycle lowers diagnostic Peirce from 142.74 to
56.02 ms/f and Peirce lattice from 364.76 to 92.09 ms/f, which makes closure
layout and instruction-fetch behavior a first-class hypothesis.

Global O3 is not directly shippable: its ShaderBall profile consumes 25,472
more ITCM bytes. The attested full shipping roster currently uses 196,232 of
196,608 bytes, leaving only **376 bytes**. Any ITCM promotion must be paired
with a measured reclaim or stay within that limit.

### Lens and noise controls

Presets 0 and 14 have identical authored parameters. Their principal topology
difference is the surface lens: preset 0 uses the six-sector kaleidoscope and
preset 14 uses the liquid-stereo glitch lens. Their fixed peaks differ by
5.89 ms, pricing the kaleidoscope path relative to the glitch path.

The current six-sector fold computes `sqrt`, `fast_atan2`, `fmodf`,
`fast_cosf`, and `fast_sinf` for every pixel. Its disassembly is 0x194 bytes and
contains the inlined fast-trig polynomials plus a call to `fmodf`.

Preset 28 has a different risk: a dodecahedral reflection loop followed by a
mirror-tile stage and legacy stereo noise. Its 59.96-68.28 ms shader range is
large enough that the wrapped-noise double-sample interval and reflection
count must be measured separately before changing either kernel.

## Implementation plan

Each step is an independent A/B. Keep a change only if it improves the fixed
peak, preserves the framebuffer contract, and passes the full-roster size
gate.

### 1. Add low-distortion stage accounting

Add a ShaderBall-only deep-profile mode with accumulated DWT buckets for lens,
projection, planar warp, source, material, and color. Do not use an
`HS_PROFILE` registry scope around each pixel: its registry and hierarchy cost
would distort the very leaves being compared. Read CYCCNT at stage boundaries,
accumulate plain counters, and print/reset them once per profile window.

Capture fixed presets 0, 14, 19, 20, 26, 27, and 28. The controls make the
stage totals auditable: 0 versus 14 prices the kaleidoscope; 19 versus 20 prices
Peirce; 19 versus 26 prices the remaining lattice/lens/color stack.

### 2. Reproduce the small global-O3 dispatcher without a global promotion

Run two code-layout experiments:

1. Force `sample_source` out of `shade_projected` as a noinline, no-clone,
   selective-O3 flash leaf. This should shrink the caller without spending
   scarce ITCM.
2. If flash outlining is slower, place only that 888-byte leaf in ITCM after
   reclaiming at least one kilobyte elsewhere. Do not rely on the current
   376-byte margin.

Inspect the resulting ELF after each build. The intended shape is a compact
`shade_projected`, one source-dispatch symbol, no IPA clone that silently lands
in ITCM, and no duplication of inactive source arms. Profile all fixed presets,
because this layout affects every ShaderBall path.

This is the lowest-risk route to bring preset 19 near its global-O3 result and
should precede mathematical approximation work.

### 3. Replace the six-sector lens with an exact reflection fold

Map `(x,z)` into the 30-degree dihedral chamber with sign normalization,
constant 60-degree rotations/reflections, comparisons, and multiply-adds.
Reflection preserves radius, so the replacement needs no square root, inverse
trig, remainder, or forward trig.

Validate the new kernel against the old lens over a dense sphere sample and at
all sector boundaries, then run framebuffer comparisons for presets 0, 22,
23, 26, and 27. Presets 0 and 14 provide the on-device timing control. This
step is expected to recover most of preset 0's 5.89 ms lens delta and also
helps presets 26 and 27.

### 4. Build a measured Peirce fast path

Work from the DWT split rather than replacing the whole projection at once:

- Replace `atan2f` and the sine/cosine pair with the engine's validated fast
  angle kernels or a joint fast-sincos implementation.
- Express the two `asinf(sqrt(...))` terms with one measured approximation
  family; keep the existing Chebyshev elliptic recurrence, which is already a
  short multiply-add chain.
- For edge fade, avoid computing an exact angle when the consumer only applies
  `smooth_ramp(0, edge_width, distance)`. Compare a sine/cosine-domain fade to
  the current angular fade and retain exact seam classification separately.
- Hoist square-layout and `edge_distance_required=true` decisions out of the
  generic per-pixel switch for the two authored Peirce configurations. Do not
  duplicate the entire shader closure.

The current projection tests use tight cartographic oracles and seam checks.
Keep the exact kernel as the reference and add explicit maximum coordinate,
fade, and seam-error budgets for a fast renderer path; do not simply loosen
the existing tests. If arithmetic approximation cannot put preset 26 below 59
ms, prototype a seam-aware flash LUT as the fallback. A LUT must report its
flash cost, use no persistent RAM, and retain exact region/edge metadata.

### 5. Optimize the preset-28 reflection/noise path

Use the deep counters to separate:

- dodecahedral chamber reflections, including maximum and average reflection
  count;
- mirror-tile coordinate work;
- one-sample versus wrapped two-sample legacy noise;
- source/material/color cost.

Then apply only the measured lever. Candidate changes are an unrolled
three-mirror chamber fold with a bounded observed iteration count, hoisted
frame-constant noise-wrap selection, and a paired wrapped-noise kernel that
shares coordinate preparation. Preserve the periodic time seam and the exact
mirror chamber; changing either is a visual change, not a performance cleanup.

### 6. Close the transition and roster gates

After all fixed peaks are below 59 ms:

1. Re-run every affected fixed preset for 70 s with 32-frame windows.
2. Run the complete 29-preset shipping cycle with 16-frame windows and confirm
   the 28-to-0 wrap, zero spills, and no transition-owned red buckets.
3. Capture the global-O3 twin as a regression ceiling.
4. Build `pio run -e phantasm`; require the full roster to remain below the
   196,608-byte ITCM ceiling.
5. Run the native suite and ShaderBall framebuffer/oracle tests. Any approved
   fast-math path also needs dense host error scans and seam-boundary cases.
6. Update the two per-config profile reports and all three ranked profile
   READMEs from the new logs.

## Acceptance checklist

- Every fixed preset in `{0, 14, 19, 26, 27, 28}` has peak render `< 59.00 ms`
  and `0/N` spills in a 70 s shipping capture.
- Presets 1 and 20 do not regress above 59 ms.
- The complete 29-preset shipping cycle has no red or yellow bucket and wraps
  to preset 0.
- Full-roster Phantasm stays within the existing ITCM/DTCM bank allocation;
  no result is accepted from a single-effect image alone.
- Exact paths remain available as test oracles for any approximate lens or
  projection renderer.
- Framebuffer, projection seam, topology, and native tests pass.
- ELF symbol placement and byte deltas are recorded beside each device A/B.

## Implementation ledger

All timing entries are shipping-image device measurements on COM3. Positive
gain means a lower peak render time. ITCM deltas use the full Phantasm roster,
not the single-effect profile image. Instrumentation images are excluded from
timing acceptance.

### Accepted changes

| Change | Evidence | Frame gain | ITCM cost / benefit |
|---|---|---:|---:|
| Raw ShaderBall stage/detail counters | DWT buckets for lens, projection, warp, source, material, color; preset-28 detail adds mirror, legacy-noise, reflection, and one/two-sample counts | diagnostic only | 0 B shipping; compiled out |
| Flash outline of `sample_source` | `shade_projected` shrank 980→224 B; one 772 B flash source symbol, no IPA clone | preset 19: **+6.22 ms** before Peirce math | 0 B promotion |
| Exact six-sector reflection fold | Dense 129×1440 oracle and every sector boundary pass; no sqrt/inverse trig/remainder/forward trig | preset 0 early A/B: **+4.43 ms** | no retained ITCM promotion |
| Square zero-meridian Peirce renderer | 129×512 exact-reference scan, explicit seam offsets, coordinate error <0.0012, fade error <0.0002, exact metadata | preset 19: **+116.916 ms** final; preset 26: **+336.987 ms** final | flash kernels: fast 974 B, exact 990 B |
| LUT gamut clip and parabolic joint hue sin/cos | Dense 47,385-color/angle scan: max channel error 5,516/65,535, mean ≤256; exact path remains the oracle | liquid stack contributes preset 14 **+8.336 ms**, preset 27 **+16.949 ms** with the lens/layout changes | included in final net **-80 B** ITCM |
| One-step reciprocal sqrt and scale-invariant LUT hue bin | Same dense gamut oracle passes | preset 27 short A/B: about **+0.10 ms** | 0 B measured delta |
| Selective O3 polyhedral fold | Preset-28 reflection count averages 3.2–3.8, observed max 9; exact chamber unchanged | preset 28 short A/B: **+1.33 ms** | about +96 B, retained inside net saving |
| Frame-prepared legacy-noise wrap phase | Exact bitwise oracle across the blend boundary; both noise samples and the periodic seam remain unchanged | preset 28 full-phase: **+1.158 ms** | retained inside final net saving |

The final full-roster image uses **196,152/196,608 ITCM bytes**, versus the
196,232-byte baseline: **80 bytes saved**, with 456 bytes of headroom. After
rebasing onto two concurrently added presets, FLASH code is 400,488 bytes.
`pio run -e phantasm` passes every region gate.

### Rejected experiments

| Experiment | Decision evidence | Frame gain / loss | ITCM / memory effect |
|---|---|---:|---:|
| Conditional liquid `wrap_t` | Difference was measurement noise | +0.01 ms | **+112 B** ITCM; rejected |
| ITCM `warp_stage_lookup` | Preset 28 slowed 66.13→66.60 ms | **-0.47 ms** | roughly **+464 B**, only 48 B headroom; rejected |
| ITCM gamut clip helper | Preset 28 slowed 66.13→66.69 ms | **-0.56 ms** | roughly **+352 B**; rejected |
| Inline gamut clip | Preset 28 slowed 66.13→67.35 ms | **-1.22 ms** | roughly **+144 B**; rejected |
| Paired XY-offset noise preparation | 66.13→66.09 ms, below run noise | +0.04 ms | no useful roster benefit; rejected |
| RAM2 gradient-table cache | Full-roster gate failed | not accepted | **+1,024 B RAM2**, leaving 3,200 B below the 4,096 B floor; rejected |
| Unclipped hue rotation | Max channel error 51,807 and mean 1,740 | not profiled further | visual contract failure; rejected |
| RGB-line gamut clip | Max channel error 16,622 | not accepted | visual contract failure; rejected |
| Global O3 shipping topology | Twin cycle has 96/2272 spills versus shipping 66/2304 and makes preset 20 red | **30 more spills** | global promotion remains unsuitable; rejected |

### Final fixed-preset ledger

Each final row is a validated 70 s shipping capture with 32-frame windows.

| Preset | Baseline peak, ms | Final peak, ms | Gain/loss, ms | Result |
|---:|---:|---:|---:|---|
| 0 | 65.45 | **53.680** | **+11.770** | accepted, pass |
| 1 | 45.02 | **44.985** | +0.035 | control pass |
| 14 | 59.56 | **51.224** | **+8.336** | accepted, pass |
| 19 | 166.78 | **49.864** | **+116.916** | accepted, pass |
| 20 | 47.55 | **47.877** | -0.327 | control pass |
| 26 | 391.55 | **54.563** | **+336.987** | accepted, pass |
| 27 | 76.22 | **59.271** | **+16.949** | **misses target by 0.271 ms** |
| 28 | 73.92 | **68.284** | **+5.636** | **misses target by 9.284 ms** |

The shipping 29-preset cycle validates, visits all presets, and wraps 28→0.
It has 27 green buckets and two red buckets, with **66/2304 spills**. Preset 28
accounts for 64; preset 0 owns two transition-attribution spills from preset
28. Therefore the zero-spill roster acceptance gate is **not met**.

The native suite passes, including dense kaleidoscope and Peirce reference
scans, seam cases, gamut error budgets, and the exact legacy-noise wrap oracle.
The updated shipping and global-O3 reports record the two complete cycle
captures.

Two presets landed concurrently after the device matrix. The recorded shipping
and O3 cycles are the plan's 29-preset corpus; post-rebase host tests visit all
31 presets and the 31-preset full-roster size gate passes, but the two additions
do not have device timing claims in this ledger.

## Post-curation reprofile (2026-08-11)

The authored bank was subsequently curated to 22 presets. A new matched pair
at `f80495c4` replaces the 29-preset cycle as the current roster measurement.
Both captures ran on COM3, visited all 22 indices, wrapped 21-to-0, had no
epoch reset, and passed cycle-counter exactness validation.

| Measurement / experiment | Decision | Frame gain / loss | ITCM cost / benefit |
|---|---|---:|---:|
| 22-preset shipping cycle | accepted baseline: 72.95 ms peak, 198/2016 spills, 17 green + 5 red buckets | reference | no code change; full roster remains 196,152/196,608 B ITCM |
| Matched global-O3 cycle | rejected as shipping direction: 67.22 ms peak, 135/2080 spills, same 5 red buckets | +5.73 ms worst-cycle peak, but holdouts remain | +23,296 B ITCM and +27,928 B FLASH in the single-effect image |
| Fixed preset 19 | accepted holdout baseline: 74.56 ms peak, 69.65 ms worst-window shader, 544/544 spills | reference | no code change |
| Generic `HS_PROFILE_DEEP` image | rejected for ShaderBall attribution: validates but exposes no child scope on the full-frame shader/private gamut path | diagnostic only | 0 B shipping; compiled out |
| ShaderBall DWT stage image | accepted diagnostic: projection 29.781 ms/f, color 12.123, lens 7.560, source 5.662, warp 2.760, material 2.710 | diagnostic only | 0 B shipping; compiled out |

Preset 19 is the shipping holdout. Its Peirce projection consumes 49.1% of
measured stage time and about 1,677 cycles per participating pixel. The
dodecahedral lens averages 4.84 reflections per pixel, reaches 14, and ranges
from 6.03 to 10.82 ms/frame as the animated orientation changes. Optimize the
Peirce kernel first, then the iterative polyhedral fold; color is the stable
third target. Broader ITCM promotion cannot close the gap because global O3
still leaves preset 19 at 60.81 ms and preset 21 at 61.61 ms of clean shader
work before non-shader frame overhead.

## Final 23-preset closure (2026-08-12)

Preset curation and one subsequent addition produced the current 23-preset
bank. The matched shipping baseline and final captures both use 16-frame
windows, the fast-cycle choreography, a 1,400-revolution epoch, COM3, and the
real segmented 288x144 driver. Both validate all 23 indices and the 22-to-0
wrap with no epoch reset.

| Cycle | Green / red | Worst peak, ms | Spilled frames | Full-roster ITCM | FLASH code |
|---|---:|---:|---:|---:|---:|
| Matched baseline | 6 / 17 | 88.07 | 849/1680 | 196,344 B | 401,808 B |
| Final `dad52e09` | **23 / 0** | **49.52** | **0/1728** | **195,960 B** | 411,112 B |
| Net | +17 / -17 | **+38.55 ms gain** | **849 spills removed** | **384 B saved** | +9,304 B |

The final Phantasm image retains 648 bytes of ITCM padding, 12,864 bytes for
local state, and 4,224 bytes of RAM2 allocator headroom. Moving the 6 KiB hue
field from `FrameState` into persistent state reduced ShaderBall's host stack
high-water mark from 18,271 to 7,007 bytes.

### Accepted implementation ledger

Positive frame gain means lower peak render time. Component A/B rows use the
holdout capture that exercised the changed kernel; the first and last rows are
matched full-cycle measurements.

| Change | Acceptance evidence | Frame gain | ITCM cost / benefit |
|---|---|---:|---:|
| Exact bounded dodecahedral chamber fold | Dense exact-reference and unit-vector scans; observed reflection telemetry bounds the unrolled fold | preset 18 early holdout: **68.16 -> 58.35 ms (+9.81)** | retained without a standalone promotion |
| Square Peirce sector renderer | Dense coordinate/fade scan and exact seam metadata | preset 19 early holdout: **74.56 -> 58.68 ms (+15.88)** | flash-resident exact oracle retained |
| Static inverse-pipeline specializations | Exact generic-vs-specialized shade oracle; generic GUI fallback retained | preset 0: **72.28 -> 58.58 ms (+13.70)** | full all-topology expansion was trimmed to fit; retained kernels are selective |
| Prepared surface-noise phase/direction and half-radian exponential map | Dense geometry error below 2e-7 and stage capture | preset 0: about **+5.11 ms**; surface stage about **+3.37 ms** | frame constants replace per-pixel trig; no separate promotion |
| Single-traversal vector simplex surface noise | Vector-field version is explicitly `DIRECT_VECTOR_V2`; scalar V1 remains available; native noise tests pass | preset 18: **65.09 -> 56.76 ms (+8.33)** | raw vector kernel stays in flash to protect ITCM |
| Skip liquid hue conversion when deformation is zero | Exact semantic shortcut | included in the final liquid-path gain | negligible measured ITCM effect |
| Persistent 64x16 value/hue field with bilinear lookup | Added error oracle: max introduced channel error 5,278/65,535, mean 202; 911,825 ShaderBall assertions pass | preset 21: **82.28 -> 44.56 ms (+37.72)**; final worst bucket **49.52 ms** | final roster **195,960 B**, 384 B below baseline; 6 KiB persistent storage, no stack penalty |

The formerly reported holdouts are all green in the final cycle:

| Preset | Baseline peak, ms | Final peak, ms | Gain, ms |
|---:|---:|---:|---:|
| 0 | 72.28 | 49.15 | **23.13** |
| 18 | 88.07 | 43.35 | **44.72** |
| 19 | 86.99 | 48.02 | **38.97** |
| 20 | 81.99 | 47.63 | **34.36** |
| 21 | 87.44 | 44.56 | **42.88** |
| 22 | 86.39 | 49.52 | **36.87** |

### Rejected and superseded ledger

| Experiment | Decision evidence | Frame gain / loss | ITCM / memory effect |
|---|---|---:|---:|
| Compact FastNoise gradient representation | Extra decode work dominated the smaller table | preset 18 **-3.7 ms** | 48 B ITCM and 828 B table storage saved; rejected |
| Full static specialization of every topology | Full-roster image exceeded the bank | not profiled | **+2,624 B ITCM** over the ceiling; rejected |
| Static transform shortcut | Slower holdout | preset 18 **-0.92 ms** | **+16 B ITCM**; rejected |
| Flash/noinline surface leaf | Instruction fetch outweighed bank recovery | preset 18 **-3.75 ms** | **592 B ITCM saved**; rejected |
| Shared topology helper promotion | Full-roster build failed before timing | not accepted | **+3,288 B ITCM** over the ceiling |
| Ungated vector-noise sampler in ITCM | Full-roster build failed; raw flash sampler retained instead | not accepted | 1,144-1,480 B over the ceiling |
| Endpoint-only hue bases | Useful partial result but did not close the animated holdout | preset 21 **82.28 -> 74.81 ms (+7.47)** | **+176 B ITCM**; superseded |
| Endpoint 128-entry interpolated hue LUT | Interior softened values still used exact OKLab conversion | preset 21 **74.81 -> 65.44 ms (+9.37)** | another **+272 B ITCM**; superseded |
| Nearest-bin endpoint hue LUT | Quantization and cache behavior regressed interpolation | preset 21 **65.44 -> 71.64 ms (-6.20)** | 80 B ITCM saved; rejected |
| 16-value 2D hue field | Fast, but the value grid exceeded the accepted color-error envelope | 23/23 green, 49.29 ms worst | 3 KiB frame storage; max channel error 12,148; rejected |
| 6 KiB hue field inside `FrameState` | Hardware timing passed, host stack gate failed | 23/23 green, 49.74 ms worst | stack **18,271 B > 12,288 B**; storage placement rejected |
| Global O3 final twin | Compiler layout catastrophically regressed the Peirce/dodecahedral topology | shipping 49.52 vs O3 **207.37 ms** worst; preset 19 48.02 -> 205.81 ms | single-effect image **+11,904 B ITCM**, **+16,144 B FLASH**; rejected |

Final validation comprises the on-device shipping cycle, the Phantasm roster
size gate, all 68 native CTest targets, the 7,007-byte stack watermark, and the
release WASM smoke roster at 96x20 and 288x144.

## Variadic inverse-pipeline follow-on

ShaderBall can benefit from the compile-time techniques used by `Filter`, but
its inverse sampler should be a sibling abstraction rather than a direct reuse
of the forward buffer-filter API. The proposed stage order is:

```text
OuterCamera -> SurfaceProject -> PlanarWarp -> Source -> Material -> Color
```

A variadic `InversePipeline<Stages...>` would fuse those point-sampling stages
into one callable shader. Authored runtime slots resolve once per frame to a
finite set of measured pipeline instantiations; arbitrary GUI combinations
continue through the current generic fallback. Keep a single `Scan::Shader`
driver, specialize only topology clusters that win on hardware, and leave
large or cold kernels flash-resident. This preserves exact generic oracles and
avoids the ITCM failure already demonstrated by instantiating every topology.

The useful shared surface with `Filter` is policy/concept machinery--stage
composition, compile-time traits, and optional fusion--not its forward
source/destination buffer contract. This remains a design candidate; none of
the final timing claims above depend on it.
