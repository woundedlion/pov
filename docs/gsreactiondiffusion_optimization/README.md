# GSReactionDiffusion optimization plan

Status: complete. Baseline captured 2026-07-30 at `dd25b7d0`; final
code captured at `1deec495`.

The goal is a sustained 16 fps shipping result on Teensy 4.0 under the real
segmented Phantasm driver, not a host-only speedup and not a global-`-O3`
single-effect result. The hard acceptance target is:

- shipping peak render at or below 58 ms over a full 130 s reaction cycle;
- zero frames above the 62.5 ms display window;
- no regression in the native suite, arena gate, Phantasm size/layout gate,
  parameter behavior, dissolve/reseed lifecycle, or visible pole/seam quality;
- at least 4 KiB of headroom inside the active ITCM bank after every kept
  change.

The 58 ms target leaves about 4.5 ms for workload drift. A result that merely
lands at 62.4 ms is not finished.

## Measured baseline

Both captures used `POVSegmented<288,4,480>`, DMA LEDs, live flywheel/DMA
interrupts, a 32-frame window, and a 150 s epoch. Each capture ran for 130 s.
Both logs validate with no epoch reset and cycle-counter/wall agreement within
1.6 ppm.

| configuration | peak render | spilled | pass mean | cheapest window |
|---|---:|---:|---:|---:|
| shipping selective-O3 | 94.83 ms | 1011/1024 (98.7%) | 88.10 ms | 70.49 ms |
| global O3 reference | 89.73 ms | 995/1024 (97.2%) | 83.89 ms | 68.20 ms |

Global O3 saves only 5.10 ms at the peak, or 5.4%. Shipping still needs a
32.33 ms cut merely to touch 62.5 ms and a 36.83 ms cut to reach the 58 ms
acceptance target. This is primarily an algorithm/workload problem.

Peak shipping window, frames 961-992:

```text
frame                 124.98 ms  74.99 Mcyc  100%
  grd_render           94.35 ms  56.61 Mcyc   75%
    grd_rasterize      62.78 ms  37.67 Mcyc   66%
      grd_shader_draw  59.04 ms  35.42 Mcyc   94%
    grd_simulate       31.56 ms  18.94 Mcyc   33%
  canvas_buffer_wait   30.52 ms  18.31 Mcyc   24%
```

The raster is coverage-dependent, but simulation is a fixed 31.6 ms floor.
The cheapest shipping window still spends 31.6 ms simulating and 34.8 ms in
the shader, so coverage culling alone cannot cross the cadence tier.

### Exact compiler and artifact provenance

The flashed shipping and O3 ELFs and the full Phantasm ELF all contain:

```text
GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
```

| artifact | SHA-256 |
|---|---|
| profile shipping ELF | `8866821075631143B99F80E40543FB5A4CC644C59A6B3E7B0081400713CB6209` |
| profile O3 ELF | `E14143D4E33F4FC685E090E187250A72729215D040C7E2E65EE82331A5B8CAB4` |
| full Phantasm ELF | `D1D38F1C2B20A6048C760CBF6635872D100AEDB9FBDF264EDF05C97EE841ECFA` |

Stale `build/prof/sizes/` files initially made the current capture look like a
GCC 11 build. They were not produced from the flashed artifacts. Phase 0 makes
that ambiguity impossible.

### Image and FlexRAM budget

The exact profile ELFs have these sizes:

| configuration | FLASH code | RAM1 variables | ITCM code | ITCM delta |
|---|---:|---:|---:|---:|
| shipping | 54,508 B | 314,464 B | 31,320 B | - |
| global O3 | 71,524 B | 314,496 B | 47,256 B | +15,936 B |

Spending 15,936 B buys only 5.10 ms. The real full-roster Phantasm image is
tighter:

```text
RAM1 variables: 314,176 B
ITCM code:       192,152 B
ITCM ceiling:    196,608 B
headroom:          4,456 B
```

No experiment may assume global-O3-sized code fits the shipping image.

## Disassembly findings

The shipping GCC 15.2.1 ELF gives concrete codegen targets.

### Interpolation

The GS specialization of `refine_and_accumulate` is 2,384 B and 627
instructions. Its generated body has:

- a 60-byte stack frame;
- 143 `vldr` and 33 `vstr` instructions;
- 42 distance-coordinate subtracts, 42 multiplies, and 42 FMAs;
- 26 floating comparisons followed by 26 `vmrs` flag transfers;
- stack-resident `d2s[7]` and `ids[7]`;
- one type-erased fragment call per sub-sample, 41,472 calls per frame.

The split `Scan::Shader::draw<288,144,4>` body is another 892 B. GS is still
using this type-erased path. BZ already uses the templated `draw_grid` path,
refines and gathers once per pixel, and owns all four samples inline.

### Physics

The shipping `step_physics` body is only 248 B, but the inner node loop:

- reloads all five parameter values from `this` for every node update;
- branches around a six-iteration neighbor loop;
- computes two indexed addresses for every neighbor, one for A and one for B;
- is called 16 times per frame;
- performs 122,880 node updates and roughly 1.47 million gathered neighbor
  float loads per frame.

The Q16 input loop is scalar. The Q16 output loop calls the 48-byte `to_q16`
function 15,360 times per frame in shipping instead of inlining it.

Global O3 improves simulation only from 31.56 to 30.75 ms in the peak window.
Inlining/unrolling alone is not the complete physics answer, but the shipping
assembly has several cheap misses worth fixing before changing the numerical
scheme.

## Execution rules

Every lever is a separate commit and a separate device A/B. Do not bundle a
codegen annotation, an algorithm change, and a quality change into one result.

For every candidate:

1. Build and capture the shipping image and its global-O3 reference with the
   exact compiler used by the real Phantasm build.
2. Archive the compiler fingerprint, complete build flags, ELF/map hashes,
   section sizes, relevant symbol sizes, and disassembly.
3. Run a 130 s capture with `HS_PROFILE_EPOCH_REVS=1200`; validate the log
   before reading it.
4. Compare peak render, spill fraction, phase times, ISR share, ITCM, FLASH,
   stack frame, and arena high-water against the immediately preceding commit.
5. Run `pio run -e phantasm`, the native suite, and the arena/stack gates.
6. Keep the change only if it passes its stated performance and quality gate.

Host benchmarks are useful for iteration, never for a keep decision.

### Correctness classes

- **Class A, bit-identical:** pointer qualification, data-layout-only changes,
  redundant-load removal, and codegen restructuring should preserve Q16 state
  and framebuffer hashes.
- **Class B, numerically bounded:** stencil reuse and harmless reassociation
  may move near-tie pixels. Compare Q16 state, per-frame framebuffer error,
  temporal error, pole/seam crops, and transition timing over at least 640
  substeps.
- **Class C, visual redesign:** sample-count, integrator, fixed-point, and
  prefiltered-field experiments require recorded side-by-side video over the
  default reaction plus the feed/kill/diffusion/speed corner matrix. They do
  not land on timing alone.

For Class B/C, record mean and maximum channel error, changed-pixel fraction,
SSIM, lit coverage, reaction lifetime, and stabilization/dissolve frame. A
single summary hash is insufficient.

## Phase 0: make profiles self-proving

Do this before accepting any optimization result.

1. Extend `profile_one.sh` to extract the compiler string from `.comment`, the
   fully expanded compile flags, ELF SHA-256, git source state, and `teensy_size`
   output from the exact ELF it flashes.
2. Build the real `phantasm` environment in the same run and reject the profile
   if compiler identity, CPU/FPU ABI, fast-math flags, framework version, or
   core version differ.
3. Save the exact ELF and map beside the raw log, or save content-addressed
   copies and put their hashes in the log header.
4. Stop using independently generated `build/prof/sizes` files as report
   inputs. Size the archived flashed ELF.
5. Add a validation failure for a profile whose provenance block is missing or
   whose artifact hash does not match.

Acceptance: deliberately point the profile at the old GCC 11 package and prove
the run fails before flashing.

## Phase 1: deepen the baseline without distorting it

Add whole-loop scopes and counters, not per-pixel RAII scopes:

- `grd_q16_in`, `grd_substeps`, `grd_q16_out`;
- `grd_orient`, `grd_hot_flags`, `grd_grid`, and `grd_pixel_body`;
- counters for hot2-accepted pixels, refined-center changes, fallback second
  walks, supported stencil nodes, and palette hits.

Add a fixture-driven on-device microbench for `step_physics`, the pixel body,
and the Q16 conversion loops. The final result still comes from the live driver;
the microbench exists to explain the result and inspect cycles/update or
cycles/pixel without coverage drift.

The new scopes must compile out of shipping and must not change an unprofiled
Phantasm ELF's code or data size.

## Phase 2: remove raster architecture waste

This is the highest-priority performance phase.

### R1: typed pixel owner, otherwise exact

Introduce a GS `shade_pixel` and call the existing templated
`Scan::Shader::draw_grid`. Initially keep the current per-sub-sample
`refine_and_accumulate` semantics.

Purpose:

- remove 41,472 type-erased fragment calls per frame;
- eliminate `Fragment` copy/default traffic from the generic split path;
- let the compiler see cull, interpolation, palette, and accumulation as one
  body;
- establish the cost of dispatch removal independently of stencil reuse.

Gate: Class A framebuffer/state identity. Keep for at least 2 ms peak win with
no more than 2 KiB net ITCM growth, or for a net ITCM reduction.

### R2: port BZ's per-pixel stencil gather to GS

For each hot pixel:

1. refine the center once;
2. gather the seven world-space node positions and seven Q16 B values once;
3. evaluate four sub-sample weight sets over that fixed stencil;
4. accumulate the premultiplied pixel inline.

Preserve the existing exact two-ring cold-pixel cull before gathering.

This removes repeated neighbor-table and B-state gathers, the 60-byte
`d2s`/`ids` construction per sample, most float-compare/`vmrs` traffic, rare
fallback walks, and the standalone 2,384-byte GS interpolation specialization.
It follows a measured BZ precedent: the corresponding typed grid/stencil/O3
bundle moved BZ from 103 ms to 58.4 ms and a hard 16 fps lock.

Gate: Class B. Keep for at least 8 ms peak shader win, no visible boundary
sparkle, and a full Phantasm size pass. Record the exact rate at which the
center stencil differs from a sub-sample's independently refined stencil.

### R3: compare three center strategies

Measure these separately:

- center refine plus shared seven-node stencil;
- cubemap seed directly as the shared stencil center;
- a union stencil containing the center plus neighbors needed by all four
  independently refined sub-samples.

The direct seed uses only the four weight walks and may beat the BZ pattern,
but it has the largest interpolation error. The union preserves more exactness
but can erase the cycle win. Choose from device cycles and image error, not
intuition.

### R4: sample-count experiments

Only if R1-R3 do not put peak render below about 75 ms:

- two diagonal samples, alternating diagonal by frame;
- one center sample in smooth interiors and four near the B threshold or a
  large local B gradient;
- temporal 2x sampling with a deterministic four-frame pattern;
- analytic center value plus first-order gradient estimates for the four
  sub-samples.

Adaptive 1x/4x was a measured pessimization for BZ because most dense spiral
fronts selected 4x. GS has an exact cold-pixel cull and different coverage, so
one instrumented retry is justified, not assumed.

Gate: Class C and at least 8 ms peak win. Reject temporal shimmer, threshold
sparkle, and pole/seam asymmetry even if still-frame metrics pass.

### R5: prefiltered spherical field, moonshot

If direct kernel sampling remains dominant, prototype a frame-local scalar
field:

- precompute filtered B at either the 7,680 lattice nodes or the 24,576 cubemap
  cells once per frame;
- raster samples then perform a cheap nearest/bilinear lookup and palette map;
- alternatively forward-splat active nodes into the screen quadrant and
  normalize accumulated weights.

This changes work from roughly `pixels * samples * stencil` toward
`nodes * stencil + pixels * cheap_sample`. It is allowed to replace the current
renderer only if Class C video shows no crawling holes, cubemap-face seams, or
orientation-dependent blur.

## Phase 3: cut the fixed simulation floor

### S1: surgical shipping-codegen cleanup

Apply and measure one at a time:

1. snapshot `feed`, `k`, `d_a`, `d_b`, and `dt` outside the node loop;
2. add proven non-aliasing (`restrict`) to the four float buffers;
3. force-inline `to_q16` and inspect that all 15,360 calls disappear;
4. unroll the fixed six-neighbor gather;
5. software-pipeline the A/B loads and FMAs, with optional `pld` of the next
   neighbor row;
6. compare an O2 function region with O3. O2 may keep the useful scheduling
   without O3's code expansion.

Do not combine coefficients until the bit-identical variants are exhausted.
Then test precomputed `d_a*dt`, `d_b*dt`, `feed*dt`, and `(k+feed)*dt` as a
Class B reassociation experiment.

Gate: each kept change saves at least 0.5 ms simulation or removes ITCM while
remaining neutral. The complete S1 set should target 3 ms.

### S2: interleaved A/B float scratch

Replace four SoA float buffers with two ping-pong arrays of:

```cpp
struct alignas(8) AB {
  float a;
  float b;
};
```

Each neighbor index then produces one address and one adjacent A/B load pair.
This attacks address generation and DTCM traffic without changing the
algorithm. Preserve arithmetic order first for a Class A attempt.

Inspect for paired `vldmia`/`ldrd`, fewer address calculations, no new spills,
and improved dual-issue scheduling. Gate: at least 2 ms simulation win.

### S3: normalize work to wall-clock cadence

The current effect performs 16 substeps at about 8 fps: 128 substeps/s.
At the target 16 fps, eight substeps/frame also perform 128 substeps/s. A fixed
eight-step final configuration therefore preserves today's real-time reaction
speed while halving the per-frame simulation workload.

Test eight substeps only after the raster path is fast enough to approach
16 fps. Convert lifecycle constants that are meant to represent time from
rendered-frame counts to substep or elapsed-revolution counts, so a cadence
tier crossing does not silently halve grow/dissolve duration.

Expected first-order saving: about 15.8 ms/frame. Gate: Class C, 130 s capture,
same substeps/s and matched wall-clock morphology video against the 16-step
8 fps baseline.

### S4: persistent float state

Keeping A/B as float removes both full-lattice Q16 conversion passes and needs
only one alternate float generation. It also changes the quantization contract.
Prototype it only after recording the cost of `grd_q16_in/out`.

The memory version must pair with compact hot flags and the Phase 4 node-table
move. Otherwise persistent growth plus the raster scratch peak exceeds the
298 KiB device arena.

Gate: Class C. If the conversion scopes are under 1.5 ms combined, reject this
complexity unless it also enables a better data layout.

### S5: packed fixed-point/DSP kernel

Prototype packed Q15 A/B state in one 32-bit word and use Cortex-M7 DSP
instructions for diffusion sums and saturation. The intended shape is:

- `__SMLAD`/`__SMUAD` or paired halfword operations for neighbor accumulation;
- Q-format coefficients prepared once per frame;
- 32- or 64-bit intermediates for `a*b*b`;
- `__SSAT16`/packed clamp at the output.

This can eliminate float ping-pong, conversion passes, and much of the scratch
arena. It is high-risk: signed Q15 loses one concentration bit and the cubic
term needs a carefully proven scale.

Gate: Class C, exhaustive white-box parameter corners, explicit overflow proof,
and at least a 2x simulation speedup. Do not land a fixed-point rewrite for a
small win.

### S6: fewer stable stages per simulated second

If eight Euler steps still dominate, compare:

- operator-split reaction and diffusion;
- a stabilized Runge-Kutta-Chebyshev diffusion step;
- semi-implicit/Jacobi diffusion with an explicit reaction update;
- two Euler steps fused into one tiled pass where graph dependencies permit.

The aim is the same wall-clock evolution with fewer full graph gathers, not a
larger unstable Euler `dt`. Gate this as a numerical-method change with a
stability sweep across every parameter corner.

## Phase 4: unlock a FlexRAM bank for selective codegen

The full image has only 4,456 B of ITCM headroom. The arena is 298 KiB because
GS is the binding tenant at 297,984 B; BZ is second at 279,552 B. Both keep the
same immutable 7,680-entry float node table in the per-effect arena.

### M1: place exact immutable nodes in flash

Generate one shared `PROGMEM` table containing the exact float bit patterns
currently produced by `ReactionGraph::node`. Both RD effects read it
sequentially while producing their frame-local oriented DTCM copy.

This removes 92,160 B from both GS and BZ persistent footprints without
changing the hot random-access table: all per-pixel walks still read the
oriented DTCM copy. It costs 92,160 B of flash data, which fits the current
479,240 B file headroom.

Benchmark the sequential flash read during `grd_orient`; the keep gate is
Class A output and no more than 0.5 ms/frame regression. If exact floats are
too expensive in flash, compare a 46,080 B snorm16 table as Class B.

With GS/BZ no longer binding, shrink `DEVICE_GLOBAL_ARENA_SIZE` by one 32 KiB
bank to 266 KiB. This moves one FlexRAM bank from DTCM to ITCM and raises the
code ceiling from 196,608 B to 229,376 B while keeping the 12 KiB stack floor.

### M2: compact the hot flags

`hot1` and `hot2` consume 15,360 B as byte arrays. Convert them to bitsets
(1,920 B total) or eliminate `hot1` with a two-pass packed construction.

Gate: Class A cull decisions and at least 12 KiB lower raster scratch. Speed may
be neutral or slightly worse; keep it if it enables persistent float state or
an arena-bank transition.

### M3: precompute the cubemap LUT in flash

The 49,152 B cubemap LUT is immutable after init. A generated flash table removes
its build and persistent allocation. Its per-pixel-center accesses are hotter
and less sequential than the node-table read, so test it separately and reject
cache-thrash regressions.

This does not unlock another bank by itself because DisplacementField becomes
the arena binding tenant at about 242,404 B. Reaching a second smaller arena
bank would also require roughly 3 KiB of safe DisplacementField footprint
reduction. Keep that as an optional cross-effect project, not hidden GS scope.

### M4: spend the bank deliberately

After M1 and the arena shrink, apply selective O2/O3 only to the measured GS
typed pixel driver and, if profitable, the physics loop. Diff the map after
each function. Preserve at least 4 KiB intra-bank padding.

The bank is an enabler, not a reason to promote code. A region that does not
buy at least 0.5 ms/KiB or a cadence-tier crossing is reverted.

## Integration order and decision points

Recommended landing order:

1. profile provenance enforcement;
2. deep whole-loop counters;
3. R1 typed exact pixel path;
4. R2/R3 shared-stencil winner;
5. S1 codegen cleanup and S2 A/B layout;
6. M1 shared flash nodes, arena shrink, and renewed ITCM baseline;
7. narrow selective O2/O3 on the winning typed loops;
8. S3 eight substeps at the achieved 16 fps cadence;
9. M2/S4 only if conversion or memory remains material;
10. R4/R5 or S5/S6 only if the peak is still above 58 ms.

After R2 + S1 + S2, decide:

- **Peak at or below 74 ms:** S3 should plausibly cross 58 ms; pursue the
  cadence-normalized eight-step path.
- **Peak 75-82 ms:** add the FlexRAM/typed-O3 work before changing samples.
- **Peak above 82 ms:** R2 under-delivered; instrument stencil behavior and
  test R3/R4 before touching the numerical method.
- **Peak below 58 ms with 16 steps:** keep 16 steps unless wall-clock morphology
  is undesirably fast at the new frame rate. Performance does not itself
  justify changing dynamics.

## Final validation matrix

The final candidate must pass:

- shipping and O3 130 s GS captures, both provenance-valid;
- a second shipping capture on the other attached Teensy to expose board drift;
- `pio run -e phantasm` with size/layout gate;
- native tests and 120-frame arena high-water sweep;
- stack budget measurement;
- deterministic default and parameter-corner physics tests;
- framebuffer/image metrics over growth, saturation, dissolve, and reseed;
- pole/seam crops and side-by-side video;
- full Phantasm roster smoke test to catch a shared RD/Scan regression;
- final symbol-size and disassembly audit confirming the intended inlining,
  stack-frame, load, branch, and spill changes actually reached the ELF.

The final report must state the exact compiler, flags, commit, ELF hash, map
hash, device port, capture timestamps, peak/spill numbers, ITCM/FLASH deltas,
arena high-water, and any accepted visual/numerical divergence.

## Definition of done

This effort is complete only when the shipping full-roster image holds
GSReactionDiffusion at a peak below the owner-set 59 ms ceiling, or when every
ranked architectural lever above has measured evidence and the remaining gap
is explicitly traded against an owner-approved visual/numerical change.

“Global O3 is faster,” a host benchmark, or one favorable short capture is not
a completion condition.

## Execution report

The implemented campaign uses the shipping Arm GNU 15.2.1 compiler for every
device measurement. `profile_one.sh` builds Phantasm and the profile image,
compares compiler comments, ARM attributes, PlatformIO packages, and build
flags, then archives and hashes both exact ELFs before capture.

### Final result

| Metric | Original shipping | Final shipping | Final O3 |
|---|---:|---:|---:|
| Lifecycle capture | 130 s | 130 s | 130 s |
| Windows | 48 | 64 | 64 |
| Worst render | 94.830 ms | **58.783 ms** | **58.818 ms** |
| Dense-window shader average | 59.04 ms | 42.560 ms | 42.569 ms |
| Dense-window simulation average | 31.56 ms | 11.964 ms | 11.930 ms |
| Phantasm ITCM code | 192,296 B | 192,216 B | 192,216 B |
| ITCM bank headroom | 4,312 B | 4,392 B | 4,392 B |

Shipping worst-case render fell by 36.047 ms (38.0%). No captured render frame
exceeded the owner-set 59 ms ceiling. The final Phantasm image
uses 331,340 B of flash code and 1,209,676 B of flash data, leaving 482,312 B
for files. The final shader specialization is 2,452 B and its hot body contains
673 instructions, 76 `vldr`, and 22 `vstr`; the O3 physics kernel is 356 B and
90 instructions. No indirect shader call remains.

The final 130-second captures are:

- `build/prof/gsreactiondiffusion_threshold_final_ship.log`
- `build/prof/gsreactiondiffusion_threshold_final_o3.log`

Their provenance sidecars identify
`GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`.
Both captures used COM3, `POVSegmented<288,4,480>`, a 32-frame window,
`HS_PROFILE_EPOCH_REVS=1200`, and the real build's Cortex-M7/FPU/fast-math
contract. The shipping capture completed 2026-07-30 23:02:07 PDT and the O3
capture completed 2026-07-30 23:06:14 PDT.

| Exact artifact | ELF SHA-256 | map SHA-256 |
|---|---|---|
| shipping profile | `7C260B1D1F6EDA4B885DE3C9C31F87DC58B9FC7D68FADDC5F4B23717C73D8792` | `D15DFAA367CF243697C54E3E528F59F30E80A2F7C617C4B97295AF5D00F52B62` |
| O3 profile | `D449CC80DB0F5EA813C4F39D64950FEE7DCB72A133125DD284CAA7E8B60F0729` | `5E2A553CF7F927A816BF1AF80B36D70ACF2DB45A5BCD552B1FB69642809C94DC` |
| shipping Phantasm | `2CAE8174E17C02FDDBFBD66406FBA097DCE8DE6BFC2E581CA37EB1B12D583D58` | `96A343588117AE0F27F958225849E287CAA6D8DDF53CD7FC2EE3AEF1C1CB56C0` |
| O3 Phantasm | `0B5070540C0E2BC001D4AAB73E11318B7D7E516549B65752310C833FD10DC13A` | `01C44A5513E960DEFFA01F8851FBC82F52286536EF21E1FE4EC887206185CDCD` |

### Landed design

- The typed four-sample renderer removes type erasure and indirect calls.
- Four SSAA samples share one pixel-center stencil.
- A triangle-inequality radius accepts provably nearest cubemap seeds without
  checking six neighbors; production probes report zero center mismatches.
- Physics buffers are restricted, frame parameters are hoisted, Q16 output is
  forced inline, stabilization motion accumulates directly in Q16 units, and
  the graph kernel alone receives O3 codegen.
- Six 5/3-sized Euler integrations cover the same simulated interval as ten
  smaller integrations. Lifecycle counters advance in original-size substep
  equivalents so transitions track the faster evolution.
- Nodes enter an eight-frame dissolve band that eases A/B toward rest before
  the frontier pins them, eliminating the previous hard disappearance.

The temporal block is an accepted numerical-method change: it preserves the
equations and simulated interval but is not bit-identical to eight smaller
Euler stages. Scalar-equation parity, worst-corner bounded evolution,
production framebuffer bounds, full dissolve/reseed lifecycle, arena, stack,
and all 58 native tests pass.

### Rejected experiments

- A fused refine/gather loop increased the dense peak to 64.57 ms.
- Reusing the raw 64² seed reached about 58.5 ms but exceeded hard visual-error
  and peak-channel gates.
- Checking three or four fixed neighbors also missed those gates; five passed
  but peaked at 63.45 ms.
- A compressed 256² flash seed table passed aggregate visual gates only with
  boundary fallback and regressed the peak to 69.90 ms from flash-cache traffic.
- Interleaved state, forced unrolling, and broad O3 codegen produced only small
  wins relative to their footprint and were reverted.
