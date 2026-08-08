# ShapeShifter deep optimization analysis — 2026-07-29

## Outcome

ShapeShifter’s pre-lever shipping baseline at `95cc886a` is **104.99 ms**.
Reaching **<59 ms** requires at least **45.99 ms (43.8%)** of peak-render
reduction. Global O3 reaches 83.93 ms, proving that code generation is worth
21.06 ms, but it still needs another 24.93 ms.

Levers 2–4 and the ITCM recovery landed and reduced the measured shipping peak
to **41.81 ms**: **63.18 ms (60.2%) faster** than the recaptured baseline and
**17.19 ms under the 59 ms goal**. All four shapes now sustain the display
cadence with zero spilled frames in both shipping and global-O3 captures.

The target is unlikely to come from cache tuning, tighter bounds, blending, or
ISR work. It needs a ShapeShifter-specific fused Scan driver that shares
per-pixel spherical coordinates across the seven nested layers, removes
unneeded UV/Fragment work, and keeps layer blending typed and local.

Expected combined result:

- fused/common-frame Scan + minimal typed shading: **55–68 ms**;
- plus shape-math cleanup and a selective-O3 fused driver: **42–56 ms**.

These ranges overlap; they are not additive promises. Each lever needs an
on-device A/B over the same four-shape cycle.

## Render path and measured cost

```text
draw_frame
  Canvas acquire/clear                         0.09 ms
  timeline.step                             101.60 ms
    Sprite -> draw_all                      101.54 ms
      shared camera/bases/palette loop        0.01 ms
      7 × Plot dispatch                       0.59 ms
      7 × Scan dispatch                     100.94 ms
        SDF distance                         72.7 ms  (scaled diagnostic)
        PipelineRef::plot                    14.5 ms
          terminal filter_blend               3.84 ms
        bounds/interval/pixel-loop residual   6.9 ms
        coverage + Fragment setup             4.1 ms
        type-erased fragment shader           2.8 ms
  display-sync buffer wait                   22.35 ms  (not render)
```

The corrected DWT component pass measured 718,988 tests and 677,537 shade
candidates in a 16-frame heavy SphericalPolygon window:

| Component | Cycles/event | Share of Scan | Shipping scaled ms | O3 ms |
|---|--:|--:|--:|--:|
| SDF distance | 967/test (ship), 806/test (O3) | 72–74% | 72.7 | 59.8 |
| Pipeline plot | 205/candidate, 154 O3 | 13–14% | 14.5 | 10.8 |
| Bounds/interval/loop residual | — | 6–7% | 6.9 | 5.2 |
| Coverage + Fragment setup | 54/test, 27 O3 | 3–4% | 4.1 | 2.0 |
| Fragment shader | 39/candidate, 43 O3 | 3–4% | 2.8 | 3.0 |

The heavy shipping pose tests **46,974 pixels/frame** and shades **44,026**.
Only 6.3% of tested pixels are rejected, so perfect bounds could remove at
most about 4.6 ms of the distance bucket. Bounds are not the primary lever.

The same frame blends 44,087 pixels, 4.25× the 10,368-pixel quadrant.
`filter_blend` is already 52 cycles/blend and costs 3.84 ms. Even deleting it
entirely would not close one tenth of the gap.

## Algorithm audit by shape

All seven Scan layers share one basis. Their radii, phase, color, and shape
parameters differ, but each current per-ring raster independently reconstructs
the same pixel vector and repeats most spherical-coordinate work.

### Shared work repeated per layer

For a pixel visited by several nested layers, the generic path repeats:

- `dot(p, basis.v/u/w)`;
- `fast_acos`/`angle_between` for the polar angle;
- `fast_atan2` for azimuth;
- sector width and general `wrap`/`fmodf`;
- Fragment initialization and two indirect calls.

The heavy pose averages 4.53 SDF tests and 4.25 blends per quadrant pixel.
Computing the basis frame once per pixel is therefore the central Amdahl lever.

### Redundant unit-vector normalization

Star, Flower, and SphericalPolygon call `angle_between(p, axis)`. Both operands
are unit vectors by construction, but the helper recomputes both magnitudes,
executes two square roots plus a divide, and only then evaluates fast acos.
PlanarPolygon already uses the cheaper `fast_acos(clamp(dot(...)))` form.

Using the unit-vector contract removes high-latency work from every test of
three shapes. A debug assertion can validate the contract outside the hot loop.

### UV work that ShapeShifter never reads

The fragment shader ignores `p`, `v0`–`v3`, `size`, and `age`; it writes one
precomputed color and alpha. Scan nevertheless instantiates `ComputeUVs=true`:

- Planar, Flower, and SphericalPolygon divide by radius/thickness;
- Star normalizes azimuth with another wrap;
- generic `process_pixel` fills the entire Fragment register file.

A minimal no-UV path is output-equivalent for this effect.

### Per-probe invariant sector math

Every distance call recomputes `2π / sides`, folds with a general `wrap`, and
global O3 still calls the `fmodf` veneer. `sides` is frame-invariant and shared
by all seven layers. Precompute sector and reciprocal-sector once, then fold
with multiply + floor/subtract over the bounded azimuth/phase domain.

### SphericalPolygon edge distance

SphericalPolygon is the worst shape. Each test computes:

- normalized `angle_between`;
- fast atan2 and general wrap;
- fast sin and cos of the polar angle;
- fast cos of the local sector;
- flash-routed `asinf` for great-circle edge distance;
- UV division.

The basis dot already gives `cos(polar)`; `sin(polar)` can come from
`sqrt(max(1-dot², 0))`. Near the AA edge, `asin(dp)` and `dp` differ by
O(dp³); at one 288-wide pixel the error is far below a pixel. A monotonic
`sin(pixel_width)` threshold or a near-edge-only conversion can preserve AA
while removing the per-test flash libm call.

## Cortex-M7 disassembly audit

### Code placement and caches

The hot generic Scan loops and all four shipping SDF distance bodies are in
`.text.itcm` (addresses below 0x00010000). They do not suffer instruction-cache
misses. Shape data, stack, and trig tables are in tightly coupled RAM.

Two hot helpers are flash-routed:

| Helper | Flash code | Hot use |
|---|--:|---|
| `fmodf` + `__ieee754_fmodf` | 480 B | sector fold per SDF test |
| `asinf` + `__ieee754_asinf` | 620 B | SphericalPolygon edge per test |

Those calls can incur I-cache fills after eviction and always pay call/return
and flash-path execution. The harness has no PMU event sampling, so a miss rate
was not measured; the report treats this as exposure, not a measured miss.

The two frame buffers are `DMAMEM` in RAM2. A quadrant is larger than the
32 KiB M7 D-cache, and seven separate nested passes revisit it. D-cache churn is
plausible in `PipelineRef::plot`, but that entire measured bucket is only about
14.5 ms. A fused per-pixel layer loop improves locality automatically; a
cache-only project cannot deliver 44 ms.

### Register pressure and spills

Shipping keeps SDF bodies out of line:

| SDF body | Code | Calls inside | Saved stack per probe |
|---|--:|---|---:|
| Flower | 432 B | `angle_between`, `wrap` | 48 B each direction |
| Star | 544 B | `angle_between`, `wrap` | 56 B each direction |
| PlanarPolygon | 540 B | `wrap` | 40 B each direction |
| SphericalPolygon | 600 B | `angle_between`, `wrap`, `asinf` | 56 B each direction |

At 46,974 SphericalPolygon tests/frame, its GPR/VFP push+pop traffic alone is
about **5.3 MB/frame** in DTCM. DTCM has no cache-miss penalty, but those
loads/stores occupy issue bandwidth and serialize the call boundary.

Global O3 inlines each distance body into its pixel-run lambda. The spherical
lambda grows to 1,344 B and uses:

- nine saved GPRs, six saved VFP doubles, and a 44-byte local frame;
- 31 explicit stack references, including nine VFP references;
- two remaining indirect `blx` calls per shaded pixel (shader and pipeline);
- direct veneers for `fmodf` and `asinf`.

O3 therefore trades larger code and some spills for fewer per-probe call
frames and better scheduling. It wins 18.4%, but the spills and indirect calls
show why global O3 is not the algorithmic finish line.

### Pipeline stalls

The shipping spherical body contains seven `vdiv.f32` instructions, then calls
`angle_between`, which adds `vsqrt`/`vdiv`, and finally calls `asinf`. Several
divide results have only one or two independent instructions before their
consumer. These are dependency stalls on the M7’s single-precision FPU.

The O3 body schedules more independent work but still has dependent
divide/square-root chains, flash libm calls, and two indirect branches per
shaded pixel. The strongest scheduling optimization is to delete invariant
divides and repeated transcendental work, not to hand-reorder the surviving
assembly.

## Optimization levers

Forecast ranges below came from the original 102.85 ms capture. Actual wins
use the recaptured 104.99 ms shipping baseline at `95cc886a`; they overlap and
must not be summed naively.

| Priority | Lever | Expected peak win | Actual peak win | Expected peak after lever | Actual peak after lever | Result | Risk |
|---:|---|--:|--:|--:|--:|:---:|---|
| 1 | Generic typed `Scan::MultiShape`: one quadrant walk, shared pixel frame, ordered layer pack | 22–32 ms | **4.77 ms** | 71–81 ms | **38.35 ms candidate** | ❌ | +16,448 B ITCM, Plot-order delta, O3 regression |
| 2 | Unit-vector dot-angle path, precomputed sector/reciprocal, bounded fold, `ComputeUVs=false` | 9–16 ms | **32.77 ms** | 87–94 ms | **72.22 ms** | ✅ | Low–medium; golden-image seams/AA |
| 3 | Minimal typed shader/pipeline; keep a pixel in registers across its layer loop | 7–12 ms | **9.84 ms** | 91–96 ms | **62.38 ms** | ✅ | Medium; avoid PipelineRef/Fragment semantic drift |
| 4 | Spherical edge: derive sin/cos from the axis dot and remove/limit `asinf` | 4–9 ms | **19.26 ms** | 94–99 ms | **43.12 ms** | ✅ | Medium; quantify edge error |
| 5 | Selective O3 on the fused driver, kept under the current 9,976 B ITCM-granule headroom | 5–12 ms after fusion | — | Depends on 1–4 | — | — | Medium; size gate |
| 6 | ISR wake/pack tuning | 1–2 ms realistic | — | 101–102 ms alone | — | — | Cross-effect; not ShapeShifter-specific |
| 7 | Tighter bounds/culling | 0–4.6 ms theoretical max | — | ≥98 ms | — | — | Low payoff: 93.7% of tests shade |
| 8 | Further terminal blend tuning | 0–3.8 ms theoretical max | — | ≥99 ms | — | — | Already 52 cyc/blend |

### Combined implementation forecast

The recommended implementation combines priorities 1–4 in one driver; doing
them separately would duplicate work and misstate Amdahl gains.

| Build | Forecast peak | Margin to 59 ms |
|---|--:|--:|
| Fused/common-frame only | 71–81 ms | misses by 12–22 ms |
| Fused + minimal typed path + math cleanup | 48–61 ms | center estimate clears |
| Above + selective O3 fused driver | 41–55 ms | likely clears with margin |

The present global-O3 reference is a useful ceiling check: **83.95 ms** with
`+26,112 B` RAM1 code. A new fused driver should target less than 9,976 B of
additional ITCM so the shipping image does not cross the next 32 KiB FlexRAM
allocation granule.

### Immediate quality levers

If a visual compromise is acceptable before architecture work:

- `Count=3` is estimated at **54–61 ms** from nested-area scaling and should be
  captured directly; `Count=4` is more likely 65–72 ms.
- disabling Scan for heavy poses leaves the sub-millisecond Plot path, but
  materially changes the effect;
- a pose-aware layer cap can hold cadence while retaining all seven layers in
  light orientations.

## Recommended A/B sequence

1. Add a no-UV ShapeShifter scan dispatch and verify native golden pixels.
2. Replace `angle_between` with the unit-dot contract and precompute the sector
   fold; profile all four shapes.
3. Prototype the fused driver with identical layer order and direct color
   shading; compare frame buffers byte-for-byte where the math is unchanged.
4. Add SphericalPolygon’s dot-derived sin/cos and edge approximation behind a
   visual-error test.
5. Apply `HS_O3` only to the fused driver, check the map/ITCM granule, and run
   the shipping/O3 70-second twins again.

Acceptance should require: all four shapes visited and wrapped, peak render
below 59 ms, zero spilled frames, no full-roster ITCM regression, and native
render tests passing.

## Implementation ledger

The pre-lever captures at `95cc886a` are the comparison point for lever 2:

| Build | Peak render | Spilled | Capture validation |
|---|--:|--:|---|
| Shipping selective O3 | 104.99 ms | 249/816 | PASS; all four shapes and wrap |
| Global O3 | 83.93 ms | 94/992 | PASS; all four shapes and wrap |

### ✅ Lever 2 — landed `e328591b`

| Build | Before | After | Actual peak win | Spilled after |
|---|--:|--:|--:|--:|
| Shipping selective O3 | 104.99 ms | 72.22 ms | **32.77 ms (31.2%)** | 87/1008 |
| Global O3 | 83.93 ms | 60.27 ms | **23.66 ms (28.2%)** | 0/1088 |

Both captures passed all four-shape, wrap, epoch, render-source, and cycle/wall
exactness checks. Native CTest passed 57/57 and the focused SDF test passed.
The full-roster firmware linked, but the pre-existing RAM1 code overflow grew
from 197,256 B to 200,104 B (**+2,848 B ITCM**). The measured 32.77 ms shipping
win justified landing the candidate; later levers must recover code size where
possible.

### ✅ Lever 3 — landed `baec9544`

| Build | Before | After | Actual peak win | Spilled after |
|---|--:|--:|--:|--:|
| Shipping selective O3 | 72.22 ms | 62.38 ms | **9.84 ms (13.6%)** | 0/1088 |
| Global O3 | 60.27 ms | 56.13 ms | **4.14 ms (6.9%)** | 0/1088 |

Both captures passed all four-shape, wrap, epoch, render-source, and cycle/wall
exactness checks. The candidate was captured on COM3 while lever 2 was captured
on COM4; the 9.84 ms shipping delta is well above plausible board variance.
Native CTest passed 57/57, and exact framebuffer parity across all four shapes
passed for normal and debug paths (189,869 focused assertions). Full-roster
code grew from 200,104 B to 204,328 B (**+4,224 B ITCM**).

The current one-ring-per-scan traversal cannot keep a pixel resident across the
layer loop; that part of the forecast belongs to the still-unimplemented fused
driver. This landing captures the independently measurable typed
distance/coverage/color path and removes `FragmentShaderFn`, `PipelineRef`, and
full `Fragment` construction from ShapeShifter's scan hot path.

### ✅ Lever 4 — landed `a1186b90`

| Build | Before | After | Actual peak win | Spilled after |
|---|--:|--:|--:|--:|
| Shipping selective O3 | 62.38 ms | 43.12 ms | **19.26 ms (30.9%)** | 0/1088 |
| Global O3 | 56.13 ms | 36.59 ms | **19.54 ms (34.8%)** | 0/1088 |

Both COM3 captures passed all four-shape, wrap, epoch, render-source, and
cycle/wall exactness checks. The formerly dominant spherical phase fell from
62.38 ms to 38.25 ms shipping; PlanarPolygon now sets the 43.12 ms aggregate
peak. Native CTest passed 57/57.

The sine-domain edge approximation was tested at 2,345 AA-band samples with a
maximum angular error of `1.7285347e-6` rad. Across 165,888 framebuffer samples,
415 pixels changed and the maximum RGB16 channel delta was one count. Generic
SphericalPolygon callers retain exact `asin`; only ShapeShifter's typed solid
path uses the approximation. Full-roster code grew from 204,328 B to 204,680 B
(**+352 B ITCM**).

### ❌ Finding 1 fused `Scan::MultiShape` — candidate `edfffe3f`, not landed

| Build | Before | Candidate | Actual peak win | Spilled candidate |
|---|--:|--:|--:|--:|
| Shipping selective O3 | 43.12 ms | 38.35 ms | **4.77 ms (11.1%)** | 0/1088 |
| Global O3 | 36.59 ms | 37.02 ms | **−0.43 ms (1.2% regression)** | 0/1088 |

Both COM3 captures passed all four-shape, wrap, epoch, render-source, and
cycle/wall exactness checks. The generic fused driver passed full CTest 57/57,
including stack/arena budgets, but full-roster ITCM grew from 204,872 B to
221,320 B (**+16,448 B**). The roster's existing 196,608 B derived ceiling
would be exceeded by 24,712 B.

Normal rendering also had to move all Plot layers before the fused Scan pack
instead of preserving `Plot_i`/`Scan_i` interleaving. Across the deterministic
count-4 matrix, 427–559 of 3,072 pixels changed, the largest RGB16 channel
delta was 10,032, and total RGB delta was 0.077–0.150% of full-frame range.
Debug mode retained pixel-exact legacy ordering.

The 4.77 ms shipping gain is far below the 22–32 ms forecast, global O3 became
slower, and the code-size/ordering costs are material. The candidate is
therefore rejected and was not fast-forwarded to `master`.

### ✅ ITCM recovery — landed `4bab2e8e`

| Result | Before | After | Actual peak after lever | Status |
|---|--:|--:|--:|:---:|
| Full-roster RAM1 code | 205,432 B | 194,344 B | — | ✅ |
| ITCM headroom | −8,824 B | **2,264 B** | — | ✅ |
| Shipping selective O3 | 43.12 ms | **41.81 ms** | **41.81 ms** | ✅ |
| Global O3 | 36.59 ms | **34.87 ms** | **34.87 ms** | ✅ |

The recovery reuses one coalesced run list per bounded row, folds solid debug
rendering into the typed kernel, compiles out disabled firmware pole-LOD
branches, and moves only setup/transition code to flash. The final live-master
build reclaims **11,088 B**, restores **2,264 B** of headroom (192 B more than
the historical 2,072 B state), and keeps the typed ShapeShifter dispatch.
Both device captures visited and wrapped all four shapes with zero spills.

### ❌ Generic no-UV fallback for three typed shapes — rejected

| Build | Accepted typed peak | Candidate peak | Actual peak after lever | Status |
|---|--:|--:|--:|:---:|
| Shipping selective O3 | 43.12 ms | 43.31 ms | **43.31 ms** | ❌ |
| Global O3 | 36.59 ms | 39.00 ms | **39.00 ms** | ❌ |

The fallback saved code but regressed the active render path: O3 PlanarPolygon
grew by 2.41 ms, Flower by 0.99 ms, and Star by 1.24 ms. It was removed before
landing.
