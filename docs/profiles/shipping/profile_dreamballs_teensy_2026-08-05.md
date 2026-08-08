# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-08-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DreamBalls`).
Raw capture: `build/prof/dreamballs_ship.log` in the `C:/work/hs-rc-db`
worktree, captured 2026-08-05 13:10.

Captured at master tip **`e68a4e37`** ("plot: take RasterOptions by value so its
fields stay constants"). **Supersedes the `446f5bd0` capture** that the previous
revision of this file was written from, on two counts: that capture predates
`ead9687b`, so it profiled a **4-preset** cycle where the shipping effect now
runs **5**, and it predates the `rasterize` signature change that this effect's
raster path crosses. The snubCube preset is measured here for the first time.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, board on **COM3** |
| Image | `profile` env: `-Os` + `--specs=nano.specs`, `USE_DMA_LEDS`, `HS_PROFILE_ENABLE`, arm-gcc 15.2.1. Hot path crosses the `HS_O3` regions `Plot::rasterize`, `rasterize_geodesic_strategy`, `edge_fits_one_dot` (core/render/plot.h) and filter.h's blend sink |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master, 72 px/arm, HD107S @ 24 MHz SPI |
| Effect | DreamBalls 288×144, single-entry playlist, tip `e68a4e37` (clean tree) |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 215 s capture, `HS_PROFILE_EPOCH_REVS=2000`; 213 windows / 3,408 frames, `validate` 8/8 |
| Reproduce | `just profile DreamBalls` (window 16, `HS_PROFILE_EPOCH_REVS=2000`, `--seconds 215`) |

**Why 215 s and a 2,000-rev epoch.** The cycle is now five presets at a
320-frame dwell = 1,600 frames, ≈ 100 s at the 16 fps every preset holds. 215 s
covers **2.1 full cycles**, so the wrap to preset 1 is demonstrated twice rather
than inferred from a single pass; `HS_PROFILE_EPOCH_REVS=2000` puts the epoch at
250 s so the capture never crosses an epoch boundary. The skill's tabled
DreamBalls budget (170 s / 1,600 revs) was sized for the 4-preset cycle and is
now short of one wrap plus margin.

Image size (`profile` env, the image measured): `FLASH: code:105656,
data:179940, headers:8292` / `RAM1: variables:315168, code:48408,
padding:17128, free:143584` / `RAM2: variables:520064, free:4224`.

The shipping `phantasm` image built from the same tip: `FLASH: code:347968,
data:295820, headers:8500` / `RAM1: variables:313760, code:193704,
padding:2904` — ITCM headroom **2,904 B**, size gate PASS.

Exactness cross-check: window frames 3329–3344 root counter 601,191,654 cyc ÷
600 MHz = 1,001,986.09 µs against the measured wall sum of 1,001,987 µs —
**0.9 ppm** (`tools/parse_profile.py ... validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `db_mesh_plot` avg
16.53 ms/f, worst window 41.91 ms/f (frames 3329–3344), peak frame render
**44.88 ms** (frames 193–208), spilled **0/3,408 frames (0%)**.

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px
(DreamBalls declares neither `needs_full_frame()` nor `persists_pixels()`, so
the per-frame pixel counts below are quadrant counts — the 5× coverage is
overdraw from stacked alpha wireframe copies, not a full-canvas render). Every
window of the pass holds 16 fps: wall avg 62.43 ms/f across the capture, and the
worst single frame renders in 44.88 ms — 72% of one display window, 17.62 ms of
margin (1.39× headroom). The quietest window renders in 1.17 ms/f against the
busiest hold's 42.34 ms/f, so the pass spans a 36× dynamic range inside one 16
fps tier and never leaves it. The `canvas_buffer_wait` scope
is the round-up idle to the next display flip, by design, and it absorbs the
whole spread: 22.94 ms/f at the peak, 61.28 ms/f in the trough.

## Phase-by-phase readout

Phase schedule: DreamBalls cycles **five** Archimedean-solid presets, 320 frames
each, in order 1 → 2 → 3 → 4 → 5; the capture covers two full cycles plus a
third entry into preset 1 (wrap confirmed twice, at frames 1601 and 3201, from
the `Preset: 1/5` markers). Within a preset the cost breathes by up to 25× as
the orbiting copies and the global tumble sweep the cluster in and out of the
rendered band, so the hold blocks below are the pass peak, a mid preset, the new
snubCube preset and the pass trough; the fifth is the one non-hold regime, the
sprite hand-off.

### Preset 1 hold — peak (window frames 193–208, worst of the capture)

```
frame                    61.76 ms  37.05 Mcyc  100%
  db_timeline_step       38.73 ms  23.24 Mcyc   63%
    db_draw              38.63 ms  23.18 Mcyc   63%
      db_draw_scene      38.63 ms  23.18 Mcyc   63%  x1.0
        db_mesh_plot     38.37 ms  23.02 Mcyc   62%  x18.0  2132 us/call
          filter_blend    5.19 ms   3.12 Mcyc    8%  x52221  60 cyc/blend
        db_warp_orient   114.5 us   68.7 kcyc    0%  x18.0   6.4 us/call
        db_displace      139.7 us   83.8 kcyc    0%  x18.0   7.8 us/call
      db_mesh_copy        0.8 us      500 cyc    0%  x1.0
  canvas_clear           87.7 us    52.6 kcyc    0%  x1.0
  canvas_buffer_wait     22.94 ms  13.76 Mcyc   37%  x1.0
```

Wall min/avg/max = 51.01/61.76/74.01 ms. This is the rhombicuboctahedron preset
at 18 copies — 18 × 48 = 864 edges plotted per frame, the largest edge budget in
the cycle — and where the pass's worst frame lands (44.88 ms render).
`db_mesh_plot` is effectively the entire frame: 62% against 0.4% for
`db_displace` + `db_warp_orient` combined and 0.14% for `canvas_clear`. Of
`db_mesh_plot`'s 38.37 ms, only 5.19 ms is `filter_blend`; the remaining 33.18
ms is the geodesic walk itself — arc setup, stepping, clipping and shading. Even
here 37% of the frame is idle waiting for the flip.

### Preset 4 hold — mid (window frames 1105–1120)

```
frame                    62.70 ms  37.62 Mcyc  100%
  db_timeline_step       24.50 ms  14.70 Mcyc   39%
    db_draw              24.40 ms  14.64 Mcyc   39%
      db_draw_scene      24.40 ms  14.64 Mcyc   39%  x1.0
        db_mesh_plot     24.23 ms  14.54 Mcyc   39%  x10.0  2423 us/call
          filter_blend    3.24 ms   1.94 Mcyc    5%  x32598  60 cyc/blend
        db_warp_orient    79.8 us   47.9 kcyc    0%  x10.0   8.0 us/call
        db_displace       94.3 us   56.6 kcyc    0%  x10.0   9.4 us/call
      db_mesh_copy        0.8 us      482 cyc    0%  x1.0
  canvas_clear           85.8 us    51.5 kcyc    0%  x1.0
  canvas_buffer_wait     38.11 ms  22.86 Mcyc   61%  x1.0
```

Wall min/avg/max = 44.98/62.70/79.90 ms. Icosidodecahedron at 10 copies: 600
edges/frame against preset 1's 864, and 32,598 blended px against 52,221. The
per-call `db_mesh_plot` cost is *higher* than preset 1's (2,423 vs 2,132 µs)
because each call plots a 60-edge solid rather than a 48-edge one; the shape of
the frame is otherwise identical. `db_warp_orient` and `db_displace` scale
exactly with vertex count (30 vs 24 vertices ⇒ 8.0 vs 6.4 µs/call, 1.25× for a
1.25× vertex ratio), confirming they are pure per-vertex passes and not a
scaling risk at any preset.

### Preset 5 hold — snubCube, the preset new since `446f5bd0` (window frames 1489–1504)

```
frame                    62.53 ms  37.52 Mcyc  100%
  db_timeline_step        9.69 ms   5.82 Mcyc   16%
    db_draw               9.59 ms   5.76 Mcyc   15%
      db_draw_scene       9.59 ms   5.76 Mcyc   15%  x1.0
        db_mesh_plot      9.54 ms   5.72 Mcyc   15%  x4.0  2384 us/call
          filter_blend    1.27 ms  761.0 kcyc    2%  x12754  60 cyc/blend
        db_warp_orient    25.3 us   15.2 kcyc    0%  x4.0    6.3 us/call
        db_displace       29.7 us   17.8 kcyc    0%  x4.0    7.4 us/call
      db_mesh_copy        0.8 us      484 cyc    0%  x1.0
  canvas_clear           85.1 us    51.1 kcyc    0%  x1.0
  canvas_buffer_wait     52.75 ms  31.65 Mcyc   84%  x1.0
```

Wall min/avg/max = 61.31/62.53/63.62 ms. The snub cube (V/E/F = 24/60/38) is the
**cheapest preset in the cycle** and cannot touch the peak: its preset entry asks
for 4.534 copies, which floors to **4**, so it plots 4 × 60 = 240 edges/frame
against preset 1's 864. Per *call* it is mid-pack (2,384 µs, between preset 4's
2,423 and preset 1's 2,132, as its 60-edge topology predicts) — the low total is
entirely the copy count. Its worst frame of the whole pass renders in 10.53 ms,
17% of a display window, and 84% of every frame here is display-sync idle.
Adding this preset therefore lowers the pass-wide mix without changing the
cadence-setting phase.

### Preset 2 trough (window frames 2033–2048)

```
frame                    62.46 ms  37.47 Mcyc  100%
  db_timeline_step        1.09 ms  652.6 kcyc    2%
    db_draw               1.01 ms  604.4 kcyc    2%
      db_draw_scene       1.01 ms  603.6 kcyc    2%  x1.0
        db_mesh_plot     798.8 us  479.3 kcyc    1%  x6.0  133.1 us/call
          filter_blend     0.7 us      444 cyc    0%  x7.8   57 cyc/blend
        db_warp_orient    96.6 us   57.9 kcyc    0%  x6.0   16.1 us/call
        db_displace      109.4 us   65.7 kcyc    0%  x6.0   18.2 us/call
      db_mesh_copy        1.1 us      712 cyc    0%  x1.0
  canvas_clear           84.9 us    50.9 kcyc    0%  x1.0
  canvas_buffer_wait     61.28 ms  36.77 Mcyc   98%  x1.0
```

Wall min/avg/max = 62.23/62.46/62.57 ms. The rhombicosidodecahedron cluster has
tumbled almost entirely out of the rendered quadrant: **7.8 blended px/frame**
against the peak's 52,221, a 6,700× swing, and `db_mesh_plot` collapses to 1% of
the frame. What survives is the fixed per-copy vertex work — `db_displace` and
`db_warp_orient` together now cost 206 µs, four orders of magnitude more than
`filter_blend`'s 0.7 µs, because they run over all 60 vertices whether or not a
single pixel is drawn. Wall is the flattest of any window in the pass (0.34 ms
spread) since the whole frame is display-sync idle.

### Sprite hand-off (window frames 305–320, preset 1 → 2)

```
frame                    60.33 ms  36.20 Mcyc  100%
  db_timeline_step       30.31 ms  18.19 Mcyc   50%
    db_draw              30.21 ms  18.13 Mcyc   50%
      db_draw_scene      30.21 ms  18.12 Mcyc   50%  x0.94
        db_mesh_plot     29.97 ms  17.98 Mcyc   50%  x16.9  1776 us/call
          filter_blend    4.06 ms   2.44 Mcyc    7%  x41028  59 cyc/blend
        db_warp_orient   106.3 us   63.8 kcyc    0%  x16.9   6.3 us/call
        db_displace      128.0 us   76.8 kcyc    0%  x16.9   7.6 us/call
      db_mesh_copy        0.9 us      589 cyc    0%  x0.94
  canvas_clear           86.6 us    51.9 kcyc    0%  x1.0
  canvas_buffer_wait     29.94 ms  17.96 Mcyc   50%  x1.0
```

Wall min/avg/max = 31.61/60.33/65.26 ms. The hand-off is a *cheaper* frame, not
a more expensive one: `db_draw_scene` is entered 15 times in 16 frames — one
frame draws no scene at all while the outgoing sprite retires and the incoming
one is selected — so there is no double-draw regime to budget for. The capture
contains exactly ten such windows, one per preset advance, each the last window
of a 320-frame segment (305–320, 625–640, 945–960, 1265–1280, 1585–1600,
1905–1920, 2225–2240, 2545–2560, 2865–2880, 3185–3200). The cost shows up in
wall, not render: the skipped frame lands on the far side of a display flip,
stretching a single frame's wall to 65.26 ms here. Because `spilled` counts
frames whose *render* exceeds 62.5 ms, no hand-off frame spills.

### Per-preset table

All **five** presets visited; the cycle wrapped to preset 1 **twice**, at frames
1601 and 3201 (`Preset: 1/5` markers, ten advance markers in total). Windows per
preset: 53 for preset 1 (three visits, the third truncated by the end of
capture), 40 each for presets 2/3/4/5. **Method note on the clean-hold gate:**
`parse_profile.py`'s default gate is `db_timeline_step` calls/frame, which is 1.0
in every window and therefore excludes nothing here. The load-bearing criterion
for DreamBalls is `db_draw_scene` calls == window frames; that excludes exactly
the ten hand-off windows, leaving 51 clean holds for preset 1 and 38 each for
presets 2/3/4/5. Rows are the worst clean hold of each preset, ranked by
`db_mesh_plot`. The `peak frame render` column is the bucket figure — the worst
single frame the preset owns, hand-off frames included — so it can come from a
different window than the rest of the row.

| # | Preset (solid) | V/E/F | copies | worst hold | blended px/f | `db_mesh_plot` ms | render ms | peak frame render ms | fps |
|---|---|---|--:|---|--:|--:|--:|--:|--:|
| 1 | rhombicuboctahedron | 24/48/26 | 18 | 3329–3344 | 57,162 | 41.91 | 42.34 | 44.88 (193–208) | 16 |
| 4 | icosidodecahedron | 30/60/32 | 10 | 1105–1120 | 32,598 | 24.23 | 24.59 | 32.82 (1105–1120) | 16 |
| 2 | rhombicosidodecahedron | 60/120/62 | 6 | 1937–1952 | 24,918 | 18.97 | 19.37 | 29.26 (2113–2128) | 16 |
| 3 | truncatedCuboctahedron | 48/72/26 | 6 | 673–688 | 15,751 | 12.19 | 12.54 | 15.51 (689–704) | 16 |
| 5 | snubCube | 24/60/38 | 4 | 3121–3136 | 12,810 | 9.61 | 9.85 | 10.53 (1489–1504) | 16 |

Cost tracks edges-plotted × visible coverage, not solid size: preset 2 carries
the largest solid (120 edges) but only 6 copies and a tight 0.05 orbit radius,
so it plots 720 edges against preset 1's 864 and costs less than half as much;
preset 5's snub cube is a mid-sized solid (60 edges) whose 4-copy budget makes it
the cheapest row by a wide margin. Each preset's own floor is far below its
ceiling — the minimum clean-hold `db_timeline_step` is 25.65 / 1.09 / 1.14 /
11.59 / 9.00 ms/f for presets 1/2/3/4/5 — which is why a window mean would
mis-state any of these rows. Presets 2 and 3 are the extreme cases, spanning
17.7× and 10.9× between their own floor and ceiling as the cluster tumbles out
of the rendered band — their clean-hold coverage ranges 8–24,918 and
540–15,751 blended px/f respectively. snubCube is the opposite: it sits in a
9.00–9.76 ms band on 12,204–12,904 blended px/f, a 1.06× coverage spread, so its
4-copy cluster stays essentially wholly inside the quadrant for the entire hold.
Presets 1 and 4 are intermediate at roughly 2× coverage spread each.

`parse_profile.py buckets` prints **six** rows, not five: bucket 0 is the initial
hold of preset 1, which runs from frame 1 and so is never preceded by a
`Preset:` marker, and bucket 1 is preset 1's two later visits. Both are
rhombicuboctahedron, and taken together they are the 848-frame preset-1 row
above. The README cell counts **five** presets.

### Per-pixel figures

At the pass peak (window 193–208) DreamBalls blends **52,221 px/frame** against
the 10,368 px quadrant = **5.04× coverage**, the overdraw of 18 alpha-blended
wireframe copies stacked on the same band. `filter_blend` costs **59.67
cyc/blend**. `db_mesh_plot` spends 368.35 Mcyc on those 835,533 blends (16
frames) = **440.8 cyc per blended pixel**, of which 59.7 is the blend itself and
**381.1 is the plot walk** — geodesic setup, arc stepping, clip splitting and
shading. Pass-wide the figures are 21,875 blended px/frame (2.11× the quadrant)
at 59.52 cyc/blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1154.3/frame  min/avg/max 0.68/1.91/20.97 us  cpu 3.53%
isr_pack         143.9/frame  min/avg/max 6.23/6.82/9.31 us   cpu 1.57%
isr_dma_submit   143.9/frame  min/avg/max 0.60/0.93/3.44 us   cpu 0.21%
```

- Submitting a column to the DMA engine costs 0.93 µs against 6.82 µs to pack
  it — the hand-off is about 1/7 of the packing, so the marshaling cost is
  dominated by assembling the pixel bytes, not by touching eDMA.
- 143.9 packs/frame is exactly one per displayed column (288 columns/rev, 144
  per 62.5 ms half-rev), a 434 µs column period. At 24 MHz SPI, one HD107S frame
  for 72 LEDs is (1 + 72) × 32 bits plus the end frame ≈ 2.4 kbit ≈ **99 µs on
  the wire**, all asynchronous; the CPU pays only the 7.8 µs of pack + submit out
  of that 434 µs.
- Total ISR CPU share is **5.31%** pass-wide (5.47% inside the peak window,
  where `isr_wake` runs slightly hotter). It is driver-side and cadence-fixed.
  Every render scope runs with interrupts live and therefore already contains it,
  so the peak render of 44.88 ms is directly comparable to the 62.5 ms window:
  **17.62 ms of margin, 1.39× headroom, no phase needs a speedup.**

## Summary ranking

1. `db_mesh_plot` — **62% of the frame, 38.37 ms** at the peak: the wireframe
   raster is the effect. Splits into the geodesic walk (54%, 33.18 ms) and
   `filter_blend` (8%, 5.19 ms).
2. `canvas_buffer_wait` — 37%, 22.94 ms at the peak and 45.54 ms/f pass-wide:
   pure display-sync idle.
3. `db_displace` — 0.23%, 139.7 µs: per-vertex tangent-frame displacement,
   linear in vertex count.
4. `db_warp_orient` — 0.19%, 114.5 µs: per-vertex Möbius warp + orientation.
5. `canvas_clear` — 0.14%, 87.7 µs.
6. `db_mesh_copy` — 0.001%, 0.8 µs.

No WASM/native counterpart figures exist for DreamBalls;
`docs/dreamballs_perf_analysis.md` is a device-only ledger (COMPLETE 2026-07-20)
and ends at peak frame render 73.19 ms with 6.8% spill. The trajectory since is
73.19 → 57.98 (2026-07-25 report) → 44.44 (`446f5bd0`) → **44.88** here.

### Change vs the `446f5bd0` capture

Eighty commits separate `446f5bd0` from `e68a4e37`, of which two were flagged as
plausibly moving this effect. Both are settled below. The comparison is
window-for-window valid because DreamBalls is frame-deterministic and the two
captures share their first 1,280 frames of preset schedule (presets 1–4 at 320
frames each, snubCube entering only at frame 1281): the blended-pixel counts of
the shared windows are **identical to the unit**, so any timing delta is timing
and nothing else.

| window (same animation state) | `446f5bd0` | `e68a4e37` | Δ |
|---|--:|--:|--:|
| 193–208 blended px/f | 52,221 | 52,221 | **identical** |
| 193–208 peak frame render | 44.44 ms | 44.88 ms | +1.0% |
| 113–128 blended px/f | 55,321 | 55,321 | **identical** |
| 113–128 `db_mesh_plot` | 40.13 ms/f | 40.33 ms/f | +0.5% |
| 113–128 render | 40.56 ms/f | 40.66 ms/f | +0.2% |
| 1105–1120 `db_mesh_plot` | 24.08 ms/f | 24.23 ms/f | +0.6% |
| 1105–1120 peak frame render | 32.49 ms | 32.82 ms | +1.0% |
| 673–688 `db_mesh_plot` | 12.14 ms/f | 12.19 ms/f | +0.4% |
| ISR CPU share | 5.32% | 5.31% | — |
| spilled | 0/2,688 (0%) | 0/3,408 (0%) | — |

**`e68a4e37` did not move DreamBalls measurably.** The commit is reached from
this effect — `Plot::Mesh::draw_edge` calls `rasterize<W, H>` twice
(core/render/plot.h), once with `{.edge_visible = bits}` and once with the
default `{}` — so the by-value `RasterOptions` parameter is in DreamBalls' hot
path by construction, not by accident. But every comparable window lands within
±1%, which is run-to-run variation at this scale, and the sign is if anything
adverse. What the commit bought here is ITCM, not time: phantasm RAM1 code
193,704 B with **2,904 B** of padding headroom, size gate PASS. DreamBalls'
`db_mesh_plot` is walk-bound, not call-overhead-bound — 381 of its 441 cycles per
blended pixel are inside the geodesic walk — so removing a 32-byte argument
materialization from a call made 288 times a frame is below this effect's noise
floor even when it saves real instructions elsewhere in the image.

**`fbaf81e2` cannot reach DreamBalls**, confirmed rather than assumed: neither
`effects/DreamBalls.h` nor `core/render/plot.h` includes
`core/math/spherical_field.h`, and DreamBalls runs no feedback compositor — the
`sample_bilinear_rgb` inlining that commit restored is a MeshFeedback flush-path
concern only.

**What did change is the mix, not the cadence.** snubCube is the cheapest preset
in the cycle (9.61 ms `db_mesh_plot`, 10.53 ms peak render), so extending the
cycle from four presets to five pulls every pass-wide average down without
touching the phase that sets cadence:

| pass-wide | `446f5bd0` (4 presets) | `e68a4e37` (5 presets) |
|---|--:|--:|
| `db_mesh_plot` | 17.83 ms/f | 16.53 ms/f |
| render | 18.22 ms/f | 16.89 ms/f |
| blended px/f | 23,598 | 21,875 |
| `canvas_buffer_wait` | 44.21 ms/f | 45.54 ms/f |
| wall | 62.43 ms/f | 62.43 ms/f |
| peak frame render | 44.44 ms | 44.88 ms |

These averages are **not** a like-for-like comparison and must not be read as an
improvement: they moved because a cheap fifth preset now occupies a fifth of the
cycle. The peak, the spill count and the per-preset rows are the figures that
survive the preset-count change, and all three are flat.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs) — that is the point of profiling
  under the live driver. The 5.31% ISR share is inside every figure above, so
  render times are directly comparable to the 62.5 ms window.
- `filter_blend` parents under `db_mesh_plot`; its subtree would be hidden in a
  window where that parent had 0 calls, which never happens in this capture
  (every window plots at least one mesh). Its calls ≈ blended pixels.
- Selective-O3: `db_mesh_plot`'s inner walk crosses the `HS_O3` regions
  `Plot::rasterize`, `rasterize_geodesic_strategy` and `edge_fits_one_dot` in
  core/render/plot.h, plus filter.h's blend sink. `Plot::Mesh::draw_edge`, which
  builds the `RasterOptions` argument `e68a4e37` changed, sits **outside** every
  region and so compiles at `-Os` on the shipping image.
- No dwell-compression knobs were used (no `HS_PROFILE_ORDERED_CYCLE`, no
  `HS_PROFILE_TRANS_SPEED`). `HS_PROFILE_EPOCH_REVS=2000` stretches the epoch to
  250 s so the 215 s capture never crosses an epoch boundary; it changes how long
  one effect instance lives, never its per-frame cost.
- Working tree was clean at `e68a4e37` — the run's artifact bundle recorded an
  empty `source.diff` and an empty `source_status.txt`, and the provenance file
  records `source_sha=e68a4e37ccdd5e8dad9c585032b68e26e2be0e22`.
- The capture ran on **COM3** while a peer session captured MeshFeedback on
  COM4. The two boards hold independent locks; this report's board is COM3
  throughout, and `validate` confirms a fresh boot (first frame 1, `reset cause:
  0x001 por`) with no window-name or marker mismatch.
- The four presets common to both captures keep their identity and their frame
  ranges across the preset-count change because snubCube was appended at the end
  of the table, not inserted. Only frames past 1280 diverge between the two
  captures.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DreamBalls`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=2000`; `just profile DreamBalls`
= build + flash + capture. Parsing needs no device:
`python tools/parse_profile.py <log> windows|presets|buckets|validate`.
