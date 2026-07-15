# Comets on-device profile — Teensy 4.0, segmented mode, -O3 (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_comets_o3.log` (64-frame windows, single ~75 s pass, one epoch).
The column-ISR/DMA accumulators are dumped every window. This is the **-O3**
twin of the shipping `-Os` profile; the only variable
between the two is the optimization level.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm shipping flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but built at **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` — the -O3 twin of the shipping `-Os` profile |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Comets<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 77,700 B; ITCM (RAM1 code) 58,840 B; RAM2 free 4,736 B.
(-Os was 49,460 / 33,320 / 4,736.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-2-ppm
across the steady windows: frames 449–512 root counter = 2,403,292,339 cyc =
4,005,487.23 µs vs measured `micros()` window sum 4,005,491 µs (Δ ≈ 0.9 ppm).
Two more clean windows confirm it: frames 65–128 (Δ ≈ 0.67 ppm) and frames
577–640 (Δ ≈ 1.31 ppm) — so at least three clean steady windows check out.

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Wall time snaps
up to a whole number of 62.5 ms windows.

Comets is a **light effect: it holds 16 fps (62.5 ms/frame, 1-window) across
the entire pass** — every window's wall avg sits at 62.4–62.6 ms and the render
never overflows one window, at -O3 with even more headroom than at -Os. There
are no discrete phases and no cadence flips; the render cost **breathes** with
how many comet points are live and in this segment's band. Over the run it ramps
from a near-idle start, peaks around the middle (mid-run frames 577–640), and
settles back as trails fade:

- **Peak trail density** (many comet points in band): render `cm_draw_trail`
  ≈ 21.6 ms → only ~34% of the frame, deep inside one window.
- **Light windows** (few points live): render ≈ 6.5 ms → ~10% of the frame.

`cm_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is the round-up idle and is timed separately — it fills the
balance of every 62.5 ms window (65% when render peaks, up to 90% in the light
windows), so the render numbers below are clean.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Cadence is fixed at 16 fps; only the render load moves. The two windows below
bound the breathe, with a monotonic buildup/decay between them. Columns:
time/frame, cycles/frame, % of frame; leaves add xN calls/frame and per-call
cost (cycles = µs × 600).

### Peak trail density — representative mid window (frames 577–640)

```
frame                62.50 ms  37.50 Mcyc  100%
  cm_draw_trail      21.56 ms  12.94 Mcyc   34%
    cm_point_scan    21.14 ms  12.68 Mcyc   34%  x1713  12.3 us/call
      filter_blend    1.17 ms   0.70 Mcyc    2%  x8198  85.5 cyc/blend
  cm_buffer_wait     40.85 ms  24.51 Mcyc   65%
  cm_timeline_step    0.08 ms  50.2 kcyc    0%
  cm_wipe_rebake      0.03 us     21 cyc    0%
```

Wall: min 57.1 / avg 62.5 / max 67.5 ms. `cm_point_scan`
(`Scan::Point::draw`, one call per point-subsample) is 98% of `cm_draw_trail`
and the entire algorithmic cost — 1,713 scan calls/frame writing 8,198 blended
pixels. `cm_draw_trail` self time (the trail `deep_tween` walk that dispatches
the scans) is ~0.42 ms/frame. This is the busiest window in the pass; at -O3 the
render sits at ~34% of the display window, so cadence stays 16 fps with wide
headroom.

### Light window — few comets live (frames 65–128)

```
frame                62.55 ms  37.53 Mcyc  100%
  cm_buffer_wait     56.03 ms  33.62 Mcyc   90%
  cm_draw_trail       6.47 ms   3.88 Mcyc   10%
    cm_point_scan     6.35 ms   3.81 Mcyc   10%  x463  13.7 us/call
      filter_blend    0.32 ms   0.19 Mcyc    1%  x2209  86.0 cyc/blend
  cm_timeline_step    0.05 ms  28.1 kcyc    0%
  cm_wipe_rebake      0.03 us     21 cyc    0%
```

Wall: min 58.9 / avg 62.5 / max 66.0 ms. Same 16 fps window, but only ~463
scan calls/frame (2,209 blended px), so `cm_buffer_wait` absorbs 90% of the
frame as idle. Per-scan-call cost is close to the peak window (13.7 vs 12.3 µs;
the light window has a slightly higher per-call cost from fixed setup amortized
over fewer subsamples) — the effect still scales essentially by live-point
count, not by a per-point cost change.

`cm_wipe_rebake` (palette-LUT rebake while a `ColorWipe` transition is in
flight) is 0 in most windows and peaks to only **0.42 ms/frame (~1%)** in the
windows that catch a wipe (e.g. frames 321–384, 641–704) — negligible either
way. `cm_timeline_step` (animation advance) is <0.12 ms/frame throughout.

### Per-pixel figures

`filter_blend` (the pre-existing per-pixel counter in filter.h) parents
correctly under `cm_point_scan` here — Comets' only blended-pixel path is the
point scan, so the subtree is valid, not an artifact. Blended coverage tracks
live-point density: **8,198 blended px/frame (79% of the 10,368-px quadrant)**
at the peak, **2,209 px (21%)** in the light window. The blend itself is a
rock-steady **85.4–86.0 cyc (0.143 µs) per blend** in every window (down from
226 cyc at -Os — the palette blend math is what -O3 + `-ffast-math` speeds up
most). Against the scan, `cm_point_scan` ≈ 12.68 Mcyc over 8,198 blended px ⇒
**~1,547 cyc (2.58 µs) per *blended* pixel** — so the actual palette blend is
only ~5.5% of the per-blended-pixel cost; the rest is the point SDF distance +
AA-subsample setup (~4.8 blended px per scan call).

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative peak window
(frames 577–640, 4.0 s). The pack and submit run inside the flywheel wake, on
the 1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max,
CPU share:

```
isr_wake         18433/s  0.44 / 1.49 / ~122 us  cpu 2.74%
  isr_pack        2304/s  4.64 / 5.65 / ~121 us  cpu 1.30%
  isr_dma_submit  2304/s  0.60 / 0.90 / 3.0 us   cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~0.6 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **~5.3 µs of the 434 µs column
  period (~1.2%)**. The wire transfer runs asynchronously on eDMA and costs no
  CPU.
- The `min` floors (wake 0.44 µs / 266 cyc, pack 4.64 µs / 2782 cyc, submit
  0.6 µs) are the stable per-call cost; the avg/max inflate in the busier
  windows (wake avg ~1.49 µs, max ~122 µs) from serial-dump preemption landing
  inside the scope, not real ISR growth. (The pack `min` floor is slightly
  higher at -O3, 2782 vs 2551 cyc — the ISR code is compiled at the image's opt
  level too, and here -O3 inlines the marshal marginally larger.)
- Net: the ISR machinery steals **~2.5–2.7% of the chip** (isr_wake cpu% is the
  total, pack+submit nested). Against the 62.5 ms window that is ~1.7 ms,
  leaving **~60.8 ms of render budget per window**. The peak render (21.6 ms)
  needs only ~0.35× of that budget — which is why Comets never leaves 16 fps.

## Summary ranking (peak window, share of the 62.5 ms frame)

1. `Scan::Point::draw` (`cm_point_scan` — per-point SDF rasterize + shade +
   blend) — **34%** (21.1 ms) — the entire algorithmic cost.
2. display-window sync (`cm_buffer_wait`) — **65%** (40.9 ms, idle by design;
   grows to 90% in the light windows).
3. `cm_draw_trail` self time (trail `deep_tween` walk / scan dispatch) — ~1%
   (0.42 ms).
4. everything else (`cm_timeline_step` animation advance, `cm_wipe_rebake`
   palette rebake) — <1% each (wipe peaks to ~1% only mid-transition).

Comets is point-scan-bound whenever it is doing real work: the only render
lever is the per-point SDF/AA cost (12.3 µs/scan-call at -O3) or the number of
subsample scans per point. Because the whole pass fits one display window with
~3× headroom at -O3, there is no cadence benefit to optimizing it — the win
would be CPU headroom, not frame rate. No WASM/native Comets figures are
recorded in the perf ledger for comparison.

## -O3 vs -Os

Both runs hold 16 fps (62.5 ms/frame) and never overflow a display window, so
**cadence is unchanged** — -O3 buys CPU headroom, not frame rate. The peak
window (frames 577–640) has an identical 1,713 scan calls/frame in both builds,
so the per-call comparison is apples-to-apples:

| metric (peak window 577–640) | -Os | -O3 | Δ |
|---|---|---|---|
| render `cm_draw_trail` | 33.32 ms | 21.56 ms | **1.55× faster** (−11.8 ms) |
| `cm_point_scan` per-call | 19.0 µs (x1713) | 12.3 µs (x1713) | 1.54× |
| per *blended* pixel | ~2,223 cyc | ~1,547 cyc | 1.44× |
| `filter_blend` per blend | 226 cyc | 85.5 cyc | **2.64×** |
| frame cadence | 16 fps (62.5 ms) | 16 fps (62.5 ms) | unchanged |
| render share of frame (peak) | 53% | 34% | −19 pts headroom |
| FLASH code | 49,460 B | 77,700 B | +28,240 (+57%) |
| ITCM (RAM1 code) | 33,320 B | 58,840 B | +25,520 (+77%) |

The point SDF + AA scan is ~1.55× faster and the per-pixel palette blend is the
biggest single winner at **2.64×** (the blend math is exactly what -O3 +
`-ffast-math` simplifies). The cost is size: **+28 KB FLASH and +26 KB ITCM**.
**-O3 is not the shipping config** — the full 26-effect Phantasm image overflows
FlexRAM/ITCM at -O3 (that is why the ship build is -Os); -O3 is only viable here
as a single-effect profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: here it parents under `cm_point_scan` (Comets'
  only blend path), so its percentage is meaningful; its per-pixel scope
  overhead inflates the scan counter by ~0.2% of frame at -O3.
- Epoch-start window (`frames 1-64`, `frame` min 1,310 µs) mixes the reconstruct
  ramp into its first window — skipped when picking representatives; per-call
  averages there stay valid, per-frame sums don't.
- **-O3 build**: this is the `profile_o3` env (base `-O3`, `-ffast-math
  -fno-finite-math-only`), which does **not** ship — the shipping Phantasm image
  is `-Os` because the full 26-effect roster overflows ITCM at -O3. Use these
  numbers only as a single-effect ceiling, not as the on-ship cost.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Comets.h`
  (`cm_buffer_wait`, `cm_timeline_step`, `cm_wipe_rebake`, `cm_draw_trail`,
  `cm_point_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Comets`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (the
  `-Os` twin uses `just profile Comets [seconds]`).
