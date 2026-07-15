# Comets on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_comets.log` (64-frame windows, single ~75 s pass, one epoch).
The column-ISR/DMA accumulators are dumped every window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Comets<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile Comets` |

Image size: FLASH code 49,460 B; ITCM (RAM1 code) 33,320 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
steady window frames 449–512 root counter = 2,407,046,074 cyc = 4,011,743.46 µs
vs measured `micros()` window sum 4,011,747 µs (Δ ≈ 0.9 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Wall time snaps
up to a whole number of 62.5 ms windows.

Comets is a **light effect: it holds 16 fps (62.5 ms/frame, 1-window) across
the entire pass** — every window's wall avg sits at 62.4–62.7 ms and the render
never overflows one window. There are no discrete phases and no cadence flips;
instead the render cost **breathes** with how many comet points are live and in
this segment's band. Over the run it ramps from a near-idle start, peaks around
the middle (mid-run frames 577–704), and settles back as trails fade:

- **Peak trail density** (many comet points in band): render `cm_draw_trail`
  ≈ 33 ms → still ~53% of the frame, comfortably inside one window.
- **Light windows** (few points live): render ≈ 9 ms → ~15% of the frame.

`cm_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is the round-up idle and is timed separately — it fills the
balance of every 62.5 ms window (46% when render peaks, up to 84–93% in the
light windows), so the render numbers below are clean.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Cadence is fixed at 16 fps; only the render load moves. The two windows below
bound the breathe, with a monotonic buildup/decay between them. Columns:
time/frame, cycles/frame, % of frame; leaves add xN calls/frame and per-call
cost (cycles = µs × 600).

### Peak trail density — representative mid window (frames 577–640)

```
frame                62.59 ms  37.55 Mcyc  100%
  cm_draw_trail      33.32 ms  19.99 Mcyc   53%
    cm_point_scan    32.54 ms  19.52 Mcyc   52%  x1713  19.0 us/call
      filter_blend    3.31 ms   1.98 Mcyc    5%  x8781  226 cyc/blend
  cm_buffer_wait     29.15 ms  17.49 Mcyc   47%
  cm_timeline_step    0.11 ms  68.6 kcyc    0%
  cm_wipe_rebake      0.07 us     41 cyc    0%
```

Wall: min 50.5 / avg 62.6 / max 74.1 ms. `cm_point_scan`
(`Scan::Point::draw`, one call per point-subsample) is 97% of `cm_draw_trail`
and the entire algorithmic cost — 1,713 scan calls/frame writing 8,781 blended
pixels. `cm_draw_trail` self time (the trail `deep_tween` walk that dispatches
the scans) is ~0.8 ms/frame. This is the busiest window in the pass; render
still fits one display window, so cadence stays 16 fps.

### Light window — few comets live (frames 65–128)

```
frame                62.57 ms  37.54 Mcyc  100%
  cm_buffer_wait     53.06 ms  31.84 Mcyc   85%
  cm_draw_trail       9.44 ms   5.67 Mcyc   15%
    cm_point_scan     9.22 ms   5.53 Mcyc   15%  x463  19.9 us/call
      filter_blend    0.83 ms   0.50 Mcyc    1%  x2209  226 cyc/blend
  cm_timeline_step    0.07 ms  42.3 kcyc    0%
  cm_wipe_rebake      0.07 us     41 cyc    0%
```

Wall: min 56.7 / avg 62.6 / max 68.2 ms. Same 16 fps window, but only ~463
scan calls/frame (2,209 blended px), so `cm_buffer_wait` absorbs 85% of the
frame as idle. Per-scan-call cost is essentially identical to the peak window
(19.9 vs 19.0 µs) — the effect scales purely by live-point count, not by any
per-point cost change.

`cm_wipe_rebake` (palette-LUT rebake while a `ColorWipe` transition is in
flight) is 0 in most windows and peaks to only **0.84 ms/frame (~1%)** in the
windows that catch a wipe (e.g. frames 321–384, 641–704) — negligible either
way. `cm_timeline_step` (animation advance) is <0.2 ms/frame throughout.

### Per-pixel figures

`filter_blend` (the pre-existing per-pixel counter in filter.h) parents
correctly under `cm_point_scan` here — Comets' only blended-pixel path is the
point scan, so the subtree is valid, not an artifact. Blended coverage tracks
live-point density: **8,781 blended px/frame (85% of the 10,368-px quadrant)**
at the peak, **2,209 px (21%)** in the light window. The blend itself is a
rock-steady **225.5–225.8 cyc (0.376 µs) per blend** in every window. Against
the scan, `cm_point_scan` ≈ 19.52 Mcyc over 8,781 blended px ⇒ **~2,223 cyc
(3.7 µs) per *blended* pixel** — so the actual palette blend is only ~10% of
the per-blended-pixel cost; the rest is the point SDF distance + AA-subsample
setup (~5.1 blended px per scan call).

## Column-ISR / DMA marshaling cost

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative peak window
(frames 577–640, 4.0 s). The pack and submit run inside the flywheel wake, on
the 1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max,
CPU share:

```
isr_wake         18432/s  0.62 / 1.55 / ~149 us  cpu 2.86%
  isr_pack        2304/s  4.25 / 5.09 / ~147 us  cpu 1.17%
  isr_dma_submit  2304/s  0.63 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **~6 µs of the 434 µs column
  period (~1.3%)**. The wire transfer runs asynchronously on eDMA and costs no
  CPU.
- The `min` floors (wake 0.62 µs / 371 cyc, pack 4.25 µs / 2551 cyc, submit
  0.6 µs) are the stable per-call cost; the avg/max inflate in the busier
  windows (wake avg up to ~1.55 µs, max ~149 µs) from serial-dump preemption
  landing inside the scope, not real ISR growth.
- Net: the ISR machinery steals **~2.5–2.9% of the chip** (isr_wake cpu% is the
  total, pack+submit nested). Against the 62.5 ms window that is ~1.7 ms,
  leaving **~60.7 ms of render budget per window**. Even the peak render
  (33 ms) needs only ~0.55× of that budget — which is why Comets never leaves
  16 fps.

## Summary ranking (peak window, share of the 62.5 ms frame)

1. `Scan::Point::draw` (`cm_point_scan` — per-point SDF rasterize + shade +
   blend) — **52%** (32.5 ms) — the entire algorithmic cost.
2. display-window sync (`cm_buffer_wait`) — **47%** (29.2 ms, idle by design;
   grows to 84–93% in the light windows).
3. `cm_draw_trail` self time (trail `deep_tween` walk / scan dispatch) — ~1%
   (0.8 ms).
4. everything else (`cm_timeline_step` animation advance, `cm_wipe_rebake`
   palette rebake) — <1% each (wipe peaks to ~1% only mid-transition).

Comets is point-scan-bound whenever it is doing real work: the only render
lever is the per-point SDF/AA cost (19 µs/scan-call) or the number of
subsample scans per point. Because the whole pass fits one display window with
~2× headroom, there is no cadence benefit to optimizing it — the win would be
CPU headroom, not frame rate. No WASM/native Comets figures are recorded in
the perf ledger for comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: here it parents under `cm_point_scan` (Comets'
  only blend path), so its percentage is meaningful; its per-pixel scope
  overhead inflates the scan counter by ~0.5% of frame.
- Epoch-start window (`frames 1-64`, `frame` min 1,867 µs) mixes the reconstruct
  ramp into its first window — skipped when picking representatives; per-call
  averages there stay valid, per-frame sums don't.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Comets.h`
  (`cm_buffer_wait`, `cm_timeline_step`, `cm_wipe_rebake`, `cm_draw_trail`,
  `cm_point_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Comets`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile Comets [seconds]`.
