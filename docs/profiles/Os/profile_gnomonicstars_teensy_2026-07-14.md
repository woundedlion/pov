# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_gnomonicstars.log` (64-frame windows, single ~75 s pass).
The column-ISR/DMA accumulators are dumped every window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `GnomonicStars<288, 144>` only (single-entry playlist), current working-tree state (default presets) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile GnomonicStars` |

Image size: FLASH code 41,524 B; ITCM (RAM1 code) 27,400 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to a few
ppm: window frames 641–704 root counter = 2,401,468,221 cyc = 4,002,447.0 µs
vs measured `micros()` window sum 4,002,459 µs (Δ ≈ 3.0 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). GnomonicStars
scatters ~600 fib-spiral point stars into that quadrant and rasterizes each as
a tiny `Scan::Star` sprite, so it touches only the pixels the sprites cover
(here ~3.4 k blended px/frame, ~33 % of the quadrant) — not a full-frame
shader. Wall time snaps up to a whole number of 62.5 ms windows.

The effect is LIGHT–MEDIUM: render (`gn_draw_stars` + `gn_timeline_step`) is
~17 ms/frame — it breathes between ~14 and ~21 ms across the capture as the
star field drifts and the per-sprite coverage changes — but it never comes near
the 62.5 ms boundary, so the cadence is a **steady 16 fps (62.5 ms/frame)** for
the whole run. `gn_buffer_wait` (the `buffer_free()` spin inside the `Canvas`
ctor waiting for the display flip) is timed separately and is exactly the
round-up idle: ~45 ms/frame, 64–76 % of the window.

## Phase-by-phase readout (64-frame windows, per-frame averages)

GnomonicStars has no discrete phases — it holds one steady regime and only
breathes gently with the drifting point field. One representative steady window
below (frames 641–704); every window in the capture is the same 16 fps shape,
`gn_draw_stars` 23–30 % of frame (34 % in the two densest windows). Each parent
includes its children. Columns: time/frame, cycles/frame, % of frame; leaves
add calls/frame and per-call cost (cycles = µs × 600):

### Steady star field (representative window: frames 641–704)

```
frame                62.54 ms  37.52 Mcyc  100%
  gn_draw_stars      17.07 ms  10.24 Mcyc   27%
    gn_star_scan     15.97 ms   9.58 Mcyc   26%  x600   26.6 us/star
      filter_blend    1.31 ms   0.79 Mcyc    2%  x3406  231 cyc/blend
  gn_timeline_step     38 us    22.8 kcyc    0%
  gn_buffer_wait     45.43 ms  27.26 Mcyc   73%
```

`gn_star_scan` is the per-star `Scan::Star::draw` rasterize (render leaf);
`gn_draw_stars` is the parent over the ~600-star loop, its **~1.10 ms/frame
self-time** being the `transformer.transform` + `make_basis` per star.
`gn_buffer_wait` is the display-window sync idle.

Wall: min 50.2 / avg 62.5 / max 75.0 ms — a clean single-window cadence. The
render is entirely the point-sprite rasterize: 600 stars/frame at 26.6 µs each,
of which the blend composite (`filter_blend`, 1.31 ms) is only ~8 % — the rest
is the per-star basis + `Scan::Star` SDF distance evaluation over each sprite's
small bounding box. `gn_spiral_build` (the fib-spiral cache rebuild) fires
**exactly once in the whole capture** (window 1: 1 call, 885 µs) — by design it
only runs on a "Points" slider change, so it is ~0 % of steady cost and does
not appear in the tree above.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h): **~3,406 blended
px/frame** here vs the 10,368-px quadrant — the star field covers only ~33 % of
the segment (sparse sprites, not a full stack), and this is the figure that
breathes (2,665 px/frame in the sparsest window, 4,946 in the densest). Cost is
**231 cyc (0.38 µs) per blend** (includes its own scope overhead). Against the
scan work, `gn_star_scan` ≈ 9.58 Mcyc/frame over ~3,406 blended px ⇒ **~2,814
cyc per *blended* pixel**, i.e. most of the scan cost is the per-star basis
build + `Scan::Star` distance/coverage evaluation, not the blend itself
(~5.7 blended px land per star sprite).

## Column-ISR / DMA marshaling cost (`build/profile_gnomonicstars.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative window (frames
641–704, 4.00 s). The pack and submit run inside the flywheel wake, on the
1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max, CPU
share:

```
isr_wake         18434/s  0.52 / 1.53 / ~91 us   cpu 2.82%
  isr_pack        2304/s  4.37 / 5.18 / ~89 us   cpu 1.19%
  isr_dma_submit  2304/s  0.61 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, and
includes pack + submit nested); `isr_pack` = the 72× `packPixel` marshal, once
per column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush +
eDMA kick).

- **The DMA kick itself is ~1 µs of CPU** (submit avg 0.96 µs); the pixel-pack
  marshal around it is the real per-column cost. Whole per-column CPU ≈
  **~12 µs of the 434 µs column period (~2.8 %)**. The wire transfer runs
  asynchronously on eDMA and costs no CPU.
- The avg/max figures inflate from serial-dump preemption (higher-priority
  USB/DMA-completion ISRs landing inside the scope, plus the dumps themselves)
  — the ~89–91 µs maxima are artifacts, not real ISR growth. The `min` floor is
  the stable per-call cost: pack 4.37 µs (2621 cyc), submit 0.61 µs.
- Net: the ISR machinery steals **~2.82 % of the chip**. Against the 62.5 ms
  window that's ~1.76 ms, leaving **~60.7 ms of render budget per window**. The
  ~17 ms render fits one window with ~3.5× headroom — hence the rock-solid
  16 fps.

## Summary ranking (steady window, share of the 62.5 ms frame)

1. display-window sync (`gn_buffer_wait`) — **73%** (45.4 ms, idle by design —
   the render finishes in ~27 % of the window and spins out the rest).
2. `Scan::Star::draw` rasterize (`gn_star_scan`, incl. `filter_blend`) —
   **26%** (16.0 ms) — the entire render cost: 600 sprites/frame.
3. per-star `transformer.transform` + `make_basis` (`gn_draw_stars` self) —
   **2%** (1.1 ms).
4. everything else (`gn_timeline_step`; `gn_spiral_build` fires ~once/run) —
   <1%.

GnomonicStars is a light point-sprite rasterizer: with render at ~17 ms against
a ~61 ms budget it is nowhere near window-bound, and the only real work is the
600-star `Scan::Star` loop. The levers (all low-value at this headroom) are the
per-star basis/transform and the sprite distance evaluation; there is no LUT,
cull, or geometry stage to attack. No WASM/native GnomonicStars figures are
recorded in the perf ledger for comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it parents under `gn_star_scan` (the first scope
  to enter it) and counts every blended pixel across the star loop; its printed
  % is meaningful here because that is the only blend path, but the subtree
  would move if another scope entered blend first. Per-pixel `filter_blend`
  scope overhead inflates the scan counter by a fraction of a percent.
- Epoch/reset windows (a `frame` count > 64 or a very low first-frame `min`)
  mix an undumped tail into the next window — the first window (min 21,779 µs,
  epoch-start) is skipped when picking the representative; per-call averages
  stay valid there, per-frame sums do not.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/GnomonicStars.h` (`gn_buffer_wait`, `gn_timeline_step`,
  `gn_spiral_build`, `gn_draw_stars`, `gn_star_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=GnomonicStars`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel
  pack, DMA submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile GnomonicStars [seconds]`.
