# GnomonicStars on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_gnomonicstars_o3.log` (64-frame windows, single ~75 s pass).
The column-ISR/DMA accumulators are dumped every window. This is the **-O3**
twin of the shipping `-Os` report; the only variable
between the two runs is the optimization level.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but keeping base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `GnomonicStars<288, 144>` only (single-entry playlist), current working-tree state (default presets) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 59,428 B; ITCM (RAM1 code) 44,728 B; RAM2 free 4,736 B.
(The shipping `-Os` build was 41,524 / 27,400 / 4,736 — see `## -O3 vs -Os`.)

**Exactness cross-check** — cycle counter and wall clock agree to a few ppm
across every steady window; three verified: frames 513–576 root counter =
2,400,689,863 cyc = 4,001,149.8 µs vs `micros()` sum 4,001,150 µs (Δ ≈ 0.06 ppm);
frames 321–384 = 4,022,222.8 µs vs 4,022,228 (Δ ≈ 1.3 ppm); frames 641–704 =
4,001,038.3 µs vs 4,001,049 (Δ ≈ 2.7 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). GnomonicStars
scatters ~600 fib-spiral point stars into that quadrant and rasterizes each as
a tiny `Scan::Star` sprite, so it touches only the pixels the sprites cover
(here ~3.4 k blended px/frame, ~33 % of the quadrant) — not a full-frame
shader. Wall time snaps up to a whole number of 62.5 ms windows.

The effect is LIGHT: at -O3 the render (`gn_draw_stars` + `gn_timeline_step`) is
~12–19 ms/frame — it breathes as the star field drifts and the per-sprite
coverage changes — and never comes near the 62.5 ms boundary, so the cadence is
a **steady 16 fps (62.5 ms/frame)** for the whole run. `gn_buffer_wait` (the
`buffer_free()` spin inside the `Canvas` ctor waiting for the display flip) is
timed separately and is exactly the round-up idle: ~46–48 ms/frame, 70–78 % of
the window.

## Phase-by-phase readout (64-frame windows, per-frame averages)

GnomonicStars has no discrete phases — it holds one steady regime and only
breathes gently with the drifting point field. One representative steady window
below (frames 513–576, deliberately chosen because its ~3.39 k blended px/frame
matches the -Os report's representative ~3.41 k, for a like-for-like read).
Every window in the capture is the same 16 fps shape, `gn_draw_stars` 19–29 %
of frame. Each parent includes its children. Columns: time/frame, cycles/frame,
% of frame; leaves add calls/frame and per-call cost (cycles = µs × 600):

### Steady star field (representative window: frames 513–576)

```
frame                62.52 ms  37.51 Mcyc  100%
  gn_draw_stars      14.11 ms   8.47 Mcyc   22%
    gn_star_scan     13.45 ms   8.07 Mcyc   22%  x600   22.4 us/star
      filter_blend    0.54 ms   0.33 Mcyc    1%  x3391  96 cyc/blend
  gn_timeline_step     29 us    17.7 kcyc    0%
  gn_buffer_wait     48.38 ms  29.03 Mcyc   77%
```

`gn_star_scan` is the per-star `Scan::Star::draw` rasterize (render leaf);
`gn_draw_stars` is the parent over the ~600-star loop, its **~0.66 ms/frame
self-time** being the `transformer.transform` + `make_basis` per star (down
from 1.10 ms at -Os). `gn_buffer_wait` is the display-window sync idle.

Wall: min 53.4 / avg 62.5 / max 71.8 ms — a clean single-window cadence. The
render is entirely the point-sprite rasterize: 600 stars/frame at 22.4 µs each,
of which the blend composite (`filter_blend`, 0.54 ms) is only ~4 % — the rest
is the per-star basis + `Scan::Star` SDF distance evaluation over each sprite's
small bounding box. `gn_spiral_build` (the fib-spiral cache rebuild) fires
**exactly once in the whole capture** (window 1: 1 call, 875 µs) — by design it
only runs on a "Points" slider change, so it is ~0 % of steady cost and does
not appear in the tree above.

### Per-pixel figures

`filter_blend` (pre-existing per-pixel counter in filter.h): **~3,391 blended
px/frame** here vs the 10,368-px quadrant — the star field covers only ~33 % of
the segment (sparse sprites, not a full stack), and this is the figure that
breathes (2,665 px/frame in the sparsest window, 5,204 in the densest). Cost is
**96 cyc (0.16 µs) per blend** at -O3 (includes its own scope overhead) — versus
231 cyc at -Os, a 2.4× win on the composite. Against the scan work,
`gn_star_scan` ≈ 8.07 Mcyc/frame over ~3,391 blended px ⇒ **~2,381 cyc per
*blended* pixel** (was ~2,814 at -Os), i.e. most of the scan cost is still the
per-star basis build + `Scan::Star` distance/coverage evaluation, not the blend
itself (~5.7 blended px land per star sprite).

## Column-ISR / DMA marshaling cost (`build/profile_gnomonicstars_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative window (frames
513–576, 4.00 s). The pack and submit run inside the flywheel wake, on the
1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max, CPU
share:

```
isr_wake         18434/s  0.52 / 1.51 / ~98 us   cpu 2.77%
  isr_pack        2304/s  4.65 / 5.59 / ~97 us   cpu 1.28%
  isr_dma_submit  2304/s  0.60 / 0.91 / 2.1 us   cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, and
includes pack + submit nested); `isr_pack` = the 72× `packPixel` marshal, once
per column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush +
eDMA kick).

- **The DMA kick itself is ~0.9 µs of CPU** (submit avg 0.91 µs); the pixel-pack
  marshal around it is the real per-column cost. Whole per-column CPU ≈
  **~12 µs of the 434 µs column period (~2.8 %)**. The wire transfer runs
  asynchronously on eDMA and costs no CPU.
- The ISR is the same driver code as the -Os run and is essentially unchanged
  (2.77 % here vs 2.82 % at -Os); -O3 did not meaningfully move the marshal —
  the pack `min` floor even sits slightly higher (2789 cyc / 4.65 µs vs the -Os
  2621 cyc / 4.37 µs), within run-to-run and code-layout noise.
- The avg/max figures inflate from serial-dump preemption (higher-priority
  USB/DMA-completion ISRs landing inside the scope, plus the dumps themselves)
  — the ~97–98 µs maxima are artifacts, not real ISR growth. The `min` floor is
  the stable per-call cost: pack 4.65 µs (2789 cyc), submit 0.60 µs.
- Net: the ISR machinery steals **~2.77 % of the chip**. Against the 62.5 ms
  window that's ~1.73 ms, leaving **~60.8 ms of render budget per window**. The
  ~14 ms render fits one window with ~4.3× headroom — hence the rock-solid
  16 fps.

## Summary ranking (steady window, share of the 62.5 ms frame)

1. display-window sync (`gn_buffer_wait`) — **77%** (48.4 ms, idle by design —
   the render finishes in ~23 % of the window and spins out the rest).
2. `Scan::Star::draw` rasterize (`gn_star_scan`, incl. `filter_blend`) —
   **22%** (13.5 ms) — the entire render cost: 600 sprites/frame.
3. per-star `transformer.transform` + `make_basis` (`gn_draw_stars` self) —
   **~1%** (0.66 ms).
4. everything else (`gn_timeline_step`; `gn_spiral_build` fires ~once/run) —
   <1%.

GnomonicStars is a light point-sprite rasterizer: at -O3 with render at ~14 ms
against a ~61 ms budget it is even further from window-bound than at -Os, and
the only real work is the 600-star `Scan::Star` loop. The levers (all low-value
at this headroom) are the per-star basis/transform and the sprite distance
evaluation; there is no LUT, cull, or geometry stage to attack. No WASM/native
GnomonicStars figures are recorded in the perf ledger for comparison.

## -O3 vs -Os

Same driver, same instrumentation, same effect — only the optimization level
differs. Comparing the coverage-matched representative windows (~3.4 k blended
px/frame in both) plus the coverage-independent per-call costs:

| metric | -Os | -O3 | delta |
|---|---|---|---|
| `gn_star_scan` per star (~5.66 blends/star) | 26.6 µs | 22.4 µs | **1.19× faster** |
| scan per *blended* pixel | 2,814 cyc | 2,381 cyc | 1.18× faster |
| `filter_blend` per blend | 231 cyc (0.38 µs) | 96 cyc (0.16 µs) | **2.4× faster** |
| `gn_draw_stars` self (transform + basis) | 1.10 ms | 0.66 ms | ~1.7× faster |
| render (`gn_draw_stars`, run range) | ~14–21 ms | ~12–19 ms | ~1.1–1.2× |
| frame cadence | 16 fps (62.5 ms) | 16 fps (62.5 ms) | **unchanged** |
| FLASH code | 41,524 B | 59,428 B | +17,904 B (+43 %) |
| ITCM (RAM1 code) | 27,400 B | 44,728 B | +17,328 B (+63 %) |
| RAM2 free | 4,736 B | 4,736 B | unchanged |

The scan work is largely per-pixel SDF distance/coverage evaluation that -O3
speeds up only modestly (~1.2×); the biggest single win is the blend composite
(2.4×), but it is such a small slice of the scan (~4 %) that the overall render
only drops ~15–20 %. Because the effect was already ~4× under the display-window
budget at -Os, that render win buys **no cadence change** — it stays a
rock-solid 16 fps, now with even more idle headroom.

**-O3 is NOT the shipping config.** It costs +18 KB FLASH and +17 KB ITCM for
this single-effect image; the full 26-effect Phantasm roster overflows FlexRAM
(ITCM) at -O3, which is exactly why the ship build is `-Os`. -O3 is only viable
here as a single-effect profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it parents under `gn_star_scan` (the first scope
  to enter it) and counts every blended pixel across the star loop; its printed
  % is meaningful here because that is the only blend path, but the subtree
  would move if another scope entered blend first. Per-pixel `filter_blend`
  scope overhead inflates the scan counter by a fraction of a percent.
- **Coverage/phase differs between the two runs**: the blended-px/frame breathes
  with the drifting field, so raw per-frame ms are only comparable at matched
  coverage. The -O3-vs-Os table normalizes on per-star / per-blend cost (the
  representative -O3 window was picked to match the -Os coverage) rather than on
  raw window ms.
- Epoch/reset windows (a `frame` count > 64 or a very low first-frame `min`)
  mix an undumped tail into the next window — the first window (min 17,794 µs,
  epoch-start, and the sole `gn_spiral_build` call) is skipped when picking the
  representative; per-call averages stay valid there, per-frame sums do not.
- **-O3 build** (`profile_o3` env): this does **not** ship. The shipping
  Phantasm image is `-Os` because the full 26-effect roster overflows ITCM at
  -O3; this single-effect -O3 image is a profiling artifact only.
- Run with the uncommitted `HS_PROFILE` instrumentation in
  `effects/GnomonicStars.h` (`gn_buffer_wait`, `gn_timeline_step`,
  `gn_spiral_build`, `gn_draw_stars`, `gn_star_scan`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=GnomonicStars`, `-D HS_PROFILE_WINDOW=<frames>`); also
  dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel
  pack, DMA submit) each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (the
  -O3 twin of the shipping `-Os` `profile` env; no `just` recipe for -O3).
