# Raymarch on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_raymarch_o3.log` (64-frame windows, ~72 s single pass, 9 windows /
576 frames). The column-ISR/DMA accumulators are dumped every window. This is the
**-O3** twin of the shipping `-Os` report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm shipping flags (**`-O3`** `-ffast-math -fno-finite-math-only`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` — the -O3 twin of the shipping `-Os` profile (only the optimization level differs) |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Raymarch<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 93,124 B; ITCM (RAM1 code) 47,960 B; RAM2 free 4,736 B.
(The -Os shipping image is 65,388 / 36,056 B — see [`## -O3 vs -Os`](#-o3-vs--os).)

**Exactness cross-check** — the cycle counter and wall clock agree to ~1 ppm:
window frames 449–512 root counter = 4,794,768,118 cyc = 7,991,280.2 µs vs
measured `micros()` window sum 7,991,291 µs (Δ ≈ 1.4 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Raymarch is a
pure volume ray-marcher — `Scan::Volume::draw` walks every one of those pixels
with a per-pixel ray-march loop (no cull, no early-out band), so the frame cost
is essentially `pixels × ray-march-steps-per-pixel`. Wall time snaps up to a
whole number of 62.5 ms windows.

The render is **steady** (no discrete phases): at -O3 the shader costs **~92–102 ms
per frame** in every window, comfortably between one and two windows. So the
cadence is a **solid 2-window — 125 ms/frame (8 fps)** in every frame of the
capture. Per-frame wall stays in a tight **120–129 ms** band across all nine
windows (avg ~124–125 ms), and the per-frame max never approaches the third
window (~187.5 ms). This is the visible -O3 win: the shader drops well clear of
the 125 ms boundary, so the occasional 3-window (5 fps) spills seen at -Os are
gone — 8 fps is now rock-solid.

`rm_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle; with the
faster render it absorbs more of the fixed 125 ms frame and ranges **18–26 %
(≈22.6–33.4 ms/frame)**, so the render numbers below are clean.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Raymarch has no phases; one representative window bounds it. Columns:
time/frame, cycles/frame, % of frame; the shader leaf adds calls/frame and
per-pixel cost (cycles = µs × 600).

### Steady ray-march (representative: frames 449–512)

```
frame                124.86 ms  74.92 Mcyc  100%
  rm_timeline_step    98.11 ms  58.87 Mcyc   79%
    rm_shader_draw    98.08 ms  58.85 Mcyc   79%  x1  10368 px  5676 cyc/px
  rm_buffer_wait      26.75 ms  16.05 Mcyc   21%
```

Wall: min 122.6 / avg 124.9 / max 128.6 ms. The single `rm_shader_draw` leaf
(the whole per-pixel `Scan::Volume::draw` ray-march loop) is the entire render
cost; the pixel work self-parents here because the shader runs via a Sprite
callback inside `timeline.step`, so `rm_timeline_step` self-time (animation
advance outside the shader) is only ~32 µs/frame. This holds across the whole
capture: `rm_shader_draw` is **73–82 %** of the frame in every window, render
91.6–102.1 ms, and the only thing that moves between windows is small
frame-to-frame ray-march variation (which shifts the `rm_buffer_wait` idle
share, not the cadence — every frame is 2-window).

### Per-pixel figures

The ray-marcher shades the full quadrant unconditionally: **10,368 px/frame**,
no cull. `Scan::Volume::draw` composites in place and does **not** route through
`filter_blend`, so there is no separate blended-pixel counter here — the
per-pixel cost is derived from the fixed 10,368-px quadrant coverage:
**5,676 cyc (9.46 µs) per pixel** in this window, ranging **5,300–5,907 cyc
(8.83–9.85 µs)** across the capture. Even at -O3 this remains the **highest
per-pixel cost of any effect in the roster**, and it is entirely the ray-march
step count per pixel (SDF evaluations along each ray until hit/exit), not a
blend or palette stage.

## Column-ISR / DMA marshaling cost (`build/profile_raymarch_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative window (frames
449–512, 7.99 s). The pack and submit run inside the flywheel wake, on the
1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max, CPU
share:

```
isr_wake         18432/s  0.56 / 2.03 / ~161 us  cpu 3.74%
  isr_pack        2304/s  4.78 / 9.25 / ~159 us  cpu 2.13%
  isr_dma_submit  2304/s  0.68 / 0.95 / 1.26 us  cpu 0.21%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid,
includes pack+submit nested); `isr_pack` = the 72× `packPixel` marshal, once
per column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush +
eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **16 µs of the 434 µs column
  period (~3.7 %)**. The wire transfer runs asynchronously on eDMA and costs no
  CPU. This ISR machinery is driver code and is essentially unchanged from the
  -Os run (pack min-floor 2865 vs 2551 cyc is layout noise, not an effect win).
- The `min` floor (pack 4.78 µs / 2865 cyc, submit 0.68 µs / 408 cyc) is the
  stable per-call cost; the avg/max inflate from serial-dump preemption in the
  long ~125 ms frames, not real ISR growth.
- Net: the ISR machinery steals **~3.7 % of the chip** (`isr_wake` cpu). Against
  the 62.5 ms window that's ~2.3 ms, leaving **~60 ms of render budget per
  window**. The ~98 ms shader needs **~1.6×** that, which is why it holds a
  steady 8 fps (2-window) — but now with ~27 ms of headroom under 125 ms, so it
  no longer tips into a 3-window (5 fps) frame the way -Os did.

## Summary ranking (share of the 124.9 ms representative frame)

1. `Scan::Volume::draw` (per-pixel volume ray-march loop) — **79%** (98.1 ms)
   — the entire algorithmic cost, the most expensive shader in the roster.
2. display-window sync (`rm_buffer_wait`) — **21%** (26.8 ms, idle by design).
3. everything else (`rm_timeline_step` self-time: animation advance) — <1%.

Raymarch is **100 % per-pixel-shader-bound**: a flat volume ray-march over the
whole 10,368-px quadrant with no cull, LUT, or geometry stage to attack. The
only levers are the per-pixel ray-march itself (fewer/larger march steps,
earlier ray termination, cheaper SDF, or a coarser supersample). -O3 shaves the
per-pixel step cost but does not change the algorithm's structure. No
WASM/native Raymarch figures are recorded in the perf ledger for comparison.

## -O3 vs -Os

Only the base optimization level differs between this run and the shipping
`-Os` report; the driver, image geometry, and instrumentation are identical.

| | -Os (ship) | -O3 | Δ |
|---|---|---|---|
| `rm_shader_draw` render (rep window) | 121.05 ms | 98.08 ms | **−22.97 ms, 1.23× faster** |
| Per-pixel cost | 7,006 cyc/px (11.68 µs) | 5,676 cyc/px (9.46 µs) | −1,330 cyc/px, 1.23× |
| Render range (capture) | 115.9–123.8 ms | 91.6–102.1 ms | ~1.22× faster |
| Frame wall (rep avg) | 140.4 ms | 124.9 ms | −15.5 ms |
| Cadence | 8 fps, occasional 5 fps (3-window) spills (max ~198 ms) | **steady 8 fps**, no spills (max ~129 ms) | spills eliminated |
| FLASH code | 65,388 B | 93,124 B | +27,736 B (+42%) |
| ITCM (RAM1) | 36,056 B | 47,960 B | +11,904 B (+33%) |
| RAM2 free | 4,736 B | 4,736 B | — |

**Headline: -O3 makes the ray-march ~1.23× faster (−23 ms/frame, 7,006 → 5,676
cyc/px).** The frame stays at 8 fps because the shader is still >62.5 ms (one
window) and cannot reach 16 fps, but the ~23 ms saving pulls render clear of the
125 ms boundary, converting the flaky -Os 8/5-fps mix into a **rock-solid 8 fps**.

-O3 is **not** the shipping config: it costs +27.7 KB FLASH and +11.9 KB ITCM
here, and the full 26-effect Phantasm roster overflows FlexRAM/ITCM at -O3
(exactly why the ship build is -Os). -O3 is only viable as this single-effect
profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Volume::draw` composites in place, so the
  per-pixel blended-pixel counter never fires; the per-pixel figures are derived
  from the fixed 10,368-px quadrant coverage.
- The pixel work self-parents under `rm_timeline_step` (the shader runs via a
  Sprite callback inside `timeline.step`), so `rm_timeline_step` ≈
  `rm_shader_draw`; the ~32 µs gap between them is the animation-advance
  self-time.
- Frames 1–64 are the epoch/startup window (min wall 105.2 ms) — skipped when
  choosing the representative window.
- **-O3 build** (the `profile_o3` env): this does **not** ship — the shipping
  Phantasm image is `-Os` because the full 26-effect roster overflows ITCM at
  -O3. This image profiles Raymarch alone.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Raymarch.h`
  (three scopes: `rm_buffer_wait`, `rm_timeline_step`, `rm_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Raymarch`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` for the
  -O3 image; the -Os twin uses `just profile Raymarch [seconds]`.
</content>
</invoke>
