# Liquid2D on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_liquid2d.log` (64-frame windows, ~75 s single pass). The
column-ISR/DMA accumulators are dumped every window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Liquid2D<288, 144>` only (single-entry playlist), current working-tree state (default presets) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile Liquid2D` |

Image size: FLASH code 42,220 B; ITCM (RAM1 code) 26,248 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 449–512 root counter = 2,399,923,668 cyc = 3,999,872.78 µs vs
measured `micros()` window sum 3,999,879 µs (Δ ≈ 1.6 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Liquid2D is the
stereographic sibling of Flyby — a pure full-quadrant per-pixel shader:
`Scan::Shader::draw` evaluates every one of those 10,368 pixels (no cull, no
early-out), so the frame cost is essentially `pixels × per-pixel-shader`. Wall
time snaps up to a whole number of 62.5 ms windows.

Unlike Flyby (whose preset lerp breathes the shader across the window boundary
into an 8 fps regime), Liquid2D is **steady 16 fps — 62.5 ms/frame — for the
whole run**. The shader holds ~28–30 ms/frame, comfortably under one 62.5 ms
window, so it never crosses into the 2-window (8 fps) cadence. Its per-pixel
cost is roughly half Flyby's: ~1,720 cyc/px vs Flyby's cheap-regime 2,876 and
expensive-regime 4,453 cyc/px, which is exactly why it stays pinned at 16 fps.

`lq_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle — with
the shader at ~30 ms it is the ~33 ms remainder of the 62.5 ms window.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Liquid2D has no discrete phases and no cadence flips — every steady window is
near-identical (shader 44–47 % of frame across the whole capture, a slight
upward breathe from ~28.1 ms early to a ~29.8 ms plateau by frame ~700). The
first window (frames 1–64) is the startup ramp (min 27.6 ms wall) and is
skipped. One representative settled window below. Columns: time/frame,
cycles/frame, % of frame, calls/frame, per-pixel cost.

### Steady 16 fps (representative: frames 705–768)

```
frame                 62.49 ms  37.50 Mcyc  100%
  lq_buffer_wait      32.70 ms  19.62 Mcyc   52%
  lq_shader_draw      29.76 ms  17.85 Mcyc   48%  x1  10368 px  1722 cyc/px
  lq_timeline_step     41 us     24.7 kcyc    0%
```

Wall: min 62.10 / avg 62.49 / max 62.79 ms. The single `lq_shader_draw` leaf
(stereographic project → liquid warp → grid sample → palette, per pixel; the
per-pixel lambda is intentionally not scoped) is the entire render cost;
`lq_timeline_step` only advances the animation and is negligible.
`lq_buffer_wait` (52 %) is the larger slice but is pure display-sync idle — the
render (shader + timeline ≈ 29.8 ms) fits one window with ~33 ms to spare, so
the buffer wait absorbs the round-up. Numbers are stable across every steady
window: shader 28.1 ms (frames 65–128) rising to a 29.8 ms plateau
(frames 705–896), never approaching the 62.5 ms window ceiling.

### Per-pixel figures

The shader shades the full quadrant unconditionally: **10,368 px/frame**, no
cull. `Scan::Shader::draw` writes/composites in place and does **not** route
through `filter_blend`, so there is no separate blended-pixel counter here — the
per-pixel cost is derived from the fixed 10,368-px quadrant coverage. Per-pixel
scan cost: **1,722 cyc (2.87 µs) per pixel** in the steady state (17.85 Mcyc ÷
10,368 px) — about 0.60× of Flyby's cheap-regime 2,876 cyc/px and well under
half its expensive 4,453, which keeps the whole shader inside one window.

## Column-ISR / DMA marshaling cost (`build/profile_liquid2d.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative steady window (frames 705–768, ~4.0 s). The
pack and submit run inside the flywheel wake, on the 1-in-8 wakes that render a
column. Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.62 / 1.66 / ~192 us  cpu 3.06%
  isr_pack        2304/s  4.26 / 6.03 / ~190 us  cpu 1.38%
  isr_dma_submit  2304/s  0.64 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, and
includes the nested pack+submit); `isr_pack` = the 72× `packPixel` marshal, once
per column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush +
eDMA kick).

- **The DMA kick itself is ~1 µs of CPU** (submit min 0.64 µs); the pixel-pack
  marshal around it is the real per-column cost. Whole per-column CPU ≈ **13 µs
  of the 434 µs column period (~3.06 %)**. The wire transfer runs asynchronously
  on eDMA and costs no CPU.
- The `min` floor (wake 0.62 µs, pack 4.26 µs / 2,558 cyc, submit 0.64 µs) is
  the stable per-call cost; the avg/max figures (pack avg 6.0 µs, ~190 µs
  maxima) inflate from serial-dump preemption inside the scope during long
  frames, not real ISR growth — every window's floor is identical.
- Net: the ISR machinery steals **~3.1 % of the chip**. Against the 62.5 ms
  window that's ~1.9 ms, leaving **~60.6 ms of render budget per window**. The
  ~29.8 ms shader needs only ~0.5× of that, so 16 fps holds with wide margin.

## Summary ranking (steady 16 fps, share of the 62.5 ms frame)

1. `Scan::Shader::draw` (per-pixel stereo project + liquid warp + grid sample +
   palette) — **48%** (29.8 ms) — the entire algorithmic render cost.
2. display-window sync (`lq_buffer_wait`) — 52% (32.7 ms, idle by design; the
   larger slice but pure round-up, not work).
3. everything else (`lq_timeline_step`: animation advance) — <1% (41 µs).

Liquid2D is 100 % rasterizer/shader-bound: a single flat per-pixel loop with no
cull, LUT, or geometry stage to attack. Because its per-pixel math is ~0.6× of
Flyby's cheap regime, the whole quadrant clears one 62.5 ms window and the effect
never leaves 16 fps. The only levers are the per-pixel math (cheaper warp, fewer
trig taps) or supersample factor. No WASM/native Liquid2D figures are recorded
in the perf ledger for comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Shader::draw` composites in place, so the
  per-pixel blended-pixel counter never fires; the per-pixel figures are derived
  from the fixed 10,368-px quadrant coverage. The per-pixel shader lambda is
  intentionally not scoped, so no sub-shader breakdown is available.
- Startup window (frames 1–64, `frame` avg 61.4 ms with a 27.6 ms first-frame
  ramp) mixes the instance warm-up into the window — skipped when picking the
  representative; every later window is steady.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Liquid2D.h`
  (three scopes: `lq_buffer_wait`, `lq_timeline_step`, `lq_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Liquid2D`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel pack, DMA
  submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile Liquid2D [seconds]`.
