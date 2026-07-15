# Liquid2D on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_liquid2d_o3.log` (64-frame windows, ~75 s single pass). The
column-ISR/DMA accumulators are dumped every window. This is the **-O3** twin of
the shipping `-Os` profile;
the only variable between the two runs is the optimization level.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm shipping flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but keeping the base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of forcing `-Os` — the -O3 twin of the shipping `-Os` profile |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Liquid2D<288, 144>` only (single-entry playlist), current working-tree state (default presets) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 73,788 B; ITCM (RAM1 code) 54,728 B; RAM2 free 4,736 B.
(vs `-Os`: FLASH 42,220 B, ITCM 26,248 B — see the `## -O3 vs -Os` section.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 705–768 root counter = 2,399,697,744 cyc = 3,999,496.24 µs vs
measured `micros()` window sum 3,999,500 µs (Δ ≈ 0.94 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Liquid2D is the
stereographic sibling of Flyby — a pure full-quadrant per-pixel shader:
`Scan::Shader::draw` evaluates every one of those 10,368 pixels (no cull, no
early-out), so the frame cost is essentially `pixels × per-pixel-shader`. Wall
time snaps up to a whole number of 62.5 ms windows.

Liquid2D is **steady 16 fps — 62.5 ms/frame — for the whole run** at -O3 exactly
as at -Os. The shader now holds ~24–26 ms/frame (down from ~28–30 ms at -Os),
still comfortably under one 62.5 ms window, so it never crosses into the 2-window
(8 fps) cadence. It was already pinned at 16 fps at -Os with ~33 ms of slack, so
-O3 does not change the visible cadence at all — it only **widens the idle**: the
~3.7 ms it shaves off the shader is absorbed one-for-one by `lq_buffer_wait`.

`lq_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle — with
the shader at ~26 ms it is the ~36 ms remainder of the 62.5 ms window (vs ~33 ms
at -Os).

## Phase-by-phase readout (64-frame windows, per-frame averages)

Liquid2D has no discrete phases and no cadence flips — every steady window is
near-identical (shader 39–42 % of frame across the whole capture, a slight
upward breathe from ~24.1 ms early to a ~26.1 ms plateau by frame ~640). The
first window (frames 1–64) is the startup ramp (min 23.8 ms wall first frame)
and is skipped. One representative settled window below. Columns: time/frame,
cycles/frame, % of frame, calls/frame, per-pixel cost.

### Steady 16 fps (representative: frames 705–768)

```
frame                 62.49 ms  37.50 Mcyc  100%
  lq_buffer_wait      36.40 ms  21.84 Mcyc   58%
  lq_shader_draw      26.07 ms  15.64 Mcyc   42%  x1  10368 px  1509 cyc/px
  lq_timeline_step     27 us     16.1 kcyc    0%
```

Wall: min 61.98 / avg 62.49 / max 63.05 ms. The single `lq_shader_draw` leaf
(stereographic project → liquid warp → grid sample → palette, per pixel; the
per-pixel lambda is intentionally not scoped) is the entire render cost;
`lq_timeline_step` only advances the animation and is negligible.
`lq_buffer_wait` (58 %) is the larger slice but is pure display-sync idle — the
render (shader + timeline ≈ 26.1 ms) fits one window with ~36 ms to spare, so
the buffer wait absorbs the round-up. Numbers are stable across every steady
window: shader 24.1 ms (frames 65–128) rising to a ~26.1 ms plateau
(frames 641–896), never approaching the 62.5 ms window ceiling. Three settled
windows agree tightly: 705–768 = 26.07 ms / 1509 cyc/px, 641–704 = 26.17 ms /
1514 cyc/px, 833–896 = 26.09 ms / 1510 cyc/px.

### Per-pixel figures

The shader shades the full quadrant unconditionally: **10,368 px/frame**, no
cull. `Scan::Shader::draw` writes/composites in place and does **not** route
through `filter_blend`, so there is no separate blended-pixel counter here — the
per-pixel cost is derived from the fixed 10,368-px quadrant coverage. Per-pixel
scan cost: **1,509 cyc (2.51 µs) per pixel** in the steady state (15.64 Mcyc ÷
10,368 px), down from 1,722 cyc (2.87 µs) at -Os — a 213 cyc/px (12.4 %)
reduction that keeps the whole shader well inside one window.

## Column-ISR / DMA marshaling cost (`build/profile_liquid2d_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators (the CycleCounter tree is
main-loop-only), representative steady window (frames 705–768, ~4.0 s). The
pack and submit run inside the flywheel wake, on the 1-in-8 wakes that render a
column. Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18432/s  0.50 / 1.61 / ~187 us  cpu 2.96%
  isr_pack        2304/s  4.65 / 6.52 / ~185 us  cpu 1.50%
  isr_dma_submit  2304/s  0.65 / 0.91 / 1.2 us   cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid, and
includes the nested pack+submit); `isr_pack` = the 72× `packPixel` marshal, once
per column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush +
eDMA kick).

- **The DMA kick itself is ~1 µs of CPU** (submit min 0.65 µs); the pixel-pack
  marshal around it is the real per-column cost. Whole per-column CPU ≈ **12.9 µs
  of the 434 µs column period (~2.96 %)**. The wire transfer runs asynchronously
  on eDMA and costs no CPU.
- The `min` floor (wake 0.50 µs, pack 4.65 µs / 2,789 cyc, submit 0.65 µs) is
  the stable per-call cost; the avg/max figures (pack avg 6.5 µs, ~185 µs
  maxima) inflate from serial-dump preemption inside the scope during long
  frames, not real ISR growth — every window's floor is identical. The ISR path
  is essentially unchanged by -O3 (it is dominated by the fixed packPixel
  marshal, not by shader math).
- Net: the ISR machinery steals **~3.0 % of the chip**. Against the 62.5 ms
  window that's ~1.9 ms, leaving **~60.6 ms of render budget per window**. The
  ~26.1 ms shader needs only ~0.43× of that, so 16 fps holds with even wider
  margin than at -Os.

## Summary ranking (steady 16 fps, share of the 62.5 ms frame)

1. `Scan::Shader::draw` (per-pixel stereo project + liquid warp + grid sample +
   palette) — **42%** (26.1 ms) — the entire algorithmic render cost.
2. display-window sync (`lq_buffer_wait`) — 58% (36.4 ms, idle by design; the
   larger slice but pure round-up, not work — it grew by exactly the shader's
   -O3 saving).
3. everything else (`lq_timeline_step`: animation advance) — <1% (27 µs).

Liquid2D is 100 % rasterizer/shader-bound: a single flat per-pixel loop with no
cull, LUT, or geometry stage to attack. -O3 buys a straight ~14 % speedup on the
per-pixel math (1,722 → 1,509 cyc/px) but does not move the cadence — the effect
was already pinned at 16 fps at -Os with wide slack, so the only visible effect
is more idle in `lq_buffer_wait`. The only cadence levers remain the per-pixel
math (cheaper warp, fewer trig taps) or the supersample factor. No WASM/native
Liquid2D figures are recorded in the perf ledger for comparison.

## -O3 vs -Os

Same board, same driver, same effect and window selection (frames 705–768) — the
only variable is the base optimization level.

| Metric | -Os (ship) | -O3 | Δ |
|---|---|---|---|
| Shader `lq_shader_draw` per frame | 29.76 ms | 26.07 ms | **−3.69 ms (1.14× faster)** |
| Shader cyc/px | 1,722 | 1,509 | −213 cyc/px (−12.4 %) |
| Shader % of frame | 48 % | 42 % | −6 pts |
| `lq_buffer_wait` (idle round-up) | 32.70 ms | 36.40 ms | +3.70 ms (absorbs the saving) |
| Frame wall / cadence | 62.5 ms / **16 fps** | 62.5 ms / **16 fps** | unchanged |
| FLASH code | 42,220 B | 73,788 B | +31,568 B (+74.8 %) |
| ITCM (RAM1 code) | 26,248 B | 54,728 B | +28,480 B (+108.5 %) |

Headline: **-O3 is ~1.14× (14 %) faster on the shader, −3.69 ms/frame**, but the
effect was already at 16 fps at -Os with ~33 ms of slack, so the win lands
entirely in the idle `lq_buffer_wait` — **no cadence change**. The cost is a
~75 % larger FLASH image and a **>2× larger ITCM footprint**.

**-O3 is NOT the shipping config.** The full 26-effect Phantasm image overflows
FlexRAM/ITCM at -O3 (the ITCM more than doubling for a single shader here is
exactly why) — that is precisely why the ship build is `-Os`. -O3 is only viable
as this single-effect profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Shader::draw` composites in place, so the
  per-pixel blended-pixel counter never fires; the per-pixel figures are derived
  from the fixed 10,368-px quadrant coverage. The per-pixel shader lambda is
  intentionally not scoped, so no sub-shader breakdown is available.
- Startup window (frames 1–64, `frame` avg 61.3 ms with a 23.8 ms first-frame
  ramp) mixes the instance warm-up into the window — skipped when picking the
  representative; every later window is steady.
- **`-O3` build** (the `profile_o3` env), which does **NOT** ship: the shipping
  Phantasm image is `-Os` because the full 26-effect roster overflows ITCM at
  -O3. This report exists only to quantify the optimization-level delta.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Liquid2D.h`
  (three scopes: `lq_buffer_wait`, `lq_timeline_step`, `lq_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Liquid2D`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel pack, DMA
  submit) each window.
- `pio run -e profile_o3 -t upload` + `tools/profile_capture.py` (the -O3 env has
  no `just` recipe; the -Os twin uses `just profile Liquid2D [seconds]`).
</content>
</invoke>
