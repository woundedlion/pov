# Flyby on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw captures:
`build/profile_flyby.log` (128-frame windows, ~2 epochs),
`build/profile_flyby_w32.log` (32-frame windows, phase resolution). The
column-ISR/DMA accumulators are dumped every window in both.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Flyby<288, 144>` only (single-entry playlist), current working-tree state (default presets, drift 0.7) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 128 (resp. 32) frames then reset |
| Reproduce | `just profile Flyby` |

Image size: FLASH code 40,044 B; ITCM (RAM1 code) 24,888 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
coarse window frames 385–512 (2nd epoch) root counter = 9,601,869,962 cyc =
16,003,116.6 µs vs measured `micros()` window sum 16,003,130 µs (Δ ≈ 0.8 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Flyby is a pure
full-quadrant shader — `Scan::Shader::draw` evaluates every one of those pixels
(no cull, no early-out), so the frame cost is essentially `pixels ×
per-pixel-shader`. Wall time snaps up to a whole number of 62.5 ms windows.

Flyby has no discrete phases; instead the per-pixel shader cost **breathes**
with the preset lerp (`LERP_FRAMES = 480`, 5 presets, `ease_in_out_sin`), and
it swings right across the 62.5 ms window boundary, so the cadence flips:

- **Expensive presets** (high warp-scale / pole spread): shader ~77 ms/frame →
  2-window cadence, **125 ms/frame (8 fps)**.
- **Cheap presets**: shader ~50 ms/frame → 1-window cadence, **62.5 ms/frame
  (16 fps)**.

`fly_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor) is timed
separately and is exactly the round-up idle — ~13 ms when the shader fits one
window, ~48 ms when it overflows into a second.

## Phase-by-phase readout (32-frame windows, per-frame averages)

Per the preset cycle: the shader cost oscillates continuously; the two regimes
below bound it, with a monotonic transition between (w32 frames 289–320 caught
mid-swing at avg 84.6 ms). Epoch (960 revolutions = 120 s) reconstructs the
effect. Columns: time/frame, cycles/frame, % of frame, calls/frame,
per-pixel cost.

### Expensive-preset regime — 8 fps (representative: w32 frames 385–416)

```
frame                125.00 ms  75.00 Mcyc  100%
  fly_shader_draw     76.95 ms  46.17 Mcyc   62%  x1  10368 px  4453 cyc/px
  fly_buffer_wait     48.03 ms  28.82 Mcyc   38%
  fly_timeline_step    13 us     7.9 kcyc     0%
```

Wall: min 122.6 / avg 125.0 / max 127.3 ms. The single `fly_shader_draw` leaf
(stereographic project → noise warp → grid sample → pole normalize → palette →
hue-rotate, per pixel) is the entire render cost; `fly_timeline_step` only
advances the rotation + preset-lerp animations and is negligible. This regime
holds across the whole expensive half of the cycle (coarse frames 385–640,
1153–1280): shader 61–62 % of the 125 ms frame in every window.

### Cheap-preset regime — 16 fps (representative: w32 frames 97–128)

```
frame                 62.49 ms  37.49 Mcyc  100%
  fly_shader_draw     49.71 ms  29.82 Mcyc   80%  x1  10368 px  2876 cyc/px
  fly_buffer_wait     12.77 ms   7.67 Mcyc   20%
  fly_timeline_step    11 us     6.6 kcyc     0%
```

Wall: min 56.2 / avg 62.5 / max 67.0 ms. Same 10,368-pixel shader, but the
active preset pair costs ~1.55× less per pixel (2876 vs 4453 cyc/px), so it
clears the 62.5 ms window and the cadence doubles to 16 fps. The per-pixel
delta is the whole story — the pole attenuation and warp displacement of the
cheap presets push more of the sphere into the near-transparent pole band and
keep the pattern-frequency trig arguments out of the expensive range-reduction
bands.

### Per-pixel figures

The shader shades the full quadrant unconditionally: **10,368 px/frame**, no
cull. `Scan::Shader::draw` writes/composites in place and does not route through
`filter_blend`, so there is no separate blended-pixel counter here. Per-pixel
scan cost: **2,876 cyc (4.79 µs)** in the cheap regime, **4,453 cyc (7.42 µs)**
in the expensive regime — dominated by the OpenSimplex2 noise warp and the two
`fast_sinf`/`fast_cosf` pattern taps, not the palette blend.

## Column-ISR / DMA marshaling cost (`build/profile_flyby.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative cheap-regime window
(coarse frames 897–1024, 8.0 s). The pack and submit run inside the flywheel
wake, on the 1-in-8 wakes that render a column. Columns: rate, per-call
min/avg/max, CPU share:

```
isr_wake         18432/s  0.60 / 1.77 / ~245 us  cpu 3.3%
  isr_pack        2304/s  4.25 / 6.49 / ~243 us  cpu 1.5%
  isr_dma_submit  2304/s  0.15 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **14 µs of the 434 µs column
  period (~3.3 %)**. The wire transfer runs asynchronously on eDMA and costs no
  CPU.
- The avg/max figures inflate in the 8 fps windows (isr_wake avg climbs to
  ~2.0 µs, cpu to ~6 %) — this is preemption by the serial dumps themselves,
  not real ISR growth; the `min` floor (wake 0.58 µs, pack 4.25 µs) is stable
  across every window.
- Net: the ISR machinery steals **~3.3 % of the chip**. Against the 62.5 ms
  window that's ~2.1 ms, leaving **~60 ms of render budget per window**. The
  expensive shader (77 ms) needs **~1.28×** to hold 16 fps; the cheap shader
  (50 ms) already fits.

## Summary ranking (expensive regime, share of the 125 ms frame)

1. `Scan::Shader::draw` (per-pixel stereo project + noise warp + trig + palette
   + hue-rotate) — **62%** (77.0 ms) — the entire algorithmic cost.
2. display-window sync (`fly_buffer_wait`) — 38% (48.0 ms, idle by design).
3. everything else (`fly_timeline_step`: rotation + preset lerp) — <1%.

Flyby is 100 % rasterizer/shader-bound: a single flat per-pixel loop with no
cull, LUT, or geometry stage to attack. The only levers are the per-pixel math
(cheaper noise, fewer trig taps, earlier pole/alpha reject) or supersample
factor. No WASM/native Flyby figures are recorded in the perf ledger for
comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Shader::draw` composites in place, so
  the per-pixel blended-pixel counter never fires; the per-pixel figures are
  derived from the fixed 10,368-px quadrant coverage.
- Epoch-boundary windows (coarse `frame` count > 128, e.g. 131/135 calls) mix
  the torn-down instance's undumped tail into the new instance's first window —
  per-call averages stay valid, per-frame sums don't; those windows are skipped
  when picking representatives.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Flyby.h`
  (three scopes: `fly_buffer_wait`, `fly_timeline_step`, `fly_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Flyby`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile Flyby [seconds]`.
