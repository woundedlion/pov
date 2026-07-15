# Flyby on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). This is the **-O3** twin of the
shipping `-Os` report at `docs/profiles/Os/profile_flyby_teensy_2026-07-14.md`.
Raw capture: `build/profile_flyby_o3.log` (64-frame windows, ~72 s single pass,
18 windows). The column-ISR/DMA accumulators are dumped every window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm shipping flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, built at **`-O3`** (`-ffast-math -fno-finite-math-only`) — the -O3 twin of the shipping `-Os` profile |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Flyby<288, 144>` only (single-entry playlist), current working-tree state (default presets, drift 0.7) |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` |

Image size: FLASH code 65,764 B; ITCM (RAM1 code) 48,280 B; RAM2 free 4,736 B.
(-Os shipped 40,044 / 24,888 / 4,736.)

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm
across three steady windows: frames 961–1024 root counter = 2,400,386,954 cyc =
4,000,644.9 µs vs measured `micros()` window sum 4,000,643 µs (Δ ≈ 0.5 ppm);
frames 65–128 Δ ≈ 1.1 ppm; frames 385–448 Δ ≈ 2.2 ppm.

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Flyby is a pure
full-quadrant shader — `Scan::Shader::draw` evaluates every one of those pixels
(no cull, no early-out), so the frame cost is essentially `pixels ×
per-pixel-shader`. Wall time snaps up to a whole number of 62.5 ms windows.

Flyby has no discrete phases; the per-pixel shader cost **breathes** with the
preset lerp (`LERP_FRAMES = 480`, 5 presets, `ease_in_out_sin`). At **-O3 the
whole swing now fits inside one 62.5 ms window**:

- **Expensive presets** (high warp-scale / pole spread): shader **~42.4 ms/frame**
  (peak windows 385–448 / 449–512) — still one-window cadence, **62.5 ms/frame
  (16 fps)**.
- **Cheap presets**: shader **~32.0 ms/frame** (windows 897–960 / 961–1024) →
  **62.5 ms/frame (16 fps)**.

Every one of the 18 windows measured a frame wall of ~62.5 ms (avg 62.43–62.60
ms) — **there is no 8 fps regime at -O3**. The breathing is entirely absorbed
by `fly_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor), the
round-up idle: ~30.5 ms when the shader is cheap, ~20.1 ms when expensive, both
rounding the frame up to a single 62.5 ms window.

## Phase-by-phase readout (64-frame windows, per-frame averages)

The shader cost oscillates continuously; the two regimes below bound it, with a
monotonic transition between (windows 321–384 / 641–704 catch mid-swing at
~38–40 ms). Unlike the -Os run, both regimes now clear the 62.5 ms window, so
the cadence stays at 16 fps throughout. Columns: time/frame, cycles/frame,
% of frame, calls/frame, per-pixel cost.

### Expensive-preset regime — 16 fps (representative: frames 385–448)

```
frame                 62.53 ms  37.52 Mcyc  100%
  fly_shader_draw     42.41 ms  25.44 Mcyc   68%  x1  10368 px  2454 cyc/px
  fly_buffer_wait     20.11 ms  12.07 Mcyc   32%
  fly_timeline_step   11 us      6.6 kcyc     0%
```

Wall: min 59.3 / avg 62.5 / max 64.9 ms. The single `fly_shader_draw` leaf
(stereographic project → noise warp → grid sample → pole normalize → palette →
hue-rotate, per pixel) is the entire render cost; `fly_timeline_step` only
advances the rotation + preset-lerp animations and is negligible. This is the
peak of the cycle (windows 385–448 and 449–512, 42.4 / 42.2 ms) and it still
sits **20 ms under** the 62.5 ms window — the cadence never drops.

### Cheap-preset regime — 16 fps (representative: frames 961–1024)

```
frame                 62.51 ms  37.51 Mcyc  100%
  fly_shader_draw     32.00 ms  19.20 Mcyc   51%  x1  10368 px  1852 cyc/px
  fly_buffer_wait     30.50 ms  18.30 Mcyc   49%
  fly_timeline_step    9 us      5.5 kcyc     0%
```

Wall: min 61.2 / avg 62.5 / max 63.6 ms. Same 10,368-pixel shader; the active
preset pair costs ~1.33× less per pixel (1852 vs 2454 cyc/px). The per-pixel
delta is the whole story — the pole attenuation and warp displacement of the
cheap presets push more of the sphere into the near-transparent pole band and
keep the pattern-frequency trig arguments out of the expensive range-reduction
bands.

### Per-pixel figures

The shader shades the full quadrant unconditionally: **10,368 px/frame**, no
cull. `Scan::Shader::draw` writes/composites in place and does not route through
`filter_blend`, so there is no separate blended-pixel counter here. Per-pixel
scan cost: **1,852 cyc (3.09 µs)** in the cheap regime, **2,454 cyc (4.09 µs)**
in the expensive regime — dominated by the OpenSimplex2 noise warp and the two
`fast_sinf`/`fast_cosf` pattern taps, not the palette blend.

## Column-ISR / DMA marshaling cost (`build/profile_flyby_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative steady window (frames
897–960, 4.0 s). The pack and submit run inside the flywheel wake, on the 1-in-8
wakes that render a column. Columns: rate, per-call min/avg/max, CPU share:

```
isr_wake         18412/s  0.52 / 1.67 / ~196 us  cpu 3.07%
  isr_pack        2304/s  4.65 / 6.77 / ~194 us  cpu 1.55%
  isr_dma_submit  2304/s  0.66 / 0.91 / 1.19 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **13 µs of the 434 µs column
  period (~3.1 %)**. The wire transfer runs asynchronously on eDMA and costs no
  CPU.
- The avg/max figures inflate from preemption by the serial dumps themselves,
  not real ISR growth; the `min` floor (wake 0.52 µs, pack 4.65 µs / 2551 cyc→
  here 2789 cyc, submit ~0.66 µs) is stable across every window.
- Net: the ISR machinery steals **~3.1 % of the chip**. Against the 62.5 ms
  window that's ~1.9 ms, leaving **~60 ms of render budget per window**. At -O3
  even the expensive shader (42.4 ms) fits with ~18 ms to spare, so 16 fps holds
  across the entire preset cycle.

## Summary ranking (expensive regime, share of the 62.5 ms frame)

1. `Scan::Shader::draw` (per-pixel stereo project + noise warp + trig + palette
   + hue-rotate) — **68%** (42.4 ms) — the entire algorithmic cost.
2. display-window sync (`fly_buffer_wait`) — 32% (20.1 ms, idle by design; grows
   to ~30.5 ms / 49% in the cheap regime as the shader shrinks).
3. everything else (`fly_timeline_step`: rotation + preset lerp) — <1%.

Flyby is 100 % rasterizer/shader-bound: a single flat per-pixel loop with no
cull, LUT, or geometry stage to attack. The only levers are the per-pixel math
(cheaper noise, fewer trig taps, earlier pole/alpha reject) or supersample
factor. No WASM/native Flyby figures are recorded in the perf ledger for
comparison.

## -O3 vs -Os

-O3 buys a large, uniform per-pixel speedup that collapses the -Os 8↔16 fps
breathing into a **steady 16 fps**:

| Regime | -Os shader | -O3 shader | speedup | cyc/px -Os→-O3 |
|---|---|---|---|---|
| Expensive | 77.0 ms (4453 cyc/px) | **42.4 ms** (2454 cyc/px) | **1.82×** (−34.6 ms) | −45% |
| Cheap | 49.7 ms (2876 cyc/px) | **32.0 ms** (1852 cyc/px) | **1.55×** (−17.7 ms) | −36% |

- **Cadence**: at -Os the expensive presets overran the 62.5 ms window (~77 ms
  shader → 125 ms/frame, **8 fps**) and only the cheap presets held 16 fps — the
  frame audibly breathed 8↔16 fps. At -O3 the peak shader is 42.4 ms, ~20 ms
  under budget, so **every window renders at 16 fps** (frame wall 62.4–62.6 ms in
  all 18 windows); the breathing now only moves the `fly_buffer_wait` idle
  (20↔30 ms), not the displayed cadence.
- **Size cost**: FLASH +25,720 B (40,044 → 65,764, +64 %); ITCM +23,392 B
  (24,888 → 48,280, +94 %); RAM2 free unchanged (4,736 B).
- **Not the shipping config**: -O3 is viable here only because this is a
  single-effect profile image. The full 26-effect Phantasm roster overflows
  FlexRAM/ITCM at -O3 — that is exactly why the ship build is `-Os`.

## Caveats

- **-O3 build, does not ship**: this is the `profile_o3` env (base `-O3`, else
  identical to the shipping `profile` env). The shipping Phantasm image is `-Os`
  because the full 26-effect roster overflows ITCM at -O3; the only variable vs
  the -Os report is the optimization level.
- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Shader::draw` composites in place, so
  the per-pixel blended-pixel counter never fires; the per-pixel figures are
  derived from the fixed 10,368-px quadrant coverage.
- The first window (frames 1–64) straddles instance start (wall min 31.8 ms) —
  its per-frame sums are skewed and it is skipped when picking representatives;
  per-call averages stay valid.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Flyby.h`
  (three scopes: `fly_buffer_wait`, `fly_timeline_step`, `fly_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Flyby`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (there
  is no `just` recipe for the -O3 env; the -Os twin is `just profile Flyby`).
