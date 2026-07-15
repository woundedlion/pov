# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_sphericalharmonics.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `SphericalHarmonics<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile SphericalHarmonics` |

Image size: FLASH code 33,508 B; ITCM (RAM1 code) 21,400 B; RAM2 free 4,736 B —
one of the **smallest images** in the roster.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 513–576 root counter = 2,400,367,319 cyc = 4,000,612.20 µs vs
measured `micros()` window sum 4,000,619 µs (Δ ≈ 1.7 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). SphericalHarmonics
is a full-quadrant mesh/field rasterizer: a single `Scan::rasterize` sweeps the
whole quadrant and evaluates the harmonic field per pixel in the shader.
`filter_blend` confirms it shades essentially the **entire** quadrant every
frame (a constant 10,658 blends/frame). Wall time snaps up to a whole number of
62.5 ms windows.

The render (rasterize + timeline + the un-scoped per-frame field build) costs
only **~20 ms/frame**, comfortably inside one window, so the cadence is a
**rock-steady 16 fps (62.5 ms/frame)** across the entire pass — the render cost
breathes (18.4–22.4 ms as the field animates) but never comes near the window
boundary, so it never flips to 8 fps.

`sh_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle —
~42 ms/frame, the bulk of the window the cheap render leaves empty.

## Phase-by-phase readout (64-frame windows, per-frame averages)

SphericalHarmonics has **no discrete phases**: every steady window looks alike
(`sh_rasterize` 29–35 % of frame, `sh_buffer_wait` 64–70 %). One representative
window suffices. Columns: time/frame, cycles/frame, % of frame; leaves add
xN calls/frame and per-call cost. The first window (`frames 1-64`, min ≈ 18.5 ms
epoch-start) is skipped.

### Steady render — 16 fps (representative: frames 513–576)

```
frame                 62.51 ms  37.51 Mcyc  100%
  sh_rasterize        20.14 ms  12.08 Mcyc   32%
    filter_blend       3.02 ms   1.81 Mcyc        x10658  170 cyc/blend
  sh_buffer_wait      42.35 ms  25.41 Mcyc   68%
  sh_timeline_step      13.9 us   8.3 kcyc    0%
```

Wall: min 59.6 / avg 62.5 / max 66.4 ms. The single `sh_rasterize` leaf (the
full-sphere `Scan::rasterize`, per-pixel harmonic-field eval in the shader) is
the whole render cost; `sh_timeline_step` only advances animations and is
negligible. The per-frame field build (`decode_lm` ×2 + `HarmonicField` ctor) is
cheap un-scoped setup — the un-instrumented remainder is **~0.9 µs/frame**.
This regime holds across the whole pass: `sh_rasterize` swings 18.4–22.4 ms
(29–35 % of frame) as the harmonic field animates, but coverage is fixed, so
the cadence never leaves 16 fps.

### Per-pixel figures

`filter_blend` (the pre-existing per-pixel counter in filter.h): a **constant**
10,658 blended pixels/frame in every window — the harmonic field covers
essentially the **entire 10,368-px quadrant** every frame (the swing in
`sh_rasterize` is per-pixel cost, not coverage) — at **170 cyc (0.28 µs) per
blend** (includes its own scope overhead). `sh_rasterize` ≈ 12.08 Mcyc over
10,658 blended px ⇒ **~1,134 cyc per blended pixel**, of which the blend is only
~170 cyc; the remaining **~964 cyc/px is the spherical-harmonic field
evaluation** in the shader (the dominant cost). Across windows the field eval
breathes from ~1,034 to ~1,260 cyc/blended px.

## Column-ISR / DMA marshaling cost (`build/profile_sphericalharmonics.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative steady window
(frames 513–576, 4.0 s). The pack and submit run inside the flywheel wake, on
the 1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max,
CPU share:

```
isr_wake         18417/s  0.56 / 1.68 / ~202 us  cpu 3.1%
  isr_pack        2304/s  4.37 / 6.25 / ~201 us  cpu 1.4%
  isr_dma_submit  2304/s  0.65 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick).

- **The DMA kick itself is ~1 µs of CPU** (submit floor 0.63 µs); the pixel-pack
  marshal around it is the real per-column cost (floor 4.37 µs). Whole
  per-column ISR CPU ≈ **13 µs of the 434 µs column period (~3.1 %)**. The wire
  transfer runs asynchronously on eDMA and costs no CPU.
- The avg/max figures inflate (pack avg 6.25 µs, max ~201 µs) — this is
  preemption by the serial dumps themselves, not real ISR growth; the `min`
  floor (pack 4.37 µs / 2621 cyc, submit 0.63 µs) is stable across every window.
- Net: the ISR machinery steals **~3.1 % of the chip**. Against the 62.5 ms
  window that's ~1.9 ms, leaving **~60 ms of render budget per window** — the
  ~20 ms render already fits ~3× over, hence the unwavering 16 fps.

## Summary ranking (steady window, share of the 62.5 ms frame)

1. display-window sync (`sh_buffer_wait`) — **68%** (42.4 ms), **idle by
   design**: the cheap render leaves most of the window empty, and wall time
   rounds up to the whole 62.5 ms window.
2. `Scan::rasterize` (full-quadrant harmonic-field shade + blend) — **32%**
   (20.1 ms) — the **entire algorithmic render cost**; `filter_blend` (~3.0 ms)
   is just the composite, ~964 cyc/px is the field eval.
3. everything else (`sh_timeline_step` + un-scoped field build) — <1% (~15 µs).

SphericalHarmonics is a cheap, fully render-bound effect that never stresses the
window: the only real lever is the per-pixel harmonic-field math (`sh_rasterize`
at ~1,134 cyc/blended px), and even halving it would not change the 16 fps
cadence — it is already idle-bound, not compute-bound.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it nests under `sh_rasterize` (the only scope
  that enters it) and its printed % is parent-relative (~15 % of rasterize,
  ~4.8 % of frame); its calls/cycles are correct and constant every window.
- Per-pixel `filter_blend` scope overhead is folded into the ~170 cyc/blend and
  slightly inflates `sh_rasterize`; the coarse `sh_*` scopes are negligible.
- The first window (`frames 1-64`, `min ≈ 18.5 ms` epoch-start) mixes the
  instance's cold first frames and is skipped when picking the representative.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with uncommitted `HS_PROFILE` instrumentation in
  `effects/SphericalHarmonics.h` (three scopes: `sh_buffer_wait`,
  `sh_timeline_step`, `sh_rasterize`); the per-frame field build
  (`decode_lm` ×2 + `HarmonicField` ctor) is left un-scoped.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=SphericalHarmonics`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel
  pack, DMA submit) each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile SphericalHarmonics [seconds]`.
