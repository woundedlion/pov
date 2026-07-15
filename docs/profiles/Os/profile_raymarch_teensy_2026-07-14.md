# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_raymarch.log` (64-frame windows, ~75 s single pass, 8 windows /
512 frames). The column-ISR/DMA accumulators are dumped every window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` PlatformIO env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `Raymarch<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `just profile Raymarch` |

Image size: FLASH code 65,388 B; ITCM (RAM1 code) 36,056 B; RAM2 free 4,736 B.

**Exactness cross-check** — the cycle counter and wall clock agree to ~1 ppm:
window frames 449–512 root counter = 5,390,792,559 cyc = 8,984,654.3 µs vs
measured `micros()` window sum 8,984,665 µs (Δ ≈ 1.2 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). Raymarch is a
pure volume ray-marcher — `Scan::Volume::draw` walks every one of those pixels
with a per-pixel ray-march loop (no cull, no early-out band), so the frame cost
is essentially `pixels × ray-march-steps-per-pixel`. Wall time snaps up to a
whole number of 62.5 ms windows.

The render is **steady** (no discrete phases): the shader costs **~117–124 ms
per frame** in every window, which is just under two windows. So the cadence is
**mostly 2-window — 125 ms/frame (8 fps)** — with occasional frames spilling
into a third window (**~187.5 ms, 5 fps**) when render + jitter crosses the
125 ms boundary; per-frame max reaches ~198 ms. Window-average wall therefore
drifts 127–148 ms as the 2-window / 3-window mix shifts.

`rm_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle — it
ranges **8–17 % (≈10–19 ms/frame)** with the cadence mix, so the render numbers
below are clean.

## Phase-by-phase readout (64-frame windows, per-frame averages)

Raymarch has no phases; one representative window bounds it. Columns:
time/frame, cycles/frame, % of frame; the shader leaf adds calls/frame and
per-pixel cost (cycles = µs × 600).

### Steady ray-march (representative: frames 449–512)

```
frame                140.39 ms  84.23 Mcyc  100%
  rm_timeline_step   121.09 ms  72.65 Mcyc   86%
    rm_shader_draw   121.05 ms  72.63 Mcyc   86%  x1  10368 px  7006 cyc/px
  rm_buffer_wait      19.29 ms  11.58 Mcyc   14%
```

Wall: min 121.9 / avg 140.4 / max 196.6 ms. The single `rm_shader_draw` leaf
(the whole per-pixel `Scan::Volume::draw` ray-march loop) is the entire render
cost; the pixel work self-parents here because the shader runs via a Sprite
callback inside `timeline.step`, so `rm_timeline_step` self-time (animation
advance outside the shader) is only ~40 µs/frame. This holds across the whole
capture: `rm_shader_draw` is **82–91 %** of the frame in every window, render
115.9–123.8 ms, and the only thing that moves between windows is how many
frames spill into a third display window (which inflates `rm_buffer_wait` and
the wall average, not the render).

### Per-pixel figures

The ray-marcher shades the full quadrant unconditionally: **10,368 px/frame**,
no cull. `Scan::Volume::draw` composites in place and does **not** route through
`filter_blend`, so there is no separate blended-pixel counter here — the
per-pixel cost is derived from the fixed 10,368-px quadrant coverage:
**7,006 cyc (11.68 µs) per pixel** in this window, ranging **6,704–7,162 cyc
(11.2–11.9 µs)** across the capture. This is the **highest per-pixel cost of any
effect in the roster** — ~1.6× Flyby's expensive-preset shader (4,453 cyc/px)
— and it is entirely the ray-march step count per pixel (SDF evaluations along
each ray until hit/exit), not a blend or palette stage.

## Column-ISR / DMA marshaling cost (`build/profile_raymarch.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative window (frames
449–512, 8.98 s). The pack and submit run inside the flywheel wake, on the
1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max, CPU
share:

```
isr_wake         18432/s  0.59 / 2.01 / ~135 us  cpu 3.70%
  isr_pack        2304/s  4.25 / 8.50 / ~133 us  cpu 1.95%
  isr_dma_submit  2304/s  0.65 / 0.96 / 1.2 us   cpu 0.22%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid,
includes pack+submit nested); `isr_pack` = the 72× `packPixel` marshal, once
per column; `isr_dma_submit` = `submitFrame` (overrun check + dcache flush +
eDMA kick).

- **The DMA kick itself is ~1 µs of CPU**; the pixel-pack marshal around it is
  the real per-column cost. Whole per-column CPU ≈ **13 µs of the 434 µs column
  period (~3.0 %)**. The wire transfer runs asynchronously on eDMA and costs no
  CPU.
- The `min` floor (pack 4.25 µs / 2551 cyc, submit 0.65 µs) is the stable
  per-call cost; the avg/max inflate from serial-dump preemption in the long
  ~140–198 ms frames, not real ISR growth.
- Net: the ISR machinery steals **~3.7 % of the chip** (`isr_wake` cpu). Against
  the 62.5 ms window that's ~2.3 ms, leaving **~60 ms of render budget per
  window**. The 121 ms shader needs **~2×** that, which is why it holds 8 fps
  (2-window) and tips to 5 fps whenever render + jitter overruns 125 ms.

## Summary ranking (share of the 140 ms representative frame)

1. `Scan::Volume::draw` (per-pixel volume ray-march loop) — **86%** (121.1 ms)
   — the entire algorithmic cost, the most expensive shader in the roster.
2. display-window sync (`rm_buffer_wait`) — **14%** (19.3 ms, idle by design).
3. everything else (`rm_timeline_step` self-time: animation advance) — <1%.

Raymarch is **100 % per-pixel-shader-bound**: a flat volume ray-march over the
whole 10,368-px quadrant with no cull, LUT, or geometry stage to attack. The
only levers are the per-pixel ray-march itself (fewer/larger march steps,
earlier ray termination, cheaper SDF, or a coarser supersample). No WASM/native
Raymarch figures are recorded in the perf ledger for comparison.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **No `filter_blend` subtree**: `Scan::Volume::draw` composites in place, so the
  per-pixel blended-pixel counter never fires; the per-pixel figures are derived
  from the fixed 10,368-px quadrant coverage.
- The pixel work self-parents under `rm_timeline_step` (the shader runs via a
  Sprite callback inside `timeline.step`), so `rm_timeline_step` ≈
  `rm_shader_draw`; the ~40 µs gap between them is the animation-advance
  self-time.
- Per-frame max ~198 ms frames are the occasional 3-window (5 fps) frames; the
  bulk are 2-window (8 fps). Windows are otherwise all similar — one
  representative bounds the whole capture.
- `-Os` build: the shipping Phantasm config; an `-O3` profile would differ.
- Run with the uncommitted `HS_PROFILE` instrumentation in `effects/Raymarch.h`
  (three scopes: `rm_buffer_wait`, `rm_timeline_step`, `rm_shader_draw`).

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=Raymarch`, `-D HS_PROFILE_WINDOW=<frames>`); also dumps
  the `HS_ISR_PROFILE` column-ISR accumulators each window.
- `pio run -e profile -t upload` + `tools/profile_capture.py`, or just:
  `just profile Raymarch [seconds]`.
</content>
</invoke>
