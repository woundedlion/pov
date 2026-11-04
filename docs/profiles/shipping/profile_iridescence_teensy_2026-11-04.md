# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-04, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Iridescence`).
Raw capture: `build/prof/iridescence_ship.log` in the capture worktree,
captured **2026-11-04 08:49** America/Los_Angeles on COM4.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Iridescence 288×144: potential-flow displacement, equirectangular grid, authored palette, two independent baked-cube noise layers; source `e314aed39` plus the playlist commit |
| Method | `HS_PROFILE`, 32-frame windows, 70 seconds, default epoch; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=<tree> bash tools/profile_one.sh Iridescence profile 70 32` |

Full Phantasm image with all three in the playlist: `FLASH: code:1183368,
data:793184`; `RAM1: variables:313760, code:193240`; `RAM2: variables:520064`.
The size and layout gate passes under the 32 KiB flash reserve.

| Region | Master `e314aed39` | Final | Delta |
|---|--:|--:|--:|
| RAM1 code | 191,960 | 193,240 | +1,280 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH code | 1,130,216 | 1,183,368 | +53,152 |
| FLASH data | 788,920 | 793,184 | +4,264 |

Exactness cross-check: root cycles/600 MHz differ from the measured wall sum
by **7.7 ppm** in frames 1-32, outside the validator's 5 ppm bound. Per-frame
render telemetry is unaffected; the residual is the display-phase `micros()`
truncation bias that follows the image layout.

## Frame cadence

**Steady peak: 49.785 ms**, frame 27; mean **48.111 ms**,
**0/1087 spilled (0.00%)**, and **0/1087 at or above 58 ms**.
Frames 2-1088 in complete windows are included. Setup frame 1 is
**92.959 ms** and spills separately.

README cells: peak 🟢 49.785, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. Iridescence shades one 10,368-pixel
quadrant per frame and holds the 16 fps cadence. `canvas_buffer_wait` is the
display-sync idle and is subtracted from render.

## Phase-by-phase readout

The fixed preset runs continuously: the potential field's terms advance, the
grid displaces, and palettes and noise clocks advance. There are no preset
markers.

### Worst steady window (frames 257-288)

```
frame                     62.383 ms 37.430 Mcyc 100.0%
  pov_preserve_half        0.138 ms  0.083 Mcyc   0.2%
  fx_shader_draw          43.464 ms 26.078 Mcyc  69.7%
    fx_prepare_column      1.179 ms  0.707 Mcyc   1.9% x144  8.19 us/col
  fx_palette_base_bake     0.070 ms  0.042 Mcyc   0.1%
  fx_prepare_frame         3.023 ms  1.814 Mcyc   4.8%
  fx_advance               0.927 ms  0.556 Mcyc   1.5%
    palette_cycle_step     0.902 ms  0.541 Mcyc   1.4%
  fx_timeline_step         0.018 ms  0.011 Mcyc   0.0%
  canvas_clear             0.087 ms  0.052 Mcyc   0.1%
  canvas_buffer_wait      14.320 ms  8.592 Mcyc  23.0%
```

Wall min/mean/max: **59.707/62.382/65.045 ms**; render mean/max
48.062/49.672 ms.

### Per-pixel figures

At 10,368 shaded pixels per frame, the 43.5 ms shader draw costs about
2,515 cycles per pixel including column preparation.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.61/2.01/36.99 us  cpu 3.70%
isr_pack         144/frame  min/avg/max 8.61/9.35/11.94 us  cpu 2.15%
isr_dma_submit   144/frame  min/avg/max 0.70/0.94/1.04 us  cpu 0.21%
```

- The wire transfer itself runs asynchronously on DMA.
- The inclusive `isr_wake` share is 3.70%, leaving about 60.2 ms of each
  62.5 ms window for render. Pack and submit run nested inside wake and are not
  added to it. DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 69.7% of the frame, 43.5 ms: per-pixel displacement,
   grid and color.
2. `fx_prepare_column` — 1.9%, 1.2 ms: 144 interval samples and stage
   preparations per frame.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- The per-pixel potential-field kernels and the column phase tail run from flash.
- No dwell-compression knobs were used.
- Captured from the landing worktree: master `e314aed39` plus the playlist commit, with no uncommitted changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Iridescence`,
`HS_PROFILE_WINDOW=32`; `just profile Iridescence` = build + flash + capture.
