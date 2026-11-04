# LiquidWindows on-device profile — Teensy 4.0, segmented mode (2026-11-04, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile LiquidWindows`).
Raw capture: `build/prof/liquidwindows_ship.log` in the capture worktree,
captured **2026-11-04 08:45** America/Los_Angeles on COM4.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | LiquidWindows 288×144: potential-flow displacement, equirectangular grid, generated analogous palette; source `e314aed39` plus the playlist commit |
| Method | `HS_PROFILE`, 32-frame windows, 70 seconds, default epoch; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=<tree> bash tools/profile_one.sh LiquidWindows profile 70 32` |

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
by **2.2 ppm in frames 321-352**, inside the validator's 5 ppm bound.

## Frame cadence

**Steady peak: 36.719 ms**, frame 259; mean **31.942 ms**,
**0/1087 spilled (0.00%)**, and **0/1087 at or above 58 ms**.
Frames 2-1088 in complete windows are included. Setup frame 1 is
**59.340 ms** and spills separately.

README cells: peak 🟢 36.719, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. LiquidWindows shades one 10,368-pixel
quadrant per frame and holds the 16 fps cadence. `canvas_buffer_wait` is the
display-sync idle and is subtracted from render.

## Phase-by-phase readout

The fixed preset runs continuously: the potential field's terms advance, the
grid displaces, and palettes and noise clocks advance. There are no preset
markers.

### Worst steady window (frames 257-288)

```
frame                     62.409 ms 37.445 Mcyc 100.0%
  pov_preserve_half        0.142 ms  0.085 Mcyc   0.2%
  fx_shader_draw          29.464 ms 17.679 Mcyc  47.2%
    fx_prepare_column      1.040 ms  0.624 Mcyc   1.7% x144  7.22 us/col
  fx_palette_base_bake     1.092 ms  0.655 Mcyc   1.7%
    fx_palette_rgb_bake    0.957 ms  0.574 Mcyc   1.5%
  fx_prepare_frame         0.464 ms  0.279 Mcyc   0.7%
  fx_advance               0.963 ms  0.578 Mcyc   1.5%
    palette_cycle_step     0.945 ms  0.567 Mcyc   1.5%
  fx_timeline_step         0.018 ms  0.011 Mcyc   0.0%
  canvas_clear             0.093 ms  0.056 Mcyc   0.1%
  canvas_buffer_wait      29.971 ms 17.983 Mcyc  48.0%
```

Wall min/mean/max: **56.475/62.409/68.344 ms**; render mean/max
32.438/36.719 ms.

### Per-pixel figures

At 10,368 shaded pixels per frame, the 29.5 ms shader draw costs about
1,705 cycles per pixel including column preparation.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.61/2.00/28.07 us  cpu 3.69%
isr_pack         144/frame  min/avg/max 8.61/9.31/11.99 us  cpu 2.14%
isr_dma_submit   144/frame  min/avg/max 0.61/0.94/7.62 us  cpu 0.21%
```

- The wire transfer itself runs asynchronously on DMA.
- The inclusive `isr_wake` share is 3.69%, leaving about 60.2 ms of each
  62.5 ms window for render. Pack and submit run nested inside wake and are not
  added to it. DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 47.2% of the frame, 29.5 ms: per-pixel displacement,
   grid and color.
2. `fx_prepare_column` — 1.7%, 1.0 ms: 144 interval samples and stage
   preparations per frame.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- The per-pixel potential-field kernels and the column phase tail run from flash.
- No dwell-compression knobs were used.
- Captured from the landing worktree: master `e314aed39` plus the playlist commit, with no uncommitted changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=LiquidWindows`,
`HS_PROFILE_WINDOW=32`; `just profile LiquidWindows` = build + flash + capture.
