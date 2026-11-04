# GridWarp on-device profile — Teensy 4.0, segmented mode (2026-11-04, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GridWarp`).
Raw capture: `build/prof/gridwarp_ship.log` in the capture worktree,
captured **2026-11-04 08:47** America/Los_Angeles on COM4.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridWarp 288×144: outer camera, potential-flow displacement, equirectangular grid, direct displacement and baked-cube noise modulation; source `e314aed39` plus the playlist commit |
| Method | `HS_PROFILE`, 32-frame windows, 70 seconds, default epoch; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=<tree> bash tools/profile_one.sh GridWarp profile 70 32` |

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
by **1.0 ppm in frames 1-32**, inside the validator's 5 ppm bound.

## Frame cadence

**Steady peak: 44.038 ms**, frame 257; mean **39.225 ms**,
**0/1087 spilled (0.00%)**, and **0/1087 at or above 58 ms**.
Frames 2-1088 in complete windows are included. Setup frame 1 is
**69.816 ms** and spills separately.

README cells: peak 🟢 44.038, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. GridWarp shades one 10,368-pixel
quadrant per frame and holds the 16 fps cadence. `canvas_buffer_wait` is the
display-sync idle and is subtracted from render.

## Phase-by-phase readout

The fixed preset runs continuously: the potential field's terms advance, the
grid displaces, and palettes and noise clocks advance. There are no preset
markers.

### Worst steady window (frames 257-288)

```
frame                     62.356 ms 37.414 Mcyc 100.0%
  pov_preserve_half        0.141 ms  0.085 Mcyc   0.2%
  fx_shader_draw          33.676 ms 20.205 Mcyc  54.0%
    fx_prepare_column      1.637 ms  0.982 Mcyc   2.6% x144 11.37 us/col
  fx_palette_base_bake     1.064 ms  0.638 Mcyc   1.7%
    fx_palette_rgb_bake    0.936 ms  0.562 Mcyc   1.5%
  fx_prepare_frame         3.439 ms  2.064 Mcyc   5.5%
  fx_advance               0.980 ms  0.588 Mcyc   1.6%
    palette_cycle_step     0.951 ms  0.570 Mcyc   1.5%
  fx_timeline_step         0.083 ms  0.050 Mcyc   0.1%
  canvas_clear             0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait      22.634 ms 13.580 Mcyc  36.3%
```

Wall min/mean/max: **56.298/62.356/68.332 ms**; render mean/max
39.722/44.038 ms.

### Per-pixel figures

At 10,368 shaded pixels per frame, the 33.7 ms shader draw costs about
1,950 cycles per pixel including column preparation.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.52/1.98/25.11 us  cpu 3.65%
isr_pack         144/frame  min/avg/max 8.60/9.29/11.85 us  cpu 2.13%
isr_dma_submit   144/frame  min/avg/max 0.59/0.93/6.19 us  cpu 0.21%
```

- The wire transfer itself runs asynchronously on DMA.
- The inclusive `isr_wake` share is 3.65%, leaving about 60.2 ms of each
  62.5 ms window for render. Pack and submit run nested inside wake and are not
  added to it. DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 54.0% of the frame, 33.7 ms: per-pixel displacement,
   grid and color.
2. `fx_prepare_column` — 2.6%, 1.6 ms: 144 interval samples and stage
   preparations per frame.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- The per-pixel potential-field kernels and the column phase tail run from flash.
- No dwell-compression knobs were used.
- Captured from the landing worktree: master `e314aed39` plus the playlist commit, with no uncommitted changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GridWarp`,
`HS_PROFILE_WINDOW=32`; `just profile GridWarp` = build + flash + capture.
