# Involution on-device profile — Teensy 4.0, segmented mode (2026-11-04, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Involution`).
Raw capture: `build/prof/involution_ship.log` in the capture worktree,
captured **2026-11-04 08:42** America/Los_Angeles on COM4.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM shader, noise and color kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Involution 288×144, authored `jewel-melt` preset, cubic atlas at Area Correction 0.7, source `662de1922` |
| Method | `HS_PROFILE`, 32-frame windows, 70 seconds, default epoch; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=<tree> bash tools/profile_one.sh Involution profile 70 32` |

Profile image: `FLASH: code:119744, data:159840, headers:9184`;
`RAM1: variables:314144, code:29448, padding:3320, free:177376`;
`RAM2: variables:520064, free:4224`.

Full Phantasm image with Involution in the playlist: `FLASH: code:1130104,
data:788924`; `RAM1: variables:313760, code:191960`; `RAM2: variables:520064`.
The size and layout gate passes.

| Region | Master `e314aed39` | Final | Delta |
|---|--:|--:|--:|
| RAM1 code | 191,960 | 191,960 | 0 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH code | 1,130,216 | 1,130,104 | −112 |
| FLASH data | 788,920 | 788,924 | +4 |

Exactness cross-check: root cycles/600 MHz differ from the measured wall sum
by **5.0 ppm** in frames 1-32, at the validator's 5 ppm bound.

## Frame cadence

**Steady peak: 44.986 ms**, frame window 33-64; mean **43.375 ms**,
**0/1087 spilled (0.00%)**, and **0/1087 at or above 58 ms**.
Frames 2-1088 in complete windows are included. Setup frame 1 is
**82.342 ms** and spills separately.

README cells: peak 🟢 44.986, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. Involution shades one 10,368-pixel
quadrant per frame and holds the 16 fps cadence. Setup frame 1 renders the
full canvas without a preserved half. `canvas_buffer_wait` is the
display-sync idle and is subtracted from render.

The atlas kernel picks a face from half the selectors on centrally symmetric
solids, tests a fixed side count, and derives Area Correction from the face
cosine with one divide and one square root. Before that change the atlas
measured a 53.399 ms steady peak (`15d4750f8`, COM3), and S2 measured
48.380 ms (`18601fa2d`, COM4). The earlier kernel at `225c1678a` sampled at
48.5–49.9 ms per window under the PC sampler.

## Phase-by-phase readout

The single authored preset runs continuously: the camera and projection walk,
the periodic affine lattice turns, and palettes and noise clocks advance. There
are no preset markers.

### Worst steady window (frames 33-64)

```
frame                     62.392 ms 37.435 Mcyc 100.0%
  pov_preserve_half        0.137 ms  0.082 Mcyc   0.2%
  fx_shader_draw          41.751 ms 25.051 Mcyc  66.9%
    fx_prepare_column      1.964 ms  1.179 Mcyc   3.1% x144 13.64 us/col
  fx_palette_base_bake     0.089 ms  0.053 Mcyc   0.1%
  fx_prepare_frame         0.024 ms  0.014 Mcyc   0.0%
  fx_advance               1.043 ms  0.626 Mcyc   1.7%
    palette_cycle_step     0.975 ms  0.585 Mcyc   1.6%
  fx_timeline_step         0.052 ms  0.031 Mcyc   0.1%
  canvas_clear             0.095 ms  0.057 Mcyc   0.2%
  canvas_buffer_wait      18.930 ms 11.358 Mcyc  30.3%
```

Wall min/mean/max: **59.779/62.391/64.889 ms**; render mean/max
43.462/44.986 ms.

### Per-pixel figures

At 10,368 shaded pixels per frame, the 41.8 ms shader draw costs about
2,420 cycles per pixel including column preparation.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.50/1.98/23.83 us  cpu 3.64%
isr_pack         144/frame  min/avg/max 8.23/9.33/11.71 us  cpu 2.14%
isr_dma_submit   144/frame  min/avg/max 0.60/0.94/5.36 us  cpu 0.21%
```

- Packing costs about ten times the DMA submission per column.
- The wire transfer itself runs asynchronously on DMA.
- The inclusive `isr_wake` share is 3.64%, leaving about 60.2 ms of each
  62.5 ms window for render; the 45.0 ms steady peak needs no further speedup.
  Pack and submit run nested inside wake and are not added to it.
  DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 66.9% of the frame, 41.8 ms: per-pixel shading.
2. `fx_prepare_column` — 3.1%, 2.0 ms: 144 interval samples and stage
   preparations per frame.
3. `fx_advance` / `palette_cycle_step` — 1.7%, 1.0 ms: palette rebakes on
   alternate frames.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- No `filter_blend` scope is active for this shader.
- The shader, simplex noise, lattice and color finishing run in the selective
  `HS_O3`/ITCM regions; the shader body and atlas projection run from flash.
- No dwell-compression knobs were used.
- Captured from the landing worktree at `662de1922` with no uncommitted
  changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Involution`,
`HS_PROFILE_WINDOW=32`; `just profile Involution` = build + flash + capture.
