# Involution on-device profile — Teensy 4.0, segmented mode (2026-11-03, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Involution`).
Raw capture: `build/prof/involution_ship.log`, captured **2026-11-03 14:44**
America/Los_Angeles on COM4. First shipping report for Involution.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM shader, noise and color kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Involution 288×144, authored `jewel-melt` preset, source `18601fa2d` |
| Method | `HS_PROFILE`, 32-frame windows, 70 seconds, default epoch; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere bash tools/profile_one.sh Involution profile 70 32` |

Profile image: `FLASH: code:116904, data:159020, headers:8748`;
`RAM1: variables:314112, code:29592, padding:3176, free:177408`;
`RAM2: variables:520064, free:4224`.

Full Phantasm image with Involution in the playlist: `FLASH: code:1087456,
data:788692`; `RAM1: variables:313760, code:195240`; `RAM2: variables:520064`.
The size and layout gate passes with 1,368 B to the RAM1 code ceiling.

| Region | Master `c261d0b25` | Final | Delta |
|---|--:|--:|--:|
| RAM1 code | 195,160 | 195,240 | +80 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH code | 1,063,048 | 1,087,456 | +24,408 |
| FLASH data | 780,944 | 788,692 | +7,748 |

Exactness cross-check: root cycles/600 MHz differ from the measured wall sum
by **10.8 ppm** in frames 1-32 (8–12 ppm in every window), above the
validator's 5 ppm bound. Frame starts are phase-locked to the display window,
so `micros()` truncation at each frame end biases wall time by up to 1 µs per
frame; the offset moved with code placement and appears identically on COM3
and COM4. An earlier pass of the same shader code (`83f2548fd`, before the
cold-code placement commit) validated at 3.4 ppm with a 48.306 ms steady peak.

## Frame cadence

**Steady peak: 48.380 ms**, frame 934; mean **46.368 ms**,
**0/1087 spilled (0.00%)**, and **0/1087 at or above 58 ms**.
Frames 2-1088 in complete windows are included. Setup frame 1 is
**84.427 ms** and spills separately.

README cells: peak 🟢 48.380, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. Involution shades one 10,368-pixel
quadrant per frame and holds the 16 fps cadence. Setup frame 1 renders the
full canvas without a preserved half, so it costs about twice a steady frame.
`canvas_buffer_wait` is the display-sync idle and is subtracted from render.

Before this change the same preset measured a 79.774 ms steady peak and a
77.370 ms mean, with 351/351 steady frames spilled (`c261d0b25`, COM3).

## Phase-by-phase readout

The single authored preset runs continuously: the camera and projection walk,
the periodic affine lattice turns, and palettes and noise clocks advance. There
are no preset markers.

### Startup-containing window (frames 1-32)

```
frame                     62.578 ms 37.547 Mcyc 100.0%
  pov_preserve_half        0.134 ms  0.080 Mcyc   0.2%
  fx_shader_draw          45.722 ms 27.433 Mcyc  73.1%
    fx_prepare_column      2.942 ms  1.765 Mcyc   4.7% x144 20.43 us/col
  fx_palette_base_bake     0.093 ms  0.056 Mcyc   0.1%
  fx_prepare_frame         0.121 ms  0.073 Mcyc   0.2%
  fx_advance               1.065 ms  0.639 Mcyc   1.7%
    palette_cycle_step     0.928 ms  0.557 Mcyc   1.5%
  fx_timeline_step         0.088 ms  0.053 Mcyc   0.1%
  canvas_clear             0.095 ms  0.057 Mcyc   0.2%
  canvas_buffer_wait      14.952 ms  8.971 Mcyc  23.9%
```

Wall min/mean/max: **44.346/62.577/84.427 ms**. Scopes are inclusive; child
rows overlap their parent. The window holds the full-canvas setup frame.

### Worst steady window (frames 641-672)

```
frame                     62.417 ms 37.450 Mcyc 100.0%
  pov_preserve_half        0.137 ms  0.082 Mcyc   0.2%
  fx_shader_draw          44.902 ms 26.941 Mcyc  71.9%
    fx_prepare_column      3.102 ms  1.861 Mcyc   5.0% x144 21.54 us/col
  fx_palette_base_bake     0.092 ms  0.055 Mcyc   0.1%
  fx_prepare_frame         0.031 ms  0.019 Mcyc   0.0%
  fx_advance               1.020 ms  0.612 Mcyc   1.6%
    palette_cycle_step     0.940 ms  0.564 Mcyc   1.5%
  fx_timeline_step         0.082 ms  0.049 Mcyc   0.1%
  canvas_clear             0.091 ms  0.055 Mcyc   0.1%
  canvas_buffer_wait      15.753 ms  9.452 Mcyc  25.2%
```

Wall min/mean/max: **59.737/62.416/65.072 ms**; render mean/max
46.664/48.063 ms. Shader evaluation dominates. Per-column preparation copies
the interval frame state and prepares every stage 144 times per frame.
`palette_cycle_step` alternates a display rebake every second frame.

### Per-pixel figures

At 10,368 shaded pixels per frame, the 44.9 ms shader draw costs about
2,600 cycles per pixel including column preparation. A PC-sampling pass of the
same shader ranks the per-pixel work as the shader body (rotation, S2
projection, periodic warp and coverage, ~24%), the direct simplex lift noise
(~12%), S2 projection (~8%), color finishing with gamut mapping (~12%) and the
lattice primitive (~5%).

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.58/2.00/23.26 us  cpu 3.69%
isr_pack         144/frame  min/avg/max 8.19/9.38/12.18 us  cpu 2.16%
isr_dma_submit   144/frame  min/avg/max 0.65/0.94/6.12 us  cpu 0.21%
```

- Packing costs about ten times the DMA submission per column.
- The wire transfer itself runs asynchronously on DMA.
- The inclusive `isr_wake` share is 3.69%, leaving about 60.2 ms of each
  62.5 ms window for render; the 48.4 ms steady peak needs no further speedup.
  Pack and submit run nested inside wake and are not added to it.
  DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 71.9% of the frame, 44.9 ms: per-pixel shading.
2. `fx_prepare_column` — 5.0%, 3.1 ms: 144 interval samples and stage
   preparations per frame.
3. `fx_advance` / `palette_cycle_step` — 1.6%, 1.0 ms: palette rebakes on
   alternate frames.

The speedup over the 79.8 ms baseline came from four changes. Modulation is
evaluated once per pixel at the column phase instead of at both temporal
endpoints. Palette endpoints blend through an interval OKLab cache, which also
serves frames whose color endpoints coincide. The periodic warp's displacement
uses centered remainders instead of libm `remainderf`/`hypotf`. The affine and
S2 kernels multiply by prepared reciprocals.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- No `filter_blend` scope is active for this shader.
- The shader, simplex noise, lattice and color finishing run in the selective
  `HS_O3`/ITCM regions; the shader body and S2 projection run from flash.
- No dwell-compression knobs were used.
- The capture tree carried the uncommitted playlist and documentation edits
  that accompany this report; the profile image does not include them.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Involution`,
`HS_PROFILE_WINDOW=32`; `just profile Involution` = build + flash + capture.
