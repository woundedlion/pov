# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-08.md).
Raw capture: `build/prof/shaderball_o3.log`. The smooth two-sample lens
crossfade is restored while direct gamut-direction clipping remains active.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `b9e121de` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:92456, data:147668, headers:8708` / `RAM1:
variables:314976, code:56920, padding:8616, free:143776` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 4625-4640 root counter cycles divided by
600 MHz match measured wall time within **1.3 ppm**.

## Frame cadence

`sb_shader_draw` averages **46.96 ms/frame** and total render averages
**53.94 ms/frame**. Peak frame render is **86.69 ms** and **1437/4928 frames
(29.2%)** spill. Presets 1, 5, and 6 spill; the other nine presets remain at
16 fps and peak at 56.93 ms or less.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.46 ms  37.48 Mcyc  100%
  sb_shader_draw       28.53 ms  17.12 Mcyc   45%
  sb_timeline_step      63 us    37.8 kcyc     0%
  canvas_clear          86 us    51.9 kcyc     0%
  canvas_buffer_wait   27.04 ms  16.22 Mcyc   43%
```

Render averages 35.42 ms and peaks at 35.71 ms.

### Expensive lens crossfade (frames 4625-4640, worst of capture)

```
frame                 124.99 ms  75.00 Mcyc  100%
  sb_shader_draw       79.34 ms  47.61 Mcyc   63%
  sb_timeline_step      50 us    30.3 kcyc     0%
  canvas_clear          85 us    51.0 kcyc     0%
  canvas_buffer_wait   38.61 ms  23.17 Mcyc   30%
```

Wall min/avg/max is 124.78/124.99/125.12 ms. Render averages 86.38 ms and
peaks at 86.69 ms.

### Per-preset table

Initial unlabeled frames are folded into preset 1.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 1 | 86.69 | 8-16 fps | 479/545 |
| 6 | 85.87 | 8 fps | 480/480 |
| 5 | 69.60 | 8 fps | 478/480 |
| 4 | 56.93 | 16 fps | 0/136 |
| 12 | 55.68 | 16 fps | 0/480 |
| 7 | 53.99 | 16 fps | 0/480 |
| 10 | 53.85 | 16 fps | 0/480 |
| 11 | 50.92 | 16 fps | 0/480 |
| 9 | 46.12 | 16 fps | 0/480 |
| 8 | 46.09 | 16 fps | 0/480 |
| 2 | 36.35 | 16 fps | 0/174 |
| 3 | 35.74 | 16 fps | 0/145 |

### Per-pixel figures

The worst shader window costs about **7.65 us/pixel**, or **4,592
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        2305/frame  min/avg/max 0.64/1.87/11.22 us  cpu 3.45%
isr_pack         288/frame  min/avg/max 6.32/6.99/9.40 us   cpu 1.61%
isr_dma_submit   288/frame  min/avg/max 0.75/0.91/1.02 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. Global O3 is a
single-effect compiler ceiling, not a shippable roster image.

## Global O3 vs selective O3

Global O3 lowers peak render from 94.38 to 86.69 ms (**1.09x**) and aggregate
render average from 61.43 to 53.94 ms. It adds **29,448 B FLASH** and **26,768
B ITCM** to the single-effect image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
