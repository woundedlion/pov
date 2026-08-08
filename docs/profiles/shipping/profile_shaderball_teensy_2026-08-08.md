# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the one-sample lens capture;
the smooth two-sample output crossfade is restored while direct gamut-direction
clipping remains active.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; the closure `Scan::Shader::draw` selective-O3 region is active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `b9e121de` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:63008, data:147624, headers:8504` / `RAM1:
variables:314976, code:30152, padding:2616, free:176544` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 1233-1248 root counter cycles divided by
600 MHz match measured wall time within **1.7 ppm**.

## Frame cadence

`sb_shader_draw` averages **54.44 ms/frame** and total render averages
**61.43 ms/frame**. Peak frame render is **94.38 ms** and **1456/4912 frames
(29.6%)** spill. Presets 1, 4, 5, and 6 spill; the other eight presets spill
nothing and peak at 61.80 ms or less.

A display window is 62.5 ms. ShaderBall shades one 72x144 quadrant, or 10,368
pixels. The `canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.45 ms  37.47 Mcyc  100%
  sb_shader_draw       32.85 ms  19.71 Mcyc   52%
  sb_timeline_step      69 us    41.5 kcyc     0%
  canvas_clear          87 us    52.1 kcyc     0%
  canvas_buffer_wait   22.70 ms  13.62 Mcyc   36%
```

Render averages 39.76 ms and peaks at 40.30 ms in this window, leaving ample
16 fps margin.

### Expensive lens crossfade (frames 1265-1280, worst of capture)

```
frame                 125.06 ms  75.04 Mcyc  100%
  sb_shader_draw       86.04 ms  51.63 Mcyc   68%
  sb_timeline_step      72 us    43.4 kcyc     0%
  canvas_clear          85 us    51.1 kcyc     0%
  canvas_buffer_wait   31.94 ms  19.16 Mcyc   25%
```

Wall min/avg/max is 123.29/125.06/126.34 ms. Render averages 93.12 ms and
peaks at 94.38 ms. Fractional lens morphs evaluate both directions and blend
their pattern outputs, preserving the smooth authored transition.

### Per-preset table

Initial unlabeled frames are folded into preset 1. Peak render uses exact
per-frame ownership from `parse_profile.py ... buckets`.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 6 | 94.38 | 8 fps | 480/480 |
| 1 | 93.36 | 8-16 fps | 479/545 |
| 5 | 78.61 | 8 fps | 479/480 |
| 4 | 63.05 | 8-16 fps | 18/136 |
| 10 | 61.80 | 16 fps | 0/480 |
| 7 | 61.72 | 16 fps | 0/480 |
| 12 | 61.70 | 16 fps | 0/480 |
| 11 | 57.86 | 16 fps | 0/480 |
| 9 | 53.15 | 16 fps | 0/480 |
| 8 | 52.63 | 16 fps | 0/480 |
| 3 | 40.30 | 16 fps | 0/145 |
| 2 | 40.22 | 16 fps | 0/158 |

### Per-pixel figures

The worst shader window costs about **8.30 us/pixel**, or **4,979
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        2306/frame  min/avg/max 0.67/1.94/11.46 us  cpu 3.58%
isr_pack         288/frame  min/avg/max 6.23/6.99/9.28 us   cpu 1.61%
isr_dma_submit   288/frame  min/avg/max 0.65/0.94/1.11 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. LED transfer continues
asynchronously, and `filter_blend` is not exposed on this closure-shader path.

## Summary ranking

1. `sb_shader_draw` - 86.04 ms/frame in the worst lens crossfade.
2. Palette-cycle and frame work outside the shader - about 7.1 ms/frame.
3. Timeline and clear - below 0.2 ms/frame combined.

The direct gamut-direction optimization lowers this peak from the 97.93 ms
2026-08-07 baseline without changing the lens transition.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
