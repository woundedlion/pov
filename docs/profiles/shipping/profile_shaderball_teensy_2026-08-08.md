# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the smooth two-sample lens
capture: lens transitions now morph projected coordinates and evaluate the
warp/pattern once, and ShaderBall uses one cycling liquid palette.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; the closure `Scan::Shader::draw` selective-O3 region is active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, capture tip `e5992197` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:62632, data:147628, headers:8876` / `RAM1:
variables:314976, code:29992, padding:2776, free:176544` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 1249-1264 root counter cycles divided by
600 MHz match measured wall time within **1.4 ppm**.

## Frame cadence

`sb_shader_draw` averages **45.07 ms/frame** and total render averages
**47.15 ms/frame**. Peak frame render is **63.05 ms** and **3/6368 frames
(0.05%)** spill. Preset 6 owns all three spills; the other 11 presets spill
nothing and peak at 61.81 ms or less.

A display window is 62.5 ms. ShaderBall shades one 72x144 quadrant, or 10,368
pixels. The `canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.46 ms  37.48 Mcyc  100%
  sb_shader_draw       30.15 ms  18.09 Mcyc   48%
  sb_timeline_step      64 us    38.3 kcyc     0%
  canvas_clear          86 us    51.6 kcyc     0%
  canvas_buffer_wait   30.25 ms  18.15 Mcyc   48%
```

Render averages 32.21 ms and peaks at 32.76 ms in this window.

### Expensive lens morph (frames 6017-6032, worst of capture)

```
frame                  74.01 ms  44.40 Mcyc  100%
  sb_shader_draw       57.72 ms  34.63 Mcyc   77%
  sb_timeline_step      60 us    35.9 kcyc     0%
  canvas_clear          89 us    53.2 kcyc     0%
  canvas_buffer_wait   14.18 ms   8.51 Mcyc   19%
```

Wall min/avg/max is 56.73/74.01/125.23 ms. Render averages 59.83 ms and peaks
at 63.05 ms. Fractional lens morphs compute both inexpensive projections,
interpolate their stereographic coordinates, and run the expensive
warp/pattern/palette path exactly once.

### Per-preset table

Initial unlabeled frames are folded into preset 1. Peak render uses exact
per-frame ownership from `parse_profile.py ... buckets`.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 6 | 63.05 | 8-16 fps | 3/960 |
| 1 | 61.81 | 16 fps | 0/633 |
| 7 | 55.46 | 16 fps | 0/657 |
| 10 | 54.88 | 16 fps | 0/480 |
| 12 | 54.56 | 16 fps | 0/480 |
| 11 | 51.37 | 16 fps | 0/480 |
| 5 | 48.64 | 16 fps | 0/960 |
| 9 | 46.05 | 16 fps | 0/480 |
| 8 | 45.32 | 16 fps | 0/480 |
| 4 | 33.94 | 16 fps | 0/286 |
| 2 | 32.87 | 16 fps | 0/232 |
| 3 | 32.76 | 16 fps | 0/240 |

### Per-pixel figures

The worst shader window costs about **5.57 us/pixel**, or **3,340
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        2306/frame  min/avg/max 0.68/1.91/11.17 us  cpu 3.52%
isr_pack         288/frame  min/avg/max 6.23/6.77/8.95 us   cpu 1.56%
isr_dma_submit   288/frame  min/avg/max 0.70/0.94/1.06 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. LED transfer continues
asynchronously, and `filter_blend` is not exposed on this closure-shader path.

## Summary ranking

1. `sb_shader_draw` - 57.72 ms/frame in the worst lens-transition window.
2. Palette-cycle and frame work outside the shader - about 2.08 ms/frame.
3. Timeline and clear - below 0.2 ms/frame combined.

Generated code contains one call to the projected warp/pattern sampler from
the pixel shader. The hot scan body shrank from 1,232 to 512 bytes and total
ITCM fell by 240 bytes versus the two-sample image. Peak render fell from
94.38 to 63.05 ms (**1.50x**).

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- The profile image compresses the epoch only; it does not change per-frame
  effect work.
- The capture ran from the clean isolated implementation commit shown above.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
