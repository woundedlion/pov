# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-08.md).
Raw capture: `build/prof/shaderball_o3.log`. Lens transitions morph projected
coordinates and evaluate the warp/pattern once, and ShaderBall uses one
cycling liquid palette.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, capture tip `e5992197` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:91400, data:147672, headers:8736` / `RAM1:
variables:315008, code:56600, padding:8936, free:143744` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 4625-4640 root counter cycles divided by
600 MHz match measured wall time within **0.4 ppm**.

## Frame cadence

`sb_shader_draw` averages **39.57 ms/frame** and total render averages
**41.64 ms/frame**. Peak frame render is **56.42 ms** and **0/6368 frames**
spill. Every preset holds 16 fps.

A display window is 62.5 ms. ShaderBall shades one 72x144 quadrant, or 10,368
pixels. The `canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.46 ms  37.47 Mcyc  100%
  sb_shader_draw       27.21 ms  16.32 Mcyc   43%
  sb_timeline_step      56 us    33.5 kcyc     0%
  canvas_clear          86 us    51.4 kcyc     0%
  canvas_buffer_wait   33.20 ms  19.92 Mcyc   53%
```

Render averages 29.25 ms and peaks at 29.55 ms in this window.

### Worst lens morph (frames 4561-4576)

```
frame                  62.08 ms  37.25 Mcyc  100%
  sb_shader_draw       50.17 ms  30.10 Mcyc   80%
  sb_timeline_step      56 us    33.5 kcyc     0%
  canvas_clear          88 us    52.8 kcyc     0%
  canvas_buffer_wait    9.81 ms   5.89 Mcyc   15%
```

Wall min/avg/max is 60.79/62.08/63.88 ms. Render averages 52.27 ms and peaks
at 56.42 ms, leaving 6.08 ms of peak margin below the 62.5 ms display window.

### Per-preset table

Initial unlabeled frames are folded into preset 1. Peak render uses exact
per-frame ownership from `parse_profile.py ... buckets`.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 1 | 56.42 | 16 fps | 0/633 |
| 6 | 55.90 | 16 fps | 0/960 |
| 7 | 49.61 | 16 fps | 0/657 |
| 12 | 49.43 | 16 fps | 0/480 |
| 10 | 48.31 | 16 fps | 0/480 |
| 11 | 46.32 | 16 fps | 0/480 |
| 5 | 43.12 | 16 fps | 0/960 |
| 8 | 40.25 | 16 fps | 0/480 |
| 9 | 40.10 | 16 fps | 0/480 |
| 4 | 31.39 | 16 fps | 0/286 |
| 2 | 29.91 | 16 fps | 0/232 |
| 3 | 29.55 | 16 fps | 0/240 |

### Per-pixel figures

The worst shader window costs about **4.84 us/pixel**, or **2,903
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        2290/frame  min/avg/max 0.68/1.87/11.44 us  cpu 3.44%
isr_pack         286/frame  min/avg/max 6.32/7.05/9.36 us   cpu 1.62%
isr_dma_submit   286/frame  min/avg/max 0.75/0.91/1.02 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. Global O3 is a
single-effect compiler ceiling, not a shippable roster image.

## Summary ranking

1. `sb_shader_draw` - 50.17 ms/frame in the worst lens-transition window.
2. Palette-cycle and frame work outside the shader - about 2.07 ms/frame.
3. Timeline and clear - below 0.2 ms/frame combined.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- The profile image compresses the epoch only; it does not change per-frame
  effect work.
- The capture ran from the clean isolated implementation commit shown above.

## Global O3 vs selective O3

Global O3 lowers peak render from 63.05 to 56.42 ms (**1.12x**) and aggregate
render average from 47.15 to 41.64 ms. It adds **28,768 B FLASH** and **26,608
B ITCM** to the single-effect image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
