# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-08.md).
Raw capture: `build/prof/shaderball_o3.log`. It covers the same 13-preset
roster and exact implementation commit as the shipping pass.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, capture tip `17fbf726` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3600`; all 13 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 400 16 "-D HS_PROFILE_EPOCH_REVS=3600"` |

Image size: `FLASH: code:91736, data:147692, headers:8380` / `RAM1:
variables:315008, code:56776, padding:8760, free:143744` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 6289-6304 root counter cycles divided by
600 MHz match measured wall time within **0.3 ppm**.

## Frame cadence

`sb_shader_draw` averages **36.17 ms/frame** and total render averages
**38.25 ms/frame**. Peak frame render is **50.20 ms** and **0/6368 frames**
spill. Every preset holds 16 fps, with 12.30 ms of worst-frame margin.

ShaderBall shades one 72x144 quadrant, or 10,368 pixels. The
`canvas_buffer_wait` scope is cadence-quantized display-sync idle.

## Phase-by-phase readout

### Held liquid regime (frames 321-336)

```
frame                  62.45 ms  37.47 Mcyc  100%
  sb_shader_draw       27.13 ms  16.28 Mcyc   43%
  sb_timeline_step      63 us    37.8 kcyc     0%
  canvas_clear          86 us    51.6 kcyc     0%
  canvas_buffer_wait   33.27 ms  19.96 Mcyc   53%
```

Render averages 29.18 ms and peaks at 29.40 ms in this window.

### Fractional-pattern preset (frames 6289-6304, worst of capture)

```
frame                  62.56 ms  37.53 Mcyc  100%
  sb_shader_draw       47.05 ms  28.23 Mcyc   75%
  sb_timeline_step      52 us    31.1 kcyc     0%
  canvas_clear          88 us    52.6 kcyc     0%
  canvas_buffer_wait   13.42 ms   8.05 Mcyc   21%
```

Wall min/avg/max is 61.11/62.56/64.18 ms. Render averages 49.14 ms and
peaks at 50.20 ms.

### Per-preset table

Initial unlabeled frames are folded into preset 1. Peak render uses exact
per-frame ownership from `parse_profile.py ... buckets`.

| # | peak render ms | cadence | spilled/frames |
|---:|--:|---:|---:|
| 6 | 50.20 | 16 fps | 0/657 |
| 7 | 50.18 | 16 fps | 0/480 |
| 1 | 48.73 | 16 fps | 0/545 |
| 13 | 45.30 | 16 fps | 0/480 |
| 12 | 44.74 | 16 fps | 0/480 |
| 8 | 43.15 | 16 fps | 0/480 |
| 5 | 41.99 | 16 fps | 0/960 |
| 11 | 41.19 | 16 fps | 0/480 |
| 10 | 39.88 | 16 fps | 0/480 |
| 9 | 38.38 | 16 fps | 0/480 |
| 4 | 32.39 | 16 fps | 0/286 |
| 3 | 30.88 | 16 fps | 0/240 |
| 2 | 30.72 | 16 fps | 0/232 |

### Per-pixel figures

The worst shader window costs about **4.54 us/pixel**, or **2,723
cycles/pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        1154/frame  min/avg/max 0.67/1.87/11.45 us  cpu 3.44%
isr_pack         144/frame  min/avg/max 6.34/7.07/9.27 us   cpu 1.63%
isr_dma_submit   144/frame  min/avg/max 0.81/0.91/1.02 us   cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. Global O3 is a
single-effect compiler ceiling, not a shippable roster image.

## Summary ranking

1. `sb_shader_draw` - 47.05 ms/frame in the worst fractional-pattern window.
2. Palette-cycle and frame work outside the shader - about 2.09 ms/frame.
3. Timeline and clear - below 0.2 ms/frame combined.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- The profile image compresses the epoch only; it does not change per-frame
  effect work.
- The capture ran from the clean isolated implementation commit shown above.

## Global O3 vs selective O3

Global O3 lowers peak render from 53.99 to 50.20 ms (**1.08x**) and aggregate
render average from 41.72 to 38.25 ms. It adds **28,480 B FLASH** and **26,320
B ITCM** to the single-effect image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
