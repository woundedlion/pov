# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-08.md).
Raw capture: `build/prof/shaderball_o3.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `85f58fc7` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:87880, data:147668, headers:9188` / `RAM1:
variables:314976, code:52344, padding:13192, free:143776` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 1201-1216 root counter cycles divided by
600 MHz match measured wall time within **0.1 ppm**.

## Frame cadence

`sb_shader_draw` averages **40.24 ms/frame**; total render averages
**47.20 ms/frame**, peaks at **61.28 ms**, and spills **0/6368 frames**. All
12 presets remain at 16 fps. The previous O3 capture peaked at 89.21 ms and
spilled 1437/4928 frames.

The optimized O3 image is cadence-clean, but its 61.28 ms peak is still above
the stricter 59 ms render target.

## Representative windows

### Held liquid regime (frames 321-336)

```
frame                  62.46 ms  37.47 Mcyc  100%
  render               35.78 ms  21.47 Mcyc   57%
  sb_shader_draw       28.91 ms  17.34 Mcyc   46%
  canvas_buffer_wait   26.68 ms  16.01 Mcyc   43%
```

### Worst render window (frames 1201-1216)

```
frame                  62.46 ms  37.48 Mcyc  100%
  render               59.27 ms  35.56 Mcyc   95%
  sb_shader_draw       52.64 ms  31.58 Mcyc   84%
  sb_timeline_step      48 us    28.6 kcyc     0%
  canvas_clear          90 us    53.9 kcyc     0%
  canvas_buffer_wait    3.19 ms   1.92 Mcyc     5%
```

## Per-preset cadence

Initial unlabeled frames are folded into preset 1.

| # | peak render ms | spilled/frames |
|---:|--:|---:|
| 6 | 61.28 | 0/960 |
| 1 | 60.23 | 0/545 |
| 12 | 56.24 | 0/480 |
| 7 | 55.58 | 0/657 |
| 10 | 54.38 | 0/480 |
| 11 | 51.23 | 0/480 |
| 5 | 47.97 | 0/960 |
| 9 | 46.99 | 0/480 |
| 8 | 46.94 | 0/480 |
| 4 | 36.27 | 0/286 |
| 2 | 36.26 | 0/232 |
| 3 | 36.10 | 0/240 |

## ISR figures

```
isr_wake        1152/frame  avg 1.85 us  cpu 3.41%
isr_pack         144/frame  avg 6.91 us  cpu 1.59%
isr_dma_submit   144/frame  avg 0.91 us  cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. Global O3 is a
single-effect compiler ceiling, not a shippable roster image.

## Global O3 vs selective O3

Global O3 lowers the new worst-frame peak from 67.33 to 61.28 ms (**1.10x**)
and aggregate render average from 52.75 to 47.20 ms. It adds **24,872 B FLASH**
and **22,192 B ITCM** to the single-effect image.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
