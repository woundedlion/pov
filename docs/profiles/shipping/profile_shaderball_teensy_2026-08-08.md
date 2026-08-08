# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-08, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShaderBall`). Raw capture:
`build/prof/shaderball_ship.log`. This replaces the 2026-08-07 capture after
the one-sample lens morph and direct gamut-direction optimizations.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile`: `-Os`, newlib-nano, DMA LEDs, `HS_PROFILE_ENABLE`; the closure `Scan::Shader::draw` selective-O3 region is active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `85f58fc7` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 400 s capture, `HS_PROFILE_EPOCH_REVS=3300`; all 12 presets and the wrap to preset 1 are present |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 400 16 "-D HS_PROFILE_EPOCH_REVS=3300"` |

Image size: `FLASH: code:63008, data:147624, headers:8504` / `RAM1:
variables:314976, code:30152, padding:2616, free:176544` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 1233-1248 root counter cycles divided by
600 MHz match measured wall time within **0.8 ppm**
(`tools/parse_profile.py build/prof/shaderball_ship.log validate`).

## Frame cadence

`sb_shader_draw` averages **45.75 ms/frame**; total render averages
**52.75 ms/frame**, peaks at **67.33 ms**, and spills **536/5824 frames
(9.2%)**. Presets 1 and 6 spill; the other ten presets spill nothing and peak
at 61.73 ms or less.

The previous capture peaked at 97.93 ms and spilled 1501/4864 frames. The two
changes recover 30.60 ms from the peak and reduce the aggregate spill rate from
30.9% to 9.2%, but the shipping image still misses the 59 ms target.

## Representative windows

### Held liquid regime (frames 321-336)

```
frame                  62.46 ms  37.47 Mcyc  100%
  sb_shader_draw       32.67 ms  19.60 Mcyc   52%
  sb_timeline_step      68 us    40.7 kcyc     0%
  canvas_clear          87 us    52.3 kcyc     0%
  canvas_buffer_wait   22.88 ms  13.73 Mcyc   37%
```

### Worst render window (frames 4561-4576)

```
frame                 116.85 ms  70.11 Mcyc  100%
  render               64.11 ms  38.47 Mcyc   55%
  sb_shader_draw       57.02 ms  34.21 Mcyc   49%
  sb_timeline_step      59 us    35.3 kcyc     0%
  canvas_clear          86 us    51.4 kcyc     0%
  canvas_buffer_wait   52.74 ms  31.65 Mcyc   45%
```

The approximately 7.1 ms between the shader scope and total render is mostly
palette-cycle maintenance outside `sb_shader_draw`. That is now the clearest
remaining route toward the 59 ms target.

## Per-preset cadence

Initial unlabeled frames are folded into preset 1. Peak render is based on
exact per-frame ownership from `parse_profile.py ... buckets`.

| # | family | peak render ms | cadence | spilled/frames |
|---:|---|--:|---:|---:|
| 1 | liquid mild / wrap blend | 67.33 | 8-16 fps | 302/545 |
| 6 | grid orbit A | 67.14 | 8-16 fps | 234/593 |
| 12 | grid/liquid B | 61.73 | 16 fps | 0/480 |
| 10 | grid orbit E | 61.71 | 16 fps | 0/480 |
| 7 | grid orbit B | 61.44 | 16 fps | 0/480 |
| 11 | grid/liquid A | 57.81 | 16 fps | 0/480 |
| 5 | liquid-to-grid bridge | 53.94 | 16 fps | 0/960 |
| 9 | grid orbit D | 53.31 | 16 fps | 0/480 |
| 8 | grid orbit C | 52.32 | 16 fps | 0/480 |
| 3 | liquid fine | 39.91 | 16 fps | 0/240 |
| 2 | liquid deep | 39.80 | 16 fps | 0/232 |
| 4 | liquid direct | 39.72 | 16 fps | 0/286 |

## Per-pixel and ISR figures

The worst shader window shades 10,368 pixels/frame at about **5.50 us/pixel**,
or **3,300 cycles/pixel**.

```
isr_wake        2155/frame  avg 1.93 us  cpu 3.55%
isr_pack         269/frame  avg 7.08 us  cpu 1.63%
isr_dma_submit   269/frame  avg 0.94 us  cpu 0.21%
```

All scopes absorb ISR time because CYCCNT free-runs. LED transfer continues
asynchronously, and `filter_blend` is not exposed on this closure-shader path.

## Code-generation check

The shipping disassembly contains one call to ShaderBall's expensive
`sample_direction` lambda in the scan loop. `gamut_clip_preserve_chroma` calls
the normalized-direction `gamut_max_chroma(float, float, float)` overload
directly, with no `atan2f`, `cosf`, or `sinf` round trip in the clip path.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.
