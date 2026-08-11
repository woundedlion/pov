# ShaderBall on-device profile - Teensy 4.0, segmented mode (2026-08-10, **-O3**)

Global-O3 twin of the [shipping capture](../shipping/profile_shaderball_teensy_2026-08-10.md).
Raw capture: `build/prof/shaderball_o3.log`. It covers the same 26 authored
presets, dwell compression, hardware, and source commit as the shipping pass.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect compiler ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShaderBall 288x144, single-entry playlist, tip `ffbd4807` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 140 s capture, two-revolution profile-only preset holds, `HS_PROFILE_EPOCH_REVS=1400` |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 140 16 "-D HS_PROFILE_SHADERBALL_FAST_CYCLE -D HS_PROFILE_EPOCH_REVS=1400"` |

Image size: `FLASH: code:144328, data:169972, headers:8260` / `RAM1:
variables:323232, code:97592, padding:712, free:102752` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 865-880 root counter cycles divided by
600 MHz match the measured wall sum within **1.2 ppm**. The capture has 138
windows, monotonically increasing frame numbers, all 26 preset indices, and a
confirmed wrap from preset 25 to preset 0.

## Frame cadence

Peak frame render is **60.99 ms** and **0/2208 frames** spill. All **26/26
presets are green** at the 62.5 ms display cadence. The worst clean hold is
preset 25 at **52.70 ms/frame** in `sb_shader_draw`; preset 13 owns the 60.99
ms worst transition frame.

## Per-preset results

| # | authored regime | shader ms | peak render ms | cadence | spilled/frames |
|---:|---|--:|--:|---:|---:|
| 0 | liquid stereo 0 | 37.91 | 58.06 | green | 0/99 |
| 1 | liquid stereo 1 | 37.96 | 43.37 | green | 0/102 |
| 2 | liquid stereo 2 | 37.98 | 43.34 | green | 0/102 |
| 3 | liquid stereo 3 | 33.65 | 47.75 | green | 0/102 |
| 4 | liquid stereo 4 | 45.14 | 53.25 | green | 0/102 |
| 5 | liquid stereo 5 | 46.15 | 50.87 | green | 0/102 |
| 6 | liquid stereo 6 | 50.03 | 57.45 | green | 0/102 |
| 7 | liquid stereo 7 | 43.83 | 54.78 | green | 0/102 |
| 8 | liquid stereo 8 | 42.37 | 48.25 | green | 0/102 |
| 9 | liquid stereo 9 | 43.76 | 47.17 | green | 0/102 |
| 10 | liquid stereo 10 | 50.21 | 55.82 | green | 0/102 |
| 11 | liquid stereo 11 | 25.07 | 55.75 | green | 0/102 |
| 12 | liquid stereo 12 | 35.39 | 35.72 | green | 0/102 |
| 13 | liquid stereo 13 | 49.66 | 60.99 | green | 0/69 |
| 14 | liquid stereo 14 | 49.52 | 55.37 | green | 0/68 |
| 15 | diagnostic 0 | 49.51 | 55.45 | green | 0/68 |
| 16 | diagnostic 1 | 39.37 | 55.67 | green | 0/68 |
| 17 | diagnostic 2 | 47.58 | 52.80 | green | 0/68 |
| 18 | diagnostic 3 | 46.01 | 53.09 | green | 0/68 |
| 19 | diagnostic 4 | 37.42 | 52.17 | green | 0/68 |
| 20 | wave shear liquid | 25.33 | 45.39 | green | 0/68 |
| 21 | kaleidoscope mirror | 26.36 | 31.75 | green | 0/68 |
| 22 | gnomonic grid + kaleidoscope | 26.30 | 31.93 | green | 0/68 |
| 23 | gnomonic grid + glitch | 19.47 | 31.69 | green | 0/68 |
| 24 | Bonne primitive lattice | 45.41 | 50.36 | green | 0/68 |
| 25 | Peirce primitive lattice | 52.70 | 58.20 | green | 0/68 |

## Global O3 vs selective O3

Global O3 lowers the worst clean-hold shader cost from 57.37 to 52.70 ms and
the pass peak render from 67.71 to 60.99 ms. It changes the result from 23/26
green with 31 spills to 26/26 green with zero spills. In this single-effect
image it adds **46,584 B FLASH** and **43,776 B ITCM**.

Global O3 is a compiler ceiling, not a shippable full-roster image. All scopes
absorb ISR time because CYCCNT free-runs, and the compressed dwell changes only
how quickly the capture visits presets.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShaderBall`,
`HS_PROFILE_WINDOW=16`.


